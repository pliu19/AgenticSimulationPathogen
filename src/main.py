#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Main simulation loop for the ESKAPE multi-pathogen ICU ABM.

Two trial arms are supported:
  ADE     — Antimicrobial De-Escalation: at 72 h switch to the narrowest
             agent that covers the identified pathogen.
  control — Standard care: stay on broad-spectrum (level 4) after 72-h
             lab results unless Reserve (level 5) is strictly required.

Usage examples:
    python main.py --arm ADE     --num_runs 1000
    python main.py --arm control --num_runs 1000
"""

import os
import copy
import argparse
import math
import time
import pickle

import scipy.integrate as integrate
import scipy.stats as stats

from abm import HealthCareWorker
from config import PATHOGENS, DRUG_LEVELS
from scheme_func import (
    draw, initial_patient, random_schedule, bulk_initialization,
    ph_interaction, hcw_cleanUp,
    infection_development, increase_treatment_icu_time,
    drug_change_control, drug_change_ADE,
    treatment_completion, discharge_admission, death_event,
    intrinsic_mutation,
)


# ─── Utility functions ────────────────────────────────────────────────────────

def truncnorm_func(mean, std_rate):
    """Sample from a truncated normal distribution centred on mean."""
    if std_rate == 0 or mean == 0:
        return mean
    bounded = abs(mean) * std_rate
    std     = abs(mean) * std_rate
    lower   = (mean - bounded - mean) / std   # = -1
    upper   = (mean + bounded - mean) / std   # = +1
    X = stats.truncnorm(lower, upper, loc=mean, scale=std)
    return max(1e-5, X.rvs())


def _constant(x, lam):
    return lam


def read_csv(file_path):
    """
    Read a day → probability CSV and return a dict of survival probabilities
    computed via the hazard-integral method used in the original model.
    """
    record = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                day_idx = int(parts[0])
                prob    = float(parts[1])
                record[day_idx] = prob

    result = {}
    for key, value in record.items():
        I = integrate.quad(_constant, 0, key, args=(value,))[0]
        result[key] = value * math.exp(-I)
    return result


def recordEverything(data, arm, file_path):
    os.makedirs(file_path, exist_ok=True)
    timestr  = time.strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(file_path, f'{arm}_{timestr}.pkl')
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


# ─── Metrics initialisation ───────────────────────────────────────────────────

def make_reference_results():
    """
    Build the flat counter dict that accumulates simulation events.
    Keys are generated from the ESKAPE pathogen registry so adding a new
    species/phenotype to config.py automatically extends tracking.
    """
    r = {
        'admission': 0, 'discharge': 0, 'death': 0,
        'cumsuperinfection': 0, 'misempiric': 0, 'tempempiric': 0,
        'deescalation': 0, 'escalation': 0,
        'cuminfection_none': 0,
    }
    for level in DRUG_LEVELS:
        r[f'drug_level_{level}_use'] = 0

    for sp, cfg in PATHOGENS.items():
        for ph in cfg['phenotypes']:
            r[f'cuminfection_{sp}_{ph}']  = 0
            r[f'labresult_{sp}_{ph}']     = 0
            r[f'colonization_{sp}_{ph}']  = 0
            r[f'transmission_{sp}_{ph}']  = 0
            r[f'mutation_{sp}_{ph}']      = 0

    return r


# ─── Main simulation loop ─────────────────────────────────────────────────────

def main(args, drug_change_func):
    """
    Run one complete simulation of args.days days.
    Returns final_results: a dict of daily time-series and cumulative counters.
    """
    reference_results = make_reference_results()

    # Daily census arrays
    final_results = {'superinfection': [0] * args.days}
    for sp, cfg in PATHOGENS.items():
        for ph in cfg['phenotypes']:
            final_results[f'infection_{sp}_{ph}'] = [0] * args.days
    for key in reference_results:
        final_results[key] = []

    assert args.num_patient % args.num_hcw == 0, \
        "num_patient must be divisible by num_hcw"

    death_probs     = read_csv(args.death_prob_file)
    discharge_probs = read_csv(args.discharge_prob_file)

    patient_list    = bulk_initialization(args)
    patient_idx     = args.num_patient - 1
    current_p_names = list(range(args.num_patient))

    hcw_list = {f'h{i}': HealthCareWorker() for i in range(args.num_hcw)}

    current_hour = 0
    current_day  = 0

    for _ in range(args.days):

        # ── Three 8-hour shifts ───────────────────────────────────────────
        for _shift in range(3):
            schedule  = random_schedule(current_p_names, args.num_hcw)
            visit_idx = 0

            for _ in range(args.num_patient // args.num_hcw):
                for hcw_idx in range(args.num_hcw):
                    hcw_id    = f'h{hcw_idx}'
                    pt_key    = schedule[hcw_id][visit_idx]
                    hcw_list[hcw_id], patient_list[pt_key] = ph_interaction(
                        patient_list[pt_key], hcw_list[hcw_id],
                        args, current_day, current_hour, reference_results)

                visit_idx    += 1
                current_hour += args.time_interval

            # End-of-shift: guaranteed HCW decontamination
            for hid in hcw_list:
                hcw_list[hid] = hcw_cleanUp(
                    hcw_list[hid], current_hour, args.eta, end_of_shift=True)

        # ── End-of-day events ────────────────────────────────────────────
        patient_list = infection_development(
            patient_list, current_p_names, current_day, reference_results)

        patient_list = intrinsic_mutation(
            patient_list, current_p_names, args, current_hour, reference_results)

        patient_list = increase_treatment_icu_time(
            patient_list, current_p_names)

        patient_list = drug_change_func(
            reference_results, patient_list, current_p_names, current_hour)

        patient_list = treatment_completion(
            reference_results, patient_list, current_p_names, current_hour)

        patient_list, current_p_names, patient_idx = discharge_admission(
            patient_list, current_hour, args, discharge_probs,
            current_p_names, patient_idx, reference_results)

        patient_list, current_p_names, patient_idx = death_event(
            patient_list, current_p_names, patient_idx, current_hour,
            args, death_probs, reference_results)

        # ── Daily snapshot ───────────────────────────────────────────────
        for key in reference_results:
            final_results[key].append(reference_results[key])

        for pt_key in current_p_names:
            patient = patient_list[pt_key]
            if patient.status == 'I' and patient.species != 'none':
                ikey = f'infection_{patient.species}_{patient.phenotype}'
                final_results[ikey][current_day] += 1
            if patient.super_infe:
                final_results['superinfection'][current_day] += 1

        current_day += 1

    return final_results


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='ESKAPE ICU ABM — multi-pathogen antimicrobial stewardship simulation')

    # ── Trial design ─────────────────────────────────────────────────────
    parser.add_argument('--arm', default='ADE', choices=['ADE', 'control'],
                        help='Stewardship arm (default: ADE)')
    parser.add_argument('--days', type=int, default=1460,
                        help='Simulation duration in days (default: 1460 = 4 years)')
    parser.add_argument('--num_runs', type=int, default=1000,
                        help='Monte Carlo replicates (default: 1000)')

    # ── Ward configuration ────────────────────────────────────────────────
    parser.add_argument('--num_patient', type=int, default=64,
                        help='ICU bed capacity (default: 64)')
    parser.add_argument('--num_hcw', type=int, default=16,
                        help='Number of healthcare workers (default: 16)')
    parser.add_argument('--time_interval', type=int, default=2,
                        help='Hours between HCW visits (default: 2)')

    # ── Transmission parameters ───────────────────────────────────────────
    parser.add_argument('--p', type=float, default=0.30,
                        help='HCW→patient transmission probability (default: 0.30)')
    parser.add_argument('--q', type=float, default=0.30,
                        help='Patient→HCW contamination probability (default: 0.30)')
    parser.add_argument('--r', type=float, default=0.30,
                        help='Super-infection acquisition probability, early phase (default: 0.30)')
    parser.add_argument('--s', type=float, default=0.015,
                        help='Super-infection probability, late phase (default: 0.015)')
    parser.add_argument('--eta', type=float, default=0.50,
                        help='Hand-hygiene compliance (default: 0.50)')

    # ── Pathogen dynamics ────────────────────────────────────────────────
    parser.add_argument('--epsilon', type=float, default=0.03,
                        help='Intrinsic mutation hazard (default: 0.03)')
    # sigmax / sigmac kept for backward compatibility;
    # per-pathogen sigma values are now defined in config.PATHOGENS
    parser.add_argument('--sigmax', type=float, default=0.16,
                        help='Fallback sigma for unregistered susceptible strains')
    parser.add_argument('--sigmac', type=float, default=0.45,
                        help='Fallback sigma for unregistered resistant strains')

    # ── Admission parameters ─────────────────────────────────────────────
    parser.add_argument('--a', type=float, default=0.10,
                        help='Probability of ESKAPE colonization on admission (default: 0.10)')
    parser.add_argument('--m', type=float, default=0.60,
                        help='Prior pathogen exposure probability (default: 0.60)')

    # ── Clinical outcome hazard ratios ────────────────────────────────────
    parser.add_argument('--kappa_mu', type=float, default=0.74,
                        help='Discharge hazard ratio for infected patients (default: 0.74)')
    parser.add_argument('--kappa_nu', type=float, default=1.04,
                        help='Death hazard ratio for inadequately treated patients (default: 1.04)')

    # ── Probability data files ────────────────────────────────────────────
    parser.add_argument('--death_prob_file', default='../deathprob.csv')
    parser.add_argument('--discharge_prob_file', default='../dischargeprob.csv')

    # ── Stochastic noise ─────────────────────────────────────────────────
    parser.add_argument('--std', type=float, default=0.0,
                        help='Noise level for Monte Carlo parameter sampling (0 = deterministic)')

    args = parser.parse_args()

    # Select drug-change function based on trial arm
    drug_change_func = drug_change_ADE if args.arm == 'ADE' else drug_change_control

    base_dir = f'./log/{args.arm}_p{args.num_patient}_h{args.num_hcw}/'
    os.makedirs(base_dir, exist_ok=True)

    # Parameters perturbed per run when std > 0
    noisy_params = ['a', 'p', 'q', 'r', 's', 'epsilon',
                    'sigmax', 'sigmac', 'eta', 'kappa_mu', 'kappa_nu']

    whole_res = {}

    for i in range(args.num_runs):
        run_args = copy.deepcopy(args)
        if args.std > 0.0:
            for param in noisy_params:
                val = getattr(args, param)
                setattr(run_args, param, truncnorm_func(val, args.std))

        whole_res[i] = main(run_args, drug_change_func)

        if i % 50 == 0:
            print(f'[{args.arm}] run {i}/{args.num_runs}')

    recordEverything(whole_res, args.arm, base_dir)
    print(f'Results saved to {base_dir}')
