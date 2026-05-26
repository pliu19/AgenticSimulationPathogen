#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Main simulation loop for the ESKAPE multi-pathogen ICU ABM.

Drug change rule: at 72 h keep the empiric agent if it covers the identified
phenotype; otherwise escalate to the minimum covering agent (last-resort
only if no lower-level agent covers the pathogen).

Usage examples:
    python main.py --num_runs 1000
    python main.py --num_runs 10 --days 30
"""

import os
import copy
import argparse
import math
import time
import pickle

import numpy as np

from abm import HealthCareWorker
from config import PATHOGENS, AGENT_LEVEL
from scheme_func import (
    draw, initial_patient, random_schedule, bulk_initialization,
    ph_interaction, hcw_cleanUp,
    infection_development, increase_treatment_icu_time,
    drug_change,
    treatment_completion, discharge_admission, death_event,
    intrinsic_mutation,
)


# ─── Utility functions ────────────────────────────────────────────────────────

def truncnorm_func(mean, std_rate):
    """Sample from a truncated normal distribution centred on mean (±1 std)."""
    if std_rate == 0 or mean == 0:
        return mean
    std = abs(mean) * std_rate
    while True:
        sample = np.random.normal(mean, std)
        if mean - std <= sample <= mean + std:
            return max(1e-5, sample)


def read_csv(file_path):
    """
    Read a day → probability CSV and return a dict of survival probabilities
    computed via the hazard-integral method used in the original model.
    Integrating a constant lam from 0 to t gives lam * t.
    """
    record = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                day_idx = int(parts[0])
                prob    = float(parts[1])
                record[day_idx] = prob

    return {key: value * math.exp(-value * key) for key, value in record.items()}


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
        'cumsuperinfection': 0, 'misempiric': 0,
        'cuminfection_none': 0,
    }
    for agent in AGENT_LEVEL:
        r[f'drug_agent_{agent}_use'] = 0

    for sp, cfg in PATHOGENS.items():
        for ph in cfg['phenotypes']:
            r[f'cuminfection_{sp}_{ph}']  = 0
            r[f'labresult_{sp}_{ph}']     = 0
            r[f'colonization_{sp}_{ph}']  = 0
            r[f'transmission_{sp}_{ph}']  = 0
            r[f'mutation_{sp}_{ph}']      = 0

    return r


# ─── Main simulation loop ─────────────────────────────────────────────────────

def main(args):
    """
    Run one complete simulation of args.days days.
    Returns final_results: a dict of daily time-series and cumulative counters.
    """
    reference_results = make_reference_results()

    # Per-patient lists (one entry per discharge/death) — variable length
    treat_lengths = []
    icu_los       = []

    # Per-run first-mutation tracking: (species, phenotype) -> simulation day
    first_mutation_times = {
        (sp, ph): None
        for sp, cfg in PATHOGENS.items()
        for ph in cfg['phenotypes']
    }

    # Daily reset counters (reset to 0 after each day's snapshot)
    daily_counters = {
        'empiric_group_1': 0,
        'empiric_group_2': 0,
        'misempiric_daily': 0,
    }

    # Daily census arrays — pre-allocated numpy arrays (int32)
    final_results = {'superinfection': np.zeros(args.days, dtype=np.int32)}
    for sp, cfg in PATHOGENS.items():
        for ph in cfg['phenotypes']:
            final_results[f'infection_{sp}_{ph}']        = np.zeros(args.days, dtype=np.int32)
            final_results[f'colonization_prev_{sp}_{ph}'] = np.zeros(args.days, dtype=np.int32)
    final_results['hcw_contamination_burden'] = np.zeros(args.days, dtype=np.float32)
    final_results['empiric_group_1_daily']    = np.zeros(args.days, dtype=np.int32)
    final_results['empiric_group_2_daily']    = np.zeros(args.days, dtype=np.int32)
    final_results['misempiric_daily']         = np.zeros(args.days, dtype=np.int32)
    for key in reference_results:
        final_results[key] = np.zeros(args.days, dtype=np.int32)

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
        hcw_burden_sum = 0.0
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

            # Snapshot burden before cleanup (avg strains per HCW this shift)
            hcw_burden_sum += sum(len(hcw_list[h].strain_set) for h in hcw_list) / len(hcw_list)

            # End-of-shift: guaranteed HCW decontamination
            for hid in hcw_list:
                hcw_list[hid] = hcw_cleanUp(
                    hcw_list[hid], current_hour, args.eta, end_of_shift=True)

        # ── End-of-day events ────────────────────────────────────────────
        patient_list = infection_development(
            patient_list, current_p_names, current_day, reference_results,
            args, daily_counters)

        patient_list = intrinsic_mutation(
            patient_list, current_p_names, args, current_hour,
            reference_results, first_mutation_times)

        patient_list = increase_treatment_icu_time(
            patient_list, current_p_names)

        patient_list = drug_change(
            reference_results, patient_list, current_p_names, current_hour,
            args, daily_counters)

        patient_list = treatment_completion(
            reference_results, patient_list, current_p_names, current_hour)

        patient_list, current_p_names, patient_idx = discharge_admission(
            patient_list, current_hour, args, discharge_probs,
            current_p_names, patient_idx, reference_results,
            treat_lengths, icu_los)

        patient_list, current_p_names, patient_idx = death_event(
            patient_list, current_p_names, patient_idx, current_hour,
            args, death_probs, reference_results, treat_lengths, icu_los)

        # ── Daily snapshot ───────────────────────────────────────────────
        for key in reference_results:
            final_results[key][current_day] = reference_results[key]

        for pt_key in current_p_names:
            patient = patient_list[pt_key]
            if patient.status == 'I' and patient.species != 'none':
                final_results[f'infection_{patient.species}_{patient.phenotype}'][current_day] += 1
            if patient.super_infe:
                final_results['superinfection'][current_day] += 1

        # Colonization prevalence — single pass over patients
        col_counts = {}
        for k in current_p_names:
            pt = patient_list[k]
            if pt.status == 'C' and pt.species != 'none':
                col_counts[(pt.species, pt.phenotype)] = col_counts.get((pt.species, pt.phenotype), 0) + 1
        for sp, cfg in PATHOGENS.items():
            for ph in cfg['phenotypes']:
                final_results[f'colonization_prev_{sp}_{ph}'][current_day] = col_counts.get((sp, ph), 0)

        # HCW contamination burden — daily average across 3 shifts (pre-cleanup)
        final_results['hcw_contamination_burden'][current_day] = hcw_burden_sum / 3

        # Daily counters snapshot then reset
        final_results['empiric_group_1_daily'][current_day] = daily_counters['empiric_group_1']
        final_results['empiric_group_2_daily'][current_day] = daily_counters['empiric_group_2']
        final_results['misempiric_daily'][current_day]      = daily_counters['misempiric_daily']
        for k in daily_counters:
            daily_counters[k] = 0

        current_day += 1

    final_results['treatment_length_on_discharge'] = treat_lengths
    final_results['icu_los_on_discharge']          = icu_los
    final_results['first_mutation_day'] = {
        f'{sp}_{ph}': day
        for (sp, ph), day in first_mutation_times.items()
    }
    return final_results


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='ESKAPE ICU ABM — multi-pathogen antimicrobial stewardship simulation')

    # ── Trial design ─────────────────────────────────────────────────────
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

    # ── Empiric regimen ───────────────────────────────────────────────────
    parser.add_argument('--empiric_regimen', default='fixed',
                        choices=['fixed', 'cycling', 'mixing'],
                        help='Empiric agent selection regimen (default: fixed)')
    parser.add_argument('--cycle_months', type=int, default=6,
                        help='Months per group in cycling regimen (default: 6)')
    parser.add_argument('--empiric_hours', type=int, default=72,
                        choices=[24, 48, 72],
                        help='Duration of empiric therapy before drug review (default: 72)')

    # ── Pathogen dynamics ────────────────────────────────────────────────
    parser.add_argument('--epsilon', type=float, default=0.03,
                        help='Intrinsic mutation hazard (default: 0.03)')

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

    # ── Output ────────────────────────────────────────────────────────────
    parser.add_argument('--run_name', type=str, default='',
                        help='Optional subfolder under log/ for batch organisation (e.g. run2)')

    args = parser.parse_args()

    subdir = f'{args.run_name}/' if args.run_name else ''
    base_dir = (f'./log/{subdir}{args.empiric_regimen}_h{args.empiric_hours}'
                f'_p{args.p}_q{args.q}_r{args.r}_n{args.num_patient}_w{args.num_hcw}/')
    os.makedirs(base_dir, exist_ok=True)

    # Parameters perturbed per run when std > 0
    noisy_params = ['p', 'q', 'r', 's', 'epsilon', 'eta', 'kappa_mu', 'kappa_nu']

    whole_res = {}

    for i in range(args.num_runs):
        run_args = copy.deepcopy(args)
        if args.std > 0.0:
            for param in noisy_params:
                val = getattr(args, param)
                setattr(run_args, param, truncnorm_func(val, args.std))

        whole_res[i] = main(run_args)

        if i % 50 == 0:
            print(f'run {i}/{args.num_runs}')

    recordEverything(whole_res, 'sim', base_dir)
    print(f'Results saved to {base_dir}')
