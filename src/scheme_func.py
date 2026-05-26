#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Simulation event functions for the ESKAPE multi-pathogen ICU ABM.

All event logic reads species/phenotype relationships from config.py rather
than hardcoded strain-drug if/elif chains.
"""

import random
import numpy as np

from abm import HealthCareWorker, Patient
from config import (
    PATHOGENS, INIT_SPECIES_PROBS, INIT_PHENOTYPE_PROBS,
    EMPIRIC_AGENTS_1, EMPIRIC_AGENTS_2, MUTATION_PATHWAYS,
    is_covered, min_covering_agent,
)


class SimulationError(Exception):
    def __init__(self, m):
        self.message = m

    def __str__(self):
        return self.message


def draw():
    """Draw a uniform random number in [0, 1)."""
    return np.random.random_sample()


# ─── Initialisation ───────────────────────────────────────────────────────────

def initial_patient(args, name, current_day):
    """
    Create a newly admitted patient.

    Species is drawn from INIT_SPECIES_PROBS, which includes 'none' (uncolonized).
    If a species is drawn, phenotype is drawn from INIT_PHENOTYPE_PROBS.
    """
    species_list  = list(INIT_SPECIES_PROBS.keys())
    species_probs = list(INIT_SPECIES_PROBS.values())
    species = np.random.choice(species_list, p=species_probs)

    if species == 'none':
        phenotype = 'none'
    else:
        pheno_dict      = INIT_PHENOTYPE_PROBS[species]
        phenotype_list  = list(pheno_dict.keys())
        phenotype_probs = list(pheno_dict.values())
        phenotype = np.random.choice(phenotype_list, p=phenotype_probs)

    patient = Patient(args, name, current_day, species, phenotype)
    patient.recordStatus('initialization', 0)
    return patient


def bulk_initialization(args):
    """Populate the full ICU ward at simulation start."""
    return {idx: initial_patient(args, idx, 0) for idx in range(args.num_patient)}


def random_schedule(current_patients, num_hcw=4):
    """
    Randomly partition the current patient list among HCWs for one shift.
    Returns a dict {hcw_id: [patient_indices]}.
    """
    shuffled = list(current_patients)
    np.random.shuffle(shuffled)
    gap = len(shuffled) // num_hcw
    return {f'h{i}': shuffled[i * gap:(i + 1) * gap] for i in range(num_hcw)}


# ─── Patient–HCW Contact Events ───────────────────────────────────────────────

def _acquire_from_hcw(patient, healthW, args, current_day, reference_results, prob):
    """
    Patient (species='none') may acquire one strain carried by the HCW.
    Returns modified (healthW, patient).
    """
    if len(healthW.strain_set) == 0:
        return healthW, patient
    if draw() < prob:
        sp, ph = random.choice(list(healthW.strain_set))
        patient.species   = sp
        patient.phenotype = ph
        patient.reset_infection_time(args, current_day)
        key = f'colonization_{sp}_{ph}'
        reference_results[key] = reference_results.get(key, 0) + 1
    return healthW, patient


def _contaminate_hcw(patient, healthW, prob):
    """HCW may acquire the patient's colonising strain."""
    if patient.species != 'none' and draw() < prob:
        healthW.strain_set.add((patient.species, patient.phenotype))
    return healthW


def ph_interaction(patient, healthW, args, current_day, current_hour,
                   reference_results):
    """
    Dispatch to the appropriate contact handler based on patient state.
    Mirrors the original three-case structure (colonised, early infected,
    late infected) but uses the AWaRe coverage matrix instead of string
    comparisons.
    """
    p = float(args.p)
    q = float(args.q)
    r = float(args.r)
    s = float(args.s)

    review_hours = args.empiric_hours

    if patient.status == 'C':
        healthW, patient = _contact_colonised(
            patient, healthW, args, current_day, reference_results, p, q)

    elif patient.status == 'I' and patient.treat_time is not None \
            and patient.treat_time <= review_hours:
        healthW, patient = _contact_infected_early(
            patient, healthW, args, current_day, current_hour,
            reference_results, q, r)

    elif patient.status == 'I' and patient.treat_time is not None \
            and patient.treat_time > review_hours:
        healthW, patient = _contact_infected_late(
            patient, healthW, current_hour, reference_results, q, s)

    return healthW, patient


def _contact_colonised(patient, healthW, args, current_day,
                        reference_results, p, q):
    """
    Colonised patient contact:
      • species='none'  → patient can acquire HCW strain at rate p.
      • ESKAPE carrier  → HCW can acquire patient's strain at rate q.
    """
    if patient.species == 'none':
        healthW, patient = _acquire_from_hcw(
            patient, healthW, args, current_day, reference_results, p)
    else:
        healthW = _contaminate_hcw(patient, healthW, q)
    return healthW, patient


def _contact_infected_early(patient, healthW, args, current_day, current_hour,
                             reference_results, q, r):
    """
    Infected patient, early treatment (treat_time ≤ 72 h).

    HCW acquires patient's strain at rate q.
    Patient can acquire any HCW strain that is NOT covered by the current
    empiric drug level, at rate r (super-infection / strain upgrade).
    """
    # HCW contamination
    healthW = _contaminate_hcw(patient, healthW, q)

    # Super-infection: find uncovered HCW strains
    candidates = _uncovered_hcw_strains(healthW, patient)
    if candidates and draw() < r:
        h_sp, h_ph = random.choice(candidates)
        _apply_superinfection(patient, h_sp, h_ph, current_hour, reference_results)

    return healthW, patient


def _contact_infected_late(patient, healthW, current_hour,
                            reference_results, q, s):
    """
    Infected patient, late treatment (treat_time > 72 h).

    Same as early phase but uses the super-infection probability s instead
    of the transmission probability r, reflecting targeted therapy in use.
    """
    healthW = _contaminate_hcw(patient, healthW, q)

    candidates = _uncovered_hcw_strains(healthW, patient)
    if candidates and draw() < s:
        h_sp, h_ph = random.choice(candidates)
        _apply_superinfection(patient, h_sp, h_ph, current_hour, reference_results)

    return healthW, patient


def _uncovered_hcw_strains(healthW, patient):
    """
    Return a list of (species, phenotype) carried by the HCW that are NOT
    covered by the patient's current drug level.

    Includes same-species strains with a more resistant phenotype (resistance
    upgrade within a species).
    """
    dagent = patient.drug_agent
    result = []
    for h_sp, h_ph in healthW.strain_set:
        # Skip identical strain
        if h_sp == patient.species and h_ph == patient.phenotype:
            continue
        if not is_covered(dagent, h_sp, h_ph):
            result.append((h_sp, h_ph))
    return result


def _apply_superinfection(patient, h_sp, h_ph, current_hour, reference_results):
    """Apply a super-infection: update patient's dominant pathogen."""
    patient.species   = h_sp
    patient.phenotype = h_ph
    patient.super_infe = True
    patient.convt_time = current_hour
    reference_results['cumsuperinfection'] = \
        reference_results.get('cumsuperinfection', 0) + 1
    reference_results['misempiric'] = reference_results.get('misempiric', 0) + 1
    tkey = f'transmission_{h_sp}_{h_ph}'
    reference_results[tkey] = reference_results.get(tkey, 0) + 1


# ─── HCW Decontamination ──────────────────────────────────────────────────────

def hcw_cleanUp(healthW, current_hour, eta=0.5, end_of_shift=False):
    """
    HCW decontamination.  Always clears at the end of a shift; otherwise
    clears with probability eta (hand-hygiene compliance).
    """
    if end_of_shift or draw() < eta:
        healthW.strain_set.clear()
    return healthW


# ─── Daily Events ─────────────────────────────────────────────────────────────

def get_empiric_agent(species, args, current_day):
    """
    Select the empiric agent for a species based on the active regimen:
      fixed   — always group 1
      cycling — group 1 for N months, then group 2 for N months, repeating
      mixing  — group 1 or group 2 chosen with equal probability per patient
    Falls back to 'meropenem' for unregistered species.
    """
    regimen = getattr(args, 'empiric_regimen', 'fixed')

    if regimen == 'cycling':
        current_month  = current_day // 30
        cycle_position = (current_month // args.cycle_months) % 2
        agents = EMPIRIC_AGENTS_1 if cycle_position == 0 else EMPIRIC_AGENTS_2
    elif regimen == 'mixing':
        agents = EMPIRIC_AGENTS_1 if draw() < 0.5 else EMPIRIC_AGENTS_2
    else:  # fixed
        agents = EMPIRIC_AGENTS_1

    return agents.get(species, 'meropenem')


def infection_development(patients_record, current_patients, current_day,
                           reference_results, args, daily_counters):
    """
    Colonised patients whose infct_time has arrived progress to active
    infection.  Species-specific empiric agent is assigned via get_empiric_agent().
    """
    for key in current_patients:
        patient = patients_record[key]
        if (patient.status == 'C'
                and patient.infct_flag
                and patient.infct_time == current_day):

            patient.status     = 'I'
            patient.treat_time = 0
            agent = get_empiric_agent(patient.species, args, current_day)
            patient.drug_agent = agent
            # Track which empiric group was assigned
            if agent == EMPIRIC_AGENTS_1.get(patient.species, 'meropenem'):
                daily_counters['empiric_group_1'] += 1
            else:
                daily_counters['empiric_group_2'] += 1

            if patient.species != 'none':
                patient.lab_result = (patient.species, patient.phenotype)
                ci_key = f'cuminfection_{patient.species}_{patient.phenotype}'
                reference_results[ci_key] = reference_results.get(ci_key, 0) + 1
                lr_key = f'labresult_{patient.species}_{patient.phenotype}'
                reference_results[lr_key] = reference_results.get(lr_key, 0) + 1
            else:
                reference_results['cuminfection_none'] = \
                    reference_results.get('cuminfection_none', 0) + 1

    return patients_record


def increase_treatment_icu_time(patients_record, current_patients):
    """Increment ICU and treatment times by 24 h at the end of each day."""
    for key in current_patients:
        patient = patients_record[key]
        patient.time_inICU += 24
        if isinstance(patient.treat_time, int):
            patient.treat_time       += 24
            patient.total_treat_time += 24
    return patients_record


def drug_change(reference_results, patients_record, current_patients,
                current_hour, args, daily_counters):
    """
    At empiric_hours: if the empiric agent covers the identified phenotype, keep it.
    Otherwise escalate to the lowest-level non-last-resort covering agent,
    or last-resort if no other agent covers the pathogen.
    Re-evaluates after any strain conversion.
    """
    review_hours = args.empiric_hours

    for key in current_patients:
        patient = patients_record[key]

        if patient.treat_time == review_hours and patient.lab_result is not None:
            sp, ph = patient.lab_result
            if not is_covered(patient.drug_agent, sp, ph):
                new_agent = min_covering_agent(sp, ph)
                patient.drug_agent = new_agent
                dakey = f'drug_agent_{new_agent}_use'
                reference_results[dakey] = reference_results.get(dakey, 0) + 1
                reference_results['misempiric'] = \
                    reference_results.get('misempiric', 0) + 1
                daily_counters['misempiric_daily'] += 1

        elif (patient.convt_time is not None
              and (current_hour - patient.convt_time) == review_hours):
            sp, ph = patient.species, patient.phenotype
            if not is_covered(patient.drug_agent, sp, ph):
                new_agent = min_covering_agent(sp, ph)
                patient.drug_agent = new_agent
                dakey = f'drug_agent_{new_agent}_use'
                reference_results[dakey] = reference_results.get(dakey, 0) + 1

    return patients_record


def check_treatment_completion(patient, current_hour):
    """
    Return True when the patient's treatment course is finished.
    Uses per-species/phenotype durations from the pathogen registry.
    """
    if patient.species == 'none' or patient.treat_time is None:
        return False

    treat_days = PATHOGENS[patient.species]['treatment_days'].get(
        patient.phenotype, 10)

    # Straightforward treatment course (no strain conversion)
    if patient.convt_time is None and patient.treat_time == treat_days * 24:
        return True

    # Post-conversion: 10-day course from the conversion point
    if patient.convt_time is not None \
            and (current_hour - patient.convt_time) == 10 * 24:
        return True

    return False


def treatment_completion(reference_results, patients_record, current_patients,
                          current_hour):
    """Clear infection status for patients whose treatment course is complete."""
    for key in current_patients:
        patient = patients_record[key]
        if check_treatment_completion(patient, current_hour):
            patient.status     = 'C'
            patient.treat_time = None
            patient.convt_time = None
            patient.lab_result = None
            patient.drug_agent = None
    return patients_record


def intrinsic_mutation(patients_record, current_patients, args, current_hour,
                        reference_results, first_mutation_times):
    """
    Resistance emergence under antibiotic selective pressure.

    Mutation occurs when the current drug COVERS the patient's pathogen
    (i.e. the drug kills susceptibles, selecting for resistant mutants).
    The organism advances one step along its species' resistance ladder
    as defined in MUTATION_PATHWAYS.
    """
    epsilon = float(args.epsilon)

    for key in current_patients:
        patient = patients_record[key]
        if patient.status != 'I' or patient.species == 'none':
            continue

        # Only mutate under effective drug pressure
        if not is_covered(patient.drug_agent, patient.species, patient.phenotype):
            continue

        if draw() < epsilon:
            next_ph_options = MUTATION_PATHWAYS.get((patient.species, patient.phenotype))
            if next_ph_options is not None:
                next_ph = random.choice(next_ph_options)
                patient.phenotype  = next_ph
                if patient.convt_time is None:
                    patient.convt_time = current_hour
                mkey = f'mutation_{patient.species}_{next_ph}'
                reference_results[mkey] = reference_results.get(mkey, 0) + 1
                # Record day of first mutation to this phenotype
                fkey = (patient.species, next_ph)
                if first_mutation_times.get(fkey) is None:
                    first_mutation_times[fkey] = current_hour // 24

    return patients_record


def discharge_admission(patients_record, current_hour, args, discharge_probs,
                         current_patients, patient_last_idx, reference_results,
                         treat_lengths, icu_los):
    """
    Stochastic discharge events.  Discharged patients are immediately
    replaced by a newly initialised admission.

    Infected patients have a modified discharge hazard (kappa_mu).
    Patients with super-infection are not eligible for discharge.
    Force-discharges any patient who has reached 40 days.
    """
    current_day  = current_hour // 24
    kappa_mu     = float(args.kappa_mu)
    copy_pts     = current_patients.copy()

    for key in copy_pts:
        patient = patients_record[key]

        # Mandatory discharge at 40 days
        if patient.time_inICU == 40 * 24:
            patient_last_idx, current_patients = _replace_patient(
                patients_record, current_patients, key, patient_last_idx,
                args, current_day, reference_results, treat_lengths, icu_los,
                event='discharge')
            continue

        if patient.super_infe or patient.time_inICU < 2 * 24:
            continue

        day  = patient.time_inICU // 24
        prob = discharge_probs.get(day, 0)
        if patient.status == 'I':
            prob *= kappa_mu

        if draw() < prob:
            patient_last_idx, current_patients = _replace_patient(
                patients_record, current_patients, key, patient_last_idx,
                args, current_day, reference_results, treat_lengths, icu_los,
                event='discharge')

    return patients_record, current_patients, patient_last_idx


def death_event(patients_record, current_patients, patient_last_idx,
                current_hour, args, death_probs, reference_results,
                treat_lengths, icu_los):
    """
    Stochastic death events.

    Patients on inadequate therapy (drug agent doesn't cover their pathogen)
    face an elevated mortality hazard (kappa_nu).
    """
    current_day = current_hour // 24
    kappa_nu    = float(args.kappa_nu)
    copy_pts    = current_patients.copy()

    for key in copy_pts:
        patient   = patients_record[key]
        day       = patient.time_inICU // 24
        base_prob = death_probs.get(day, 0)

        inadequate = (patient.species != 'none'
                      and not is_covered(patient.drug_agent,
                                         patient.species, patient.phenotype))
        prob = base_prob * kappa_nu if inadequate else base_prob

        if draw() < prob:
            patient_last_idx, current_patients = _replace_patient(
                patients_record, current_patients, key, patient_last_idx,
                args, current_day, reference_results, treat_lengths, icu_los,
                event='death')

    return patients_record, current_patients, patient_last_idx


def _replace_patient(patients_record, current_patients, old_key,
                      patient_last_idx, args, current_day,
                      reference_results, treat_lengths, icu_los, event):
    """
    Replace a departing patient with a newly admitted one.
    Increments admission + event (discharge/death) counters.
    Records total treatment hours and ICU length of stay for the departing patient.
    """
    treat_lengths.append(patients_record[old_key].total_treat_time)
    icu_los.append(patients_record[old_key].time_inICU)
    patient_last_idx += 1
    patients_record[patient_last_idx] = initial_patient(args, patient_last_idx,
                                                         current_day)
    current_patients.remove(old_key)
    current_patients.append(patient_last_idx)
    reference_results['admission'] = reference_results.get('admission', 0) + 1
    reference_results[event]       = reference_results.get(event, 0) + 1
    return patient_last_idx, current_patients
