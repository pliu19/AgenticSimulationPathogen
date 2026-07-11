#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Runtime config layer for the ESKAPE multi-pathogen ICU ABM.

Data lives in config_data.py (auto-generated from the two excel workbooks
by build_config.py at the repo root). This module adds the helper
functions used by scheme_func.py and main.py:

  * effective_coverage(drug, species, resistances)
  * best_agent(species, resistances)
  * mutation_target(species, drug)
  * empiric_group_1(species), empiric_group_2(species)

Patient phenotype state in this version is a FROZENSET of resistance
labels (e.g. frozenset({'Pip-Tazo-R', 'Meropenem-R'})). The empty set is
'Susceptible'. Coverage of a drug against a patient is the minimum
coverage across every label in the set.
"""

from __future__ import annotations

from config_data import (
    PATHOGENS,
    INIT_SPECIES_PROBS,
    INIT_PHENOTYPE_PROBS,
    DRUG_TIER,
    DRUG_TO_R_PHENOTYPE,
    COVERAGE_MATRIX,
)

__all__ = [
    'PATHOGENS', 'INIT_SPECIES_PROBS', 'INIT_PHENOTYPE_PROBS',
    'DRUG_TIER', 'DRUG_TO_R_PHENOTYPE', 'COVERAGE_MATRIX',
    'COVERAGE_RANK', 'SUSCEPTIBLE_LABEL',
    'effective_coverage', 'best_agent', 'mutation_target',
    'empiric_group_1', 'empiric_group_2', 'treatment_days_for',
    'all_drugs',
]


# ─── Coverage state ordering ─────────────────────────────────────────────────
# Worst → best. The min over a resistance set composes monotonically.

COVERAGE_RANK = {'none': 0, 'partial': 1, 'covers': 2}
COVERAGE_LABEL = {0: 'none', 1: 'partial', 2: 'covers'}
SUSCEPTIBLE_LABEL = 'Susceptible'


def all_drugs() -> list[str]:
    """Drugs ordered by (tier, name) for deterministic iteration."""
    return sorted(DRUG_TIER, key=lambda d: (DRUG_TIER[d], d))


def effective_coverage(drug: str, species: str, resistances) -> str:
    """Coverage of drug against a patient with the given resistance set.

    Composition rule: min over every label in `resistances` (worst case
    wins). Empty set → look up the Susceptible column. Returns one of
    'covers' | 'partial' | 'none'. Unknown drug / species returns 'none'.
    """
    if species == 'none' or species is None:
        return 'covers'
    sp_cov = COVERAGE_MATRIX.get(species)
    if sp_cov is None:
        return 'none'
    drug_cov = sp_cov.get(drug)
    if drug_cov is None:
        return 'none'

    if not resistances:
        return drug_cov.get(SUSCEPTIBLE_LABEL, 'none')

    worst = COVERAGE_RANK['covers']
    for r in resistances:
        state = drug_cov.get(r)
        if state is None:
            # Unknown resistance label for this species — be conservative.
            return 'none'
        rank = COVERAGE_RANK[state]
        if rank < worst:
            worst = rank
            if worst == 0:
                return 'none'
    return COVERAGE_LABEL[worst]


def best_agent(species: str, resistances) -> str | None:
    """Lowest-tier drug that covers the (species, resistance set).

    Prefer 'covers' (🟢) over 'partial' (🟡) at any tier. Among ties,
    pick the alphabetically first drug for determinism. Returns None if
    no drug has even partial coverage.
    """
    if species == 'none' or species is None:
        return None

    best_cov_drug = None
    best_cov_tier = 999
    best_par_drug = None
    best_par_tier = 999

    for drug in all_drugs():
        state = effective_coverage(drug, species, resistances)
        tier  = DRUG_TIER[drug]
        if state == 'covers' and tier < best_cov_tier:
            best_cov_tier = tier
            best_cov_drug = drug
        elif state == 'partial' and tier < best_par_tier:
            best_par_tier = tier
            best_par_drug = drug

    return best_cov_drug if best_cov_drug is not None else best_par_drug


def mutation_target(species: str, drug: str) -> str | None:
    """If treating `species` with `drug` can induce resistance, return the
    resulting R-phenotype label; otherwise None.

    Mutation is possible only when the corresponding 'Drug-R' column
    exists in that species' coverage sheet.
    """
    if species == 'none' or species is None or drug is None:
        return None
    r_label = DRUG_TO_R_PHENOTYPE.get(drug)
    if r_label is None:
        return None
    cfg = PATHOGENS.get(species)
    if cfg is None:
        return None
    if r_label in cfg['phenotypes']:
        return r_label
    return None


def treatment_days_for(species: str, resistances) -> int:
    """Treatment length given current resistance set.

    7d for Susceptible (empty set), 10d for single-R, 14d for multi-R.
    """
    n = len(resistances) if resistances else 0
    if n == 0:
        return 7
    if n == 1:
        return 10
    return 14


# ─── Empiric agent precomputation ────────────────────────────────────────────
# At infection onset the resistance set is the patient's current set. But the
# *empiric* group is selected by species only (the lab hasn't returned a
# phenotype yet, only Gram stain / species). Group 1 = the most stewardship-
# friendly drug that covers the Susceptible phenotype for that species.
# Group 2 = the next such drug. Override per species via EMPIRIC_OVERRIDE_*.

EMPIRIC_OVERRIDE_GROUP_1: dict[str, str] = {
    # Clinical defaults that match the prior model where applicable.
    'kpneumoniae':  'Piperacillin-Tazobactam',
    'enterobacter': 'Piperacillin-Tazobactam',
    'paeruginosa':  'Piperacillin-Tazobactam',
    'abaumannii':   'Ampicillin-Sulbactam',
    'saureus':      'Vancomycin',
    'efaecium':     'Vancomycin',
}

EMPIRIC_OVERRIDE_GROUP_2: dict[str, str] = {
    'kpneumoniae':  'Cefepime',
    'enterobacter': 'Cefepime',
    'paeruginosa':  'Cefepime',
    'abaumannii':   'Meropenem',
    'saureus':      'Tigecycline',
    'efaecium':     'Ampicillin-Sulbactam',
}


def _auto_empiric(species: str, exclude: set[str]) -> str | None:
    """Lowest-tier 🟢 drug for Susceptible, skipping `exclude`. Falls back
    to 🟡 if no 🟢 candidate is available."""
    best_cov_drug = None
    best_cov_tier = 999
    best_par_drug = None
    best_par_tier = 999
    for drug in all_drugs():
        if drug in exclude:
            continue
        state = effective_coverage(drug, species, frozenset())
        tier  = DRUG_TIER[drug]
        if state == 'covers' and tier < best_cov_tier:
            best_cov_tier = tier
            best_cov_drug = drug
        elif state == 'partial' and tier < best_par_tier:
            best_par_tier = tier
            best_par_drug = drug
    return best_cov_drug if best_cov_drug is not None else best_par_drug


def empiric_group_1(species: str) -> str | None:
    """Group-1 empiric agent for `species`."""
    if species in EMPIRIC_OVERRIDE_GROUP_1:
        return EMPIRIC_OVERRIDE_GROUP_1[species]
    return _auto_empiric(species, exclude=set())


def empiric_group_2(species: str) -> str | None:
    """Group-2 empiric agent (must differ from Group 1)."""
    if species in EMPIRIC_OVERRIDE_GROUP_2:
        return EMPIRIC_OVERRIDE_GROUP_2[species]
    g1 = empiric_group_1(species)
    return _auto_empiric(species, exclude={g1} if g1 else set())
