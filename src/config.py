#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Configuration registry for the ESKAPE multi-pathogen ICU ABM.

All pathogen biology, drug coverage, and ecological parameters live here.
Event functions in scheme_func.py read from this module rather than using
hardcoded if/elif chains.
"""

# ─── ESKAPE Pathogen Registry ──────────────────────────────────────────────────
# sigma: per-phenotype probability that colonization progresses to active infection
# treatment_days: standard treatment course length in days
PATHOGENS = {
    'efaecium': {
        'full_name': 'Enterococcus faecium',
        'phenotypes': ['susceptible', 'VRE'],
        'sigma': {'susceptible': 0.094, 'VRE': 0.094},
        'treatment_days': {'susceptible': 7, 'VRE': 14},
    },
    'saureus': {
        'full_name': 'Staphylococcus aureus',
        'phenotypes': ['MSSA', 'MRSA'],
        'sigma': {'MSSA': 0.094, 'MRSA': 0.094},
        'treatment_days': {'MSSA': 7, 'MRSA': 14},
    },
    'kpneumoniae': {
        'full_name': 'Klebsiella pneumoniae',
        'phenotypes': ['susceptible', 'ESBL', 'KPC'],
        'sigma': {'susceptible': 0.126, 'ESBL': 0.126, 'KPC': 0.126},
        'treatment_days': {'susceptible': 7, 'ESBL': 10, 'KPC': 14},
    },
    'abaumannii': {
        'full_name': 'Acinetobacter baumannii',
        'phenotypes': ['susceptible', 'CRAB'],
        'sigma': {'susceptible': 0.200, 'CRAB': 0.200},
        'treatment_days': {'susceptible': 7, 'CRAB': 14},
    },
    'paeruginosa': {
        'full_name': 'Pseudomonas aeruginosa',
        'phenotypes': ['susceptible', 'CRPA', 'DTR'],
        'sigma': {'susceptible': 0.450, 'CRPA': 0.450, 'DTR': 0.450},
        'treatment_days': {'susceptible': 7, 'CRPA': 10, 'DTR': 14},
    },
    'enterobacter': {
        'full_name': 'Enterobacter spp.',
        'phenotypes': ['susceptible', 'ESBL', 'AmpC'],
        'sigma': {'susceptible': 0.094, 'ESBL': 0.094, 'AmpC': 0.094},
        'treatment_days': {'susceptible': 7, 'ESBL': 10, 'AmpC': 10},
    },
}

# ─── WHO AWaRe Drug Classification ────────────────────────────────────────────
DRUG_LEVELS = {
    1: {'label': 'Access/Narrow',    'examples': 'Nafcillin, Penicillin G'},
    2: {'label': 'Access/Moderate',  'examples': 'Amp-Sulbactam, Gentamicin'},
    3: {'label': 'Watch/Medium',     'examples': 'Ceftriaxone, Ciprofloxacin'},
    4: {'label': 'Watch/Broad',      'examples': 'Pip-Tazo, Cefepime, Meropenem'},
    5: {'label': 'Reserve',          'examples': 'CZA, Colistin, Linezolid'},
}

# ─── Drug Coverage Matrix ─────────────────────────────────────────────────────
# Minimum AWaRe level required to adequately cover each (species, phenotype).
# Clinical basis:
#   Level 1  → MSSA only (nafcillin/penicillinase-resistant penicillins)
#   Level 2  → Susceptible Enterococcus, susceptible GN community organisms
#   Level 3  → Susceptible Enterobacteriaceae (ceftriaxone, ciprofloxacin)
#   Level 4  → Pseudomonas, ESBL-producers, A. baumannii susceptible (pip-tazo/meropenem)
#   Level 5  → KPC, CRAB, MRSA, VRE, CRPA, DTR (reserve agents only)
COVERAGE_REQUIREMENTS = {
    ('saureus',      'MSSA'):        1,
    ('saureus',      'MRSA'):        5,
    ('efaecium',     'susceptible'): 2,
    ('efaecium',     'VRE'):         5,
    ('kpneumoniae',  'susceptible'): 3,
    ('kpneumoniae',  'ESBL'):        4,
    ('kpneumoniae',  'KPC'):         5,
    ('abaumannii',   'susceptible'): 4,
    ('abaumannii',   'CRAB'):        5,
    ('paeruginosa',  'susceptible'): 4,
    ('paeruginosa',  'CRPA'):        5,
    ('paeruginosa',  'DTR'):         5,
    ('enterobacter', 'susceptible'): 3,
    ('enterobacter', 'ESBL'):        4,
    ('enterobacter', 'AmpC'):        4,
    ('none',         'none'):        0,  # uncolonized; any drug level "covers"
}


def min_drug_level(species, phenotype):
    """Return the minimum AWaRe level needed to cover (species, phenotype)."""
    return COVERAGE_REQUIREMENTS.get((species, phenotype), 4)


def is_covered(drug_level, species, phenotype):
    """Return True if drug_level provides adequate coverage for this pathogen."""
    if drug_level is None:
        return False
    return drug_level >= min_drug_level(species, phenotype)


# ─── Intrinsic Mutation Pathways ───────────────────────────────────────────────
# Resistance escalation under treatment pressure (drug covers current phenotype,
# selecting for more-resistant mutants).
# Format: (species, current_phenotype) -> next_phenotype
MUTATION_PATHWAYS = {
    ('kpneumoniae',  'susceptible'): 'ESBL',
    ('kpneumoniae',  'ESBL'):        'KPC',
    ('abaumannii',   'susceptible'): 'CRAB',
    ('paeruginosa',  'susceptible'): 'CRPA',
    ('paeruginosa',  'CRPA'):        'DTR',
    ('saureus',      'MSSA'):        'MRSA',
    ('efaecium',     'susceptible'): 'VRE',
    ('enterobacter', 'susceptible'): 'ESBL',
    ('enterobacter', 'ESBL'):        'AmpC',
}

# ─── Initial Population Distribution ─────────────────────────────────────────
# Probability of being colonized with each ESKAPE species on ICU admission.
# Derived from midpoints of reported ICU prevalence ranges, normalized to sum=1.
# (E. faecium 25.95%, S. aureus 27.5%, K. pneumoniae 32.95%, A. baumannii 16.4%,
#  P. aeruginosa 21.9%, Enterobacter 5.35%) / 130.05 total
INIT_SPECIES_PROBS = {
    'efaecium':    0.200,
    'saureus':     0.212,
    'kpneumoniae': 0.253,
    'abaumannii':  0.126,
    'paeruginosa': 0.168,
    'enterobacter': 0.041,
}

# Within-species resistance phenotype distribution for newly admitted patients
# (reflects ICU-level resistance burden in modern surveillance data)
INIT_PHENOTYPE_PROBS = {
    'efaecium':    {'susceptible': 0.40, 'VRE':  0.60},
    'saureus':     {'MSSA':        0.40, 'MRSA': 0.60},
    'kpneumoniae': {'susceptible': 0.25, 'ESBL': 0.45, 'KPC':  0.30},
    'abaumannii':  {'susceptible': 0.30, 'CRAB': 0.70},
    'paeruginosa': {'susceptible': 0.50, 'CRPA': 0.35, 'DTR':  0.15},
    'enterobacter': {'susceptible': 0.40, 'ESBL': 0.40, 'AmpC': 0.20},
}
