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
# treatment_days: susceptible=7, resistance phenotypes=10, dual=14
PATHOGENS = {
    'efaecium': {
        'full_name': 'Enterococcus faecium',
        'phenotypes': ['susceptible', 'VREfm', 'LRE', 'dual'],
        'sigma': {'susceptible': 0.08, 'VREfm': 0.08, 'LRE': 0.08, 'dual': 0.08},
        'treatment_days': {'susceptible': 7, 'VREfm': 10, 'LRE': 10, 'dual': 14},
    },
    'saureus': {
        'full_name': 'Staphylococcus aureus',
        'phenotypes': ['susceptible', 'VRSA', 'MRSA', 'dual'],
        'sigma': {'susceptible': 0.04, 'VRSA': 0.04, 'MRSA': 0.04, 'dual': 0.04},
        'treatment_days': {'susceptible': 7, 'VRSA': 10, 'MRSA': 10, 'dual': 14},
    },
    'kpneumoniae': {
        'full_name': 'Klebsiella pneumoniae',
        'phenotypes': ['susceptible', 'PTR-KP', 'CRKP', 'dual'],
        'sigma': {'susceptible': 0.22, 'PTR-KP': 0.22, 'CRKP': 0.22, 'dual': 0.22},
        'treatment_days': {'susceptible': 7, 'PTR-KP': 10, 'CRKP': 10, 'dual': 14},
    },
    'abaumannii': {
        'full_name': 'Acinetobacter baumannii',
        'phenotypes': ['susceptible', 'CRAB', 'SAMR-AB', 'dual'],
        'sigma': {'susceptible': 0.278, 'CRAB': 0.278, 'SAMR-AB': 0.278, 'dual': 0.278},
        'treatment_days': {'susceptible': 7, 'CRAB': 10, 'SAMR-AB': 10, 'dual': 14},
    },
    'paeruginosa': {
        'full_name': 'Pseudomonas aeruginosa',
        'phenotypes': ['susceptible', 'PTR-PA', 'CRPA', 'dual'],
        'sigma': {'susceptible': 0.23, 'PTR-PA': 0.23, 'CRPA': 0.23, 'dual': 0.23},
        'treatment_days': {'susceptible': 7, 'PTR-PA': 10, 'CRPA': 10, 'dual': 14},
    },
    'enterobacter': {
        'full_name': 'Enterobacter spp.',
        'phenotypes': ['susceptible', 'CRE', 'ESBL', 'dual'],
        'sigma': {'susceptible': 0.03, 'CRE': 0.03, 'ESBL': 0.03, 'dual': 0.03},
        'treatment_days': {'susceptible': 7, 'CRE': 10, 'ESBL': 10, 'dual': 14},
    },
}

# ─── WHO AWaRe Drug Classification ────────────────────────────────────────────
DRUG_LEVELS = {
    1: {'label': 'Access',  'examples': 'Ampicillin-sulbactam, Nafcillin/Oxacillin'},
    2: {'label': 'Watch',   'examples': 'Vancomycin, Pip-Tazo, Cefepime, Meropenem'},
    3: {'label': 'Reserve', 'examples': 'Linezolid, Last-resort'},
}

# ─── Agent-Level Mapping ──────────────────────────────────────────────────────
AGENT_LEVEL = {
    'nafcillin':     1,
    'amp-sulbactam': 1,
    'vancomycin':    2,
    'pip-tazo':      2,
    'cefepime':      2,
    'meropenem':     2,
    'linezolid':     3,
    'last-resort':   3,
}

# ─── Agent Coverage Matrix ────────────────────────────────────────────────────
# Maps each agent to the set of (species, phenotype) it covers.
# 'last-resort' covers all phenotypes and is handled as a special case in is_covered().
AGENT_COVERAGE = {
    'nafcillin': {
        ('saureus', 'susceptible'), ('saureus', 'VRSA'),
    },
    'amp-sulbactam': {
        ('abaumannii', 'susceptible'), ('abaumannii', 'CRAB'),
    },
    'vancomycin': {
        ('efaecium', 'susceptible'),
        ('saureus',  'susceptible'), ('saureus',  'MRSA'),
    },
    'pip-tazo': {
        ('kpneumoniae', 'susceptible'), ('kpneumoniae', 'CRKP'),
        ('paeruginosa', 'susceptible'), ('paeruginosa', 'CRPA'),
    },
    'cefepime': {
        ('kpneumoniae',  'susceptible'), ('kpneumoniae',  'PTR-KP'),
        ('enterobacter', 'susceptible'), ('enterobacter', 'CRE'),
    },
    'meropenem': {
        ('abaumannii',   'susceptible'), ('abaumannii',   'SAMR-AB'),
        ('paeruginosa',  'susceptible'), ('paeruginosa',  'PTR-PA'),
        ('enterobacter', 'susceptible'), ('enterobacter', 'ESBL'),
    },
    'linezolid': {
        ('efaecium', 'susceptible'), ('efaecium', 'LRE'), ('efaecium', 'VREfm'),
    },
}

# ─── Empiric Agents ───────────────────────────────────────────────────────────
# First agent given at infection onset, by species.
# Species is assumed known at onset; phenotype is confirmed at 72 h by lab.
# Group is selected by empiric regimen: fixed (always group 1), cycling
# (alternate every N months), or mixing (50/50 random per patient).

EMPIRIC_AGENTS_1 = {
    'efaecium':     'vancomycin',
    'saureus':      'vancomycin',
    'kpneumoniae':  'pip-tazo',
    'abaumannii':   'meropenem',
    'paeruginosa':  'meropenem',
    'enterobacter': 'meropenem',
}

EMPIRIC_AGENTS_2 = {
    'efaecium':     'linezolid',
    'saureus':      'nafcillin',
    'kpneumoniae':  'cefepime',
    'abaumannii':   'amp-sulbactam',
    'paeruginosa':  'pip-tazo',
    'enterobacter': 'cefepime',
}

# ─── Intrinsic Mutation Pathways ───────────────────────────────────────────────
# Resistance escalation under treatment pressure (drug covers current phenotype,
# selecting for more-resistant mutants).
# Format: (species, current_phenotype) -> [list of possible next phenotypes]
# When mutation fires, one next phenotype is chosen at random from the list.
MUTATION_PATHWAYS = {
    ('efaecium',     'susceptible'): ['VREfm', 'LRE'],
    ('efaecium',     'VREfm'):       ['dual'],
    ('efaecium',     'LRE'):         ['dual'],
    ('saureus',      'susceptible'): ['VRSA', 'MRSA'],
    ('saureus',      'VRSA'):        ['dual'],
    ('saureus',      'MRSA'):        ['dual'],
    ('kpneumoniae',  'susceptible'): ['PTR-KP', 'CRKP'],
    ('kpneumoniae',  'PTR-KP'):      ['dual'],
    ('kpneumoniae',  'CRKP'):        ['dual'],
    ('abaumannii',   'susceptible'): ['CRAB', 'SAMR-AB'],
    ('abaumannii',   'CRAB'):        ['dual'],
    ('abaumannii',   'SAMR-AB'):     ['dual'],
    ('paeruginosa',  'susceptible'): ['PTR-PA', 'CRPA'],
    ('paeruginosa',  'PTR-PA'):      ['dual'],
    ('paeruginosa',  'CRPA'):        ['dual'],
    ('enterobacter', 'susceptible'): ['CRE', 'ESBL'],
    ('enterobacter', 'CRE'):         ['dual'],
    ('enterobacter', 'ESBL'):        ['dual'],
}

# ─── Initial Population Distribution ─────────────────────────────────────────
# Probability of being colonized with each ESKAPE species on ICU admission.
INIT_SPECIES_PROBS = {
    'none':         0.554,
    'saureus':      0.300,
    'kpneumoniae':  0.050,
    'paeruginosa':  0.026,
    'efaecium':     0.030,
    'abaumannii':   0.030,
    'enterobacter': 0.010,
}

# Within-species resistance phenotype distribution for newly admitted patients.
# PLACEHOLDER: equal distribution — update with ICU surveillance data.
INIT_PHENOTYPE_PROBS = {
    'efaecium':     {'susceptible': 0.922, 'VREfm': 0.058, 'LRE': 0.02,  'dual': 0.0},
    'saureus':      {'susceptible': 0.799, 'VRSA':  0.001, 'MRSA': 0.2,  'dual': 0.0},
    'kpneumoniae':  {'susceptible': 0.32,  'PTR-KP': 0.28, 'CRKP': 0.4,  'dual': 0.0},
    'abaumannii':   {'susceptible': 0.49,  'CRAB':  0.06,  'SAMR-AB': 0.45, 'dual': 0.0},
    'paeruginosa':  {'susceptible': 0.70,  'PTR-PA': 0.18, 'CRPA': 0.12, 'dual': 0.0},
    'enterobacter': {'susceptible': 0.736, 'CRE':   0.014, 'ESBL': 0.25, 'dual': 0.0},
}


# ─── Coverage Helper Functions ────────────────────────────────────────────────

def is_covered(drug_agent, species, phenotype):
    """Return True if drug_agent provides adequate coverage for (species, phenotype)."""
    if drug_agent is None:
        return False
    if species == 'none':
        return True
    if drug_agent == 'last-resort':
        return True
    return (species, phenotype) in AGENT_COVERAGE.get(drug_agent, set())


def min_covering_agent(species, phenotype):
    """
    Return the lowest-level non-last-resort agent that covers (species, phenotype).
    Falls back to 'last-resort' if no other agent covers it.
    """
    best_agent = None
    best_level = 999
    for agent, coverage in AGENT_COVERAGE.items():
        if (species, phenotype) in coverage:
            level = AGENT_LEVEL[agent]
            if level < best_level:
                best_level = level
                best_agent = agent
    return best_agent if best_agent is not None else 'last-resort'
