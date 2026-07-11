#!/usr/bin/env python3
"""
Read the two ESKAPE workbooks and emit src/config_data.py.

The workbooks are the human-editable source of truth. This script extracts
the coverage matrix, tier ranking, phenotype lists per species, and the
drug→resistance-phenotype mapping into a generated Python module that
config.py imports at runtime.

Re-run this script whenever the workbooks change:

    python build_config.py
"""

from __future__ import annotations

import os
from pathlib import Path
from pprint import pformat

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent
COVERAGE_XLSX = REPO_ROOT / 'ESKAPE_6_Species_CrossResistance_Coverage.xlsx'
TIER_XLSX     = REPO_ROOT / 'Global_20_Drug_Tier_Ranking_MultiSheet.xlsx'
OUT_PATH      = REPO_ROOT / 'src' / 'config_data.py'


# Excel sheet name -> internal species key. Stable across re-runs.
SHEET_TO_SPECIES = {
    'K_pneumoniae':         'kpneumoniae',
    'Enterobacter_cloacae': 'enterobacter',
    'P_aeruginosa':         'paeruginosa',
    'A_baumannii':          'abaumannii',
    'S_aureus':             'saureus',
    'E_faecium':            'efaecium',
}

# Full species name for documentation / output.
SPECIES_FULL_NAME = {
    'efaecium':     'Enterococcus faecium',
    'saureus':      'Staphylococcus aureus',
    'kpneumoniae':  'Klebsiella pneumoniae',
    'abaumannii':   'Acinetobacter baumannii',
    'paeruginosa':  'Pseudomonas aeruginosa',
    'enterobacter': 'Enterobacter cloacae',
}

# Per-species colonization→infection probability (carried over from the
# previous config — flat across phenotypes for now).
SPECIES_SIGMA = {
    'efaecium':     0.08,
    'saureus':      0.04,
    'kpneumoniae':  0.22,
    'abaumannii':   0.278,
    'paeruginosa':  0.23,
    'enterobacter': 0.03,
}

# Admission species distribution (carried from previous config). Sums to 1.
INIT_SPECIES_PROBS = {
    'none':         0.554,
    'saureus':      0.300,
    'kpneumoniae':  0.050,
    'paeruginosa':  0.026,
    'efaecium':     0.030,
    'abaumannii':   0.030,
    'enterobacter': 0.010,
}

# Drug name (verbatim from coverage workbook) → matching resistance-phenotype
# column label. Drugs whose abbreviated label differs from the full name
# need an explicit entry.
DRUG_TO_R_PHENOTYPE = {
    'Piperacillin-Tazobactam': 'Pip-Tazo-R',
    'Cefepime':                'Cefepime-R',
    'Ceftazidime':             'Ceftazidime-R',
    'Ceftriaxone':             'Ceftriaxone-R',
    'Aztreonam':               'Aztreonam-R',
    'Meropenem':               'Meropenem-R',
    'Imipenem':                'Imipenem-R',
    'Ertapenem':               'Ertapenem-R',
    'Ceftazidime-Avibactam':   'CZA-R',
    'Meropenem-Vaborbactam':   'MVB-R',
    'Ceftolozane-Tazobactam':  'C/T-R',
    'Cefiderocol':             'Cefiderocol-R',
    'Ampicillin-Sulbactam':    'Ampicillin-Sulbactam-R',
    'Sulbactam-Durlobactam':   'Sulbactam-Durlobactam-R',
    'Amikacin':                'Amikacin-R',
    'Ciprofloxacin':           'Ciprofloxacin-R',
    'Tigecycline':             'Tigecycline-R',
    'Eravacycline':            'Eravacycline-R',
    'Colistin':                'Colistin-R',
    'Vancomycin':              'Vancomycin-R',
}

# Map workbook glyph → coverage state string. The simulation works with
# 'covers' / 'partial' / 'none' rather than emoji to keep the code grep-friendly.
GLYPH_TO_STATE = {
    '🟢': 'covers',
    '🟡': 'partial',
    '❌': 'none',
}


# ─── Extractors ──────────────────────────────────────────────────────────────

def load_coverage(path: Path) -> tuple[dict, dict]:
    """Return (PATHOGENS, COVERAGE_MATRIX) from the coverage workbook.

    PATHOGENS[species] = {full_name, phenotypes, sigma, treatment_days}
    COVERAGE_MATRIX[species][drug][phenotype] = 'covers' | 'partial' | 'none'
    """
    xl = pd.ExcelFile(path)
    pathogens: dict = {}
    coverage:  dict = {}

    for sheet, sp_key in SHEET_TO_SPECIES.items():
        if sheet not in xl.sheet_names:
            raise ValueError(f'coverage workbook missing sheet: {sheet}')
        df = xl.parse(sheet)
        phenotypes = list(df.columns[1:])  # column 0 is 'Drug'

        pathogens[sp_key] = {
            'full_name':      SPECIES_FULL_NAME[sp_key],
            'phenotypes':     phenotypes,
            'sigma':          SPECIES_SIGMA[sp_key],
            'treatment_days': _treatment_days_table(phenotypes),
        }

        per_drug: dict = {}
        for _, row in df.iterrows():
            drug = row['Drug']
            per_phen = {}
            for ph in phenotypes:
                glyph = row[ph]
                if glyph not in GLYPH_TO_STATE:
                    raise ValueError(
                        f'unknown coverage glyph {glyph!r} in {sheet} '
                        f'drug={drug} phenotype={ph}')
                per_phen[ph] = GLYPH_TO_STATE[glyph]
            per_drug[drug] = per_phen
        coverage[sp_key] = per_drug

    return pathogens, coverage


def _treatment_days_table(phenotypes: list[str]) -> dict[str, int]:
    """Treatment days by phenotype. Susceptible=7; any -R phenotype=10.

    Multi-drug-resistance state (resistance set with >=2 entries) uses
    14 days; that's computed at simulation time, not stored per-phenotype.
    """
    days = {}
    for ph in phenotypes:
        if ph in ('Susceptible', 'Std'):
            days[ph] = 7
        else:
            days[ph] = 10
    return days


def load_tiers(path: Path) -> dict[str, int]:
    """Return DRUG_TIER[drug] = 1 | 2 | 3."""
    xl = pd.ExcelFile(path)
    tier_map: dict = {}
    for tier_num, sheet in [(1, 'Tier 1 Drugs'), (2, 'Tier 2 Drugs'), (3, 'Tier 3 Drugs')]:
        df = xl.parse(sheet)
        for drug in df['Drug']:
            # Strip parenthetical alias e.g. "Ceftazidime-Avibactam (CZA)"
            name = drug.split('(')[0].strip()
            tier_map[name] = tier_num
    return tier_map


def scaffold_init_phenotype_probs(pathogens: dict) -> dict:
    """Placeholder admission phenotype distribution: 90% Susceptible, 10%
    spread uniformly across all -R phenotypes per species. TODO replace with
    surveillance data."""
    probs: dict = {}
    for sp_key, cfg in pathogens.items():
        phenotypes = cfg['phenotypes']
        # Susceptible label is the first column ('Susceptible' for most,
        # 'Std' for older sheets — only 'Susceptible' here).
        susc_label = phenotypes[0]
        r_phenotypes = phenotypes[1:]
        per_species = {susc_label: 0.90}
        if r_phenotypes:
            share = 0.10 / len(r_phenotypes)
            for ph in r_phenotypes:
                per_species[ph] = share
        probs[sp_key] = per_species
    return probs


# ─── Emitter ─────────────────────────────────────────────────────────────────

HEADER = '''"""
Auto-generated by build_config.py — DO NOT EDIT BY HAND.

Run `python build_config.py` after editing either ESKAPE workbook to
regenerate this module.

Source workbooks:
    {coverage_src}
    {tier_src}
"""

# fmt: off
'''


def emit(path: Path, pathogens, coverage, tier, init_phenotype_probs):
    parts = [HEADER.format(coverage_src=COVERAGE_XLSX.name, tier_src=TIER_XLSX.name)]

    parts.append('PATHOGENS = ' + pformat(pathogens, sort_dicts=False, width=100) + '\n\n')

    parts.append(
        'INIT_SPECIES_PROBS = '
        + pformat(INIT_SPECIES_PROBS, sort_dicts=False, width=100) + '\n\n'
    )

    parts.append(
        '# TODO: replace placeholder distribution with ICU surveillance data.\n'
        '# Currently 90% Susceptible per species, 10% spread uniformly over R-phenotypes.\n'
        'INIT_PHENOTYPE_PROBS = '
        + pformat(init_phenotype_probs, sort_dicts=False, width=100) + '\n\n'
    )

    parts.append('DRUG_TIER = ' + pformat(tier, sort_dicts=False, width=100) + '\n\n')

    parts.append(
        'DRUG_TO_R_PHENOTYPE = '
        + pformat(DRUG_TO_R_PHENOTYPE, sort_dicts=False, width=100) + '\n\n'
    )

    parts.append(
        '# COVERAGE_MATRIX[species][drug][phenotype] -> "covers" | "partial" | "none"\n'
        'COVERAGE_MATRIX = '
        + pformat(coverage, sort_dicts=False, width=120) + '\n'
    )

    path.write_text(''.join(parts))


# ─── Validation ──────────────────────────────────────────────────────────────

def validate(pathogens, coverage, tier):
    """Cross-check that all drugs in the coverage matrix have tiers and that
    every drug→R-phenotype mapping points at a real column somewhere."""
    drugs_in_coverage = set()
    for sp_key, per_drug in coverage.items():
        for drug in per_drug:
            drugs_in_coverage.add(drug)

    missing_tier = drugs_in_coverage - set(tier)
    if missing_tier:
        raise ValueError(f'drugs in coverage with no tier assignment: {sorted(missing_tier)}')

    missing_cov = set(tier) - drugs_in_coverage
    if missing_cov:
        print(f'[WARN] tiered drugs not in any coverage sheet: {sorted(missing_cov)}')

    # Audit which (drug, species) pairs have a corresponding R-phenotype column.
    can_mutate = []
    no_mutation = []
    for drug, r_label in DRUG_TO_R_PHENOTYPE.items():
        for sp_key, cfg in pathogens.items():
            if r_label in cfg['phenotypes']:
                can_mutate.append((drug, sp_key))
            else:
                no_mutation.append((drug, sp_key))
    print(f'[info] {len(can_mutate)} (drug, species) pairs can fire mutation; '
          f'{len(no_mutation)} cannot (no matching R-phenotype column).')


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    if not COVERAGE_XLSX.exists():
        raise FileNotFoundError(f'missing coverage workbook: {COVERAGE_XLSX}')
    if not TIER_XLSX.exists():
        raise FileNotFoundError(f'missing tier workbook: {TIER_XLSX}')

    pathogens, coverage = load_coverage(COVERAGE_XLSX)
    tier                = load_tiers(TIER_XLSX)
    init_phen_probs     = scaffold_init_phenotype_probs(pathogens)

    validate(pathogens, coverage, tier)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    emit(OUT_PATH, pathogens, coverage, tier, init_phen_probs)
    print(f'[ok] wrote {OUT_PATH}  ({len(pathogens)} species, '
          f'{len(tier)} drugs, {sum(len(p["phenotypes"]) for p in pathogens.values())} '
          f'(species, phenotype) states)')


if __name__ == '__main__':
    main()
