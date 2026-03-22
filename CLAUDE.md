# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Stochastic Agent-Based Model (ABM) of multi-pathogen ESKAPE dynamics and antimicrobial stewardship in a 64-bed ICU. The simulation models transmission of the full ESKAPE pathogen group (E. faecium, S. aureus, K. pneumoniae, A. baumannii, P. aeruginosa, Enterobacter spp.) between patients and healthcare workers (HCWs), using the WHO AWaRe drug classification (levels 1–5). Two trial arms are implemented: ADE (de-escalation) and control (broad-spectrum).

## Running the Simulation

```bash
# Run from src/ directory (CSV paths are relative to src/)
cd src

# ADE arm — de-escalate at 72h to minimum effective AWaRe level
python main.py --arm ADE --num_runs 1000

# Control arm — stay on broad-spectrum unless Reserve is required
python main.py --arm control --num_runs 1000

# Quick test run (10 replicates, 30 days)
python main.py --arm ADE --num_runs 10 --days 30

# With stochastic parameter noise
python main.py --arm ADE --num_runs 1000 --std 0.1
```

Output is written to `./log/{arm}_p{num_patient}_h{num_hcw}/` as timestamped `.pkl` files.

## Key Parameters

| Arg | Default | Description |
|-----|---------|-------------|
| `--arm` | `ADE` | `ADE` or `control` trial arm |
| `--days` | 1460 | Simulation duration (4 years) |
| `--num_runs` | 1000 | Monte Carlo replicates |
| `--num_patient` | 64 | ICU bed capacity (must be divisible by `num_hcw`) |
| `--num_hcw` | 16 | Healthcare workers |
| `--p` | 0.30 | HCW→patient transmission probability |
| `--q` | 0.30 | Patient→HCW contamination probability |
| `--eta` | 0.50 | Hand-hygiene compliance |
| `--a` | 0.10 | Admission colonization probability |
| `--epsilon` | 0.03 | Intrinsic mutation hazard |
| `--std` | 0.0 | Parameter noise for Monte Carlo (0 = deterministic) |

## Architecture

### Source Files (`src/`)

**`config.py`** — Single source of truth for all biological parameters. Contains:
- `PATHOGENS` — ESKAPE species registry: phenotypes, per-phenotype sigma (colonization→infection probability), treatment durations.
- `COVERAGE_REQUIREMENTS` — Minimum AWaRe drug level needed per (species, phenotype). Drives all drug-effectiveness logic. Edit here when adding new pathogens or drugs.
- `MUTATION_PATHWAYS` — Within-species resistance escalation graph, e.g. `kpneumoniae susceptible → ESBL → KPC`.
- `INIT_SPECIES_PROBS` / `INIT_PHENOTYPE_PROBS` — Admission colonization distribution (from ICU surveillance data).
- `min_drug_level(species, phenotype)` / `is_covered(level, species, phenotype)` — helper functions used throughout.

**`abm.py`** — Entity classes:
- `Patient` — tracks `species`, `phenotype`, `drug_level` (int 1–5), `treat_time`, `convt_time`, `time_inICU`, `super_infe`, `infct_flag/infct_time`. Per-species sigma is looked up from `config.PATHOGENS` inside `_get_sigma()`.
- `HealthCareWorker` — `strain_set` holds `(species, phenotype)` tuples.

**`scheme_func.py`** — All simulation event functions:
- `ph_interaction` — dispatches contact events based on patient state (colonised / early infected / late infected).
- `infection_development` — colonised → infected when `infct_time` is reached; assigns empiric level 4.
- `drug_change_ADE` / `drug_change_control` — the only functions that differ between arms. Both fire at `treat_time == 72h` and at `convt_time + 72h`.
- `intrinsic_mutation` — resistance emerges when drug **covers** the current phenotype (selection pressure), advancing one step along `MUTATION_PATHWAYS`.
- `discharge_admission`, `death_event` — stochastic patient turnover using CSV hazard tables; inadequate therapy (drug doesn't cover pathogen) multiplies death hazard by `kappa_nu`.

**`main.py`** — Simulation orchestrator:
- `main(args, drug_change_func)` runs one replicate; `drug_change_func` is either `drug_change_ADE` or `drug_change_control`.
- Outer loop runs `num_runs` Monte Carlo replicates with optional per-run parameter noise (`truncnorm_func`).
- `make_reference_results()` builds the metric counter dict from the ESKAPE registry — adding a pathogen in `config.py` automatically extends tracking.

### Simulation Flow

Each simulated day consists of 3 × 8-hour shifts. Within each shift, HCWs visit their assigned patients. At shift end, all HCWs are decontaminated. End-of-day events in order: infection development → intrinsic mutation → increment time counters → drug change → treatment completion → discharge → death.

### AWaRe Drug Levels

| Level | Category | Typical Agents | Covers |
|-------|----------|---------------|--------|
| 1 | Access/Narrow | Nafcillin, Pen-G | MSSA only |
| 2 | Access/Moderate | Amp-Sulbactam, Gent | Susceptible GN, VSE |
| 3 | Watch/Medium | Ceftriaxone, Cipro | Susceptible Enterobacteriaceae |
| 4 | Watch/Broad | Pip-Tazo, Cefepime, Meropenem | ESBL-Kp, PA susceptible, A. baumannii susceptible |
| 5 | Reserve | CZA, Colistin, Linezolid | KPC, CRAB, MRSA, VRE, CRPA, DTR |

### Trial Arms

- **ADE arm**: At 72 h, de-escalates to `min_drug_level(lab_result)` — the narrowest agent that covers the cultured pathogen.
- **Control arm**: At 72 h, stays at level 4 unless the identified pathogen requires level 5 (Reserve).

The "selective window" effect emerges naturally: ADE de-escalation may leave gaps that allow a secondary pathogen (not covered by the narrower drug) to colonise or super-infect.

### Adding New Pathogens or Drugs

1. Add the species to `config.PATHOGENS` with `phenotypes`, `sigma`, and `treatment_days`.
2. Add all `(species, phenotype)` entries to `config.COVERAGE_REQUIREMENTS`.
3. Add mutation edges to `config.MUTATION_PATHWAYS` if applicable.
4. Add initial prevalence to `config.INIT_SPECIES_PROBS` and `config.INIT_PHENOTYPE_PROBS`.
5. No changes needed to `abm.py`, `scheme_func.py`, or `main.py` — all logic is data-driven.

### Output

Results written to `./log/{arm}_p{N}_h{N}/` as `.pkl` files. Each file is a dict `{run_index: final_results}` where `final_results` contains:
- Daily arrays: `infection_{species}_{phenotype}`, `superinfection`
- Cumulative time-series: `cuminfection_*`, `labresult_*`, `colonization_*`, `transmission_*`, `mutation_*`
- Clinical outcomes: `admission`, `discharge`, `death`
- Stewardship metrics: `deescalation`, `escalation`, `misempiric`, `drug_level_{1-5}_use`

### Data Files

- `deathprob.csv` / `dischargeprob.csv` — Day-indexed probability tables. Loaded by `read_csv()` in `main.py` and integrated via `scipy.integrate.quad` to produce hazard-based CDFs.
