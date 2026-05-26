# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Stochastic Agent-Based Model (ABM) of multi-pathogen ESKAPE dynamics and antimicrobial stewardship in a 64-bed ICU. The simulation models transmission of the full ESKAPE pathogen group (E. faecium, S. aureus, K. pneumoniae, A. baumannii, P. aeruginosa, Enterobacter spp.) between patients and healthcare workers (HCWs), using a 3-level WHO AWaRe drug classification (Access / Watch / Reserve). Drug assignment is agent-specific. Empiric therapy is assigned by species at infection onset using one of two agent groups, selected via a configurable regimen (fixed, cycling, or mixing).

## Setup

See [`README.md`](./README.md#setup-on-a-new-machine) for the venv + `pip install -r requirements.txt` flow. `run_experiments.py` auto-detects `.venv/bin/python` (Unix) or `.venv\Scripts\python.exe` (Windows) at the repo root, so use that interpreter for all commands below.

## Running the Simulation

```bash
# Run from src/ directory (CSV paths are relative to src/)
cd src

# Standard run (fixed empiric regimen, 72h review)
python main.py --num_runs 1000

# Quick test run (10 replicates, 30 days)
python main.py --num_runs 10 --days 30

# Cycling regimen — alternate group 1/2 every 3 months, 48h review
python main.py --empiric_regimen cycling --cycle_months 3 --empiric_hours 48

# Mixing regimen — 50/50 random group assignment per patient
python main.py --empiric_regimen mixing

# With stochastic parameter noise
python main.py --num_runs 1000 --std 0.1
```

Output is written to `./log/{regimen}_h{empiric_hours}_p{p}_q{q}_r{r}_n{num_patient}_w{num_hcw}/` as timestamped `.pkl` files.

## Running Batch Experiments

`run_experiments.py` (project root) orchestrates the full 3 regimens × 36 (p, q) combinations = 108 simulations:

```bash
# Run from project root
python run_experiments.py
```

Key design:
- All 108 jobs are submitted as threads at startup; a semaphore (`MAX_PARALLEL=20`) limits concurrent processes to 20, acting as a rolling pool — a new job starts the moment any slot frees.
- **Skip logic**: jobs whose output directory already contains a `.pkl` file are skipped, enabling safe resume after interruption.
- **Race-condition safe**: a `claimed` set (guarded by `done_lock`) prevents two threads from claiming the same job simultaneously.
- Output lands in `src/log/{RUN_NAME}/` — set `RUN_NAME` at the top of `run_experiments.py` before each batch (e.g. `'run2'`, `'run3'`).
- `r = p` and `s = 0.1 × p` — both derived from `p` and passed explicitly.
- **Path resolution is automatic**: `REPO_ROOT`, `SRC_DIR`, `MAIN` are computed from `__file__`; the interpreter is auto-resolved as `ABM_PYTHON` env var → `<repo>/.venv/bin/python` (or `.venv\Scripts\python.exe` on Windows) → `sys.executable`. `ABM_SRC_DIR` and `ABM_MAIN` are available as one-off overrides. No per-machine edits required.
- **Monitor thread**: prints progress, rate, ETA, process count, and RAM every 2 minutes to the log file. Memory probe is cross-platform — `ps` on Linux/macOS, `tasklist` on Windows.

## Key Parameters

| Arg | Default | Description |
|-----|---------|-------------|
| `--days` | 1460 | Simulation duration (4 years) |
| `--num_runs` | 1000 | Monte Carlo replicates |
| `--num_patient` | 64 | ICU bed capacity (must be divisible by `num_hcw`) |
| `--num_hcw` | 16 | Healthcare workers |
| `--p` | 0.30 | HCW→patient transmission probability |
| `--q` | 0.30 | Patient→HCW contamination probability |
| `--r` | 0.30 | Super-infection probability, early phase; set to `p` in batch runs |
| `--s` | 0.015 | Super-infection probability, late phase; set to `0.1 × p` in batch runs |
| `--eta` | 0.50 | Hand-hygiene compliance |
| `--epsilon` | 0.03 | Intrinsic mutation hazard |
| `--time_interval` | 2 | Hours between HCW patient visits within a shift |
| `--kappa_mu` | 0.74 | Discharge hazard ratio for infected patients |
| `--kappa_nu` | 1.04 | Death hazard multiplier for inadequately treated patients |
| `--empiric_regimen` | `fixed` | Empiric group selection: `fixed`, `cycling`, `mixing` |
| `--cycle_months` | 6 | Months per group in cycling regimen |
| `--empiric_hours` | 72 | Hours before drug review (24, 48, or 72) |
| `--std` | 0.0 | Parameter noise for Monte Carlo (0 = deterministic) |

When `--std > 0`, the per-run perturbation in `main.py` covers: `p`, `q`, `r`, `s`, `epsilon`, `eta`, `kappa_mu`, `kappa_nu`.

## Architecture

### Source Files (`src/`)

**`config.py`** — Single source of truth for all biological parameters. Contains:
- `PATHOGENS` — ESKAPE species registry: phenotypes (susceptible + 2 resistance phenotypes + dual), sigma (colonization→infection probability — keyed per phenotype, but currently set to a single per-species value across all four phenotypes), treatment durations (7 / 10 / 14 days).
- `DRUG_LEVELS` — 3-level AWaRe classification: Access (1), Watch (2), Reserve (3).
- `AGENT_LEVEL` — Maps each named agent to its AWaRe level.
- `AGENT_COVERAGE` — Maps each agent to the set of `(species, phenotype)` it covers. `last-resort` covers all phenotypes.
- `EMPIRIC_AGENTS_1` — Group 1 empiric agents (Watch-level; Vancomycin for Gram+, Pip-Tazo/Meropenem for Gram−).
- `EMPIRIC_AGENTS_2` — Group 2 empiric agents (mixed levels; Linezolid, Nafcillin, Cefepime, Amp-sulbactam, Pip-Tazo, Cefepime).
- `MUTATION_PATHWAYS` — Within-species resistance escalation graph with branching paths, e.g. `kpneumoniae susceptible → [PTR-KP, CRKP]`. Each entry is a list; one next phenotype is chosen at random when mutation fires.
- `INIT_SPECIES_PROBS` — Admission species distribution including `'none'` (uncolonized). Values sum to 1.
- `INIT_PHENOTYPE_PROBS` — Within-species phenotype distribution for admitted patients.
- `is_covered(drug_agent, species, phenotype)` / `min_covering_agent(species, phenotype)` — helper functions used throughout.

**`abm.py`** — Entity classes:
- `Patient` — tracks `species`, `phenotype`, `drug_agent` (str, name of current treating agent), `treat_time`, `convt_time`, `time_inICU`, `super_infe`, `infct_flag/infct_time`. Per-species sigma is looked up from `config.PATHOGENS` inside `_get_sigma()`.
- `HealthCareWorker` — `strain_set` holds `(species, phenotype)` tuples.

**`scheme_func.py`** — All simulation event functions:
- `get_empiric_agent(species, args, current_day)` — selects empiric agent based on regimen: `fixed` always uses group 1; `cycling` alternates group 1/2 every `cycle_months` months; `mixing` randomly selects group 1 or 2 with equal probability.
- `ph_interaction` — dispatches contact events based on patient state (colonised / early / late infected). Early/late boundary is `empiric_hours`.
- `infection_development` — colonised → infected when `infct_time` is reached; calls `get_empiric_agent()` to assign drug.
- `drug_change` — fires at `treat_time == empiric_hours` and at `convt_time + empiric_hours`. Keeps empiric agent if it covers the identified phenotype; otherwise escalates to `min_covering_agent()`, using `last-resort` only if necessary.
- `intrinsic_mutation` — resistance emerges when drug covers current phenotype (selection pressure); advances one step along `MUTATION_PATHWAYS` by random choice.
- `discharge_admission`, `death_event` — stochastic patient turnover using CSV hazard tables; inadequate therapy multiplies death hazard by `kappa_nu`.

**`main.py`** — Simulation orchestrator:
- `main(args)` runs one replicate.
- Outer loop runs `num_runs` Monte Carlo replicates with optional per-run parameter noise (`truncnorm_func`).
- `make_reference_results()` builds the metric counter dict from the ESKAPE registry.
- All daily output arrays are pre-allocated `numpy.ndarray` (int32/float32) — no Python list growth.
- Colonization prevalence uses a single patient pass per day (not 24 separate scans).
- No scipy dependency: `read_csv()` uses `value * key` arithmetic; `truncnorm_func()` uses numpy rejection sampling.
- `--run_name` arg inserts a subfolder under `log/` (e.g. `run2`, `run3`) for batch organisation.
- `hcw_contamination_burden` is measured as the average strains per HCW across all 3 shifts, snapshotted **before** each end-of-shift cleanup.

### Simulation Flow

Each simulated day consists of 3 × 8-hour shifts. Within each shift, HCWs visit their assigned patients. At shift end, all HCWs are decontaminated. End-of-day events in order: infection development → intrinsic mutation → increment time counters → drug change → treatment completion → discharge → death.

### AWaRe Drug Levels and Agents

| Level | Label | Agent | Covers |
|-------|-------|-------|--------|
| 1 | Access | Nafcillin/Oxacillin | S. aureus susceptible, VRSA |
| 1 | Access | Ampicillin-sulbactam | A. baumannii susceptible, CRAB |
| 2 | Watch | Vancomycin | E. faecium susceptible, S. aureus susceptible/MRSA |
| 2 | Watch | Pip-Tazo | K. pneumoniae susceptible/CRKP, P. aeruginosa susceptible/CRPA |
| 2 | Watch | Cefepime | K. pneumoniae susceptible/PTR-KP, Enterobacter susceptible/CRE |
| 2 | Watch | Meropenem | A. baumannii susceptible/SAMR-AB, P. aeruginosa susceptible/PTR-PA, Enterobacter susceptible/ESBL |
| 3 | Reserve | Linezolid | E. faecium susceptible/LRE/VREfm |
| 3 | Reserve | Last-resort | All phenotypes |

### Empiric Agent Groups

| Species | Group 1 | Group 2 |
|---------|---------|---------|
| *E. faecium* | Vancomycin (2) | Linezolid (3) |
| *S. aureus* | Vancomycin (2) | Nafcillin (1) |
| *K. pneumoniae* | Pip-Tazo (2) | Cefepime (2) |
| *A. baumannii* | Meropenem (2) | Amp-sulbactam (1) |
| *P. aeruginosa* | Meropenem (2) | Pip-Tazo (2) |
| *Enterobacter* | Meropenem (2) | Cefepime (2) |
| Fallback | Meropenem | Meropenem |

### Empiric Regimens

| Regimen | Behavior |
|---------|----------|
| `fixed` | Always group 1 (default) |
| `cycling` | Group 1 for `cycle_months` months, then group 2, repeating |
| `mixing` | Group 1 or 2 chosen at random (50/50) per patient at infection onset |

### Drug Change Rule

At `empiric_hours` (24 / 48 / 72 h), for each infected patient:
1. If current agent **covers** the identified phenotype → **keep** agent
2. If not → escalate to **lowest-level non-last-resort** covering agent
3. If no such agent exists → give **last-resort**

Re-evaluates at `convt_time + empiric_hours` after any strain conversion.

### Phenotypes and Mutation Pathways

Each species has: `susceptible`, two resistance phenotypes, and `dual` (resistant to both).

| Species | Resistance phenotypes | Mutation paths |
|---------|----------------------|----------------|
| *E. faecium* | VREfm, LRE | susceptible → VREfm or LRE → dual |
| *S. aureus* | VRSA, MRSA | susceptible → VRSA or MRSA → dual |
| *K. pneumoniae* | PTR-KP, CRKP | susceptible → PTR-KP or CRKP → dual |
| *A. baumannii* | CRAB, SAMR-AB | susceptible → CRAB or SAMR-AB → dual |
| *P. aeruginosa* | PTR-PA, CRPA | susceptible → PTR-PA or CRPA → dual |
| *Enterobacter* | CRE, ESBL | susceptible → CRE or ESBL → dual |

Treatment days: susceptible = 7, resistance phenotypes = 10, dual = 14.

### Adding New Pathogens or Drugs

1. Add species to `config.PATHOGENS` with `phenotypes`, `sigma`, and `treatment_days`.
2. Add agent to `config.AGENT_LEVEL` and `config.AGENT_COVERAGE`.
3. Add empiric agents to `config.EMPIRIC_AGENTS_1` and `config.EMPIRIC_AGENTS_2`.
4. Add mutation edges to `config.MUTATION_PATHWAYS`.
5. Add initial prevalence to `config.INIT_SPECIES_PROBS` and `config.INIT_PHENOTYPE_PROBS`.
6. No changes needed to `abm.py`, `scheme_func.py`, or `main.py`.

### Output

Results written to `./log/{regimen}_h{empiric_hours}_p{p}_q{q}_r{r}_n{N}_w{N}/` as `.pkl` files. Batch runs are organised into `src/log/run1/`, `src/log/run2/`, etc. Each file is a dict `{run_index: final_results}` where `final_results` contains:

**Daily prevalence arrays** (length = days):
- `infection_{species}_{phenotype}` — daily infected patient count per phenotype
- `colonization_prev_{species}_{phenotype}` — daily colonized carrier count per phenotype
- `superinfection` — daily super-infection count
- `hcw_contamination_burden` — daily avg strains per HCW, averaged across 3 shifts (measured pre-cleanup each shift)

**Cumulative time-series** (daily snapshot of running total):
- `cuminfection_*`, `labresult_*`, `colonization_*`, `transmission_*`, `mutation_*`
- `admission`, `discharge`, `death`
- `misempiric`, `drug_agent_{agent}_use`

**Daily reset counters** (actual events that day, not cumulative):
- `empiric_group_1_daily`, `empiric_group_2_daily` — patients assigned each empiric group
- `misempiric_daily` — empiric failures per day

**Per-patient lists** (one entry per discharge or death):
- `treatment_length_on_discharge` — total hours on treatment
- `icu_los_on_discharge` — total ICU hours at exit

**Per-run scalar dict**:
- `first_mutation_day` — `{sp_ph: day}` of first mutation to each phenotype (`None` if never)

### Data Files

- `deathprob.csv` / `dischargeprob.csv` — Day-indexed probability tables. Loaded by `read_csv()` in `main.py`; hazard integral computed as `value * key` (integrating a constant).

## Calibration

### Metric

For each run independently, compute the **infection-based within-species phenotype fraction** over the last 365 days:

```
fraction(sp, ph) = sum(infection_{sp}_{ph}[-365:]) / sum(infection_{sp}_{all_ph}[-365:])
```

Runs where a species produces zero infections in the last 365 days are excluded (return None). The reported calibration statistic is the **% of 1000 runs where this fraction falls within the observed target range**.

### Observed Target Ranges

| Phenotype | Target range |
|-----------|-------------|
| efaecium VREfm | 19.4 – 59.6% |
| efaecium LRE | < 2% |
| saureus VRSA | < 2.9% |
| saureus MRSA | 6.5 – 81.5% |
| kpneumoniae PTR-KP | 42.6 – 92.9% |
| abaumannii CRAB | > 26.1% |
| abaumannii SAMR-AB | 48 – 92.1% |
| paeruginosa PTR-PA | 13 – 53% |
| paeruginosa CRPA | 20.4 – 50% |
| enterobacter CRE | < 37.8% |
| enterobacter ESBL | 13.1 – 70.5% |

### Calibration Status (run3)

Overall pass rate (all 11 constraints simultaneously): **0% across all 36 (p,q) pairs and all 3 regimens**.

Varying p,q has no effect on the 6 failing phenotypes — the steady-state resistance fraction is anchored by `INIT_PHENOTYPE_PROBS`, not transmission rates. Updating `INIT_PHENOTYPE_PROBS` in `config.py` is required before re-running.

### Calibration Script

Save as `analyze_calibration.py` in the project root. See memory file `session_20260406_analysis.md` for the full script.
