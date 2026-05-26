# ABM-Antimicrobial

Stochastic agent-based model of multi-pathogen ESKAPE dynamics and antimicrobial stewardship in a 64-bed ICU. Models transmission of the full ESKAPE group (*E. faecium*, *S. aureus*, *K. pneumoniae*, *A. baumannii*, *P. aeruginosa*, *Enterobacter* spp.) between patients and healthcare workers, using the 3-level WHO AWaRe drug classification (Access / Watch / Reserve). Empiric therapy is assigned by species via a configurable regimen (fixed, cycling, or mixing).

Full model design notes live in [`CLAUDE.md`](./CLAUDE.md).

## Requirements

- Python 3.10+ (3.12 recommended — matches the development environment)
- The only third-party dependency is `numpy` (see [`requirements.txt`](./requirements.txt)). Everything else is the Python standard library.

## Setup on a new machine

```bash
# 1. Clone
git clone https://github.com/pliu19/AgenticSimulationPathogen.git
cd AgenticSimulationPathogen

# 2. Create and activate a virtual environment
python -m venv .venv

# Windows (PowerShell)
.\.venv\Scripts\Activate.ps1
# Windows (cmd)
.\.venv\Scripts\activate.bat
# macOS / Linux
source .venv/bin/activate

# 3. Install dependencies
python -m pip install --upgrade pip
pip install -r requirements.txt
```

To leave the venv later: `deactivate`.

## Running a single simulation

CSV hazard tables (`deathprob.csv`, `dischargeprob.csv`) are loaded with paths relative to `src/`, so run from there:

```bash
cd src

# Standard run (fixed empiric regimen, 72h review, 1000 replicates)
python main.py --num_runs 1000

# Quick test (10 replicates, 30 days)
python main.py --num_runs 10 --days 30

# Cycling regimen — alternate group 1/2 every 3 months, 48h review
python main.py --empiric_regimen cycling --cycle_months 3 --empiric_hours 48

# Mixing regimen — 50/50 random group assignment per patient
python main.py --empiric_regimen mixing

# With stochastic per-run parameter noise
python main.py --num_runs 1000 --std 0.1
```

Output `.pkl` files land in `src/log/{regimen}_h{empiric_hours}_p{p}_q{q}_r{r}_n{num_patient}_w{num_hcw}/`.

## Running batch experiments

`run_experiments.py` (project root) orchestrates the full sweep of 3 regimens × 36 (p, q) pairs = 108 simulations:

```bash
# From the project root, with venv active
python run_experiments.py
```

Key knobs at the top of the script:

- `RUN_NAME` — output subfolder under `src/log/` (e.g. `'run2'`, `'run3'`); set this before each batch.
- `MAX_PARALLEL` — semaphore-limited concurrent process count (default 20).

Paths are derived from the script location, so no edits are needed per machine. The interpreter is resolved in this order: `ABM_PYTHON` env var → `<repo>/.venv/bin/python` (or `.venv\Scripts\python.exe` on Windows) → the current `python`. `ABM_SRC_DIR` and `ABM_MAIN` env vars are available for one-off overrides.

Jobs whose output directory already contains a `.pkl` are skipped, so the batch can resume safely after interruption. A monitor thread prints progress, rate, ETA, process count, and RAM every 2 minutes (cross-platform: `ps` on Linux/macOS, `tasklist` on Windows).

## Key parameters

| Arg | Default | Description |
|-----|---------|-------------|
| `--days` | 1460 | Simulation duration (4 years) |
| `--num_runs` | 1000 | Monte Carlo replicates |
| `--num_patient` | 64 | ICU bed capacity (must be divisible by `--num_hcw`) |
| `--num_hcw` | 16 | Healthcare workers |
| `--p` | 0.30 | HCW→patient transmission probability |
| `--q` | 0.30 | Patient→HCW contamination probability |
| `--r` | 0.30 | Super-infection probability, early phase (set to `p` in batch runs) |
| `--s` | 0.015 | Super-infection probability, late phase (set to `0.1 × p` in batch runs) |
| `--eta` | 0.50 | Hand-hygiene compliance |
| `--epsilon` | 0.03 | Intrinsic mutation hazard |
| `--empiric_regimen` | `fixed` | `fixed`, `cycling`, or `mixing` |
| `--cycle_months` | 6 | Months per group in cycling regimen |
| `--empiric_hours` | 72 | Hours before drug review (24, 48, or 72) |
| `--std` | 0.0 | Parameter noise for Monte Carlo (0 = deterministic) |
| `--run_name` | — | Optional subfolder under `src/log/` for batch organisation |

See [`CLAUDE.md`](./CLAUDE.md) for architecture, AWaRe drug levels, empiric groups, mutation pathways, output schema, and calibration notes.

## Repository layout

```
ABM-Antimicrobial/
├── README.md
├── CLAUDE.md              # model design + architecture notes
├── EXPANSION_SUMMARY.md
├── requirements.txt
├── run_experiments.py     # batch orchestrator (3 regimens × 36 (p,q))
├── deathprob.csv          # daily death-hazard table
├── dischargeprob.csv      # daily discharge-hazard table
└── src/
    ├── main.py            # simulation orchestrator (single-replicate entry point)
    ├── abm.py             # Patient + HealthCareWorker classes
    ├── scheme_func.py     # all event functions (interaction, mutation, drug change, …)
    └── config.py          # ESKAPE registry, AWaRe levels, empiric groups, mutation graph
```

Simulation outputs in `src/log/` are **not** tracked in git (they grow into hundreds of GB).
