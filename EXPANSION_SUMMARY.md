# Expansion Summary: Original vs Refactored ESKAPE ICU ABM

## 1. Pathogen Model

| | Original | Refactored |
|--|---------|------------|
| Representation | Abstract codes: `0`, `1`, `2`, `12`, `xa`, `xn` | Real ESKAPE species + phenotype: `(kpneumoniae, PTR-KP)`, `(saureus, MRSA)`, etc. |
| Number of states | 6 flat strings | 6 species × 4 phenotypes = **24 distinct pathogen states** |
| Phenotypes per species | Not modelled | susceptible, 2 resistance phenotypes, dual |
| Infection probability | Two values: `sigmax=0.16` (xa/xn), `sigmac=0.45` (resistant) | **Per-species**: Ef=0.08, Sa=0.04, Kp=0.22, Ab=0.278, Pa=0.23, Eb=0.03 |
| "Susceptible" concept | xa/xn as abstract placeholders | `species='none'` (no ESKAPE colonisation) |

---

## 2. Drug Model

| | Original | Refactored |
|--|---------|------------|
| Representation | Single letters: `A`, `B`, `C`, `L` | **Named agents** mapped to 3 AWaRe levels |
| Number of levels | 4 | **3**: Access (1), Watch (2), Reserve (3) |
| Number of agents | 4 | **8**: Nafcillin, Amp-sulbactam, Vancomycin, Pip-Tazo, Cefepime, Meropenem, Linezolid, Last-resort |
| Coverage logic | ~200 lines of hardcoded `if/elif` strain-drug combinations | `is_covered(drug_agent, species, phenotype)` lookup against `AGENT_COVERAGE` dict |
| Coverage granularity | Level-based (all level-N agents equivalent) | **Agent-specific**: Vancomycin ≠ Pip-Tazo even at same level |
| Empiric therapy | Drug `B` at infection onset | **Species-specific empiric agent** from `EMPIRIC_AGENTS` (Vancomycin for Gram+, Pip-Tazo/Meropenem for Gram−) |
| Reserve agents | No equivalent | Last-resort covers all phenotypes; used only when no Watch-level agent covers the pathogen |

---

## 3. Treatment Rule

| | Original | Refactored |
|--|---------|------------|
| Variants | `v1` vs `v2` — subtle differences in `convt_time` logic | **Single drug-change rule** — no trial arm distinction |
| Review timing | Fixed 72 h | Configurable: **24 / 48 / 72 h** via `--empiric_hours` |
| Rule | Targeted after 72h based on lab (v1) or variant logic (v2) | At review time: keep empiric agent if it covers identified phenotype; else escalate to minimum covering agent |
| Escalation | Not explicitly modelled | Prefer non-last-resort; use last-resort only if necessary |
| Empiric groups | Single fixed drug | **Two agent groups** (Group 1: Watch-level; Group 2: mixed levels) |
| Group selection | N/A | **Three regimens**: `fixed` (always group 1), `cycling` (alternate every N months), `mixing` (50/50 per patient) |
| Scientific purpose | Sensitivity/variant analysis | Models stewardship strategies: fixed vs cycling vs mixing antibiogram-guided empiric therapy |

---

## 4. Ward Scale

| | Original | Refactored |
|--|---------|------------|
| Patients | 16 | **64** |
| HCWs | 4 | **16** (same 4:1 ratio, larger absolute population) |
| Transmission p/q defaults | p=0.1, q=0.05 | **p=0.30, q=0.30** (updated from clinical evidence) |
| Hand hygiene η | 0.5 | 0.50 (unchanged) |

---

## 5. Initial Population

| | Original | Refactored |
|--|---------|------------|
| Method | Hardcoded: first 4 patients fixed as strains 0,1,2,12; rest from flat 6/16 ratios | `initial_patient()` draws species from full distribution including `none` |
| Admission colonization parameter | `a`, `m` (opaque scaling factors) | **Removed** — encoded directly in `INIT_SPECIES_PROBS` |
| Species distribution | Uniform across 6 abstract strains | Surveillance-based: Sa=0.30, Kp=0.05, Ef=0.03, Ab=0.03, Pa=0.026, Eb=0.01, none=0.554 |
| Phenotype distribution | Not modelled | Per-species from ICU surveillance (dual=0 at admission; emerges only via mutation) |

---

## 6. Resistance Mutation

| | Original | Refactored |
|--|---------|------------|
| Trigger condition | Hardcoded: specific `(strain, drug)` pairs in `if/elif` blocks | Drug **covers** current phenotype → selection pressure → mutation |
| Direction | `0→2`, `1→12` under specific abstract drugs | **Branching** per-species ladders: susceptible → [resistance_1, resistance_2] → dual |
| Branching | Not modelled | Random choice from list of next phenotypes |
| Representation | Hardcoded in `intrinsic_mutation_v1/v2` | Data-driven `MUTATION_PATHWAYS` dict (list-valued) in `config.py` |

---

## 7. Code Architecture

| | Original | Refactored |
|--|---------|------------|
| Files | 3 files: `abm.py`, `scheme_func.py`, `main.py` | **4 files**: + `config.py` as biological data layer |
| Extensibility | Adding a pathogen = rewrite `scheme_func.py` | Adding a pathogen = **edit `config.py` only** |
| Drug-pathogen logic | Inline `if/elif` chains | Agent coverage matrix: `AGENT_COVERAGE[agent]` |
| Patient drug tracking | `drug_level` (int 1–5) | `drug_agent` (str, named agent) |
| Empiric agents | Single hardcoded drug | Two groups (`EMPIRIC_AGENTS_1`, `EMPIRIC_AGENTS_2`) selected by `get_empiric_agent()` |
| Empiric regimen | Not modelled | `fixed` / `cycling` / `mixing` via `--empiric_regimen` |
| Review timing | Hardcoded 72 h | Configurable 24 / 48 / 72 h via `--empiric_hours` |
| Metric tracking | ~60 manually listed keys in `main()` | Generated dynamically from pathogen registry |
| Trial arms | Two functions: `drug_change_ADE` / `drug_change_control` | Single function: `drug_change` |
| CLI | `--arm`, `--a`, `--m` | Removed; new: `--empiric_regimen`, `--cycle_months`, `--empiric_hours` |

---

## 8. What Is Not Yet Changed

These aspects carry over from the original design and are candidates for future expansion:

- **Single dominant pathogen per patient** — no poly-microbial co-infection modelled yet
- **HCW-to-HCW transmission** — not modelled (HCWs only interact with patients)
- **Environmental reservoir** — sink/surface contamination not modelled
- **Immune status** — patients are epidemiologically homogeneous (no immunocompromised flag)
- **Lab result timing** — culture result assumed available exactly at 72 h (no delay or false-negative modelled)
- **Antibiotic pharmacokinetics** — treatment is binary effective/ineffective, no PK/PD curves
