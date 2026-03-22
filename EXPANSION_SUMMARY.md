# Expansion Summary: Original vs New ESKAPE ICU ABM

## 1. Pathogen Model

| | Original | New |
|--|---------|-----|
| Representation | Abstract codes: `0`, `1`, `2`, `12`, `xa`, `xn` | Real ESKAPE species + phenotype: `(kpneumoniae, KPC)`, `(saureus, MRSA)`, etc. |
| Number of states | 6 flat strings | 6 species × 2–3 phenotypes = **13 distinct pathogen states** |
| Infection probability | Two values: `sigmax=0.16` (xa/xn), `sigmac=0.45` (resistant) | **Per-species**: PA 45%, CRAB 20%, KPC 12.6%, MRSA/VRE 9.4% |
| "Susceptible" concept | xa/xn as abstract placeholders | `species='none'` (no ESKAPE colonisation) |

---

## 2. Drug Model

| | Original | New |
|--|---------|-----|
| Representation | Single letters: `A`, `B`, `C`, `L` | **WHO AWaRe integer levels 1–5** |
| Number of options | 4 | 5 levels covering 9 named agents |
| Coverage logic | ~200 lines of hardcoded `if/elif` strain-drug combinations | `is_covered(level, species, phenotype)` lookup against `COVERAGE_REQUIREMENTS` dict |
| Empiric therapy | Drug `B` at infection onset | Level 4 (Pip-Tazo / Cefepime / Meropenem) at infection onset |
| Reserve agents | No equivalent | Level 5: CZA (KPC), Colistin (CRAB), Linezolid (MRSA/VRE) |

---

## 3. Trial Arms

| | Original | New |
|--|---------|-----|
| Variants | `v1` vs `v2` — subtle differences in `convt_time` logic | **ADE arm vs Control arm** — clinically meaningful stewardship strategies |
| ADE logic | Not modelled | At 72h: de-escalate to `min_drug_level(lab_result)` |
| Control logic | v1: targeted after 72h based on lab | At 72h: stay on level 4 unless level 5 strictly required |
| Scientific purpose | Sensitivity/variant analysis | **Clinical trial design simulation** |

---

## 4. Ward Scale

| | Original | New |
|--|---------|-----|
| Patients | 16 | **64** |
| HCWs | 4 | **16** (same 4:1 ratio, larger absolute population) |
| Transmission p/q defaults | p=0.1, q=0.05 | **p=0.30, q=0.30** (updated from clinical evidence) |
| Hand hygiene η | 0.5 | 0.50 (unchanged) |

---

## 5. Initial Population

| | Original | New |
|--|---------|-----|
| Method | Hardcoded: first 4 patients fixed as strains 0,1,2,12; rest from flat 6/16 ratios | `initial_patient()` called for all patients using epidemiological priors |
| Parameters | `a`, `m`, `r1`, `r2` (opaque scaling factors) | `a=0.10` (admission colonisation rate); species/phenotype drawn from **ICU surveillance prevalence data** |
| Species distribution | Uniform across 6 abstract strains | Normalised from reported ICU ranges (K. pneu 25.3%, S. aureus 21.2%, E. faecium 20.0%, PA 16.8%, A. baumannii 12.6%, Enterobacter 4.1%) |
| Phenotype distribution | Not modelled (strain assigned directly) | Per-species: e.g. K. pneumoniae → 25% susceptible / 45% ESBL / 30% KPC |

---

## 6. Resistance Mutation

| | Original | New |
|--|---------|-----|
| Trigger condition | Hardcoded: specific `(strain, drug)` pairs in `if/elif` blocks | Drug **covers** current phenotype → selection pressure → mutation |
| Direction | `0→2`, `1→12` under specific abstract drugs | Per-species resistance ladder: `susceptible→ESBL→KPC`, `susceptible→CRPA→DTR`, `MSSA→MRSA`, etc. |
| Representation | Hardcoded in `intrinsic_mutation_v1/v2` | Data-driven `MUTATION_PATHWAYS` dict in `config.py` |

---

## 7. Code Architecture

| | Original | New |
|--|---------|-----|
| Files | 3 files: `abm.py`, `scheme_func.py`, `main.py` | **4 files**: + `config.py` as biological data layer |
| Extensibility | Adding a pathogen = rewrite `scheme_func.py` | Adding a pathogen = **edit `config.py` only** |
| Drug-pathogen logic | Inline `if/elif` chains | Matrix lookup: `COVERAGE_REQUIREMENTS[(species, phenotype)]` |
| Metric tracking | ~60 manually listed keys in `main()` | Generated dynamically from pathogen registry |

---

## 8. What Is Not Yet Changed

These aspects carry over from the original design and are candidates for future expansion:

- **Single dominant pathogen per patient** — no poly-microbial co-infection modelled yet
- **HCW-to-HCW transmission** — not modelled (HCWs only interact with patients)
- **Environmental reservoir** — sink/surface contamination not modelled
- **Immune status** — patients are epidemiologically homogeneous (no immunocompromised flag)
- **Lab result timing** — culture result assumed available exactly at 72 h (no delay or false-negative modelled)
- **Antibiotic pharmacokinetics** — treatment is binary effective/ineffective, no PK/PD curves
