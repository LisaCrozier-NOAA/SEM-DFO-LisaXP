# SEM-DFO-LisaXP
This is my working directory for the DFO-NMFS project now moved from Doug analyzeAKindices

## Current workflow status

This repository is in **Phase 1 (active)** of a staged workflow:

1. **Phase 1 (active):** Alaska predator/prey index development for SEM integration.
2. **Phase 2 (on hold):** Manuscript drafting from a later formal writing prompt.

Only Phase 1 work should be performed until Phase 2 is explicitly resumed.

## Phase 1 objective (mechanistic focus)

Build mechanistic Alaska predator/prey indices to explain additional marine mortality in Snake River spring/summer Chinook under warm-ocean conditions, centered on Gulf of Alaska / Bering Sea dynamics.

Priority is to implement literature-guided causal mechanisms first, then reconcile with observed indicator behavior.

## Phase 1 hypotheses to encode

1. **Steller sea lions (SSL):** In warm years with low preferred forage (e.g., capelin/sand lance), SSLs switch toward salmon predation due to prey limitation and potential thermal habitat compression.
2. **Salmon sharks:** Predation pressure rises in warm years through higher metabolic demand and potentially increased spatial overlap with salmon.

## Key files and canonical data flow

### Core SEM and mapping files

- `Doug_code/guilds.excludecol.csv`
  - Authoritative list of indicators currently in model.
  - Contains `SEMnode` and `guild` assignments.
  - Joins to `indicators.csv` by `shortName`.
- `Doug_code/indicators.csv`
  - Indicator list and source information.
- `outputs_4/sem_master_data.csv`
  - Prepped annual modeling table used by current SEM workflows.
  - Includes MARSS-derived trend indicators used in the SEM.
- `outputs_4/var_lookup_NCC_AK.SAR.csv`
  - Naming-convention crosswalk for harmonizing aliases/synonyms for the same indicator.

### Alaska raw ecosystem data (for mechanistic functional relationships)

- `Ferris-DFAIndicators-goa/data/EGOA_EcoState_Data_Jan2023.csv`
- `Ferris-DFAIndicators-goa/data/WGOA_EcoState_Data_Jan2023.csv`
- `Ferris-DFAIndicators-goa/data/EGOA_metadata.csv`
- `Ferris-DFAIndicators-goa/data/WGOA_metadata.csv`
- `Ferris-DFAIndicators-goa/MARSS results/climate_trends.csv` (climate DFA/MARSS trend outputs)

These files are annual time step (with seasonal indicator summaries embedded where relevant) and should be used to construct mechanistic functional indices before mapping into SEM latent structure.

### References and context

- `refs/` for stock assessments and supporting literature.
- `presentations/Project update 07242026.pptx` (indicator trends summary, especially slides 23–29).

## SEM structure notes (critical)

- Latent nodes include: `PreyNCC`, `PreyAK`, `PredNCC`, `PredAK`, salmon nodes (`Growth`, `Abundance`), and `SAR`.
- Salmon guild labels include `06.Cond1NCC` and `07.Cond2NCC`; `07.Cond2NCC` is expected to remain in the model.
- `guilds.excludecol.csv` is the source of truth for current node assignments.

## Recommended Phase 1 implementation workflow

1. Harmonize indicator names across `guilds.excludecol.csv`, `indicators.csv`, and `var_lookup_NCC_AK.SAR.csv` using `shortName` as canonical key.
2. Build annual Alaska mechanistic design table from Ferris GOA raw datasets + metadata.
3. Use SST as primary warm forcing (start with spring/summer SST). Run `climate_trends.csv` as an alternate/composite warm-state forcing sensitivity.
4. Construct mechanistic predation indices (SSL and shark) and map outputs into `PreyAK`/`PredAK` pathways.
5. Refit/evaluate SEM with lag/sensitivity tests emphasizing mechanistic sign consistency.

## Proposed formula scaffolding (to implement/refine in code)

### 1) SSL prey-switching (Type III functional response)

Use a Hill/sigmoid form to represent switching toward salmon when preferred forage declines under warm conditions:

```
p_salmon_t = Nsalmon_t^m / (Nsalmon_t^m + (alpha * F_t * exp(-beta * W_t))^m)
I_SSL_t    = SSL_t * p_salmon_t * C_t
```

Where:
- `F_t` = preferred forage availability index (capelin/sand lance-focused).
- `W_t` = warm-state index (SST-based initially; compare with climate DFA trend from `climate_trends.csv`).
- `C_t` = encounter/compression modifier (start simple; improve later with depth/habitat metrics).
- `m > 1` for Type III switching behavior.

### 2) Salmon shark metabolic predation index

```
M_t        = Q10^((T_t - T_ref)/10)
I_Shark_t  = Shark_t * M_t * O_t
```

Where:
- `T_t` = SST (spring/summer first pass).
- `O_t` = overlap/encounter proxy.
- `Shark_t` = abundance/biomass/CPUE proxy.

### 3) Combined Alaska predation-pressure index (optional composite)

```
I_PredAK_t = z(I_SSL_t) + z(I_Shark_t)
```

Use as `PredAK` candidate input or as a diagnostic composite during sensitivity testing.

## Sensitivity/evaluation priorities

- Test lags (0–2 years).
- Sensitivity over `m`, `Q10`, SST season choice, warm-state choice (SST vs climate DFA trend), and forage weighting.
- Prioritize theory-consistent directionality over strict short-term fit if discrepancies occur.

## Session handoff / restart checklist

Before ending a session, preserve continuity by updating this README with:

1. Any revised formulas (SSL/shark/composite).
2. Final variable mappings selected for `PreyAK` and `PredAK`.
3. Chosen SST season(s), warm-state forcing selection, lag structure, and sensitivity bounds.
4. Files/scripts added or modified for preprocessing and index generation.
5. Outstanding data gaps or assumptions to resolve next.
