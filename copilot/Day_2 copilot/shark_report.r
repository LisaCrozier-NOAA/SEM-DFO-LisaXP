report_path <- file.path(out_dir, "final_model_writeup_phase1_phase5.Rmd")
cat(
  '---
title: "SEM-DFO LisaXP: Finalized Model Specification (SSL + Shark)"
output: html_document
---

## 1. Final model structure

Response:
- `x16_sar`

Baseline covariate:
- `x07_dfa_cpue_int_spr_jun_hw`

Model:
`x16_sar ~ x07_dfa_cpue_int_spr_jun_hw + I_SSL + I_Shark`

Fit window:
- 1998-2021
- Reported separately for short (n=19) and long (n=24) model groups.

## 2. SSL submodel (fixed)

Temperature driver:
- `sst_wgoa_coastwatch_junjulaug` (scaled)

Forage for switching:
- `stka_herr_matbiom`
- `mid_il_capelin`
(composited, scaled)

Switching function:
`p_switch_ssl = 1 / (1 + exp(-(k0 + kT*T_ssl - kF*F)))`

Chosen parameters:
- `k0 = 0`
- `kT = 1.5`
- `kF = 2.0`

SSL index:
`I_SSL = SSL_scaled * p_switch_ssl`

## 3. Shark submodel (primary choice)

Species:
- Run separate models for salmon shark and Pacific sleeper shark.

Primary temperature driver:
- `pdo_djf`

Metabolic term:
`M = Q10^((T_shark_raw - Tref_shark_raw)/10)`, with `Q10=2`.

Overlap term (primary):
`O = max(0.1, 1 + slope * T_shark_scaled)`, with `slope=0.4`.

Shark index:
`I_Shark = Shark_scaled * M * O`

## 4. Why these choices

- Excluded `other_sharks` due to low biological specificity.
- Q10 sensitivity was negligible; fixed to standard `Q10=2`.
- Overlap is supported as non-constant, but linear vs logistic was similar; linear chosen for parsimony.
- PDO performed strongly and is easy to interpret in salmon context.

## 5. Sensitivity set to retain

Recommend checking:
- Shark temp driver: `enso_dj`, `swln_temp_spr_176to226m` (vs `pdo_djf`)
- Overlap form: logistic (same slope grid)
- Shark transform: `z_roll2` as a robustness check

## 6. Key output tables to archive

- `phase5_shark_sweep_ranked_no_other.csv`
- `phase5_shark_param_importance_no_other.csv`
- `phase1_fit_ranked_combos.csv`
- `phase2_forage_variable_importance_short_long.csv`

## 7. Graphs to include

- SSL switching curve and SAR vs I_SSL
- Shark overlap function families weighted by AIC weights
- SAR vs I_Shark for selected primary models
- Partial-prediction curves (varying I_SSL or I_Shark, other terms fixed)
', report_path)

message("Wrote: ", report_path)
