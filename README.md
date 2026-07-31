# SEM-DFO-LisaXP
This is my working directory for the DFO-NMFS project now moved from Doug analyzeAKindices

## Current workflow status

This repository is prepared for staged uploads of analysis files.

Planned next phases after uploads are complete:
1. Develop improved predator and prey indices for the Alaskan stage of the model, including prey-switching behavior in stellar sea lions and salmon sharks under warm/low-prey conditions.
2. Draft the manuscript using the established writing style and a later formal prompt.


## Current model status (PR #1)

We finalized the SSL prey-switching component and evaluated shark predation terms as an exploratory extension.

### Core model retained
We currently retain the core model:

\[
\text{x16\_sar}_t \sim \text{x07\_dfa\_cpue\_int\_spr\_jun\_hw}_t + I_{\text{SSL},t}
\]

where:

- \(I_{\text{SSL},t} = \text{SSL}_{z,t}\cdot p_{\text{switch},t}\)
- \(p_{\text{switch},t} = \text{logit}^{-1}(k_0 + k_T T_{\text{ssl},t} - k_F F_t)\)

Final fixed SSL settings:

- `sst_wgoa_coastwatch_junjulaug` as temperature driver
- `ssl_west_pup_pred` as SSL predictor
- forage pair: `stka_herr_matbiom + mid_il_capelin`
- \(k_0=0,\; k_T=1.5,\; k_F=2.0\)

### Shark module (exploratory; not retained in primary model)
We tested GOA and BSAI shark-species models for:

- `pacific_sleeper_shark_*`
- `salmon_shark_*`

with a fixed shark functional form using `pdo_djf`, \(Q10=2\), linear overlap, slope \(=0.4\), and z-scaled shark abundance.

In both GOA and BSAI tests, adding \(I_{\text{Shark}}\) did **not** improve model fit relative to the CPUE+SSL baseline (positive \(\Delta\)AIC and non-significant shark term). Therefore, shark effects are currently treated as hypothesis-generating sensitivity analyses rather than part of the core retained model.