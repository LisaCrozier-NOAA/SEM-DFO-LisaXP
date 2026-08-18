# Prey-switching share to salmon: p_switch_t = plogis(k0 + kW*W_t - kF*F_t)
# 
# SSL predation index: I_SSL_t = SSL_t * p_switch_t * C_t
# 
# This does what you described:
#   
#   warmer years increase switching (kW > 0)
# more forage decreases switching (kF > 0)
# no salmon abundance term in switching function
# Shark index can stay: I_Shark_t = Shark_t * Q10^((T_t - Tref)/10) * O_t
# 
# Yes. Here’s exactly what your current code is doing.

## 1) Functional form being constructed

For each year \(t\), it builds two mechanistic predator indices:
  
  ### SSL prey-switching index

In code terms:
  - `Nsalmon` = shifted/scaled salmon proxy (from `x07` + `x16` combo)
- `F_t` = shifted/scaled forage composite
- `W_t` = warm-state (`W_sst` or `W_clim`)
- `m=2`, `alpha=1`, `beta=0.8`
- `C_t=1` currently (no extra compression covariate yet)

Interpretation: when warm conditions rise and/or forage is low, the denominator term can shrink relative to salmon term, increasing `p_salmon`, so SSL predation index can increase.

---
  
  ### Shark metabolic index
with:
  - `Q10=2`
- `T_t = SST_raw`, `T_ref = mean(SST_raw)`
- `O_t = max(0.2, 1 + 0.3 W_t)` (simple overlap proxy tied to warm state)

Interpretation: warmer water increases metabolism (`M_t`), and overlap term also increases with warmth.

---
  
  ### Composite (diagnostic) index
  \[
    I_{\text{PredAK},t} = z(I_{\text{SSL},t}) + z(I_{\text{Shark},t})
    \]
Used mainly for output/diagnostic; not the main fitted predictor in your current linear fit.

---
  
  ## 2) How it is fit to salmon data
  
  For each candidate combo of:
  - one SST column
- one SSL column
- one shark column
- one forcing choice (`W_sst` vs `W_clim`)

it rebuilds indices, then fits this regression on **1998–2021 only**, lag 0:
  
  \[
    x16\_sar_t = \beta_0 + \beta_1 I_{\text{SSL},t} + \beta_2 I_{\text{Shark},t} + \varepsilon_t
    \]

In code:
  - response = `col_salmon_adult` (currently `"x16_sar"`)
- model = `lm(x16_sar ~ I_SSL + I_Shark)`
- rows filtered to years 1998–2021 and complete cases only.

Then each combo is scored by:
  - AIC
- \(R^2\)
- signs/p-values of `I_SSL`, `I_Shark`

and ranked.

---
  
  ## 3) Important subtlety about your current setup
  
  Right now, salmon is used in **two places**:
  1. inside index construction (`Nsalmon` in the SSL functional response),
2. and as fit response (`x16_sar` via `col_salmon_adult`).

That can create circularity/endogeneity risk (predictor partly built from outcome series).  
It may still be useful for exploratory screening, but for inference we should usually decouple.

---
  
  ## 4) Recommended next refinement (quick)
  
  For cleaner fitting:
  - use `x07` only (juvenile abundance proxy) inside `Nsalmon`,
- keep `x16_sar` only as response.

That preserves your smolt-year focus and reduces leakage from response into predictor.

If you want, I can provide a precise patch that makes that change and reruns ranking with identical workflow.




# A) Shark abundance candidate forms
#     For each shark column S_raw, test:
#   
#             S_z = z(S_raw)
#           logS_z = z(log1p(S_raw - min(S_raw,na.rm=T))) (if positive-ish)
#           rolling mean (2-yr) to reduce noise
#       
#       
# B) Shark temperature driver
#         Keep SSL on sst_wgoa_coastwatch_junjulaug, but for sharks test:
#         
#         surface SST candidate(s)
#       depth temperatures (e.g., 50m, 100m, 150m…)
#       maybe annual mean vs summer
#       
#       
# C) Metabolic term
#     M_t = Q10^((T_shark - Tref)/10) with Q10 ∈ {1.5, 2.0, 2.5, 3.0}
# 
# D) Overlap term alternatives
#     Current: O_t = pmax(0.2, 1 + 0.3*T_use)
# 
# Also test:
#   
#     linear clipped: pmax(0.1, 1 + a*T_shark) with a ∈ {0.1,0.2,...,0.6}
#     logistic: plogis(o0 + o1*T_shark) (more biologically bounded)
#     constant overlap (=1) as null check
# 
# 
# E) Final shark index variants
#     multiplicative: I_Shark = S * M * O
#     additive on log-scale style: I_Shark = S * (wM*M + wO*O) where wM,wO fixed simple values for sensitivity (e.g., 0.5/0.5)


#============
# Recommended final choices (primary)
# For each shark species model (run separately):
#   
#   Shark species: salmon_shark_* and pacific_sleeper_shark_* (yes, keep both)
# Shark temperature driver: pdo_djf
# Q10: 2
# Overlap form: linear
# Overlap slope: 0.4 (midpoint, since 0.2–0.8 all similar)
# Shark transform: z (simplest)
# Why:
#   
#   PDO is strong and familiar in your salmon context.
# Linear overlap is easiest to explain and avoids over-claiming logistic threshold behavior you may not be able to resolve with this sample size.
# Mid slope + z-transform are parsimonious defaults given weak sensitivity.
# SSL + Shark formulas (copy-friendly, plain text)
# I’ll give this in plain text first, then an Rmd writer block you can run.
# 
# SSL switching p_switch_ssl_t = logistic(k0 + kT * T_ssl_t - kF * F_t)
# where logistic(x) = 1 / (1 + exp(-x))
# 
# I_SSL_t = SSL_t_scaled * p_switch_ssl_t
# 
# Chosen parameters: k0 = 0
# kT = 1.5
# kF = 2.0
# T_ssl_t = scaled sst_wgoa_coastwatch_junjulaug
# F_t = scaled composite of stka_herr_matbiom + mid_il_capelin
# 
# Shark index M_t = Q10^((T_shark_t - Tref_shark) / 10)
# 
# O_t = max(0.1, 1 + slope * T_shark_t_scaled) [linear overlap]
# 
# I_Shark_t = Shark_t_scaled * M_t * O_t
# 
# Chosen parameters: Q10 = 2
# slope = 0.4
# T_shark_t = pdo_djf (primary choice)
# 
# Final fitted model x16_sar_t = b0 + b1x07_dfa_cpue_int_spr_jun_hw_t + b2*I_SSL_t + b3*I_Shark_t + error_t