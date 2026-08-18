#summary for github copilot:

Here is a concise summary designed for GitHub Copilot (or project documentation) that outlines the generated products, details the qualitative scoring framework, and explains how to integrate the integrated risk tier as a prey-switching index.

---
  
  ### **Executive Summary & Pipeline Deliverables**
  
  **Key Deliverable File:**
  `master_csl_predation_risk_table_5ocean_prey_1998_2024.csv`

This exercise established a multi-tiered, empirical predation risk model for California Sea Lions (CSLs) preying on Columbia River Chinook salmon (1998–2024). It integrates **in-river estuary predator-prey dynamics** with **regional ocean forage availability** using a 5-series alternative prey ensemble.

---
  
  ### **Qualitative Scoring & Index Logic (Columns 1–4)**
  
  The first 4 columns of the master output file establish the primary risk exposure and buffer scores:
  
  | Column Header | Data Type / Format | Qualitative Scoring & Operational Definition |
  | --- | --- | --- |
  | **`year`** | Integer (1998–2024) | Primary time-series indexing key. |
  | **`csl_during_chinook`** | Continuous Numeric | Raw count of CSLs present during the spring Chinook migration passage window. |
  | **`overall_chinook_risk`** | Categorical Factor | **Baseline Estuary Risk Tier** (`LOW RISK` to `VERY HIGH RISK`). Represents predator pressure relative to in-river alternative prey buffers (Eulachon and Shad). |
  | **`ocean_buffer_score`** | Ordinal Integer (0–5) | **Regional Forage Buffer Count**. Represents the count of the 5 primary ocean forage metrics exhibiting **Moderate** or **High Biomass** ($>33\text{rd}$ percentile) in a given year. |
  
  ---
  
  ### **Quantitative Framework for the Ocean Buffer & Integrated Risk**
  
  1. **Ocean Forage Classification (Quantile Thresholds):**
  Each of the 5 standardized ocean forage metrics (Hake, Sardine, Anchovy, NCC Herring, Sitka Herring) is classified into 3 terciles:
  * **Upper Risk / Low Biomass:** $\le 33\text{rd}\text{ percentile}$ (Score contribution = $0$)
* **Moderate Biomass:** $>33\text{rd}\text{ to } \le 66\text{th}\text{ percentile}$ (Score contribution = $+1$)
* **High Biomass (Buffer):** $>66\text{th}\text{ percentile}$ (Score contribution = $+1$)


2. **Ocean Buffer Score:**
  
  $$\text{Ocean Buffer Score} = \sum_{i=1}^{5} \mathbb{I}(\text{Prey}_i > 33\text{rd Percentile})$$
  
  
  * **5 = Maximum Ocean Buffer:** Alternative ocean prey is abundant; CSLs are likely buffered offshore, decreasing per-capita predation pressure in the estuary.
* **0 = Depleted Ocean Buffer:** Alternative ocean prey is severely constrained; CSLs are forced into estuary foraging (high prey-switching motivation).


3. **Integrated Risk Tier (`integrated_csl_risk`):**
  Combines baseline estuary risk with the ocean buffer:
  * **CRITICAL RISK:** `overall_chinook_risk` is `VERY HIGH` AND `ocean_buffer_score` $\le 1$.
* **HIGH RISK:** `overall_chinook_risk` is `HIGH` or `VERY HIGH`.
* **MODERATE RISK:** Standard baseline pressure without mitigating ocean buffers.
* **LOW-MODERATE / LOW RISK:** Mitigated by high ocean buffer scores ($\ge 4$) or low raw CSL counts.



---
  
  ### **Cross-Index Harmonization (Prey-Switching Integration)**
  
  To compare or combine this index with your other prey-switching indices:
  
  * **Ordinal Linear Scale (0–4 Risk Index):**
  Map `integrated_csl_risk` to a standard numeric scale:
  
  $$\text{Risk Rank} = \{\text{LOW: } 0, \text{ LOW-MODERATE: } 1, \text{ MODERATE: } 2, \text{ HIGH: } 3, \text{ CRITICAL/VERY HIGH: } 4\}$$
  
  
  * **Inverted Ocean Pressure Index (0–5):**
  Convert the `ocean_buffer_score` into a direct **Ocean Risk Factor**:
  
  $$\text{Ocean Risk Factor} = 5 - \text{ocean\_buffer\_score}$$
  
  
  
  This aligns the scale direction so that higher values consistently denote **higher predation/switching pressure** across all indices.