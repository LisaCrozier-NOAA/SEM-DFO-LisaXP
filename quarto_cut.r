

```{r}
#| label: tbl-obj3-peff
#| tbl-cap: 'Comparison of nominal predator biomass models versus top-ranked modulated $P_{\text{eff}}$ formulations across the 24-year series (1998--2021).'
#| tbl-pos: "H"
#| results: asis

library(tidyverse)
library(knitr)

# 1. Load pre-computed 24-yr sweep results & predator crosswalk
sweep_24_path  <- "Outputs/data/obj3_data_peff_sweep_all_models_24yr.csv"
crosswalk_path <- "metadata/master_name_crosswalk.csv"

res_24 <- read_csv(sweep_24_path, show_col_types = FALSE)

lookup_dict <- character(0)
if (file.exists(crosswalk_path)) {
  cw <- read_csv(crosswalk_path, show_col_types = FALSE)
  lookup_dict <- c(setNames(cw$LisaName, cw$dfa_cols), setNames(cw$LisaName, cw$DFAname))
}

shorten_pred_label <- function(x) {
  x <- gsub("^(X08_|X09_|X10_|X11_|X15_)", "", x)
  x <- gsub("DFA_", "", x)
  x <- gsub("_predAK", "", x)
  x <- gsub("Biomass", "", x)
  ifelse(nchar(x) > 24, paste0(substr(x, 1, 22), ".."), x)
}

# 2. Format Panel A (24-yr only) with conditional green & red text highlighting
panel_a_df <- res_24 %>%
  mutate(
    raw_mod_variant = mod_variant,
    Predator_Clean = case_when(
      predator %in% names(lookup_dict) ~ shorten_pred_label(lookup_dict[predator]),
      TRUE ~ shorten_pred_label(predator)
    ),
    # Escape underscores for LaTeX safety
    Predator_Clean = gsub("_", "\\\\_", Predator_Clean),
    mod_prey_clean = gsub("_", "\\\\_", mod_prey),
    
    # Compact mechanism labels
    mod_variant = gsub("Joint Temp \\+ AltPrey", "Joint Temp+Prey", mod_variant),
    mod_variant = gsub("AltPrey Buffered", "AltPrey Only", mod_variant),
    mod_variant = gsub("Temp Overlap Only", "Temp Only", mod_variant),
    mod_variant = gsub("Nominal Biomass \\(Base\\)", "Base (Nominal)", mod_variant),
    
    # Raw formatted strings
    base_str = paste0(sprintf("%.2f", base_beta), " (", sprintf("%.2f", base_p), ")"),
    mod_str  = paste0(sprintf("%.2f", mod_beta), " (", sprintf("%.2f", mod_p), ")"),
    daic_str = sprintf("%.2f", delta_aic),
    
    # Conditional Coloring:
    # Base beta (p): Green if Nominal selected & beta < 0 & p <= 0.05; Red if beta > 0 & p < 0.20
    Base_Display = case_when(
      (raw_mod_variant == "Nominal Biomass (Base)" | mod_variant == "Base (Nominal)") &
        base_p <= 0.05 & base_beta < 0 ~ paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", base_str, "}}"),
      base_beta > 0 & base_p < 0.20   ~ paste0("\\textcolor[HTML]{B71C1C}{\\textbf{", base_str, "}}"),
      TRUE                            ~ base_str
    ),
    # Peff beta (p): Green if modulated & beta < 0 & p <= 0.05
    Mod_Display = if_else(
      (raw_mod_variant != "Nominal Biomass (Base)" & mod_variant != "Base (Nominal)") &
        mod_p <= 0.05 & mod_beta < 0,
      paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", mod_str, "}}"),
      mod_str
    ),
    # Delta AIC: Green if negative
    Delta_Display = if_else(
      delta_aic < 0,
      paste0("\\textcolor[HTML]{1B5E20}{\\textbf{", daic_str, "}}"),
      daic_str
    ),
    Base_AIC_Display = sprintf("%.2f", base_aic),
    Mod_AIC_Display  = sprintf("%.2f", mod_aic)
  ) %>%
  arrange(mod_aic) %>%
  select(
    Predator_Clean, Base_Display, Base_AIC_Display,
    mod_variant, mod_prey_clean, Mod_Display, Mod_AIC_Display, Delta_Display
  )

# 3. Build raw LaTeX tabular
tbl_latex <- panel_a_df %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = "",
    col.names = c(
      "\\textbf{Predator}",
      "\\textbf{Base $\\beta$ ($p$)}",
      "\\textbf{Base AIC}",
      "\\textbf{Mechanism}",
      "\\textbf{Alt. Prey Pool}",
      "\\textbf{$P_{\\text{eff}}$ $\\beta$ ($p$)}",
      "\\textbf{$P_{\\text{eff}}$ AIC}",
      "\\textbf{$\\Delta$AIC}"
    ),
    align = c("l", "c", "r", "l", "l", "c", "r", "r")
  )

# Add spanned group headers (Nominal vs Effective Predation)
spanner_header <- paste0(
  "\\\\toprule\n",
  " & \\\\multicolumn{2}{c}{\\\\textbf{Nominal (Base)}} & ",
  "\\\\multicolumn{5}{c}{\\\\textbf{Effective Predation ($P_{\\\\text{eff}}$)}} \\\\\\\\\n",
  "\\\\cmidrule(lr){2-3} \\\\cmidrule(lr){4-8}\n"
)
tbl_latex <- sub("\\\\toprule", spanner_header, tbl_latex)

# 4. Print with compact spacing + source note underneath
cat(
  "\\begingroup\n",
  "\\scriptsize\n",
  "\\setlength{\\tabcolsep}{3.5pt}\n",
  "\\renewcommand{\\arraystretch}{1.2}\n",
  tbl_latex, "\n",
  "\\vspace{4pt}\n",
  "\\begin{minipage}{\\linewidth}\n",
  "\\footnotesize\\textit{Note:} $P_{\\text{eff}}$ attenuates predation via temperature gating ($P \\times \\text{SST}$) and/or alternative prey buffering ($P \\times [1 - \\text{AltPrey}]$). Bold green values indicate statistically significant negative top-down impacts ($\\beta < 0$, $p \\le 0.05$) or model improvements ($\\Delta\\text{AIC} < 0$). Bold red values in Base $\\beta$ ($p$) indicate counter-intuitive positive nominal relationships ($\\beta > 0$, $p < 0.20$).\n",
  "\\end{minipage}\n",
  "\\endgroup\n",
  sep = ""
)
```




```{r}
#| label: tbl-obj0-indicator_list
#| tbl-cap: "Summary of ecosystem indicator screening and subguild consolidation from raw environmental series to Structural Equation Model (SEM) latent nodes."
#| tbl-pos: "H"
#| eval: false

library(tidyverse)
library(knitr)
library(kableExtra)

path <- "C:/Users/Lisa.Crozier/Documents/Marine survival/SEM-DFO-LisaXP/Rcode_for_paper/Routput_for_paper/Supplements/guilds.csv"

doug_guilds <- read.csv(path, stringsAsFactors = FALSE, row.names = NULL, na.strings = c("", "NA")) %>%
  filter(!is.na(Lisa_15Oct2026)) %>%
  select(ID, SEMnode, region, subregion, Lisa_21oct26, indicator, source_short, source_raw) %>%
  arrange(Lisa_21oct26)

indicator_table_data <- doug_guilds %>%
  filter(SEMnode %in% c("PreyNCC", "PredNCC", "PreyAK", "PredAK", "Growth", "Abundance", "SAR")) %>%
  group_by(SEMnode) %>%
  summarise(
    N_Screened = sum(!is.na(Lisa_21oct26) & Lisa_21oct26 != ""),
    # Use \newline directly for fixed-width columns (avoids kableExtra makecell bug)
    Subguild_Codes = paste(
      gsub("_", "\\\\_", unique(Lisa_21oct26[Lisa_21oct26 != "" & !is.na(Lisa_21oct26)])),
      collapse = " \\newline "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    SEM_Node_Full = case_when(
      SEMnode == "PreyNCC"   ~ "Northern California Current \\newline (NCC Prey)",
      SEMnode == "PredNCC"   ~ "Northern California Current \\newline (NCC Predators)",
      SEMnode == "PreyAK"    ~ "BC Coast, Gulf of Alaska \\newline (AK Prey)",
      SEMnode == "PredAK"    ~ "BC Coast, Gulf of Alaska \\newline (AK Predators)",
      SEMnode == "Growth"    ~ "Salmon Growth / Condition",
      SEMnode == "Abundance" ~ "Salmon Early Abundance (CPUE)",
      SEMnode == "SAR"       ~ "Smolt to Adult Survival",
      TRUE                   ~ SEMnode
    ),
    Subguild_Names = case_when(
      SEMnode == "PreyNCC"   ~ "Zooplankton (Newport, WA & WCVI surveys), Krill, Forage Fish, Competitors (Juv Sablefish, Market squid, Jellies)",
      SEMnode == "PredNCC"   ~ "Piscivorous Seabirds, Nearshore/Surface-feeding Groundfish, Marine Mammals",
      SEMnode == "PreyAK"    ~ "Lower & Mid-Trophic indices (Ferris et al, EGOA/WGOA), WCVI Krill, Forage Fish & Seabirds, Pink Salmon Competitors",
      SEMnode == "PredAK"    ~ "Predatory Groundfish & Sharks (Dogfish, Sleeper & Salmon Sharks), Marine Mammals",
      SEMnode == "Growth"    ~ "Juvenile Salmon IGF-1 Condition & Stomach Fullness",
      SEMnode == "Abundance" ~ "Juvenile Salmon CPUE Indices",
      SEMnode == "SAR"       ~ "Smolt to Adult Survival",
      TRUE                   ~ ""
    ),
    order = match(SEMnode, c("PreyNCC", "PredNCC", "PreyAK", "PredAK", "Growth", "Abundance", "SAR"))
  ) %>%
  arrange(order) %>%
  select(SEM_Node_Full, Subguild_Names, N_Screened, Subguild_Codes)

indicator_table_data %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = "\\addlinespace",
    col.names = c(
      "SEM Node",
      "Target Guilds",
      "Indicator Count",
      "SEM Subguild Codes (Lisa\\_21oct26)"
    ),
    align = c("l", "l", "c", "l")
  ) %>%
  kable_styling(
    font_size = 9,
    latex_options = "HOLD_position"
  ) %>%
  column_spec(1, width = "1.5in") %>%
  column_spec(2, width = "2.4in") %>%
  column_spec(3, width = "0.7in") %>%
  column_spec(4, width = "1.4in")
```




```{r}
#| label: fig-obj5-sem-compilation
#| fig-cap: "Confirmatory SEM path diagrams. (A) Four-panel compilation of regional and trophic sub-models (M1–M4); (B) Full multi-stage confirmatory synthesis model (M5) linking coastal upwelling dynamics, coastal hake predation, juvenile growth and outmigration density, and subarctic effective ocean predation to Chinook salmon marine survival (SAR)."
#| fig-align: center
#| fig-pos: "H"
#| out-width: "95%"
#| echo: false

# 1. Grab the 4-panel compilation
comp_fig <- list.files(
  "Outputs/figures/selected_sem_plots", 
  pattern = "^fig_compilation_regional_blocks_.*\\.png$", 
  full.names = TRUE
)

# 2. Grab the full confirmatory synthesis diagram (M5)
m5_fig <- list.files(
  "Outputs/figures/selected_sem_plots", 
  pattern = "^fig_M5_Confirmatory_Synthesis\\.png$", 
  full.names = TRUE
)

# Render both primary summary figures
knitr::include_graphics(c(comp_fig, m5_fig))
```

```{r}
#| label: fig-obj5-individual-sems
#| fig-cap: "Individual path diagrams for candidate models M1 through M4 showing standardized loadings, structural coefficients, and endogenous variance explained ($R^2$)."
#| fig-align: center
#| fig-pos: "H"
#| out-width: "90%"
#| echo: false

# Automatically pull individual sub-models M1 through M4
indiv_figs <- sort(list.files(
  "Outputs/figures/selected_sem_plots", 
  pattern = "^fig_M[1-4]_.*\\.png$", 
  full.names = TRUE
))

knitr::include_graphics(indiv_figs)
```



```{r}
#| label: tbl-obj5-fit-summary
#| tbl-cap: "Goodness-of-fit metrics, AIC, and explained variance ($R^2$) across candidate and confirmatory SEM architectures."
#| echo: false
#| message: false
#| warning: false

library(tidyverse)
library(knitr)
library(kableExtra)

# Handle path whether rendering from root or Rcode_for_paper/
fit_table_path <- if (file.exists("Rcode_for_paper/Outputs/data/obj4_selected_models_fit_and_r2_summary.csv")) {
  "Rcode_for_paper/Outputs/data/obj4_selected_models_fit_and_r2_summary.csv"
} else if (file.exists("Outputs/data/obj4_selected_models_fit_and_r2_summary.csv")) {
  "Outputs/data/obj4_selected_models_fit_and_r2_summary.csv"
} else {
  NULL
}

if (!is.null(fit_table_path)) {
  df_summary <- read_csv(fit_table_path, show_col_types = FALSE)
  
  df_summary %>%
    kable(
      col.names = c("Model ID", "AIC", "CFI", "SRMR", 
                    "Growth $R^2$", "CPUE $R^2$", "SAR $R^2$"),
      digits = 3,
      align = c("l", rep("r", 6)),
      escape = FALSE
    ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
} else {
  cat("Fit summary table file not found.\n")
}
'''

df_summary <- read_csv(fit_table_path, show_col_types = FALSE)

df_summary %>%
  kable(
    col.names = c("Model ID", "AIC", "CFI", "SRMR", 
                  "Growth $R^2$", "CPUE $R^2$", "SAR $R^2$"),
    digits = 3,
    align = c("l", rep("r", 6)),
    escape = FALSE
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
} else {
  cat("Fit summary table file not found.\n")
}


```{r}
#| label: tbl-obj5-fit-summary
#| tbl-cap: "Goodness-of-fit metrics, AIC, and explained variance ($R^2$) across candidate and confirmatory SEM architectures."
#| echo: false
#| message: false
#| warning: false

library(tidyverse)
library(knitr)
library(kableExtra)

fit_table_path <- "Outputs/data/obj4_selected_models_fit_and_r2_summary.csv"

if (file.exists(fit_table_path)) {
  df_summary <- read_csv(fit_table_path, show_col_types = FALSE)
  
  # Define a lookup dictionary for clean headers
  header_dict <- c(
    "Model_ID"  = "Model ID",
    "AIC"       = "AIC",
    "p_Chi2"    = "$\\chi^2$ $p$-value",
    "CFI"       = "CFI",
    "SRMR"      = "SRMR",
    "Growth_R2" = "Growth $R^2$",
    "CPUE_R2"   = "CPUE $R^2$",
    "SAR_R2"    = "SAR $R^2$"
  )
  
  # Match names dynamically to whichever columns are actually present
  clean_col_names <- header_dict[names(df_summary)]
  clean_col_names[is.na(clean_col_names)] <- names(df_summary)[is.na(clean_col_names)]
  
  df_summary %>%
    kable(
      col.names = clean_col_names,
      digits    = 3,
      align     = c("l", rep("r", ncol(df_summary) - 1)),
      escape    = FALSE
    ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
}
```

```{r}
#| label: tbl-obj5-fit-summary
#| tbl-cap: "Goodness-of-fit metrics, AIC, and explained variance ($R^2$) across candidate and confirmatory SEM architectures."
#| eval: false
#| echo: false
#| message: false
#| warning: false

library(tidyverse)
library(knitr)
library(kableExtra)

fit_table_path <- "Outputs/data/obj4_selected_models_fit_and_r2_summary.csv"

if (file.exists(fit_table_path)) {
  read_csv(fit_table_path, show_col_types = FALSE) %>%
    kable(
      col.names = c("Model ID", "AIC", "$\\chi^2$ $p$-value", "CFI", "SRMR", 
                    "Growth $R^2$", "CPUE $R^2$", "SAR $R^2$"),
      digits = 3,
      align = c("l", rep("r", 7)),
      escape = FALSE
    ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
}
```