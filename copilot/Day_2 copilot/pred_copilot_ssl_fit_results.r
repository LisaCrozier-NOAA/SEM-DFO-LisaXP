#testing whether forage prey switching idea is supported at all over a temperature-only 
#The model comparison does not support a forage-mediated SSL switch at this stage. 
#In out-of-sample rolling CV, the temperature-only model (M1) performed best (RMSE = 0.879, MAE = 0.655) and outperformed the temperature+forage model (M2; RMSE = 0.937, MAE = 0.651). 
#M2 also collapsed to an effectively zero forage effect (kF = 0.000302), indicating forage contributes little under this formulation. 
#Both switch models improved substantially over no-switch baseline (M0), but skill remains weak overall (CV R² < 0 for all models), so results should be treated as exploratory. 
#At present, the evidence points to a temperature-dominated signal rather than a robust forage-driven mechanism.

#RESULTS
# > # 6) Console summary
#   > # ---------------------------
# > cat("\n=== In-sample metrics ===\n")
# 
# === In-sample metrics ===
#   > print(ins_metrics)
# # A tibble: 3 × 5
# model                rmse   mae   cor     r2
# <chr>               <dbl> <dbl> <dbl>  <dbl>
#   1 M0_no_switch        1.08  0.899 0.388 -0.224
# 2 M1_temp_only        0.894 0.680 0.434  0.163
# 3 M2_temp_plus_forage 0.894 0.680 0.435  0.164
# > cat("\n=== CV metrics (priority) ===\n")
# 
# === CV metrics (priority) ===
#   > print(cv_metrics)
# # A tibble: 3 × 5
# model                rmse   mae     cor     r2
# <chr>               <dbl> <dbl>   <dbl>  <dbl>
#   1 M0_no_switch        1.81  1.70  -0.278  -6.08 
# 2 M1_temp_only        0.879 0.655 -0.0607 -0.667
# 3 M2_temp_plus_forage 0.937 0.651  0.0924 -0.894
# > cat("\n=== Best-fit parameters ===\n")
# 
# === Best-fit parameters ===
#   > print(par_tbl)
# # A tibble: 7 × 3
# model               parameter    estimate
# <chr>               <chr>           <dbl>
#   1 M1_temp_only        k0         -87.4     
# 2 M1_temp_only        kT         -71.5     
# 3 M2_temp_plus_forage k0        -107.      
# 4 M2_temp_plus_forage kT         -87.4     
# 5 M2_temp_plus_forage kF           0.000302
# 6 M2_temp_plus_forage a            0.647   
# 7 M2_temp_plus_forage b            0.353   





#script================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

# ---------------------------
# 0) Inputs
# ---------------------------
path_ssl_dat <- "copilot/outputs_2/ssl.dat.csv"   # update if needed
out_dir <- "copilot/outputs_2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssl.dat <- read_csv(path_ssl_dat, show_col_types = FALSE) %>% clean_names()

# expected names (cleaned)
req <- c("year","salmon_resid","ssl_model_eric","sst_wgoa_coastwatch_junjulaug","capelin_avg","herr_avg")
miss <- setdiff(req, names(ssl.dat))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse=", "))

d0 <- ssl.dat %>%
  transmute(
    year = as.integer(year),
    y_raw = salmon_resid,
    ssl_raw = ssl_model_eric,
    sst_raw = sst_wgoa_coastwatch_junjulaug,
    f1_raw = capelin_avg,
    f2_raw = herr_avg
  ) %>%
  filter(if_all(everything(), is.finite)) %>%
  arrange(year)

if (nrow(d0) < 20) warning("Small sample size; CV may be noisy.")

# scale helper
zfit <- function(train, test) {
  mu <- mean(train, na.rm = TRUE); sdv <- sd(train, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) {
    list(train = rep(0, length(train)), test = rep(0, length(test)), mu = mu, sd = sdv)
  } else {
    list(train = (train - mu)/sdv, test = (test - mu)/sdv, mu = mu, sd = sdv)
  }
}

# ---------------------------
# 1) Model definitions
# ---------------------------

# M0: yhat = ssl
fit_m0 <- function(d) list(type="m0")
pred_m0 <- function(fit, d) d$ssl

# M1: yhat = ssl * logistic(k0 + kT*sst)
fit_m1 <- function(d) {
  obj <- function(par) {
    k0 <- par[1]; kT <- par[2]
    yhat <- d$ssl * plogis(k0 + kT*d$sst)
    sum((d$y - yhat)^2)
  }
  starts <- rbind(c(0,1), c(0,-1), c(1,1), c(-1,1), c(2,-2), matrix(runif(40, -3, 3), ncol=2))
  best <- NULL; bestv <- Inf
  for (i in seq_len(nrow(starts))) {
    fi <- try(optim(starts[i,], obj, method="BFGS", control=list(maxit=2000)), silent=TRUE)
    if (!inherits(fi, "try-error") && is.finite(fi$value) && fi$value < bestv) {
      best <- fi; bestv <- fi$value
    }
  }
  if (is.null(best)) stop("M1 optimization failed")
  list(type="m1", par = c(k0=best$par[1], kT=best$par[2]), value=best$value, conv=best$convergence)
}
pred_m1 <- function(fit, d) d$ssl * plogis(fit$par["k0"] + fit$par["kT"]*d$sst)

# M2: yhat = ssl * logistic(k0 + kT*sst - kF*(a*f1 + b*f2)),
# constraints: kF>0, a>=0,b>=0,a+b=1
theta_to_par_m2 <- function(theta) {
  k0 <- theta[1]
  kT <- theta[2]
  kF <- exp(theta[3])
  w <- theta[4:5]; w <- w - max(w); ew <- exp(w); ab <- ew/sum(ew)
  c(k0=k0, kT=kT, kF=kF, a=ab[1], b=ab[2])
}
fit_m2 <- function(d) {
  obj <- function(theta) {
    p <- theta_to_par_m2(theta)
    forage <- p["a"]*d$f1 + p["b"]*d$f2
    yhat <- d$ssl * plogis(p["k0"] + p["kT"]*d$sst - p["kF"]*forage)
    sum((d$y - yhat)^2)
  }
  starts <- rbind(
    c(0,1,log(1),0,0),
    c(0,-1,log(1),0,0),
    c(-1,2,log(2),1,-1),
    c(1,-2,log(2),-1,1),
    cbind(runif(40,-3,3), runif(40,-5,5), runif(40,log(0.05),log(20)), runif(40,-3,3), runif(40,-3,3))
  )
  best <- NULL; bestv <- Inf
  for (i in seq_len(nrow(starts))) {
    fi <- try(optim(starts[i,], obj, method="BFGS", control=list(maxit=4000)), silent=TRUE)
    if (!inherits(fi, "try-error") && is.finite(fi$value) && fi$value < bestv) {
      best <- fi; bestv <- fi$value
    }
  }
  if (is.null(best)) stop("M2 optimization failed")
  list(type="m2", theta=best$par, par=theta_to_par_m2(best$par), value=best$value, conv=best$convergence)
}
pred_m2 <- function(fit, d) {
  p <- fit$par
  forage <- p["a"]*d$f1 + p["b"]*d$f2
  d$ssl * plogis(p["k0"] + p["kT"]*d$sst - p["kF"]*forage)
}

# metrics
calc_metrics <- function(obs, pred) {
  tibble(
    rmse = sqrt(mean((obs-pred)^2)),
    mae  = mean(abs(obs-pred)),
    cor  = suppressWarnings(cor(obs,pred)),
    r2   = 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
  )
}

# ---------------------------
# 2) In-sample fits (global scaling)
# ---------------------------
z <- function(x) as.numeric(scale(x))
d <- d0 %>% mutate(
  y = z(y_raw), ssl = z(ssl_raw), sst = z(sst_raw), f1 = z(f1_raw), f2 = z(f2_raw)
)

fit0 <- fit_m0(d)
fit1 <- fit_m1(d)
fit2 <- fit_m2(d)

pred_ins <- tibble(
  year = d$year,
  y = d$y,
  m0 = pred_m0(fit0, d),
  m1 = pred_m1(fit1, d),
  m2 = pred_m2(fit2, d)
)

ins_metrics <- bind_rows(
  calc_metrics(pred_ins$y, pred_ins$m0) %>% mutate(model="M0_no_switch"),
  calc_metrics(pred_ins$y, pred_ins$m1) %>% mutate(model="M1_temp_only"),
  calc_metrics(pred_ins$y, pred_ins$m2) %>% mutate(model="M2_temp_plus_forage")
) %>% select(model, everything())

# ---------------------------
# 3) Rolling-origin CV
# ---------------------------
n <- nrow(d0)
initial <- max(12, floor(0.6*n))
if (n - initial < 5) initial <- max(8, n-5)

cv_rows <- list()
idx <- 1

for (t in seq(initial, n-1)) {
  train <- d0[1:t, ]
  test  <- d0[(t+1), , drop=FALSE]
  
  zy   <- zfit(train$y_raw,   test$y_raw)
  zssl <- zfit(train$ssl_raw, test$ssl_raw)
  zsst <- zfit(train$sst_raw, test$sst_raw)
  zf1  <- zfit(train$f1_raw,  test$f1_raw)
  zf2  <- zfit(train$f2_raw,  test$f2_raw)
  
  dtr <- tibble(year=train$year, y=zy$train, ssl=zssl$train, sst=zsst$train, f1=zf1$train, f2=zf2$train)
  dte <- tibble(year=test$year,  y=zy$test,  ssl=zssl$test,  sst=zsst$test,  f1=zf1$test,  f2=zf2$test)
  
  f0 <- fit_m0(dtr)
  f1 <- fit_m1(dtr)
  f2 <- fit_m2(dtr)
  
  cv_rows[[idx]] <- tibble(
    year = dte$year,
    y = dte$y,
    m0 = pred_m0(f0, dte),
    m1 = pred_m1(f1, dte),
    m2 = pred_m2(f2, dte)
  )
  idx <- idx + 1
}

cv_pred <- bind_rows(cv_rows)

cv_metrics <- bind_rows(
  calc_metrics(cv_pred$y, cv_pred$m0) %>% mutate(model="M0_no_switch"),
  calc_metrics(cv_pred$y, cv_pred$m1) %>% mutate(model="M1_temp_only"),
  calc_metrics(cv_pred$y, cv_pred$m2) %>% mutate(model="M2_temp_plus_forage")
) %>% select(model, everything())

# ---------------------------
# 4) Save results
# ---------------------------
par_tbl <- bind_rows(
  tibble(model="M1_temp_only", parameter=names(fit1$par), estimate=as.numeric(fit1$par)),
  tibble(model="M2_temp_plus_forage", parameter=names(fit2$par), estimate=as.numeric(fit2$par))
)

write_csv(ins_metrics, file.path(out_dir, "model_compare_in_sample_metrics.csv"))
write_csv(cv_metrics,  file.path(out_dir, "model_compare_cv_metrics.csv"))
write_csv(par_tbl,     file.path(out_dir, "model_compare_parameters.csv"))
write_csv(pred_ins,    file.path(out_dir, "model_compare_in_sample_predictions.csv"))
write_csv(cv_pred,     file.path(out_dir, "model_compare_cv_predictions.csv"))

# ---------------------------
# 5) Plots
# ---------------------------
# In-sample time series
p_ins <- pred_ins %>%
  pivot_longer(cols = c(y,m0,m1,m2), names_to = "series", values_to = "value") %>%
  ggplot(aes(year, value, color=series)) +
  geom_hline(yintercept = 0, linetype=2, color="gray55") +
  geom_line(linewidth=1) + geom_point(size=1.4) +
  theme_bw() +
  labs(title="In-sample (z): observed vs model predictions", x="Year", y="z value", color=NULL)

# CV predictions only
p_cv <- cv_pred %>%
  pivot_longer(cols = c(y,m0,m1,m2), names_to = "series", values_to = "value") %>%
  ggplot(aes(year, value, color=series)) +
  geom_hline(yintercept = 0, linetype=2, color="gray55") +
  geom_line(linewidth=1) + geom_point(size=1.6) +
  theme_bw() +
  labs(title="Rolling-origin CV (1-step ahead, z)", x="Year", y="z value", color=NULL)

# CV metric bar chart
cv_long <- cv_metrics %>%
  select(model, rmse, mae, cor, r2) %>%
  pivot_longer(-model, names_to="metric", values_to="value")

p_cvmet <- ggplot(cv_long, aes(model, value, fill=model)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=20, hjust=1)) +
  labs(title="Model comparison metrics (CV)", x=NULL, y=NULL)

ggsave(file.path(out_dir, "model_compare_in_sample_timeseries.png"), p_ins, width=10, height=5.8, dpi=170)
ggsave(file.path(out_dir, "model_compare_cv_timeseries.png"), p_cv, width=10, height=5.8, dpi=170)
ggsave(file.path(out_dir, "model_compare_cv_metrics.png"), p_cvmet, width=10, height=6, dpi=170)

# ---------------------------
# 6) Console summary
# ---------------------------
cat("\n=== In-sample metrics ===\n")
print(ins_metrics)
cat("\n=== CV metrics (priority) ===\n")
print(cv_metrics)
cat("\n=== Best-fit parameters ===\n")
print(par_tbl)

message("\nWrote outputs to: ", out_dir)