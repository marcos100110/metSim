# ==============================================================================
# examples.R — runnable examples for metSim
# ==============================================================================
# Package: metSim   Author: Marcos Filho   Year: 2026
# Version: 1.0.0-beta — package in active development
# License: GPL-3 (GNU General Public License, v3 or later) — gnu.org/licenses
# Copyright (C) 2026 Marcos Filho
#
# Run the blocks in order: later ones reuse objects (`params`, `high`) from
# earlier ones. Fitting the models uses ASReml-R and can take a while, so the
# number of Monte Carlo simulations is set in one place below (`n_sims`).
# For a quick test, lower `n_sims` (e.g. 5) and/or use fewer locations.
# ==============================================================================

library(metSim)
library(dplyr)   # for the %>% / summarise used in the examples below

n_sims <- 50     # simulations per scenario (lower it for a quick test)
# Note: n_individuals (used below) is the number of genotypes (candidates).
# Note: "location" and "environment" are synonyms (one independent trial).

# ------------------------------------------------------------------------------
# Example 1 — Basic run
# ------------------------------------------------------------------------------
params <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.60, k_axes = 3, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 12, cvg_sd = 0,
  cve_mean = 12, cve_sd = 3
)

res <- run_scenario_grid(
  params,
  n_locations   = c(5, 10, 20),
  n_reps_sim    = n_sims,
  n_individuals = 50,
  models        = c("MAIN", "CS", "FA"),
  max_k         = 4,
  seed          = 42
)

summ <- summarize_results(res)

# ------------------------------------------------------------------------------
# Example 2 — Compare low and high GxE
# ------------------------------------------------------------------------------
# Low GxE: rankings stay about the same (high cor_mean, few factors)
low <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.85, k_axes = 2, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 12, cvg_sd = 0, cve_mean = 12, cve_sd = 3
)

# High GxE: rankings change a lot (low cor_mean, more factors)
high <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.30, k_axes = 4, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 12, cvg_sd = 0, cve_mean = 12, cve_sd = 3
)

res_low <- run_scenario_grid(
  low, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 3, seed = 42
)
res_high <- run_scenario_grid(
  high, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 5, seed = 42
)

# ------------------------------------------------------------------------------
# Example 3 — Change the heritability
# ------------------------------------------------------------------------------
# More genetic CV and less error CV -> higher heritability -> better accuracy
h_low <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.55, k_axes = 3, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 7, cvg_sd = 0, cve_mean = 14, cve_sd = 3
)

h_high <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.55, k_axes = 3, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 14, cvg_sd = 0, cve_mean = 7, cve_sd = 3
)

res_hlow <- run_scenario_grid(
  h_low, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42
)
res_hhigh <- run_scenario_grid(
  h_high, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42
)

# ------------------------------------------------------------------------------
# Example 4 — Sparse testing (FA vs CS in the untested cells)
# ------------------------------------------------------------------------------
# Test each genotype in only a fraction of the locations. The gap shows up in
# the UNTESTED cells: there the model must borrow across locations. CS borrows
# through a single common covariance; FA borrows through the structured
# correlation, so it wins when GxE is strong. Here we use the high-GxE `high`
# params (from Example 2) to make the FA-over-CS contrast clear.
res_sparse <- run_scenario_grid(
  high,
  n_locations   = c(10, 20, 30),
  n_reps_sim    = n_sims,
  n_individuals = 50,
  models        = c("MAIN", "CS", "FA"),
  max_k         = 5,
  sparse_fracs  = c(1.0, 0.5),
  seed          = 42
)

# Predict the untested cells — FA should beat CS here:
res_sparse %>%
  filter(Status == "Converged", Sparse_Frac == 0.5) %>%
  group_by(Model_Label) %>%
  summarise(
    Acc_Tested   = round(mean(Acc_Loc_Tested,   na.rm = TRUE), 3),
    Acc_Untested = round(mean(Acc_Loc_Untested, na.rm = TRUE), 3),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Example 5 — Start from your own data (a real correlation matrix)
# ------------------------------------------------------------------------------
# A genetic-correlation matrix between locations (e.g. from a fitted model),
# plus a genetic CV per location. It must be square and symmetric.
cor_mat <- matrix(0.55, 6, 6)
diag(cor_mat) <- 1
cor_mat[1, 2] <- cor_mat[2, 1] <- 0.75   # two locations more alike
cor_mat[3, 4] <- cor_mat[4, 3] <- 0.35   # two locations less alike

real <- build_params(
  method     = "correlation_empirical",
  cor_matrix = cor_mat,
  cvg_vector = c(10, 12, 9, 11, 10, 13),   # one genetic CV per location
  cve_mean   = 14, cve_sd = 3,
  mu_mean = 100, cv_env = 20, n_individuals = 50,
  cov_structure = NULL
)

res_real <- run_scenario_grid(
  real, n_locations = c(6, 12), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42
)

# ------------------------------------------------------------------------------
# Example 6 — Location-specific variance (Psi) and the correlation ceiling
# ------------------------------------------------------------------------------
# psi_frac is the share of genetic variance that is specific to each location
# (it does not transfer to others). Because of it, the highest correlation you
# can ask for is 1 - psi_frac. Here 0.70 <= 1 - 0.15 = 0.85, so it is fine.
# Asking for cor_mean above 1 - psi_frac triggers a warning and is capped.
with_psi <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.70, k_axes = 3, psi_frac = 0.15,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 12, cvg_sd = 0, cve_mean = 12, cve_sd = 3
)

res_psi <- run_scenario_grid(
  with_psi, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42
)

# ------------------------------------------------------------------------------
# Example 7 — View accuracy and error per model (high GxE: biggest gaps)
# ------------------------------------------------------------------------------
res_view <- run_scenario_grid(
  high, n_locations = c(20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("MAIN", "CS", "FA"), max_k = 5, seed = 42
)

# Average each metric over the converged simulations, one row per model:
res_view %>%
  filter(Status == "Converged") %>%
  group_by(Model_Label) %>%
  summarise(
    Acc_Loc = round(mean(Acc_Loc,       na.rm = TRUE), 3),  # within location
    Acc_Glo = round(mean(Acc_Glo,       na.rm = TRUE), 3),  # across locations
    MAE     = round(mean(MAE_Glo_BLUP,  na.rm = TRUE), 3),  # mean abs error
    RMSE    = round(mean(RMSE_Glo_BLUP, na.rm = TRUE), 3),  # root mean sq error
    n       = n(),
    .groups = "drop"
  )
# Under strong GxE, MAIN drops in Acc_Loc while the higher-order FA holds up.

# ------------------------------------------------------------------------------
# Example 8 — The BLUP rescaled to four scales (BLUP, REML, Cullis, Pheno)
# ------------------------------------------------------------------------------
# These are the same BLUP prediction rescaled by one factor each, so they
# share the exact same ranking (rank accuracy is identical); only the
# scale-dependent errors (MAE / RMSE) change.
res_meth <- run_scenario_grid(
  params, n_locations = c(20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42
)

# Global error of each scale, per model:
res_meth %>%
  filter(Status == "Converged") %>%
  group_by(Model_Label) %>%
  summarise(
    MAE_BLUP   = round(mean(MAE_Glo_BLUP,   na.rm = TRUE), 3),
    MAE_REML   = round(mean(MAE_Glo_REML,   na.rm = TRUE), 3),
    MAE_Cullis = round(mean(MAE_Glo_Cullis, na.rm = TRUE), 3),
    MAE_Pheno  = round(mean(MAE_Glo_Pheno,  na.rm = TRUE), 3),
    .groups = "drop"
  )
# BLUP is the shrunk scale; REML and Cullis undo part of the shrinkage; Pheno
# goes to the phenotypic (field) scale, so it usually shows the largest error.

# ------------------------------------------------------------------------------
# Example 9 — Save and read the results
# ------------------------------------------------------------------------------
res9 <- run_scenario_grid(
  params, n_locations = c(10, 20), n_reps_sim = n_sims,
  n_individuals = 50, models = c("MAIN", "CS", "FA"), max_k = 4, seed = 42
)

saveRDS(res9, "results.rds")       # save the raw output
summ9 <- summarize_results(res9)   # average over the simulations

# ------------------------------------------------------------------------------
# Example 10 — Calibrate from your own fitted ASReml FA model
# ------------------------------------------------------------------------------
# Point `model` to your own fitted model, for example:
#   model <- asreml(fixed = Yield ~ Env, random = ~ fa(Env, K):Genotype, ...)
# The genetic term can have any name (Genotype, Hybrid, Variety, ...); only the
# fa() structure matters. This block runs only if the file below exists.
model_path <- "your_fitted_fa_model.rds"
if (file.exists(model_path)) {
  model <- readRDS(model_path)

  # Extract the loadings (!fa) and specific variances (!var) from the model
  cs <- extract_params_from_asreml(model)

  # Build parameters from it. Residual error CVe is supplied here as mean + SD;
  # or pass cve_vector = c(...) (real CVe observed per location) and metSim
  # derives the mean and SD from it, so you do not give cve_sd.
  p_fa <- build_params(
    method        = "FA",
    cov_structure = cs,
    psi_vector    = cs$psi_vector,
    cve_mean = 15, cve_sd = 3,   # or: cve_vector = c(12, 14, 11, ...)
    mu_mean = 100, cv_env = 15, n_individuals = 50
  )

  res_fa <- run_scenario_grid(
    p_fa, n_locations = c(10, 20), n_reps_sim = n_sims,
    n_individuals = 50, models = c("CS", "FA"),
    max_k = cs$n_factors, seed = 42
  )
}


