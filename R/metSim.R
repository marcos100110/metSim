# ==============================================================================
# metSim — Simulation of MET (Multi-Environment Trials) networks and model comparison
# ==============================================================================
# Package: metSim
# Author:  Marcos Filho
# Version: 1.0.0-beta — package in active development
# Year:    2026
# License: GPL-3 (GNU General Public License, version 3 or later)
#
# Copyright (C) 2026 Marcos Filho
#
# metSim is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. See the GNU General Public License
# <https://www.gnu.org/licenses/> for more details.
#
# Models: main effect, Identity, Diagonal, Compound Symmetry and FA(1..k), via ASReml-R.
# Merit reported on four scales (rescalings of the BLUP): BLUP, REML, Cullis, phenotypic.
#
# Usage: library(metSim)
# ==============================================================================

# Package dependencies are declared in DESCRIPTION (Imports) and NAMESPACE.

# ==============================================================================
# SECTION 0: UTILITIES
# ==============================================================================

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf(fmt, ...))
  cat("\n")
}

derivar_seed <- function(master, componente, sim_id) {
  as.integer((as.numeric(master) * 1e8 + as.numeric(componente) * 1e6 + as.numeric(sim_id)) %% .Machine$integer.max)
}

.empty_result <- function(model_type = "FA", model_label = "FA1", k = NA) {
  data.frame(
    Model_Type = model_type, Model_Label = model_label, K = k,
    Status = "Init", Error_Log = NA_character_,
    Var_Tot_Mean = NA_real_, Var_FA_Mean = NA_real_, Var_Psi_Mean = NA_real_,
    GxE_Est_Cor = NA_real_, Mean_GenCor = NA_real_, Perc_Var_Expl = NA_real_,
    Cor_CV_Psi = NA_real_,
    Acc_Loc = NA_real_, Acc_Glo = NA_real_,
    Spearman_Loc = NA_real_, Spearman_Glo = NA_real_,
    Acc_Loc_Tested = NA_real_, Acc_Loc_Untested = NA_real_,
    Spearman_Loc_Tested = NA_real_, Spearman_Loc_Untested = NA_real_,
    N_Tested_Obs = NA_integer_, N_Untested_Obs = NA_integer_,
    RMSE_Glo_BLUP = NA_real_, RMSE_Glo_REML = NA_real_, RMSE_Glo_Cullis = NA_real_, RMSE_Glo_Pheno = NA_real_,
    MAE_Glo_BLUP = NA_real_, MAE_Glo_REML = NA_real_, MAE_Glo_Cullis = NA_real_, MAE_Glo_Pheno = NA_real_,
    D00_Loc_BLUP = NA_real_, D00_Loc_REML = NA_real_, D00_Loc_Cullis = NA_real_, D00_Loc_Pheno = NA_real_,
    D15_Loc_BLUP = NA_real_, D15_Loc_REML = NA_real_, D15_Loc_Cullis = NA_real_, D15_Loc_Pheno = NA_real_,
    D10_Loc_BLUP = NA_real_, D10_Loc_REML = NA_real_, D10_Loc_Cullis = NA_real_, D10_Loc_Pheno = NA_real_,
    D05_Loc_BLUP = NA_real_, D05_Loc_REML = NA_real_, D05_Loc_Cullis = NA_real_, D05_Loc_Pheno = NA_real_,
    D02_Loc_BLUP = NA_real_, D02_Loc_REML = NA_real_, D02_Loc_Cullis = NA_real_, D02_Loc_Pheno = NA_real_,
    D00_Glo_BLUP = NA_real_, D00_Glo_REML = NA_real_, D00_Glo_Cullis = NA_real_, D00_Glo_Pheno = NA_real_,
    D15_Glo_BLUP = NA_real_, D15_Glo_REML = NA_real_, D15_Glo_Cullis = NA_real_, D15_Glo_Pheno = NA_real_,
    D10_Glo_BLUP = NA_real_, D10_Glo_REML = NA_real_, D10_Glo_Cullis = NA_real_, D10_Glo_Pheno = NA_real_,
    D05_Glo_BLUP = NA_real_, D05_Glo_REML = NA_real_, D05_Glo_Cullis = NA_real_, D05_Glo_Pheno = NA_real_,
    D02_Glo_BLUP = NA_real_, D02_Glo_REML = NA_real_, D02_Glo_Cullis = NA_real_, D02_Glo_Pheno = NA_real_,
    CS_10 = NA_real_, CS_20 = NA_real_,
    DG_10 = NA_real_, DG_20 = NA_real_,
    DG_Max_10 = NA_real_, DG_Max_20 = NA_real_,
    DGR_10 = NA_real_, DGR_20 = NA_real_,
    Top1_Hit = NA_real_, Top1_in_Top2 = NA_real_, Top1_in_Top3 = NA_real_,
    Top2_Hit = NA_real_, Top3_Hit = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# SECTION 0.5: INTERNAL TPE-GENERATION HELPERS
# ==============================================================================

#' Resolve CVg/CVe per virtual environment (empirical vector via spline OR parametric)
#'
#' @param cvg_params List with $vec (empirical vector) OR $mean + $sd (parametric)
#' @param n_virtual  Number of virtual environments to generate
#' @return Vector of length n_virtual, with a biological floor of 0.5%
.resolve_cvg <- function(cvg_params, n_virtual) {
  vec <- cvg_params$vec
  mu <- cvg_params$mean
  sig <- cvg_params$sd

  if (!is.null(vec)) {
    expanded <- .expand_vector_spline(vec, n_virtual)
  } else {
    if (is.null(mu) || is.null(sig)) {
      stop("CVg_params must contain 'vec' or ('mean' + 'sd').")
    }
    expanded <- rnorm(n_virtual, mean = mu, sd = max(sig, 0))
  }

  return(pmax(expanded, 0.5))
}

#' Effective K by the 80%+5% rule with a safety cap of floor(N/5)
#'
#' An eigenvalue counts as a Factor if (cumulative variance >= 80%) AND (individual >= 5%).
#' Cap: at most floor(N/5) factors (minimum of 5 observations per factor).
#'
#' @param eigenvalues Vector of eigenvalues (any order; cleaned internally)
#' @param n_envs      Number of real locations (N)
#' @return Integer K_eff >= 1
.compute_keff_80_5 <- function(eigenvalues, n_envs) {
  ev_clean <- sort(pmax(eigenvalues, 0), decreasing = TRUE)
  total_var <- sum(ev_clean)
  if (total_var < 1e-12) {
    return(1L)
  }

  prop_var <- ev_clean / total_var
  cum_var <- cumsum(prop_var)

  passed_indiv <- which(prop_var >= 0.05)
  idx_cum80 <- which(cum_var >= 0.80)[1]
  if (is.na(idx_cum80)) idx_cum80 <- length(ev_clean)

  K_eff <- length(intersect(seq_len(idx_cum80), passed_indiv))
  K_eff <- max(1L, K_eff)

  teto <- max(1L, floor(n_envs / 5L))
  K_eff <- min(K_eff, teto)

  return(as.integer(K_eff))
}

#' Univariate Skew-t fit with sanity check and Normal fallback
#'
#' Below N=25 uses Normal directly (Skew-t unstable). At N>=25 it tries Skew-t and
#' reverts to Normal if the parameters are unstable.
#'
#' Indexes the Skew-t parameters by position (xi, omega, alpha, nu) for robustness
#' against name variations across versions of the 'sn' package.
#'
#' @param col_data Numeric vector (loadings of one factor)
#' @param n_envs   Number of real locations (controls Normal vs Skew-t)
#' @return List with $dist ("normal"|"skewt") and $xi, $omega, $alpha, $nu
.fit_skewt_safe <- function(col_data, n_envs) {
  mu_col <- mean(col_data, na.rm = TRUE)
  sd_col <- sd(col_data, na.rm = TRUE)
  if (is.na(sd_col) || sd_col < 1e-10) sd_col <- 1e-6

  if (n_envs < 25L) {
    return(list(dist = "normal", xi = mu_col, omega = sd_col, alpha = 0, nu = Inf))
  }

  fit <- tryCatch(
    as.numeric(sn::st.mple(y = col_data)$dp),
    error = function(e) NULL
  )

  if (is.null(fit) || length(fit) != 4 || any(is.na(fit)) ||
    fit[4] < 2 ||
    abs(fit[3]) > 10 ||
    fit[2] < 1e-6) {
    warning("Skew-t with unstable parameters. Reverting to Normal.")
    return(list(dist = "normal", xi = mu_col, omega = sd_col, alpha = 0, nu = Inf))
  }

  return(list(
    dist = "skewt",
    xi = fit[1],
    omega = fit[2],
    alpha = fit[3],
    nu = fit[4]
  ))
}

#' Expand a small vector to n_target via a monotone spline on Hazen quantiles
#'
#' The floor (e.g. pmax(., 0.5)) is the caller's responsibility.
#'
#' @param vec      Empirical vector (CVg, CVe...)
#' @param n_target Target length (typically 1000)
#' @return Vector of length n_target
.expand_vector_spline <- function(vec, n_target) {
  n_src <- length(vec)
  if (n_src >= n_target) {
    return(sample(vec, n_target, replace = FALSE))
  }
  if (n_src == 1L) {
    return(rep(vec, n_target))
  }
  if (n_src == n_target) {
    return(vec)
  }

  sorted_vec <- sort(vec)
  probs_src <- (seq_len(n_src) - 0.5) / n_src
  probs_tgt <- (seq_len(n_target) - 0.5) / n_target

  spline_fit <- splinefun(probs_src, sorted_vec, method = "monoH.FC")
  expanded <- spline_fit(probs_tgt)

  noise_sd <- sd(vec, na.rm = TRUE) * 0.01
  expanded <- expanded + rnorm(n_target, 0, max(noise_sd, 1e-8))
  return(expanded)
}

#' Bancic engine — decomposes C = rho*J + (1-rho)*Z'Z
#'
#' Rank k_axes+1: column 1 (sqrt(rho)*1) is the global performance axis; the
#' remaining k_axes columns are the factors that generate the GxE.
#'
#' @param cor_mean      Target mean total genetic correlation in (0,1), realized in the
#'   data (what a US model would estimate, with Psi on the diagonal).
#' @param k_axes        Number of factors that generate GxE (>= 2)
#' @param n_virtual     Number of virtual environments
#' @param Means         Vector of environment means (length n_virtual)
#' @param CVg_params    List for .resolve_cvg() ($vec or $mean+$sd)
#' @param psi_frac_mean,psi_frac_sd  Mean/SD of the specific-variance fraction (Psi)
#' @return List: Lambda_TPE (n_virtual x k_axes+1), Psi_TPE, K_gen, route
.generate_bancic_tpe <- function(cor_mean, k_axes, n_virtual,
                                 Means, CVg_params,
                                 psi_frac_mean = 0.15, psi_frac_sd = 0.03) {
  Z <- matrix(rnorm(k_axes * n_virtual), nrow = k_axes, ncol = n_virtual)
  col_norms <- sqrt(colSums(Z^2))
  col_norms[col_norms < 1e-10] <- 1
  Z_norm <- scale(Z, center = FALSE, scale = col_norms)

  CVg_New <- .resolve_cvg(CVg_params, n_virtual)

  var_g_total <- pmax((CVg_New * Means / 100)^2, 1e-10)

  psi_frac_j <- rnorm(n_virtual, mean = psi_frac_mean, sd = psi_frac_sd)
  psi_frac_j <- pmin(pmax(psi_frac_j, 1e-6), 0.30)

  var_struct <- var_g_total * (1 - psi_frac_j)
  Psi_TPE <- var_g_total * psi_frac_j

  psi_frac_bar <- mean(psi_frac_j)
  cor_struct <- cor_mean / (1 - psi_frac_bar)
  if (cor_struct > 1 - 1e-7) {
    cor_struct <- 1 - 1e-7
  }

  L_corr <- cbind(
    sqrt(cor_struct) * matrix(1, n_virtual, 1),
    sqrt(1 - cor_struct) * t(Z_norm)
  )

  Lambda_TPE <- L_corr * sqrt(var_struct)
  K_gen <- k_axes + 1L

  Psi_TPE <- pmax(Psi_TPE, 1e-8)

  return(list(
    Lambda_TPE = Lambda_TPE,
    Psi_TPE    = Psi_TPE,
    K_gen      = K_gen,
    route      = "bancic"
  ))
}

#' SVD engine — extracts K_eff from a real cor_matrix and resamples n_virtual loadings
#'
#' Decomposes the genetic covariance G_cov (N x N), extracts K_eff by the 80%+5% rule,
#' fits Normal (20<=N<25) or Skew-t (N>=25) over each factor's loadings,
#' and draws n_virtual new loadings. Psi is derived from diag(G - LL') with a floor
#' relative to the data scale.
#'
#' Psi floor = 0.001 * mean(diag(G_cov)) (relative to the data scale).
#' Variance scaling via mean(Means).
#'
#' @return List: Lambda_TPE (n_virtual x K_eff), Psi_TPE, K_eff, route
.generate_svd_tpe <- function(cor_matrix, n_envs, n_virtual, Means, CVg_params) {
  CVg_real <- .resolve_cvg(CVg_params, n_envs)

  mu_escala <- mean(Means)
  var_g_real <- pmax((CVg_real * mu_escala / 100)^2, 1e-10)

  sd_g <- sqrt(var_g_real)
  G_cov <- diag(sd_g, n_envs) %*% cor_matrix %*% diag(sd_g, n_envs)

  eig <- eigen(G_cov, symmetric = TRUE)
  ev_clean <- sort(pmax(eig$values, 0), decreasing = TRUE)

  K_eff <- .compute_keff_80_5(ev_clean, n_envs)

  if (K_eff == 1L) {
    warning(paste0(
      "[K_eff=1] All variance is concentrated in a single factor. ",
      "Highly homogeneous network; selection accuracies may be inflated."
    ))
  }

  L_struct <- eig$vectors[, 1:K_eff, drop = FALSE] %*%
    diag(sqrt(ev_clean[1:K_eff]), nrow = K_eff)

  psi_raw <- diag(G_cov) - rowSums(L_struct^2)
  psi_piso <- 0.001 * mean(diag(G_cov))
  Psi_real <- pmax(psi_raw, psi_piso)

  frac_explicada <- sum(ev_clean[1:K_eff]) / sum(ev_clean)
  if (frac_explicada < 0.70) {
    warning(sprintf(
      "[SVD] K_eff=%d explains only %.1f%% of the variance of G_cov (< 70%%) for N=%d locations; the captured structure may be insufficient.",
      K_eff, frac_explicada * 100, n_envs
    ))
  }

  fit_params <- lapply(seq_len(K_eff), function(k) {
    .fit_skewt_safe(L_struct[, k], n_envs)
  })

  Lambda_TPE <- matrix(0, nrow = n_virtual, ncol = K_eff)
  for (k in seq_len(K_eff)) {
    fp <- fit_params[[k]]
    if (fp$dist == "normal") {
      Lambda_TPE[, k] <- rnorm(n_virtual, mean = fp$xi, sd = fp$omega)
    } else {
      Lambda_TPE[, k] <- sn::rst(n_virtual,
        xi = fp$xi, omega = fp$omega,
        alpha = fp$alpha, nu = fp$nu
      )
    }
  }

  psi_fit_mean <- mean(Psi_real)
  psi_fit_sd <- max(sd(Psi_real), psi_fit_mean * 0.01)
  Psi_TPE <- rnorm(n_virtual, psi_fit_mean, psi_fit_sd)
  Psi_TPE <- pmax(Psi_TPE, psi_piso)

  return(list(
    Lambda_TPE = Lambda_TPE,
    Psi_TPE    = Psi_TPE,
    K_eff      = K_eff,
    route      = if (n_envs >= 25L) "svd_skewt" else "svd_normal"
  ))
}

# ==============================================================================
# SECTION 1: INPUT PARAMETERS
# ==============================================================================

#' Extract FA parameters from a fitted ASReml model
#'
#' Reads the factor loadings (!fa) and the specific variances (!var) from the
#' model's variance components. The residual error (CVe) is not read here; it is
#' supplied later in build_params(method = "FA", ...).
#'
#' @param model asreml object (a fitted FA model, e.g. ~ fa(Env, K):Hybrid)
#' @param K     Number of factors; if NULL, it is detected from the model
#' @return List of class metSim_cov_structure
#'   (n_factors, lambda_params, psi_vector, label)
extract_params_from_asreml <- function(model, K = NULL) {
  vc <- summary(model)$varcomp
  row_names_vc <- rownames(vc)

  fa_rows <- grep("!fa\\d", row_names_vc, value = TRUE, perl = TRUE)
  if (length(fa_rows) == 0) {
    stop("The supplied model has no factor-loading components '!fa' in the variance estimates.")
  }
  if (is.null(K)) {
    K <- max(as.integer(regmatches(fa_rows, regexpr("(?<=!fa)\\d+", fa_rows, perl = TRUE))))
  }

  L_list <- lapply(1:K, function(i) {
    vc[grep(paste0("!fa", i, "(?!\\d)"), row_names_vc, perl = TRUE), "component"]
  })
  L <- do.call(cbind, L_list)

  psi_vector <- vc[grep("!var", row_names_vc), "component"]

  svd_res <- svd(L)
  min_dim <- min(nrow(L), K)
  D <- diag(svd_res$d[1:min_dim], nrow = min_dim, ncol = min_dim)

  for (f in 1:min_dim) {
    if (sum(svd_res$u[, f] < 0) / nrow(svd_res$u) > 0.5) {
      svd_res$u[, f] <- -1 * svd_res$u[, f]
    }
  }

  L_rot_scaled <- svd_res$u[, 1:min_dim] %*% D

  n_envs_asreml <- nrow(L_rot_scaled)
  Lambda_Dist_Params <- lapply(1:K, function(k) {
    .fit_skewt_safe(L_rot_scaled[, k], n_envs_asreml)
  })

  cov_structure_extracted <- list(
    n_factors     = K,
    lambda_params = Lambda_Dist_Params,
    psi_vector    = psi_vector,
    label         = paste0("FA Extracted from Model (K=", K, ")")
  )

  class(cov_structure_extracted) <- "metSim_cov_structure"
  return(cov_structure_extracted)
}

#' Build parameters manually
#'
#' @param method  Genetic simulation method: "correlation_empirical" (default) or "FA".
#' @param cve_mean,cve_sd    Mean and SD of CVe (%) — from plot-error modeling.
#' @param cve_vector         Vector of real CVe (%) observed per location.
#' @param psi_mean,psi_sd    Mean and SD of Psi (only if method = "FA").
#' @param psi_vector         Vector of observed Psi.
#' @param mu_mean            Overall mean yield of the network.
#' @param cv_env             Coefficient of variation (%) among location means (should be estimated via LMM).
#' @param n_individuals      Size of the genotype population.
#' @param cov_structure      Covariance preset (only if method = "FA", or to inherit matrix/CVg in the empirical method).
#' @param n_loc_por_candidato Number of locations per candidate genotype.
#' @param missing_plot_pct   % of missing plots.
#' @param cvg_mean,cvg_sd    Mean and SD of CVg (%) — estimated via genotypic variance (LMM).
#' @param cvg_vector         Vector of CVg (%) per location.
#' @param cor_mean,cor_sd    Mean and SD of the genetic correlation among environments (empirical-synthetic option).
#' @param cor_matrix         Symmetric real correlation matrix (empirical option).
#' @return List of class metSim_params
build_params <- function(method = "correlation_empirical",
                         cve_mean = NULL, cve_sd = NULL, cve_vector = NULL,
                         psi_mean = NULL, psi_sd = NULL, psi_vector = NULL,
                         mu_mean = 100, cv_env = 20,
                         n_individuals,
                         cov_structure = preset_moderate_gxe(),
                         n_loc_por_candidato = NULL,
                         sparse_frac = NULL,
                         missing_plot_pct = 0,
                         cvg_mean = NULL, cvg_sd = NULL, cvg_vector = NULL,
                         cor_mean = NULL, cor_sd = NULL, cor_matrix = NULL,
                         k_axes, psi_frac = 0,
                         psi_frac_mean = 0.15, psi_frac_sd = 0.03) {
  k_axes_in <- if (missing(k_axes) || is.null(k_axes)) NULL else k_axes

  if (missing(n_individuals) || is.null(n_individuals)) {
    stop("The parameter 'n_individuals' must be provided.")
  }
  if (!is.numeric(n_individuals) || n_individuals < 10) {
    stop("n_individuals must be >= 10 to ensure consistency of the ranking and selection metrics.")
  }
  n_individuals <- as.integer(n_individuals)

  if (!missing(psi_frac) && !is.null(psi_frac)) {
    psi_frac_mean <- psi_frac
    psi_frac_sd <- 0
  }
  if (psi_frac_mean < 0 || psi_frac_sd < 0) {
    stop("psi_frac_mean and psi_frac_sd must be >= 0.")
  }
  if (psi_frac_mean + 2 * psi_frac_sd > 0.30) {
    warning("psi_frac_mean + 2*psi_frac_sd > 0.30: frequent Psi samples above the 30% limit.")
  }

  .check_vector_clean <- function(vec, nome) {
    if (!is.null(vec)) {
      if (any(is.na(vec))) stop(sprintf("NA detected in '%s'.", nome))
      if (any(!is.finite(vec))) stop(sprintf("Infinite value in '%s'.", nome))
      if (any(vec < 0)) stop(sprintf("Negative value in '%s'.", nome))
    }
  }
  .check_vector_clean(cvg_vector, "cvg_vector")
  .check_vector_clean(cve_vector, "cve_vector")
  .check_vector_clean(psi_vector, "psi_vector")

  valid_methods <- c("FA", "correlation_empirical")
  if (!(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  if (missing_plot_pct < 0 || missing_plot_pct >= 100) {
    stop("missing_plot_pct must be in the interval [0, 100) (exclusive)")
  }

  if (!is.null(sparse_frac) && !is.null(n_loc_por_candidato)) {
    stop(
      "Specify 'sparse_frac' OR 'n_loc_por_candidato', not both. ",
      "sparse_frac (fraction 0-1) is recommended for paper simulations. ",
      "n_loc_por_candidato (absolute number) is recommended for practical use."
    )
  }
  if (!is.null(sparse_frac)) {
    if (!is.numeric(sparse_frac) || length(sparse_frac) != 1 ||
      sparse_frac <= 0 || sparse_frac > 1) {
      stop("sparse_frac must be a scalar in the interval (0, 1]. Use 1.0 for full testing.")
    }
  }

  mu_sd <- (cv_env * mu_mean) / 100

  if (!is.null(cve_vector)) {
    cve_mean <- mean(cve_vector, na.rm = TRUE)
    cve_sd <- sd(cve_vector, na.rm = TRUE)
    if (!is.na(cve_sd) && cve_sd > 1e-10) {
      skew <- abs(mean((cve_vector - cve_mean)^3, na.rm = TRUE) / (cve_sd^3))
      if (!is.na(skew) && skew > 1) {
        warning("Asymmetric CVe (skewness = ", round(skew, 2), "). Normal may not be ideal.")
      }
    }
  }
  if (is.null(cve_mean) || is.null(cve_sd)) stop("Provide cve_mean+cve_sd or cve_vector")

  K <- NULL
  Lambda_Dist_Params <- NULL
  Psi_Dist_Params <- list(mean = 0, sd = 0)
  rota_tpe <- NULL
  N_emp <- NULL

  if (method == "FA") {
    if (!is.null(psi_vector)) {
      psi_mean <- mean(psi_vector, na.rm = TRUE)
      psi_sd <- sd(psi_vector, na.rm = TRUE)
    }
    if (is.null(psi_mean) || is.null(psi_sd)) stop("Provide psi_mean+psi_sd or psi_vector for the FA method")
    Psi_Dist_Params <- list(mean = psi_mean, sd = psi_sd)

    if (is.null(cov_structure) || is.null(cov_structure$n_factors) || is.null(cov_structure$lambda_params)) {
      if (!is.null(cov_structure) && !is.null(cov_structure$cor_matrix)) {
        stop("Conceptual error: the 'FA' (generative) method does not accept empirical-correlation presets (such as preset_moderate_gxe()). For structured simulations from correlation matrices, set method to 'correlation_empirical'. For 'FA' simulations, use real population parameters extracted from an ASReml model (via extract_params_from_asreml()) or specify n_factors and lambda_params manually on the physical scale.")
      } else {
        stop("For the FA method, provide a valid cov_structure containing n_factors and lambda_params. Global FA presets have been removed.")
      }
    }
    K <- cov_structure$n_factors
    Lambda_Dist_Params <- cov_structure$lambda_params
    cov_label <- cov_structure$label
    if (!is.numeric(K) || K < 1L) stop("n_factors (K) must be >= 1 for the FA method.")
    rota_tpe <- "fa_bypass"
  } else if (method == "correlation_empirical") {
    user_gave_cor <- !is.null(cor_matrix) || !is.null(cor_mean)
    if (!is.null(cov_structure)) {
      if (!is.null(cov_structure$cor_matrix) && is.null(cor_matrix)) cor_matrix <- cov_structure$cor_matrix
      if (!is.null(cov_structure$cor_mean) && is.null(cor_mean)) cor_mean <- cov_structure$cor_mean
      if (!is.null(cov_structure$k_axes) && is.null(k_axes_in)) k_axes_in <- cov_structure$k_axes
      if (!is.null(cov_structure$cvg_mean) && is.null(cvg_mean)) cvg_mean <- cov_structure$cvg_mean
      if (!is.null(cov_structure$cvg_sd) && is.null(cvg_sd)) cvg_sd <- cov_structure$cvg_sd
      if (!is.null(cov_structure$cvg_vector) && is.null(cvg_vector)) cvg_vector <- cov_structure$cvg_vector
      cov_label <- cov_structure$label
    } else {
      cov_label <- if (!is.null(cor_matrix)) paste0("Empirical Correlation (", nrow(cor_matrix), " locations)") else "Empirical Correlation (Synthetic)"
    }

    if (is.null(cor_matrix) && !is.null(cor_mean)) {
      if (is.null(k_axes_in)) {
        stop("'k_axes' is required to generate synthetic TPEs (cor_mean without cor_matrix).")
      }
      cov_label <- sprintf("Bancic Synthetic (cor_mean=%.2f, k_axes=%d)", cor_mean, as.integer(k_axes_in))
      if (cor_mean > 1 - psi_frac_mean + 1e-9) {
        warning(sprintf("cor_mean (%.2f) exceeds the ceiling 1 - psi_frac (%.2f): with psi_frac = %.2f of environment-specific (non-shared) variance, the realized genetic correlation is capped at %.2f. Lower cor_mean or reduce psi_frac to reach the target.", cor_mean, 1 - psi_frac_mean, psi_frac_mean, 1 - psi_frac_mean))
      }
    }

    if (is.null(cor_matrix) && is.null(cor_mean)) {
      stop("Provide cor_matrix (real matrix) or cor_mean (+ k_axes, synthetic) for correlation_empirical.")
    }
    if (is.null(cvg_vector) && is.null(cvg_mean)) {
      stop("Provide cvg_vector or cvg_mean (+ cvg_sd) for the correlation_empirical method.")
    }

    if (!is.null(cor_mean)) {
      if (!is.numeric(cor_mean) || cor_mean <= -1 || cor_mean >= 1) {
        stop("cor_mean must be a real number in (-1, 1).")
      }
      if (cor_mean > 0.85) {
        warning("cor_mean > 0.85: accuracies may be artificially inflated.")
      }
      if (cor_mean >= 0.99) {
        warning("cor_mean >= 0.99 (GxE ~0): FA models degenerate (rank-deficient); use CS/CORH/FA1 with res_config='homo' (free residual).")
      }
    }

    if (!is.null(k_axes_in)) {
      if (!is.numeric(k_axes_in) || k_axes_in < 2 || k_axes_in > 10) {
        stop("k_axes must be an integer between 2 and 10.")
      }
      k_axes_in <- as.integer(k_axes_in)
    }

    if (!is.null(cor_matrix)) {
      if (!is.matrix(cor_matrix) || nrow(cor_matrix) != ncol(cor_matrix)) stop("cor_matrix must be a square matrix.")
      if (any(is.na(cor_matrix))) stop("cor_matrix contains NA.")
      if (any(!is.finite(cor_matrix))) stop("cor_matrix contains Inf or NaN.")
      if (!isSymmetric(cor_matrix, tol = 1e-8)) stop("cor_matrix must be symmetric.")
      if (!is.null(cvg_vector) && length(cvg_vector) != nrow(cor_matrix)) {
        stop("The length of cvg_vector must equal the number of rows of cor_matrix.")
      }
    }

    if (!is.null(cor_matrix)) {
      N_emp <- nrow(cor_matrix)
      if (N_emp < 20L) {
        warning(sprintf("Empirical matrix with N=%d locations (< 20). Using the Empirical Bancic route.", N_emp))
        rota_tpe <- "bancic_empirico"
      } else if (N_emp < 25L) {
        rota_tpe <- "svd_normal"
      } else {
        rota_tpe <- "svd_skewt"
      }
    } else {
      rota_tpe <- "bancic_sintetico"
    }

    if (!is.null(cov_structure) && !user_gave_cor) {
      message("[metSim] No correlation supplied (cor_matrix / cor_mean); using the default preset.")
      message("         -> Active structure: ", cov_label)
      message("         Pass cor_matrix or cor_mean (+ k_axes), or set cov_structure = NULL, to override.")
    }
  }

  params <- list(
    method = method,
    Mu_Dist = list(mean = mu_mean, sd = mu_sd),
    cv_env = cv_env,
    CVE_Dist = list(mean = cve_mean, sd = cve_sd),
    cve_vector = cve_vector,
    Psi_Dist = Psi_Dist_Params,
    Lambda_Dist = Lambda_Dist_Params,
    K = K,
    n_individuals = n_individuals,
    cov_label = cov_label,
    n_loc_por_candidato = n_loc_por_candidato,
    sparse_frac = sparse_frac,
    missing_plot_pct = missing_plot_pct,
    cvg_mean = cvg_mean,
    cvg_sd = cvg_sd,
    cvg_vector = cvg_vector,
    cor_mean = cor_mean,
    cor_sd = cor_sd,
    cor_matrix = cor_matrix,
    k_axes = k_axes_in,
    psi_frac = psi_frac,
    psi_frac_mean = psi_frac_mean,
    psi_frac_sd = psi_frac_sd,
    rota_tpe = rota_tpe,
    N_emp = N_emp,
    source = "manual"
  )
  class(params) <- "metSim_params"
  return(params)
}

# ==============================================================================
# SECTION 2: COVARIANCE STRUCTURES (PRESETS)
# ==============================================================================

#' Low-GxE preset — homogeneous network, highly stable ranking across locations
preset_low_gxe <- function() {
  list(
    cor_mean   = 0.80,
    k_axes     = 2,
    cvg_mean   = 10,
    cvg_sd     = 2,
    cor_matrix = NULL,
    label      = "Low GxE (cor_mean=0.80, CVg~10%, k_axes=2)"
  )
}

#' Moderate-GxE preset — partial ranking changes across locations
preset_moderate_gxe <- function() {
  list(
    cor_mean   = 0.60,
    k_axes     = 3,
    cvg_mean   = 10,
    cvg_sd     = 2,
    cor_matrix = NULL,
    label      = "Moderate GxE (cor_mean=0.60, CVg~10%, k_axes=3)"
  )
}

#' High-GxE preset — contrasting environments, frequent ranking changes
preset_high_gxe <- function() {
  list(
    cor_mean   = 0.30,
    k_axes     = 4,
    cvg_mean   = 10,
    cvg_sd     = 2,
    cor_matrix = NULL,
    label      = "High GxE (cor_mean=0.30, CVg~10%, k_axes=4)"
  )
}

# ==============================================================================
# SECTION 3: TPE GENERATOR (Target Population of Environments)
# ==============================================================================

#' Generate the virtual population of environments (TPE)
#'
#' @param params    metSim_params
#' @param n_virtual Number of virtual environments (default 1000)
generate_tpe <- function(params, n_virtual = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  rota_tpe <- params$rota_tpe
  if (is.null(rota_tpe)) {
    stop("params$rota_tpe missing. Rebuild the parameters with build_params() v3.0.")
  }
  n_virt <- as.integer(n_virtual)

  Means_New <- rnorm(n_virt, mean = params$Mu_Dist$mean, sd = params$Mu_Dist$sd)
  floor_mu <- max(params$Mu_Dist$mean - 3 * params$Mu_Dist$sd, 0.20 * params$Mu_Dist$mean)
  n_clipped <- sum(Means_New < floor_mu)
  if (n_clipped > 0) {
    warning(sprintf(
      "%d of %d environment means (%.1f%%) truncated at the floor %.2f (3-Sigma rule / 20%% of the mean).",
      n_clipped, n_virt, 100 * n_clipped / n_virt, floor_mu
    ))
  }
  Means_New <- pmax(Means_New, floor_mu)

  if (!is.null(params$cve_vector)) {
    CVEs_New <- .expand_vector_spline(params$cve_vector, n_virt)
  } else {
    CVEs_New <- rnorm(n_virt, mean = params$CVE_Dist$mean, sd = params$CVE_Dist$sd)
  }
  CVEs_New <- pmax(CVEs_New, 0.5)

  env_names <- paste0("VirtEnv_", seq_len(n_virt))

  if (rota_tpe == "fa_bypass") {
    K <- params$K
    Psi_New <- pmax(rnorm(n_virt, params$Psi_Dist$mean, params$Psi_Dist$sd), 1e-4)

    Lambda_New <- matrix(0, nrow = n_virt, ncol = K)
    for (k in seq_len(K)) {
      p_k <- params$Lambda_Dist[[k]]
      if (!is.list(p_k)) p_k <- as.list(p_k)
      is_normal <- !is.null(p_k$dist) && p_k$dist == "normal"
      if (is_normal) {
        Lambda_New[, k] <- rnorm(n_virt, mean = p_k$xi, sd = p_k$omega)
      } else {
        Lambda_New[, k] <- sn::rst(n_virt,
          xi = p_k$xi, omega = p_k$omega, alpha = p_k$alpha, nu = p_k$nu
        )
      }
    }

    return(list(
      method = "FA", route = "fa_bypass",
      Lambda_TPE = Lambda_New, Psi_TPE = Psi_New, K_gen = K,
      Means = Means_New, CVEs = CVEs_New, Env_Names = env_names
    ))
  }

  CVg_params <- list(mean = params$cvg_mean, sd = params$cvg_sd, vec = params$cvg_vector)

  if (rota_tpe == "bancic_sintetico") {
    result <- .generate_bancic_tpe(
      cor_mean      = params$cor_mean,
      k_axes        = params$k_axes,
      n_virtual     = n_virt,
      Means         = Means_New,
      CVg_params    = CVg_params,
      psi_frac_mean = params$psi_frac_mean,
      psi_frac_sd   = params$psi_frac_sd
    )
    return(modifyList(result, list(
      method    = "correlation_empirical",
      route     = "bancic_sintetico",
      Means     = Means_New,
      CVEs      = CVEs_New,
      Env_Names = env_names
    )))
  }

  if (rota_tpe == "bancic_empirico") {
    N_emp <- params$N_emp
    R_valid <- as.matrix(Matrix::nearPD(params$cor_matrix, corr = TRUE)$mat)

    cor_mean_emp <- mean(R_valid[lower.tri(R_valid)])
    eig_emp <- eigen(R_valid, symmetric = TRUE)
    K_eff_emp <- .compute_keff_80_5(pmax(eig_emp$values, 0), N_emp)

    if (K_eff_emp == 1L) {
      warning(sprintf("[K_eff=1 Empirical Bancic] %dx%d matrix with rank-1 structure.", N_emp, N_emp))
    }

    result <- .generate_bancic_tpe(
      cor_mean      = cor_mean_emp,
      k_axes        = K_eff_emp,
      n_virtual     = n_virt,
      Means         = Means_New,
      CVg_params    = CVg_params,
      psi_frac_mean = params$psi_frac_mean,
      psi_frac_sd   = params$psi_frac_sd
    )
    return(modifyList(result, list(
      method             = "correlation_empirical",
      route              = "bancic_empirico",
      K_eff              = K_eff_emp,
      K_eff_from_data    = K_eff_emp,
      cor_mean_from_data = cor_mean_emp,
      Means              = Means_New,
      CVEs               = CVEs_New,
      Env_Names          = env_names
    )))
  }

  if (rota_tpe %in% c("svd_normal", "svd_skewt")) {
    N_emp <- params$N_emp
    R_valid <- as.matrix(Matrix::nearPD(params$cor_matrix, corr = TRUE)$mat)

    result <- .generate_svd_tpe(
      cor_matrix = R_valid,
      n_envs     = N_emp,
      n_virtual  = n_virt,
      Means      = Means_New,
      CVg_params = CVg_params
    )
    return(modifyList(result, list(
      method    = "correlation_empirical",
      Means     = Means_New,
      CVEs      = CVEs_New,
      Env_Names = env_names
    )))
  }

  stop("unknown rota_tpe in generate_tpe(): ", rota_tpe)
}

# ==============================================================================
# SECTION 4: GENETIC SIMULATION
# ==============================================================================

#' Default checks (pairs with known differences)
#' Includes a pair with diff=0 for Type-I error estimation
default_checks <- function() {
  list(
    list(name_a = "Check_00A", name_b = "Check_00B", diff_pct = 0),
    list(name_a = "Check_15A", name_b = "Check_15B", diff_pct = 15),
    list(name_a = "Check_10A", name_b = "Check_10B", diff_pct = 10),
    list(name_a = "Check_05A", name_b = "Check_05B", diff_pct = 5),
    list(name_a = "Check_02A", name_b = "Check_02B", diff_pct = 2)
  )
}

#' Simulate genetic values with an FA or Correlation structure
#'
#' @param tpe                 TPE object (from generate_tpe)
#' @param n_individuals       Number of candidate genotypes
#' @param checks              List of check pairs
#' @param check_gxe_intensity 0 = stable check, 1 = normal GxE
#' @param seed                Seed
#' @return List with G_Matrix, True_Merit, Global_Ref
simulate_genetics <- function(tpe, n_individuals,
                              checks = default_checks(),
                              check_gxe_intensity = 0.5,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(!is.null(tpe$Lambda_TPE), !is.null(tpe$Psi_TPE))

  check_names <- unlist(lapply(checks, function(ch) c(ch$name_a, ch$name_b)))
  geno_names <- c(check_names, paste0("H", seq_len(n_individuals)))
  total_genos <- length(geno_names)
  n_envs_virt <- nrow(tpe$Lambda_TPE)
  K_gen <- ncol(tpe$Lambda_TPE)

  global_mean_tpe <- mean(tpe$Means)
  check_targets <- numeric(length(check_names))
  names(check_targets) <- check_names

  Scores <- matrix(rnorm(total_genos * K_gen), total_genos, K_gen)
  G_Struct <- Scores %*% t(tpe$Lambda_TPE)

  G_Specific <- matrix(
    rnorm(total_genos * n_envs_virt,
      mean = 0,
      sd = rep(sqrt(tpe$Psi_TPE), each = total_genos)
    ),
    nrow = total_genos, ncol = n_envs_virt
  )

  for (i in seq_along(checks)) {
    ch <- checks[[i]]
    idx_a <- which(geno_names == ch$name_a)
    idx_b <- which(geno_names == ch$name_b)
    diff_tgt <- global_mean_tpe * (ch$diff_pct / 100)

    dev_a <- G_Struct[idx_a, ] - mean(G_Struct[idx_a, ])
    G_Struct[idx_a, ] <- dev_a * check_gxe_intensity + diff_tgt / 2

    dev_b <- G_Struct[idx_b, ] - mean(G_Struct[idx_b, ])
    G_Struct[idx_b, ] <- dev_b * check_gxe_intensity - diff_tgt / 2

    G_Specific[idx_a, ] <- G_Specific[idx_a, ] * check_gxe_intensity
    G_Specific[idx_b, ] <- G_Specific[idx_b, ] * check_gxe_intensity

    check_targets[ch$name_a] <- diff_tgt / 2
    check_targets[ch$name_b] <- -diff_tgt / 2
  }

  G_Total <- G_Struct + G_Specific
  True_Merit <- rowMeans(G_Struct)

  rownames(G_Total) <- geno_names
  colnames(G_Total) <- tpe$Env_Names
  names(True_Merit) <- geno_names
  for (chk in check_names) True_Merit[chk] <- check_targets[chk]

  return(list(
    G_Matrix = G_Total, True_Merit = True_Merit,
    Global_Ref = global_mean_tpe, Check_Names = check_names
  ))
}

# ==============================================================================
# SECTION 5: MET REALIZATION (RCBD with blocks)
# ==============================================================================

#' Sample locations from the TPE and generate phenotypic data (RCBD)
#'
#' @param gen_obj           Object from simulate_genetics()
#' @param tpe               TPE object
#' @param n_locations       Number of sampled locations
#' @param n_reps            Reps per location (= blocks in the RCBD)
#' @param block_sigma_ratio Block variance as a fraction of sigma_e
#' @param seed              Seed
#' @return List with Data (df), Global_Mean_Real, True_Merit
realize_met <- function(gen_obj, tpe, n_locations, n_reps = 2,
                        block_sigma_ratio = 0.3, seed = NULL,
                        n_loc_por_candidato = NULL,
                        sparse_frac = NULL,
                        missing_plot_pct = 0,
                        warn_sparse = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  n_avail <- length(tpe$Env_Names)
  if (n_locations > n_avail) {
    warning(sprintf(
      "Requested number of locations (%d) is greater than available in the TPE (%d). Reducing n_locations to %d.",
      n_locations, n_avail, n_avail
    ))
    n_locations <- n_avail
  }

  idx_met <- sample(1:n_avail, n_locations)
  env_names_selected <- tpe$Env_Names[idx_met]
  df_list <- list()

  all_genos <- rownames(gen_obj$G_Matrix)
  check_names <- gen_obj$Check_Names
  candidate_names <- setdiff(all_genos, check_names)

  n_loc_efetivo <- n_locations
  sparse_frac_efetiva <- 1.0

  if (!is.null(sparse_frac) && sparse_frac < 1.0) {
    n_loc_efetivo <- max(1L, as.integer(round(n_locations * sparse_frac)))
    sparse_frac_efetiva <- sparse_frac
  } else if (!is.null(n_loc_por_candidato) && n_loc_por_candidato > 0) {
    n_loc_efetivo <- n_loc_por_candidato
    sparse_frac_efetiva <- min(1.0, n_loc_efetivo / n_locations)
  }

  if (n_loc_efetivo >= n_locations) {
    n_loc_efetivo <- n_locations
    sparse_frac_efetiva <- 1.0
  } else if (warn_sparse) {
    message(sprintf(
      "  [Sparse Testing] Each candidate in %d/%d locations (%.0f%%).",
      n_loc_efetivo, n_locations, sparse_frac_efetiva * 100
    ))
  }

  candidate_locations <- list()
  if (n_loc_efetivo < n_locations) {
    for (cand in candidate_names) {
      candidate_locations[[cand]] <- sample(env_names_selected, n_loc_efetivo)
    }
  } else {
    for (cand in candidate_names) {
      candidate_locations[[cand]] <- env_names_selected
    }
  }

  env_to_candidates <- split(
    rep(names(candidate_locations), times = lengths(candidate_locations)),
    unlist(candidate_locations)
  )

  for (j in idx_met) {
    env_name <- tpe$Env_Names[j]
    mu_j <- tpe$Means[j]
    sigma_e <- max((mu_j * tpe$CVEs[j] / 100), 1e-6)
    sigma_b <- sigma_e * block_sigma_ratio

    active_candidates <- env_to_candidates[[env_name]]
    if (is.null(active_candidates)) active_candidates <- character(0)
    active_genos <- c(check_names, active_candidates)

    g_j <- gen_obj$G_Matrix[active_genos, j]
    n_geno <- length(g_j)

    block_effects <- rnorm(n_reps, mean = 0, sd = sigma_b)

    for (r in 1:n_reps) {
      df_list[[paste0("E", j, "_B", r)]] <- data.frame(
        Hybrid = active_genos,
        Env = env_name,
        Block = paste0("B", r),
        Yield_Sim = mu_j + block_effects[r] + g_j + rnorm(n_geno, 0, sigma_e),
        G_True_Local = g_j,
        stringsAsFactors = FALSE
      )
    }
  }

  data <- do.call(rbind, df_list)
  rownames(data) <- NULL

  if (missing_plot_pct > 0 && missing_plot_pct < 100) {
    n_rows <- nrow(data)
    n_keep <- floor(n_rows * (1 - missing_plot_pct / 100))
    if (n_keep > 0 && n_keep < n_rows) {
      keep_idx <- sample(1:n_rows, n_keep)
      data <- data[keep_idx, ]
    }
  }

  data$Block <- as.factor(data$Block)
  data$Env <- as.factor(data$Env)
  data$Hybrid <- as.factor(data$Hybrid)

  return(list(
    Data = data, Global_Mean_Real = gen_obj$Global_Ref,
    True_Merit = gen_obj$True_Merit,
    Check_Names = gen_obj$Check_Names,
    G_Matrix = gen_obj$G_Matrix,
    Env_Names_Selected = env_names_selected,
    Sparse_Frac_Efetiva = sparse_frac_efetiva,
    N_Loc_Candidato = n_loc_efetivo,
    K_eff_from_tpe = if (!is.null(tpe$K_eff)) {
      tpe$K_eff
    } else if (!is.null(tpe$K_eff_from_data)) {
      tpe$K_eff_from_data
    } else {
      NULL
    },
    K_gen_from_tpe = if (!is.null(tpe$K_gen)) tpe$K_gen else NULL
  ))
}

# ==============================================================================
# SECTION 6: ASReml ENGINE (model abstraction)
# ==============================================================================

.fit_model_asreml <- function(data, model_type, k = NULL, res_config = "homo", maxiter = 1500, workspace = "6gb") {
  random_form <- switch(model_type,
    "ID"   = ~ Env:Hybrid,
    "DIAG" = ~ diag(Env):Hybrid,
    "MAIN" = ~Hybrid,
    "CS"   = ~ Hybrid + Env:Hybrid,
    "CORH" = ~ corh(Env):Hybrid,
    "FA"   = as.formula(paste0("~ fa(Env, ", k, "):Hybrid")),
    stop("invalid model_type: ", model_type)
  )

  model <- tryCatch(
    {
      if (res_config == "fixed") {
        eval(bquote(asreml(
          fixed = Yield_BLUE ~ Env, random = .(random_form),
          weights = Weight_BLUE,
          family = asr_gaussian(dispersion = 1),
          data = .(data), trace = FALSE, ai.sing = TRUE,
          maxiter = .(maxiter), workspace = .(workspace)
        )))
      } else if (res_config == "hetero") {
        res_form <- ~ dsum(~ units | Env)
        eval(bquote(asreml(
          fixed = Yield_BLUE ~ Env, random = .(random_form),
          residual = .(res_form), weights = Weight_BLUE,
          data = .(data), trace = FALSE, ai.sing = TRUE,
          maxiter = .(maxiter), workspace = .(workspace)
        )))
      } else {
        eval(bquote(asreml(
          fixed = Yield_BLUE ~ Env, random = .(random_form),
          weights = Weight_BLUE,
          data = .(data), trace = FALSE, ai.sing = TRUE,
          maxiter = .(maxiter), workspace = .(workspace)
        )))
      }
    },
    error = function(e) {
      return(e)
    }
  )

  return(model)
}

.extract_varcomp_asreml <- function(model, model_type, k = NULL, env_levels) {
  vc <- summary(model)$varcomp
  row_names_vc <- rownames(vc)
  n_env <- length(env_levels)

  if (model_type == "ID") {
    idx <- grep("Env:Hybrid", row_names_vc, fixed = TRUE)
    sigma2_ge <- if (length(idx) > 0) vc[idx[1], "component"] else NA
    G_Diag <- rep(sigma2_ge, n_env)
    names(G_Diag) <- env_levels
    return(list(
      G_Diag = G_Diag, Psi = NULL,
      GxE_Est = NA, Mean_GenCor = NA_real_, Perc_Var_Expl = NA,
      Var_Tot = sigma2_ge, Var_FA = NA, Var_Psi = NA
    ))
  } else if (model_type == "DIAG") {
    G_Diag <- numeric(n_env)
    names(G_Diag) <- env_levels
    for (e_name in env_levels) {
      idx <- grep(paste0("Env_", e_name, ":Hybrid(?!\\d)"), row_names_vc, perl = TRUE)
      if (length(idx) == 0) idx <- grep(paste0(e_name, "!(?!\\d)"), row_names_vc, perl = TRUE)
      if (length(idx) == 0) idx <- grep(paste0(e_name, "(?!\\d)"), row_names_vc, perl = TRUE)
      G_Diag[e_name] <- if (length(idx) > 0) vc[idx[1], "component"] else NA
    }
    return(list(
      G_Diag = G_Diag, Psi = NULL,
      GxE_Est = NA, Mean_GenCor = NA_real_, Perc_Var_Expl = NA,
      Var_Tot = mean(G_Diag, na.rm = TRUE), Var_FA = NA,
      Var_Psi = NA_real_
    ))
  } else if (model_type == "CS") {

    idx_G <- which(rownames(vc) == "Hybrid")
    if (length(idx_G) == 0) {
      all_hyb <- grep("Hybrid", row_names_vc, fixed = TRUE)
      all_env <- grep("Env", row_names_vc, fixed = TRUE)
      idx_G <- setdiff(all_hyb, all_env)
    }
    sigma2_G <- if (length(idx_G) > 0) vc[idx_G[1], "component"] else NA_real_

    idx_GxE <- grep("Env:Hybrid", row_names_vc, fixed = TRUE)
    if (length(idx_GxE) == 0) idx_GxE <- grep("Hybrid:Env", row_names_vc, fixed = TRUE)
    sigma2_GxE <- if (length(idx_GxE) > 0) vc[idx_GxE[1], "component"] else NA_real_

    sigma2_tot <- sigma2_G + sigma2_GxE

    gxe_pct <- if (!is.na(sigma2_tot) && sigma2_tot > 0) {
      (1 - 1 / n_env) * sigma2_GxE / sigma2_tot * 100
    } else {
      NA_real_
    }
    perc_g <- if (!is.na(sigma2_tot) && sigma2_tot > 0) {
      sigma2_G / sigma2_tot * 100
    } else {
      NA_real_
    }
    mean_cor_cs <- if (!is.na(sigma2_tot) && sigma2_tot > 0) {
      sigma2_G / sigma2_tot
    } else {
      NA_real_
    }

    G_Diag <- rep(sigma2_tot, n_env)
    names(G_Diag) <- env_levels

    return(list(
      G_Diag        = G_Diag,
      Psi           = NULL,
      GxE_Est       = gxe_pct,
      Mean_GenCor   = mean_cor_cs,
      Perc_Var_Expl = perc_g,
      Var_Tot       = sigma2_tot,
      Var_FA        = NA_real_,
      Var_Psi       = sigma2_GxE
    ))
  } else if (model_type == "CORH") {
    idx_cor <- grep("!cor$", row_names_vc)
    rho_est <- if (length(idx_cor) > 0) vc[idx_cor[1], "component"] else NA_real_

    if (!is.na(rho_est)) {
      if (rho_est < 0) rho_est <- 0
      if (rho_est > 1) rho_est <- 1
    }

    G_Diag <- numeric(n_env)
    names(G_Diag) <- env_levels
    for (e_name in env_levels) {
      idx <- grep(paste0("Env_", e_name, ":Hybrid(?!\\d)"), row_names_vc, perl = TRUE)
      if (length(idx) == 0) idx <- grep(paste0(e_name, "!(?!\\d)"), row_names_vc, perl = TRUE)
      if (length(idx) == 0) idx <- grep(paste0(e_name, "(?!\\d)"), row_names_vc, perl = TRUE)
      G_Diag[e_name] <- if (length(idx) > 0) vc[idx[1], "component"] else NA_real_
    }

    sigma2_mean <- mean(G_Diag, na.rm = TRUE)

    return(list(
      G_Diag        = G_Diag,
      Psi           = NULL,
      GxE_Est       = if (!is.na(rho_est)) (1 - rho_est) * 100 else NA_real_,
      Mean_GenCor   = rho_est,
      Perc_Var_Expl = if (!is.na(rho_est)) rho_est * 100 else NA_real_,
      Var_Tot       = sigma2_mean,
      Var_FA        = if (!is.na(rho_est)) rho_est * sigma2_mean else NA_real_,
      Var_Psi       = if (!is.na(rho_est)) (1 - rho_est) * sigma2_mean else NA_real_
    ))
  } else if (model_type == "FA") {
    lambda_mat <- matrix(0, nrow = n_env, ncol = k)
    rownames(lambda_mat) <- env_levels
    psi_vec <- numeric(n_env)
    names(psi_vec) <- env_levels

    for (e_name in env_levels) {
      idx_psi <- grep(paste0(e_name, "!var"), row_names_vc, fixed = TRUE)
      if (length(idx_psi) > 0) psi_vec[e_name] <- vc[idx_psi[1], "component"]
      for (f in 1:k) {
        idx_load <- grep(paste0(e_name, "!fa", f), row_names_vc, fixed = TRUE)
        if (length(idx_load) > 0) lambda_mat[e_name, f] <- vc[idx_load[1], "component"]
      }
    }
    psi_vec[is.na(psi_vec)] <- 0
    L_aligned <- lambda_mat

    lambdacross <- tcrossprod(L_aligned)

    VCOV <- lambdacross + diag(psi_vec)

    p_env <- nrow(VCOV)
    gxe_est <- (1 - sum(VCOV) / (p_env * sum(diag(VCOV)))) * 100
    mean_cor <- mean(cov2cor(VCOV)[upper.tri(VCOV)])

    total_factorial_variance <- sum(diag(lambdacross))
    total_specific_variance <- sum(psi_vec, na.rm = TRUE)
    total_variance <- total_factorial_variance + total_specific_variance
    prop_total_explained <- total_factorial_variance / total_variance

    return(list(
      G_Diag = diag(VCOV), Psi = psi_vec,
      GxE_Est = gxe_est, Mean_GenCor = mean_cor,
      Perc_Var_Expl = prop_total_explained * 100,
      Var_Tot = total_variance / n_env,
      Var_FA = total_factorial_variance / n_env,
      Var_Psi = total_specific_variance / n_env
    ))
  } else if (model_type == "MAIN") {
    idx <- which(rownames(vc) == "Hybrid")
    if (length(idx) == 0) idx <- grep("Hybrid", rownames(vc), fixed = TRUE)
    sigma2_h <- if (length(idx) > 0) vc[idx[1], "component"] else NA_real_
    G_Diag <- rep(sigma2_h, n_env)
    names(G_Diag) <- env_levels
    return(list(
      G_Diag = G_Diag, Psi = NULL,
      GxE_Est = NA_real_, Mean_GenCor = NA_real_, Perc_Var_Expl = NA_real_,
      Var_Tot = sigma2_h, Var_FA = NA_real_, Var_Psi = NA_real_
    ))
  }
}

# ==============================================================================
# SECTION 7: JOINT ANALYSIS (Stage 1 + Stage 2)
# ==============================================================================

.calc_all_diffs <- function(df_preds, col_val, global_mean_real) {
  vals <- df_preds[[col_val]]
  names(vals) <- as.character(df_preds$Hybrid)
  pairs <- list(
    D00 = c("Check_00A", "Check_00B"),
    D15 = c("Check_15A", "Check_15B"),
    D10 = c("Check_10A", "Check_10B"),
    D05 = c("Check_05A", "Check_05B"),
    D02 = c("Check_02A", "Check_02B")
  )
  res <- list()
  for (p_name in names(pairs)) {
    h1 <- pairs[[p_name]][1]
    h2 <- pairs[[p_name]][2]
    if (h1 %in% names(vals) && h2 %in% names(vals)) {
      res[[p_name]] <- as.numeric(((vals[h1] - vals[h2]) / global_mean_real) * 100)
    } else {
      res[[p_name]] <- NA
    }
  }
  return(res)
}

.run_stage1 <- function(sim_data) {
  sim_data <- as.data.frame(sim_data)
  ambientes <- unique(as.character(sim_data$Env))
  stage1_results <- list()
  env_cvs <- rep(NA_real_, length(ambientes))
  names(env_cvs) <- ambientes

  for (amb in ambientes) {
    dat_env <- sim_data[sim_data$Env == amb, ]
    dat_env$Hybrid <- droplevels(as.factor(dat_env$Hybrid))
    dat_env$Block <- droplevels(as.factor(dat_env$Block))
    if (nrow(dat_env) < 5) next

    mean_yield <- mean(dat_env$Yield_Sim, na.rm = TRUE)

    model_st1 <- tryCatch(
      {
        asreml(
          fixed = Yield_Sim ~ Hybrid + Block,
          data = dat_env, trace = FALSE,
          maxiter = 150,
          na.action = na.method(x = "include")
        )
      },
      error = function(e) {
        return(e)
      }
    )
    if (inherits(model_st1, "error")) next

    env_cvs[amb] <- (sqrt(model_st1$sigma2) / mean_yield) * 100

    pred <- tryCatch(
      {
        predict(model_st1, classify = "Hybrid", trace = FALSE)
      },
      error = function(e) {
        return(e)
      }
    )
    if (inherits(pred, "error") || is.null(pred$pvals)) next

    blues_df <- pred$pvals %>%
      dplyr::filter(status == "Estimable") %>%
      dplyr::select(Hybrid, predicted.value, std.error) %>%
      dplyr::rename(Yield_BLUE = predicted.value, SE = std.error) %>%
      dplyr::mutate(Env = amb, Weight_BLUE = 1 / (SE^2))

    stage1_results[[amb]] <- blues_df
  }

  stage1_data <- do.call(rbind, stage1_results)
  if (is.null(stage1_data) || nrow(stage1_data) == 0) {
    return(NULL)
  }
  stage1_data$Weight_BLUE[is.infinite(stage1_data$Weight_BLUE)] <- 1e6

  sd_blues_env <- stage1_data %>%
    group_by(Env) %>%
    summarise(
      SD_Pheno_Env = sd(Yield_BLUE[!grepl("^Check_", as.character(Hybrid))], na.rm = TRUE),
      Mean_SE2_BLUE = mean(SE[!grepl("^Check_", as.character(Hybrid))]^2, na.rm = TRUE),
      .groups = "drop"
    )

  stage1_data$Env <- as.factor(stage1_data$Env)
  stage1_data$Hybrid <- as.factor(stage1_data$Hybrid)

  return(list(
    data = stage1_data, env_cvs = env_cvs,
    sd_blues_env = sd_blues_env
  ))
}

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  sx <- sd(x[ok])
  sy <- sd(y[ok])
  if (is.na(sx) || is.na(sy) || sx < 1e-10 || sy < 1e-10) {
    return(NA_real_)
  }
  cor(x[ok], y[ok])
}

safe_cor_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  sx <- sd(x[ok])
  sy <- sd(y[ok])
  if (is.na(sx) || is.na(sy) || sx < 1e-10 || sy < 1e-10) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = "spearman")
}

.evaluate_predictions <- function(model, stage1_info, sim_data_list,
                                  varcomp_info, res_k) {
  stage1_data <- stage1_info$data
  sd_blues_env <- stage1_info$sd_blues_env
  env_cvs <- stage1_info$env_cvs
  sim_data <- as.data.frame(sim_data_list$Data)
  global_mean_real <- sim_data_list$Global_Mean_Real
  true_merit_global <- sim_data_list$True_Merit

  check_names <- sim_data_list$Check_Names
  if (is.null(check_names)) check_names <- character(0)

  G_Matrix_full <- sim_data_list$G_Matrix
  env_names_selected <- sim_data_list$Env_Names_Selected

  G_Diag <- varcomp_info$G_Diag

  if (!is.null(varcomp_info$Psi)) {
    common_envs <- intersect(names(varcomp_info$Psi), names(env_cvs))
    if (length(common_envs) > 2) {
      res_k$Cor_CV_Psi <- cor(env_cvs[common_envs], varcomp_info$Psi[common_envs],
        use = "complete.obs"
      )
    }
  }

  pred_loc <- tryCatch(
    {
      predict(model, classify = "Hybrid:Env", trace = FALSE)
    },
    error = function(e) NULL
  )

  if (!is.null(pred_loc$pvals)) {
    df_loc <- pred_loc$pvals %>%
      dplyr::filter(status == "Estimable") %>%
      dplyr::rename(BLUP_Local = predicted.value, SE_Local = std.error)

    df_loc <- df_loc %>%
      group_by(Env) %>%
      mutate(Mean_Env = mean(BLUP_Local)) %>%
      ungroup()

    df_loc$G_Local <- df_loc$BLUP_Local - df_loc$Mean_Env
    df_loc$PEV <- df_loc$SE_Local^2

    stats_env <- df_loc %>%
      group_by(Env) %>%
      summarise(
        Mean_PEV = mean(PEV, na.rm = TRUE),
        SD_BLUP = sd(BLUP_Local[!as.character(Hybrid) %in% check_names], na.rm = TRUE), .groups = "drop"
      )

    df_loc <- left_join(df_loc, stats_env, by = "Env")
    df_loc <- left_join(df_loc, sd_blues_env, by = "Env")

    df_loc$Env_Char <- as.character(df_loc$Env)
    df_loc$Var_G_Total <- G_Diag[df_loc$Env_Char]

    df_loc$T_REML <- sqrt(pmax(df_loc$Var_G_Total, 1e-6))
    h2_base <- pmax(0.01, 1 - (df_loc$Mean_PEV / (df_loc$T_REML^2)))
    df_loc$T_Cullis <- df_loc$SD_Pheno_Env * sqrt(h2_base)
    df_loc$T_Pheno <- df_loc$SD_Pheno_Env

    apply_local_method <- function(df, col_target) {
      target <- df[[col_target]]
      f_req <- target / pmax(df$SD_BLUP, 1e-6)
      if (col_target == "T_Pheno") {
        f_fin <- f_req
      } else {
        h2_rec <- pmax(1 - (df$Mean_PEV / (target^2)), 0.01)
        f_lim <- 1 / sqrt(h2_rec)
        f_fin <- pmax(1, pmin(f_req, f_lim))
      }
      return(df$Mean_Env + (df$G_Local * f_fin))
    }

    df_loc$Val_Loc_REML <- apply_local_method(df_loc, "T_REML")
    df_loc$Val_Loc_Cullis <- apply_local_method(df_loc, "T_Cullis")
    df_loc$Val_Loc_Pheno <- apply_local_method(df_loc, "T_Pheno")

    if (!is.null(G_Matrix_full)) {
      hyb_char <- as.character(df_loc$Hybrid)
      env_char <- as.character(df_loc$Env)
      valid <- hyb_char %in% rownames(G_Matrix_full) & env_char %in% colnames(G_Matrix_full)

      df_loc$G_True_Local <- NA_real_
      if (any(valid)) {
        df_loc$G_True_Local[valid] <- G_Matrix_full[cbind(hyb_char[valid], env_char[valid])]
      }

      chaves_testadas <- paste0(stage1_data$Hybrid, "|", stage1_data$Env)
      df_loc$Is_Tested <- paste0(df_loc$Hybrid, "|", df_loc$Env) %in% chaves_testadas
    } else {
      true_vals <- sim_data %>%
        dplyr::select(Hybrid, Env, G_True_Local) %>%
        distinct()
      df_loc <- inner_join(df_loc, true_vals, by = c("Hybrid", "Env"))
      df_loc$Is_Tested <- TRUE
    }

    df_val_loc_cand <- df_loc %>%
      dplyr::filter(!as.character(Hybrid) %in% check_names, !is.na(G_True_Local))

    acc_by_env <- df_val_loc_cand %>%
      group_by(Env) %>%
      summarise(c_blup = safe_cor(BLUP_Local, G_True_Local), .groups = "drop")
    res_k$Acc_Loc <- mean(acc_by_env$c_blup, na.rm = TRUE)

    acc_tested <- df_val_loc_cand %>%
      dplyr::filter(Is_Tested) %>%
      group_by(Env) %>%
      summarise(c = safe_cor(BLUP_Local, G_True_Local), .groups = "drop")
    acc_untested <- df_val_loc_cand %>%
      dplyr::filter(!Is_Tested) %>%
      group_by(Env) %>%
      summarise(c = safe_cor(BLUP_Local, G_True_Local), .groups = "drop")

    res_k$Acc_Loc_Tested <- if (nrow(acc_tested) > 0) mean(acc_tested$c, na.rm = TRUE) else NA_real_
    res_k$Acc_Loc_Untested <- if (nrow(acc_untested) > 0) mean(acc_untested$c, na.rm = TRUE) else NA_real_

    spearman_by_env <- df_val_loc_cand %>%
      group_by(Env) %>%
      summarise(s_blup = safe_cor_spearman(BLUP_Local, G_True_Local), .groups = "drop")
    res_k$Spearman_Loc <- mean(spearman_by_env$s_blup, na.rm = TRUE)

    spear_tested <- df_val_loc_cand %>%
      dplyr::filter(Is_Tested) %>%
      group_by(Env) %>%
      summarise(s = safe_cor_spearman(BLUP_Local, G_True_Local), .groups = "drop")
    spear_untested <- df_val_loc_cand %>%
      dplyr::filter(!Is_Tested) %>%
      group_by(Env) %>%
      summarise(s = safe_cor_spearman(BLUP_Local, G_True_Local), .groups = "drop")

    res_k$Spearman_Loc_Tested <- if (nrow(spear_tested) > 0) mean(spear_tested$s, na.rm = TRUE) else NA_real_
    res_k$Spearman_Loc_Untested <- if (nrow(spear_untested) > 0) mean(spear_untested$s, na.rm = TRUE) else NA_real_

    res_k$N_Tested_Obs <- sum(df_val_loc_cand$Is_Tested)
    res_k$N_Untested_Obs <- sum(!df_val_loc_cand$Is_Tested)

    diffs_loc <- df_loc %>%
      group_by(Env) %>%
      group_map(~ {
        d_b <- .calc_all_diffs(.x, "BLUP_Local", global_mean_real)
        d_re <- .calc_all_diffs(.x, "Val_Loc_REML", global_mean_real)
        d_cu <- .calc_all_diffs(.x, "Val_Loc_Cullis", global_mean_real)
        d_ph <- .calc_all_diffs(.x, "Val_Loc_Pheno", global_mean_real)
        data.frame(
          D00_B = d_b$D00, D15_B = d_b$D15, D10_B = d_b$D10, D05_B = d_b$D05, D02_B = d_b$D02,
          D00_RE = d_re$D00, D15_RE = d_re$D15, D10_RE = d_re$D10, D05_RE = d_re$D05, D02_RE = d_re$D02,
          D00_CU = d_cu$D00, D15_CU = d_cu$D15, D10_CU = d_cu$D10, D05_CU = d_cu$D05, D02_CU = d_cu$D02,
          D00_PH = d_ph$D00, D15_PH = d_ph$D15, D10_PH = d_ph$D10, D05_PH = d_ph$D05, D02_PH = d_ph$D02
        )
      }) %>%
      bind_rows()

    res_k$D00_Loc_BLUP <- mean(diffs_loc$D00_B, na.rm = TRUE)
    res_k$D15_Loc_BLUP <- mean(diffs_loc$D15_B, na.rm = TRUE)
    res_k$D10_Loc_BLUP <- mean(diffs_loc$D10_B, na.rm = TRUE)
    res_k$D05_Loc_BLUP <- mean(diffs_loc$D05_B, na.rm = TRUE)
    res_k$D02_Loc_BLUP <- mean(diffs_loc$D02_B, na.rm = TRUE)
    res_k$D00_Loc_REML <- mean(diffs_loc$D00_RE, na.rm = TRUE)
    res_k$D15_Loc_REML <- mean(diffs_loc$D15_RE, na.rm = TRUE)
    res_k$D10_Loc_REML <- mean(diffs_loc$D10_RE, na.rm = TRUE)
    res_k$D05_Loc_REML <- mean(diffs_loc$D05_RE, na.rm = TRUE)
    res_k$D02_Loc_REML <- mean(diffs_loc$D02_RE, na.rm = TRUE)
    res_k$D00_Loc_Cullis <- mean(diffs_loc$D00_CU, na.rm = TRUE)
    res_k$D15_Loc_Cullis <- mean(diffs_loc$D15_CU, na.rm = TRUE)
    res_k$D10_Loc_Cullis <- mean(diffs_loc$D10_CU, na.rm = TRUE)
    res_k$D05_Loc_Cullis <- mean(diffs_loc$D05_CU, na.rm = TRUE)
    res_k$D02_Loc_Cullis <- mean(diffs_loc$D02_CU, na.rm = TRUE)
    res_k$D00_Loc_Pheno <- mean(diffs_loc$D00_PH, na.rm = TRUE)
    res_k$D15_Loc_Pheno <- mean(diffs_loc$D15_PH, na.rm = TRUE)
    res_k$D10_Loc_Pheno <- mean(diffs_loc$D10_PH, na.rm = TRUE)
    res_k$D05_Loc_Pheno <- mean(diffs_loc$D05_PH, na.rm = TRUE)
    res_k$D02_Loc_Pheno <- mean(diffs_loc$D02_PH, na.rm = TRUE)
  }

  pred_glo <- tryCatch(
    {
      predict(model, classify = "Hybrid", trace = FALSE)
    },
    error = function(e) NULL
  )

  if (!is.null(pred_glo$pvals)) {
    df_glo <- pred_glo$pvals %>%
      dplyr::filter(status == "Estimable") %>%
      dplyr::select(Hybrid, predicted.value, std.error) %>%
      dplyr::rename(BLUP_Global = predicted.value, SE_Global = std.error)

    Mean_Glo_Pred <- mean(df_glo$BLUP_Global, na.rm = TRUE)
    df_glo$G_Hat <- df_glo$BLUP_Global - Mean_Glo_Pred
    df_glo$PEV <- df_glo$SE_Global^2
    mean_pev_glo <- mean(df_glo$PEV, na.rm = TRUE)
    sd_blup_glo <- sd(df_glo$BLUP_Global[!as.character(df_glo$Hybrid) %in% check_names], na.rm = TRUE)

    t_reml <- sqrt(if (!is.na(varcomp_info$Var_FA)) varcomp_info$Var_FA else median(G_Diag, na.rm = TRUE))
    h2_std <- pmax(1 - (mean_pev_glo / (t_reml^2)), 0.01)
    sd_pheno_global <- mean(sd_blues_env$SD_Pheno_Env, na.rm = TRUE)
    t_cullis <- sd_pheno_global * sqrt(h2_std)
    t_pheno <- local({
      vt <- stage1_data %>%
        dplyr::filter(!grepl("^Check_", as.character(Hybrid))) %>%
        dplyr::group_by(Env) %>%
        dplyr::mutate(dev_env = Yield_BLUE - mean(Yield_BLUE, na.rm = TRUE)) %>%
        dplyr::group_by(Hybrid) %>%
        dplyr::summarise(vmean = mean(dev_env, na.rm = TRUE), .groups = "drop")
      s <- sd(vt$vmean, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) sd_blup_glo else s
    })

    apply_global <- function(target, unconstrained = FALSE) {
      f_req <- target / max(sd_blup_glo, 1e-6)
      if (unconstrained) {
        f_fin <- f_req
      } else {
        h2_rec <- pmax(1 - (mean_pev_glo / (target^2)), 0.01)
        f_lim <- 1 / sqrt(h2_rec)
        f_fin <- pmax(1, pmin(f_req, f_lim))
      }
      return(list(val = Mean_Glo_Pred + (df_glo$G_Hat * f_fin), f = f_fin))
    }

    res_reml <- apply_global(t_reml)
    res_cullis <- apply_global(t_cullis)
    res_pheno <- apply_global(t_pheno, unconstrained = TRUE)

    df_glo$Val_REML <- res_reml$val
    df_glo$Val_Cullis <- res_cullis$val
    df_glo$Val_Pheno <- res_pheno$val

    df_true <- data.frame(
      Hybrid = names(true_merit_global),
      True_Merit = as.numeric(true_merit_global),
      stringsAsFactors = FALSE
    )
    df_val_glo <- inner_join(df_glo, df_true, by = "Hybrid")

    df_val_glo_cand <- df_val_glo %>%
      dplyr::filter(!as.character(Hybrid) %in% check_names)

    res_k$Acc_Glo <- cor(df_val_glo_cand$BLUP_Global, df_val_glo_cand$True_Merit, use = "complete.obs")

    cand_true <- df_val_glo_cand$True_Merit - mean(df_val_glo_cand$True_Merit, na.rm = TRUE)

    cand_blup <- df_val_glo_cand$BLUP_Global - mean(df_val_glo_cand$BLUP_Global, na.rm = TRUE)
    cand_reml <- df_val_glo_cand$Val_REML - mean(df_val_glo_cand$Val_REML, na.rm = TRUE)
    cand_cullis <- df_val_glo_cand$Val_Cullis - mean(df_val_glo_cand$Val_Cullis, na.rm = TRUE)
    cand_pheno <- df_val_glo_cand$Val_Pheno - mean(df_val_glo_cand$Val_Pheno, na.rm = TRUE)

    res_k$RMSE_Glo_BLUP <- sqrt(mean((cand_blup - cand_true)^2, na.rm = TRUE))
    res_k$RMSE_Glo_REML <- sqrt(mean((cand_reml - cand_true)^2, na.rm = TRUE))
    res_k$RMSE_Glo_Cullis <- sqrt(mean((cand_cullis - cand_true)^2, na.rm = TRUE))
    res_k$RMSE_Glo_Pheno <- sqrt(mean((cand_pheno - cand_true)^2, na.rm = TRUE))

    res_k$MAE_Glo_BLUP <- mean(abs(cand_blup - cand_true), na.rm = TRUE)
    res_k$MAE_Glo_REML <- mean(abs(cand_reml - cand_true), na.rm = TRUE)
    res_k$MAE_Glo_Cullis <- mean(abs(cand_cullis - cand_true), na.rm = TRUE)
    res_k$MAE_Glo_Pheno <- mean(abs(cand_pheno - cand_true), na.rm = TRUE)

    res_k$Spearman_Glo <- safe_cor_spearman(df_val_glo_cand$BLUP_Global, df_val_glo_cand$True_Merit)

    col_est <- "BLUP_Global"

    cand_names_orig <- setdiff(names(true_merit_global), check_names)
    n_cand_orig <- length(cand_names_orig)
    base_mean_g <- mean(true_merit_global[cand_names_orig], na.rm = TRUE)

    n_cand <- nrow(df_val_glo_cand)
    if (n_cand >= 5) {
      p10 <- max(1, round(n_cand_orig * 0.10))
      p20 <- max(1, round(n_cand_orig * 0.20))

      true_merits_orig <- sort(true_merit_global[cand_names_orig], decreasing = TRUE)
      true_top_10 <- names(true_merits_orig)[1:p10]
      true_top_20 <- names(true_merits_orig)[1:p20]

      df_est_sorted <- df_val_glo_cand %>%
        dplyr::filter(!is.na(.data[[col_est]])) %>%
        arrange(desc(.data[[col_est]]))
      est_top_10 <- na.omit(as.character(df_est_sorted$Hybrid[1:p10]))
      est_top_20 <- na.omit(as.character(df_est_sorted$Hybrid[1:p20]))

      res_k$CS_10 <- length(intersect(est_top_10, true_top_10)) / p10 * 100
      res_k$CS_20 <- length(intersect(est_top_20, true_top_20)) / p20 * 100

      dg_est_10 <- mean(true_merit_global[est_top_10], na.rm = TRUE) - base_mean_g
      dg_est_20 <- mean(true_merit_global[est_top_20], na.rm = TRUE) - base_mean_g

      dg_max_10 <- mean(true_merits_orig[1:p10], na.rm = TRUE) - base_mean_g
      dg_max_20 <- mean(true_merits_orig[1:p20], na.rm = TRUE) - base_mean_g

      res_k$DG_10 <- dg_est_10
      res_k$DG_20 <- dg_est_20
      res_k$DG_Max_10 <- dg_max_10
      res_k$DG_Max_20 <- dg_max_20

      res_k$DGR_10 <- if (abs(dg_max_10) > 1e-6) (dg_est_10 / dg_max_10) * 100 else NA_real_
      res_k$DGR_20 <- if (abs(dg_max_20) > 1e-6) (dg_est_20 / dg_max_20) * 100 else NA_real_

      true_top_1 <- names(true_merits_orig)[1]
      true_top_2 <- names(true_merits_orig)[1:2]
      true_top_3 <- names(true_merits_orig)[1:3]

      est_top_1 <- as.character(df_est_sorted$Hybrid[1])
      est_top_2 <- as.character(df_est_sorted$Hybrid[1:2])
      est_top_3 <- as.character(df_est_sorted$Hybrid[1:3])

      res_k$Top1_Hit <- ifelse(est_top_1 == true_top_1, 1, 0)
      res_k$Top1_in_Top2 <- ifelse(est_top_1 %in% true_top_2, 1, 0)
      res_k$Top1_in_Top3 <- ifelse(est_top_1 %in% true_top_3, 1, 0)
      res_k$Top2_Hit <- length(intersect(est_top_2, true_top_2)) / 2 * 100
      res_k$Top3_Hit <- length(intersect(est_top_3, true_top_3)) / 3 * 100
    }

    d_b <- .calc_all_diffs(df_val_glo, "BLUP_Global", global_mean_real)
    d_re <- .calc_all_diffs(df_val_glo, "Val_REML", global_mean_real)
    d_cu <- .calc_all_diffs(df_val_glo, "Val_Cullis", global_mean_real)
    d_ph <- .calc_all_diffs(df_val_glo, "Val_Pheno", global_mean_real)

    res_k$D00_Glo_BLUP <- d_b$D00
    res_k$D00_Glo_REML <- d_re$D00
    res_k$D00_Glo_Cullis <- d_cu$D00
    res_k$D00_Glo_Pheno <- d_ph$D00

    res_k$D15_Glo_BLUP <- d_b$D15
    res_k$D15_Glo_REML <- d_re$D15
    res_k$D15_Glo_Cullis <- d_cu$D15
    res_k$D15_Glo_Pheno <- d_ph$D15
    res_k$D10_Glo_BLUP <- d_b$D10
    res_k$D10_Glo_REML <- d_re$D10
    res_k$D10_Glo_Cullis <- d_cu$D10
    res_k$D10_Glo_Pheno <- d_ph$D10
    res_k$D05_Glo_BLUP <- d_b$D05
    res_k$D05_Glo_REML <- d_re$D05
    res_k$D05_Glo_Cullis <- d_cu$D05
    res_k$D05_Glo_Pheno <- d_ph$D05
    res_k$D02_Glo_BLUP <- d_b$D02
    res_k$D02_Glo_REML <- d_re$D02
    res_k$D02_Glo_Cullis <- d_cu$D02
    res_k$D02_Glo_Pheno <- d_ph$D02
  }

  return(res_k)
}

#' Joint analysis of a simulated MET
#'
#' @param sim_data_list  Result from realize_met()
#' @param models         Models to fit: a subset of c("MAIN", "ID", "DIAG", "CS", "FA").
#'   MAIN (~Hybrid) = main effect, cor=1 across environments.
#'   CS (~Hybrid + Env:Hybrid) = Compound Symmetry: estimates sigma2_G and sigma2_GxE;
#'   GxE_Est_Cor = sigma2_GxE/(sigma2_G+sigma2_GxE)*100.
#' @param max_k          Maximum K for FA
#' @param het_residual   TRUE = also fits versions with heterogeneous residual
#' @return data.frame with one row per model
analyze_met <- function(sim_data_list,
                        models = c("ID", "DIAG", "FA"),
                        max_k = 6,
                        res_configs = NULL,
                        het_residual = FALSE,
                        maxiter = 1500,
                        workspace = "6gb") {
  build_specs <- function(current_limit_k) {
    specs <- list()
    for (m in models) {
      m_configs <- if (is.null(res_configs)) {
        base_cfg <- if (m %in% c("DIAG", "CS", "CORH")) "fixed" else "homo"
        if (het_residual) c(base_cfg, "hetero") else base_cfg
      } else if (is.list(res_configs)) {
        if (!is.null(res_configs[[m]])) res_configs[[m]] else stop(sprintf("Config for model %s not provided in the list.", m))
      } else {
        res_configs
      }
      for (rc in m_configs) {
        get_label <- function(base) {
          if (is.null(res_configs)) {
            if (rc == "hetero" && het_residual) paste0(base, "_het") else base
          } else {
            paste0(base, "_", rc)
          }
        }
        if (m == "FA") {
          for (k_val in 1:current_limit_k) {
            specs[[length(specs) + 1]] <- list(type = "FA", k = k_val, rc = rc, label = get_label(paste0("FA", k_val)))
          }
        } else {
          specs[[length(specs) + 1]] <- list(type = m, k = NA, rc = rc, label = get_label(m))
        }
      }
    }
    return(specs)
  }

  stage1_info <- .run_stage1(sim_data_list$Data)
  if (is.null(stage1_info)) {
    n_env_raw <- length(unique(sim_data_list$Data$Env))
    K_eff_tpe <- sim_data_list$K_eff_from_tpe
    K_gen_tpe <- sim_data_list$K_gen_from_tpe
    if (!is.null(K_eff_tpe) && K_eff_tpe > 0L) {
      limit_k <- min(max_k, K_eff_tpe + 1L, n_env_raw - 1)
    } else if (!is.null(K_gen_tpe) && K_gen_tpe > 0L) {
      limit_k <- min(max_k, K_gen_tpe, n_env_raw - 1)
    } else {
      limit_k <- min(max_k, n_env_raw - 1)
    }
    if (limit_k < 1) limit_k <- 1

    model_specs_fail <- build_specs(limit_k)

    fail_rows <- lapply(model_specs_fail, function(spec) {
      .empty_result(spec$type, spec$label, spec$k) %>%
        mutate(Status = "Failed_Stage1")
    })
    return(dplyr::bind_rows(fail_rows))
  }

  stage1_data <- stage1_info$data
  sd_blues_env <- stage1_info$sd_blues_env
  env_levels <- levels(stage1_data$Env)
  n_env <- length(env_levels)

  K_eff_tpe <- sim_data_list$K_eff_from_tpe
  K_gen_tpe <- sim_data_list$K_gen_from_tpe
  if (!is.null(K_eff_tpe) && K_eff_tpe > 0L) {
    limit_k <- min(max_k, K_eff_tpe + 1L, n_env - 1)
  } else if (!is.null(K_gen_tpe) && K_gen_tpe > 0L) {
    limit_k <- min(max_k, K_gen_tpe, n_env - 1)
  } else {
    limit_k <- min(max_k, n_env - 1)
  }
  if (limit_k < 1) limit_k <- 1

  results_list <- list()

  model_specs <- build_specs(limit_k)

  for (spec in model_specs) {
    res_k <- .empty_result(spec$type, spec$label, spec$k)

    current_model <- tryCatch(
      {
        capture.output({
          mod <- .fit_model_asreml(stage1_data, spec$type, spec$k, spec$rc, maxiter = maxiter, workspace = workspace)
        })
        mod
      },
      error = function(e) {
        return(e)
      }
    )

    if (inherits(current_model, "error") || is.null(current_model) ||
      !isTRUE(current_model$converge)) {
      res_k$Status <- "Failed_Converge"
      res_k$Error_Log <- if (inherits(current_model, "error")) conditionMessage(current_model) else "No convergence"
      results_list[[length(results_list) + 1]] <- res_k
      next
    }

    varcomp_info <- tryCatch(
      {
        .extract_varcomp_asreml(current_model, spec$type, spec$k, env_levels)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (is.null(varcomp_info)) {
      res_k$Status <- "Var_Extract_Error"
      results_list[[length(results_list) + 1]] <- res_k
      next
    }

    res_k$Var_Tot_Mean <- varcomp_info$Var_Tot
    res_k$Var_FA_Mean <- varcomp_info$Var_FA
    res_k$Var_Psi_Mean <- varcomp_info$Var_Psi
    res_k$GxE_Est_Cor <- varcomp_info$GxE_Est
    res_k$Mean_GenCor <- varcomp_info$Mean_GenCor
    res_k$Perc_Var_Expl <- varcomp_info$Perc_Var_Expl

    res_k <- tryCatch(
      {
        .evaluate_predictions(
          current_model, stage1_info, sim_data_list,
          varcomp_info, res_k
        )
      },
      error = function(e) {
        res_k$Status <- "Prediction_Error"
        res_k$Error_Log <- conditionMessage(e)
        return(res_k)
      }
    )

    if (res_k$Status == "Init") {
      res_k$Status <- "Converged"
      res_k$Error_Log <- "Success"
    }

    results_list[[length(results_list) + 1]] <- res_k
  }

  final_df <- dplyr::bind_rows(results_list)
  return(final_df)
}

# ==============================================================================
# SECTION 8: SCENARIO RUNNER
# ==============================================================================

#' Run a grid of scenarios (main high-level function)
#'
#' @param params        metSim_params
#' @param n_locations   Vector of location counts: c(5, 10, 20, 30, 40)
#' @param n_reps_sim    Number of Monte Carlo simulations (default 1000)
#' @param n_individuals Number of candidate genotypes (default 30)
#' @param n_envs_tpe    Number of virtual environments in the TPE (default 1000)
#' @param n_reps_met    Reps per location in the MET (default 2)
#' @param max_k         Maximum K for FA (default 6)
#' @param checks        List of checks (default default_checks())
#' @param models        Models: a subset of c("MAIN", "ID", "DIAG", "CS", "FA").
#'   MAIN (~Hybrid) = main effect, cor=1.
#'   CS (~Hybrid + Env:Hybrid) = Compound Symmetry; estimates GxE_Est_Cor (% interaction).
#' @param het_residual  Fit versions with heterogeneous residual?
#' @param seed          Master seed
#' @param progress      Show progress?
#' @param on_error      "skip" or "stop"
#' @return data.frame with results from all simulations
run_scenario_grid <- function(params,
                              n_locations,
                              n_reps_sim = 1000,
                              n_individuals = params$n_individuals,
                              n_envs_tpe = 1000,
                              n_reps_met = 2,
                              max_k = 6,
                              checks = default_checks(),
                              models = c("ID", "DIAG", "FA"),
                              res_configs = NULL,
                              het_residual = FALSE,
                              sparse_fracs = NULL,
                              sparse_only_fa = FALSE,
                              seed = 42,
                              progress = TRUE,
                              on_error = "skip",
                              maxiter = 1500,
                              workspace = "6gb") {
  log_msg(
    ">>> Starting simulation: %d location scenarios x %d simulations <<<",
    length(n_locations), n_reps_sim
  )
  log_msg(
    "    Models: %s | Het residual: %s",
    paste(models, collapse = ", "), het_residual
  )

  if (is.null(sparse_fracs)) {
    if (!is.null(params$sparse_frac)) {
      sparse_fracs <- params$sparse_frac
    } else {
      sparse_fracs <- 1.0
    }
  }

  results_list <- list()
  counter <- 1
  total_falhas <- 0

  for (n_loc in n_locations) {
    log_msg("--- Scenario: %d Locations ---", n_loc)

    for (sfrac in sparse_fracs) {
      if (length(sparse_fracs) > 1) log_msg("    Fraction: %.2f", sfrac)

      is_sparse_scenario <- FALSE
      if (sfrac < 1.0) is_sparse_scenario <- TRUE
      if (!is.null(params$n_loc_por_candidato) && params$n_loc_por_candidato > 0 && params$n_loc_por_candidato < n_loc) {
        is_sparse_scenario <- TRUE
      }

      if (!is.null(params$n_loc_por_candidato) && params$n_loc_por_candidato > 0 && params$n_loc_por_candidato >= n_loc) {
        warning(sprintf("For the scenario with %d locations: n_loc_por_candidato (%d) >= n_locations (%d). Each candidate genotype will be tested in all locations (full testing), disabling Sparse Testing for this grid.", n_loc, params$n_loc_por_candidato, n_loc))
      }

      models_iter <- models
      if (is_sparse_scenario) {
        msg_detail <- if (sfrac < 1.0) sprintf("Fraction: %.2f", sfrac) else sprintf("Absolute: %d locations", params$n_loc_por_candidato)
        if (sparse_only_fa) {
          models_iter <- intersect(models, c("FA", "MAIN", "CS", "CORH"))
          if (length(models_iter) == 0) {
            stop(sprintf("Error: sparse_only_fa=TRUE (%s), but no 'FA'/'MAIN'/'CS' model was included.", msg_detail))
          }
          if (length(models_iter) < length(models)) {
            message(sprintf("  [!] sparse_only_fa=TRUE (%s): ID/DIAG removed; running only FA/MAIN/CS.", msg_detail))
          }
        } else if (length(intersect(models, c("ID", "DIAG"))) > 0) {
          message(sprintf("  [i] Sparse active (%s): ID/DIAG ALSO fitted. In untested cells, ID/DIAG (cor=0) predict the environment mean; MAIN (cor=1), CS (G+GxE) and FA (structured cor) replicate/borrow the genotype deviation.", msg_detail))
        }
      }

      for (sim in 1:n_reps_sim) {
        seed_tpe <- derivar_seed(seed, 1L, sim)
        seed_gen <- derivar_seed(seed, 2L, sim)
        seed_met <- derivar_seed(seed, 3L, sim) + round(sfrac * 100)

        tpe <- generate_tpe(params, n_virtual = n_envs_tpe, seed = seed_tpe)

        gen_obj <- simulate_genetics(tpe,
          n_individuals = n_individuals,
          checks = checks, seed = seed_gen
        )

        sim_res <- realize_met(gen_obj, tpe,
          n_locations = n_loc,
          n_reps = n_reps_met, seed = seed_met,
          n_loc_por_candidato = params$n_loc_por_candidato,
          sparse_frac = sfrac,
          missing_plot_pct = if (is.null(params$missing_plot_pct)) 0 else params$missing_plot_pct,
          warn_sparse = FALSE
        )

        resultado <- tryCatch(
          {
            capture.output({
              res <- analyze_met(sim_res,
                models = models_iter,
                max_k = max_k, res_configs = res_configs, het_residual = het_residual,
                maxiter = maxiter, workspace = workspace
              )
            })
            res
          },
          error = function(e) {
            if (on_error == "stop") stop(e)
            return(NULL)
          }
        )

        if (is.null(resultado) || nrow(resultado) == 0) {
          total_falhas <- total_falhas + 1
          next
        }

        resultado$N_Locais <- n_loc
        resultado$Sparse_Frac <- sim_res$Sparse_Frac_Efetiva
        resultado$N_Loc_Candidato <- sim_res$N_Loc_Candidato
        resultado$Sim_ID <- sim
        resultado$N_Attempted <- n_reps_sim

        results_list[[counter]] <- resultado
        counter <- counter + 1

        if (progress && sim %% 50 == 0) cat(".")
        if (sim %% 50 == 0) gc(verbose = FALSE)
      }
      if (progress) cat("\n")
    }
  }

  log_msg(">>> Finished! Successes: %d | Failures: %d <<<", counter - 1, total_falhas)

  df_results <- dplyr::bind_rows(results_list)
  if (is.null(df_results)) {
    warning("No simulation converged!")
    return(NULL)
  }

  df_conv <- df_results
  df_conv$Conv <- as.integer(df_conv$Status == "Converged")
  conv_rates <- aggregate(Conv ~ Model_Label + N_Locais, data = df_conv, FUN = mean)
  low_conv <- conv_rates[conv_rates$Conv < 0.70, ]
  if (nrow(low_conv) > 0) {
    warning(sprintf(
      "Models with convergence rate < 70%%: %s.",
      paste(unique(low_conv$Model_Label), collapse = ", ")
    ))
  }

  df_results <- df_results %>%
    dplyr::select(N_Locais, Sparse_Frac, N_Loc_Candidato, Sim_ID, Model_Type, Model_Label, K, Status, everything())

  class(df_results) <- c("metSim_results", class(df_results))
  return(df_results)
}

# ==============================================================================
# SECTION 9: SUMMARIZATION AND AGGREGATE METRICS
# ==============================================================================

#' Summarize results: Power, Bias, RMSE of the differences
#'
#' @param results data.frame from run_scenario_grid()
#' @return summarized data.frame
summarize_results <- function(results) {
  true_diffs <- c(D15 = 15, D10 = 10, D05 = 5, D02 = 2)

  group_cols <- c("N_Locais", "Model_Type", "Model_Label")
  if ("Scenario" %in% names(results)) {
    group_cols <- c("Scenario", group_cols)
  }
  if ("Sparse_Frac" %in% names(results)) {
    group_cols <- c(group_cols, "Sparse_Frac", "N_Loc_Candidato")
  }

  group_cols_tpe <- "N_Locais"
  if ("Scenario" %in% names(results)) {
    group_cols_tpe <- c("Scenario", group_cols_tpe)
  }
  if ("Sparse_Frac" %in% names(results)) {
    group_cols_tpe <- c(group_cols_tpe, "Sparse_Frac", "N_Loc_Candidato")
  }

  if ("N_Attempted" %in% names(results)) {
    attempts <- results %>%
      group_by(across(all_of(group_cols_tpe))) %>%
      summarise(N_Attempted = max(N_Attempted, na.rm = TRUE), .groups = "drop")
  } else {
    attempts <- results %>%
      group_by(across(all_of(group_cols_tpe))) %>%
      summarise(N_Attempted = max(Sim_ID, na.rm = TRUE), .groups = "drop")
  }

  res <- results %>% dplyr::filter(Status == "Converged")

  summary_df <- res %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      N_Sims = n(),

      Acc_Glo_se = sd(Acc_Glo, na.rm = TRUE) / sqrt(n()),
      Acc_Loc_se = sd(Acc_Loc, na.rm = TRUE) / sqrt(n()),
      GxE_Est_mean = mean(GxE_Est_Cor, na.rm = TRUE),
      GxE_Est_se = sd(GxE_Est_Cor, na.rm = TRUE) / sqrt(n()),
      Mean_GenCor_mean = mean(Mean_GenCor, na.rm = TRUE),

      Acc_Glo_mean = mean(Acc_Glo, na.rm = TRUE),
      Acc_Glo_sd = sd(Acc_Glo, na.rm = TRUE),
      Acc_Loc_mean = mean(Acc_Loc, na.rm = TRUE),
      Acc_Loc_Tested_mean = mean(Acc_Loc_Tested, na.rm = TRUE),
      Acc_Loc_Untested_mean = mean(Acc_Loc_Untested, na.rm = TRUE),
      Acc_Loc_Untested_sd = sd(Acc_Loc_Untested, na.rm = TRUE),
      Spearman_Loc_Tested_mean = mean(Spearman_Loc_Tested, na.rm = TRUE),
      Spearman_Loc_Untested_mean = mean(Spearman_Loc_Untested, na.rm = TRUE),
      N_Tested_Obs_mean = mean(N_Tested_Obs, na.rm = TRUE),
      N_Untested_Obs_mean = mean(N_Untested_Obs, na.rm = TRUE),

      RMSE_Glo_BLUP_mean = mean(RMSE_Glo_BLUP, na.rm = TRUE),
      RMSE_Glo_REML_mean = mean(RMSE_Glo_REML, na.rm = TRUE),
      RMSE_Glo_Cullis_mean = mean(RMSE_Glo_Cullis, na.rm = TRUE),
      RMSE_Glo_Pheno_mean = mean(RMSE_Glo_Pheno, na.rm = TRUE),

      MAE_Glo_BLUP_mean = mean(MAE_Glo_BLUP, na.rm = TRUE),
      MAE_Glo_REML_mean = mean(MAE_Glo_REML, na.rm = TRUE),
      MAE_Glo_Cullis_mean = mean(MAE_Glo_Cullis, na.rm = TRUE),
      MAE_Glo_Pheno_mean = mean(MAE_Glo_Pheno, na.rm = TRUE),

      Power_D15_Glo = mean(D15_Glo_BLUP > 0, na.rm = TRUE) * 100,
      Power_D10_Glo = mean(D10_Glo_BLUP > 0, na.rm = TRUE) * 100,
      Power_D05_Glo = mean(D05_Glo_BLUP > 0, na.rm = TRUE) * 100,
      Power_D02_Glo = mean(D02_Glo_BLUP > 0, na.rm = TRUE) * 100,

      Bias_D15_Glo_BLUP = mean(D15_Glo_BLUP - 15, na.rm = TRUE),
      Bias_D10_Glo_BLUP = mean(D10_Glo_BLUP - 10, na.rm = TRUE),
      Bias_D05_Glo_BLUP = mean(D05_Glo_BLUP - 5, na.rm = TRUE),
      Bias_D02_Glo_BLUP = mean(D02_Glo_BLUP - 2, na.rm = TRUE),
      Bias_D15_Glo_REML = mean(D15_Glo_REML - 15, na.rm = TRUE),
      Bias_D10_Glo_REML = mean(D10_Glo_REML - 10, na.rm = TRUE),
      Bias_D05_Glo_REML = mean(D05_Glo_REML - 5, na.rm = TRUE),
      Bias_D02_Glo_REML = mean(D02_Glo_REML - 2, na.rm = TRUE),
      Bias_D15_Glo_Cullis = mean(D15_Glo_Cullis - 15, na.rm = TRUE),
      Bias_D10_Glo_Cullis = mean(D10_Glo_Cullis - 10, na.rm = TRUE),
      Bias_D05_Glo_Cullis = mean(D05_Glo_Cullis - 5, na.rm = TRUE),
      Bias_D02_Glo_Cullis = mean(D02_Glo_Cullis - 2, na.rm = TRUE),
      Bias_D15_Glo_Pheno = mean(D15_Glo_Pheno - 15, na.rm = TRUE),
      Bias_D10_Glo_Pheno = mean(D10_Glo_Pheno - 10, na.rm = TRUE),
      Bias_D05_Glo_Pheno = mean(D05_Glo_Pheno - 5, na.rm = TRUE),
      Bias_D02_Glo_Pheno = mean(D02_Glo_Pheno - 2, na.rm = TRUE),

      RMSE_D15_Glo_BLUP = sqrt(mean((D15_Glo_BLUP - 15)^2, na.rm = TRUE)),
      RMSE_D10_Glo_BLUP = sqrt(mean((D10_Glo_BLUP - 10)^2, na.rm = TRUE)),
      RMSE_D05_Glo_BLUP = sqrt(mean((D05_Glo_BLUP - 5)^2, na.rm = TRUE)),
      RMSE_D02_Glo_BLUP = sqrt(mean((D02_Glo_BLUP - 2)^2, na.rm = TRUE)),
      RMSE_D15_Glo_REML = sqrt(mean((D15_Glo_REML - 15)^2, na.rm = TRUE)),
      RMSE_D10_Glo_REML = sqrt(mean((D10_Glo_REML - 10)^2, na.rm = TRUE)),
      RMSE_D05_Glo_REML = sqrt(mean((D05_Glo_REML - 5)^2, na.rm = TRUE)),
      RMSE_D02_Glo_REML = sqrt(mean((D02_Glo_REML - 2)^2, na.rm = TRUE)),
      RMSE_D15_Glo_Cullis = sqrt(mean((D15_Glo_Cullis - 15)^2, na.rm = TRUE)),
      RMSE_D10_Glo_Cullis = sqrt(mean((D10_Glo_Cullis - 10)^2, na.rm = TRUE)),
      RMSE_D05_Glo_Cullis = sqrt(mean((D05_Glo_Cullis - 5)^2, na.rm = TRUE)),
      RMSE_D02_Glo_Cullis = sqrt(mean((D02_Glo_Cullis - 2)^2, na.rm = TRUE)),
      RMSE_D15_Glo_Pheno = sqrt(mean((D15_Glo_Pheno - 15)^2, na.rm = TRUE)),
      RMSE_D10_Glo_Pheno = sqrt(mean((D10_Glo_Pheno - 10)^2, na.rm = TRUE)),
      RMSE_D05_Glo_Pheno = sqrt(mean((D05_Glo_Pheno - 5)^2, na.rm = TRUE)),
      RMSE_D02_Glo_Pheno = sqrt(mean((D02_Glo_Pheno - 2)^2, na.rm = TRUE)),

      MAE_D15_Glo_BLUP = mean(abs(D15_Glo_BLUP - 15), na.rm = TRUE),
      MAE_D10_Glo_BLUP = mean(abs(D10_Glo_BLUP - 10), na.rm = TRUE),
      MAE_D05_Glo_BLUP = mean(abs(D05_Glo_BLUP - 5), na.rm = TRUE),
      MAE_D02_Glo_BLUP = mean(abs(D02_Glo_BLUP - 2), na.rm = TRUE),
      MAE_D15_Glo_REML = mean(abs(D15_Glo_REML - 15), na.rm = TRUE),
      MAE_D10_Glo_REML = mean(abs(D10_Glo_REML - 10), na.rm = TRUE),
      MAE_D05_Glo_REML = mean(abs(D05_Glo_REML - 5), na.rm = TRUE),
      MAE_D02_Glo_REML = mean(abs(D02_Glo_REML - 2), na.rm = TRUE),
      MAE_D15_Glo_Cullis = mean(abs(D15_Glo_Cullis - 15), na.rm = TRUE),
      MAE_D10_Glo_Cullis = mean(abs(D10_Glo_Cullis - 10), na.rm = TRUE),
      MAE_D05_Glo_Cullis = mean(abs(D05_Glo_Cullis - 5), na.rm = TRUE),
      MAE_D02_Glo_Cullis = mean(abs(D02_Glo_Cullis - 2), na.rm = TRUE),
      MAE_D15_Glo_Pheno = mean(abs(D15_Glo_Pheno - 15), na.rm = TRUE),
      MAE_D10_Glo_Pheno = mean(abs(D10_Glo_Pheno - 10), na.rm = TRUE),
      MAE_D05_Glo_Pheno = mean(abs(D05_Glo_Pheno - 5), na.rm = TRUE),
      MAE_D02_Glo_Pheno = mean(abs(D02_Glo_Pheno - 2), na.rm = TRUE),

      Type_I_1pct_Glo_BLUP = mean(abs(D00_Glo_BLUP) > 1, na.rm = TRUE) * 100,
      Type_I_1pct_Glo_REML = mean(abs(D00_Glo_REML) > 1, na.rm = TRUE) * 100,
      Type_I_1pct_Glo_Cullis = mean(abs(D00_Glo_Cullis) > 1, na.rm = TRUE) * 100,
      Type_I_1pct_Glo_Pheno = mean(abs(D00_Glo_Pheno) > 1, na.rm = TRUE) * 100,
      Type_I_2pct_Glo_BLUP = mean(abs(D00_Glo_BLUP) > 2, na.rm = TRUE) * 100,
      Type_I_2pct_Glo_REML = mean(abs(D00_Glo_REML) > 2, na.rm = TRUE) * 100,
      Type_I_2pct_Glo_Cullis = mean(abs(D00_Glo_Cullis) > 2, na.rm = TRUE) * 100,
      Type_I_2pct_Glo_Pheno = mean(abs(D00_Glo_Pheno) > 2, na.rm = TRUE) * 100,
      D00_Glo_BLUP_mean = mean(D00_Glo_BLUP, na.rm = TRUE),
      D00_Glo_REML_mean = mean(D00_Glo_REML, na.rm = TRUE),
      D00_Glo_Cullis_mean = mean(D00_Glo_Cullis, na.rm = TRUE),
      D00_Glo_Pheno_mean = mean(D00_Glo_Pheno, na.rm = TRUE),

      Spearman_Glo_mean = mean(Spearman_Glo, na.rm = TRUE),
      Spearman_Loc_mean = mean(Spearman_Loc, na.rm = TRUE),

      CS_10_mean = mean(CS_10, na.rm = TRUE),
      CS_20_mean = mean(CS_20, na.rm = TRUE),

      DG_10_mean = mean(DG_10, na.rm = TRUE),
      DG_20_mean = mean(DG_20, na.rm = TRUE),
      DG_Max_10_mean = mean(DG_Max_10, na.rm = TRUE),
      DG_Max_20_mean = mean(DG_Max_20, na.rm = TRUE),
      DGR_10_mean = mean(DGR_10, na.rm = TRUE),
      DGR_20_mean = mean(DGR_20, na.rm = TRUE),

      P_Inv_D05_Glo = mean(D05_Glo_BLUP < 0, na.rm = TRUE) * 100,
      P_Inv_D02_Glo = mean(D02_Glo_BLUP < 0, na.rm = TRUE) * 100,

      P_Top1_Hit = mean(Top1_Hit, na.rm = TRUE) * 100,
      P_Top1_in_Top2 = mean(Top1_in_Top2, na.rm = TRUE) * 100,
      P_Top1_in_Top3 = mean(Top1_in_Top3, na.rm = TRUE) * 100,
      P_Top2_Hit = mean(Top2_Hit, na.rm = TRUE) * 100,
      P_Top3_Hit = mean(Top3_Hit, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    left_join(attempts, by = group_cols_tpe) %>%
    mutate(Conv_Rate = N_Sims / N_Attempted * 100) %>%
    dplyr::select(all_of(group_cols), N_Sims, Conv_Rate, everything(), -N_Attempted)

  return(summary_df)
}

#' Combine multiple scenarios
collect_scenarios <- function(...) {
  scenarios <- list(...)
  if (is.null(names(scenarios))) {
    names(scenarios) <- paste0("Scenario_", seq_along(scenarios))
  }
  do.call(rbind, lapply(names(scenarios), function(nm) {
    df <- scenarios[[nm]]
    df$Scenario <- nm
    df
  }))
}

# ==============================================================================
# SECTION 10: VISUALIZATIONS
# ==============================================================================

theme_metsim <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#' Plot power curves by scenario and model
#'
#' @param results  Raw result from run_scenario_grid()
#' @param level    "Glo" (global) or "Loc" (local)
plot_power_curve <- function(results, level = "Glo") {
  res <- results %>% dplyr::filter(Status == "Converged")

  diff_cols <- paste0("D", c("15", "10", "05", "02"), "_", level, "_BLUP")
  true_vals <- c(15, 10, 5, 2)

  power_data <- res %>%
    group_by(N_Locais, Model_Label) %>%
    summarise(
      Power_D15 = mean(.data[[diff_cols[1]]] > 0, na.rm = TRUE) * 100,
      Power_D10 = mean(.data[[diff_cols[2]]] > 0, na.rm = TRUE) * 100,
      Power_D05 = mean(.data[[diff_cols[3]]] > 0, na.rm = TRUE) * 100,
      Power_D02 = mean(.data[[diff_cols[4]]] > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = starts_with("Power_"),
      names_to = "Magnitude", values_to = "Power"
    ) %>%
    mutate(Magnitude = gsub("Power_", "", Magnitude))

  ggplot(power_data, aes(
    x = N_Locais, y = Power, color = Magnitude,
    linetype = Magnitude
  )) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "grey50", linewidth = 0.7) +
    annotate("text",
      x = max(power_data$N_Locais) * 0.95, y = 82,
      label = "80% poder", hjust = 1, color = "grey50", size = 3.5
    ) +
    facet_wrap(~Model_Label, scales = "fixed") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    scale_color_manual(
      values = c(
        "D02" = "#e74c3c", "D05" = "#e67e22",
        "D10" = "#2ecc71", "D15" = "#3498db"
      ),
      labels = c(
        "D02" = "2%", "D05" = "5%",
        "D10" = "10%", "D15" = "15%"
      )
    ) +
    labs(
      title = "Detection Power vs. Number of Locations",
      subtitle = "Proportion of simulations where the sign of the difference was correct (calibration-invariant)",
      x = "Number of Locations", y = "Detection Power (%)",
      color = "Diferenca Real", linetype = "Diferenca Real"
    ) +
    theme_metsim()
}

#' Accuracy boxplots by scenario and model
#'
#' @param results  Raw result from run_scenario_grid()
#' @param level    "Glo" (global) or "Loc" (local)
plot_accuracy_summary <- function(results, level = "Glo") {
  res <- results %>% dplyr::filter(Status == "Converged")
  acc_col <- paste0("Acc_", level)

  ggplot(res, aes(
    x = factor(N_Locais), y = .data[[acc_col]],
    fill = Model_Label
  )) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Prediction Accuracy (", level, ")"),
      subtitle = "Correlation between predicted and true (calibration-invariant)",
      x = "Number of Locations", y = "Acuracia (correlacao)",
      fill = "Modelo"
    ) +
    theme_metsim()
}

#' Direct comparison between models
#'
#' @param results  Raw result from run_scenario_grid()
#' @param metric   Metric to plot (default: "Acc_Glo")
plot_model_comparison <- function(results, metric = "Acc_Glo") {
  res <- results %>% dplyr::filter(Status == "Converged")

  summary_data <- res %>%
    group_by(N_Locais, Model_Label) %>%
    summarise(
      Mean = mean(.data[[metric]], na.rm = TRUE),
      SD = sd(.data[[metric]], na.rm = TRUE),
      Lower = Mean - 1.96 * SD / sqrt(n()),
      Upper = Mean + 1.96 * SD / sqrt(n()),
      .groups = "drop"
    )

  ggplot(summary_data, aes(x = N_Locais, y = Mean, color = Model_Label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Model_Label),
      alpha = 0.15, color = NA
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste0("Model Comparison: ", metric),
      subtitle = "Mean +/- 95% CI across simulations",
      x = "Number of Locations", y = metric,
      color = "Modelo", fill = "Modelo"
    ) +
    theme_metsim()
}

#' Bias and RMSE of the differences
#'
#' @param results  Raw result from run_scenario_grid()
#' @param method   Output scale (default: "Pheno")
plot_bias_rmse <- function(results, method = "Pheno") {
  res <- results %>% dplyr::filter(Status == "Converged")

  bias_data <- res %>%
    group_by(N_Locais, Model_Label) %>%
    summarise(
      Bias_D15 = mean(.data[[paste0("D15_Glo_", method)]] - 15, na.rm = TRUE),
      Bias_D10 = mean(.data[[paste0("D10_Glo_", method)]] - 10, na.rm = TRUE),
      Bias_D05 = mean(.data[[paste0("D05_Glo_", method)]] - 5, na.rm = TRUE),
      Bias_D02 = mean(.data[[paste0("D02_Glo_", method)]] - 2, na.rm = TRUE),
      RMSE_D15 = sqrt(mean((.data[[paste0("D15_Glo_", method)]] - 15)^2, na.rm = TRUE)),
      RMSE_D10 = sqrt(mean((.data[[paste0("D10_Glo_", method)]] - 10)^2, na.rm = TRUE)),
      RMSE_D05 = sqrt(mean((.data[[paste0("D05_Glo_", method)]] - 5)^2, na.rm = TRUE)),
      RMSE_D02 = sqrt(mean((.data[[paste0("D02_Glo_", method)]] - 2)^2, na.rm = TRUE)),
      .groups = "drop"
    )

  bias_long <- bias_data %>%
    pivot_longer(
      cols = starts_with("Bias_"),
      names_to = "Magnitude", values_to = "Bias"
    ) %>%
    mutate(Magnitude = gsub("Bias_", "", Magnitude))

  p_bias <- ggplot(bias_long, aes(x = N_Locais, y = Bias, color = Magnitude)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~Model_Label) +
    scale_color_manual(values = c(
      "D02" = "#e74c3c", "D05" = "#e67e22",
      "D10" = "#2ecc71", "D15" = "#3498db"
    )) +
    labs(
      title = paste0("Bias of the Estimated Differences (", method, ")"),
      subtitle = "Negative values = underestimation",
      x = "N Locations", y = "Vies (pp)", color = "Diferenca"
    ) +
    theme_metsim()

  rmse_long <- bias_data %>%
    pivot_longer(
      cols = starts_with("RMSE_"),
      names_to = "Magnitude", values_to = "RMSE"
    ) %>%
    mutate(Magnitude = gsub("RMSE_", "", Magnitude))

  p_rmse <- ggplot(rmse_long, aes(x = N_Locais, y = RMSE, color = Magnitude)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~Model_Label) +
    scale_color_manual(values = c(
      "D02" = "#e74c3c", "D05" = "#e67e22",
      "D10" = "#2ecc71", "D15" = "#3498db"
    )) +
    labs(
      title = paste0("RMSE of the Estimated Differences (", method, ")"),
      x = "N Locations", y = "RMSE (pp)", color = "Diferenca"
    ) +
    theme_metsim()

  return(list(bias = p_bias, rmse = p_rmse))
}

# ==============================================================================
# SECTION 11: POST-SIMULATION OPTIMIZATION AND DASHBOARD
# ==============================================================================

#' Optimize the number of locations for a specific metric
#'
#' @param results               Result from run_scenario_grid() or summarize_results()
#' @param metric                Name of the metric column to optimize (default: "Power_D05_Glo")
#' @param target_value          Target value of the metric (default: 85)
#' @param interpolation_method  Interpolation method: "spline" (default) or "linear"
#' @param plot                  Logical. If TRUE, plots the optimization curves
#' @return List with the optimization table (opt_table) and the plot
metSim_optimize <- function(results,
                            metric = "Power_D05_Glo",
                            target_value = 85,
                            interpolation_method = "spline",
                            plot = TRUE) {
  if ("Status" %in% names(results)) {
    summary_data <- summarize_results(results)
  } else {
    summary_data <- results
  }

  if (!(metric %in% names(summary_data))) {
    stop(paste("Metric", metric, "not found. Choose one of the summary columns (e.g. Power_D05_Glo, Acc_Glo_mean, CS_10_mean, P_Top1_in_Top2)."))
  }

  max_observed <- max(summary_data[[metric]], na.rm = TRUE)
  if (max_observed <= 1.0 && target_value > 1.0) {
    warning(paste0("The metric '", metric, "' is on the 0-1 scale, but target_value = ", target_value, ". Adjusting target_value to ", target_value / 100, "."))
    target_value <- target_value / 100
  }

  models <- unique(summary_data$Model_Label)
  opt_list <- list()
  plot_df_list <- list()

  for (m in models) {
    sub_data <- summary_data %>%
      dplyr::filter(Model_Label == m) %>%
      arrange(N_Locais)
    if (nrow(sub_data) < 2) next

    x <- sub_data$N_Locais
    y <- sub_data[[metric]]

    if (interpolation_method == "spline") {
      f_interp <- tryCatch(
        {
          splinefun(x, y, method = "monoH.FC")
        },
        error = function(e) {
          splinefun(x, y)
        }
      )
    } else {
      f_interp <- approxfun(x, y, rule = 2)
    }

    x_grid <- seq(min(x), max(x), length.out = 1000)
    y_grid <- f_interp(x_grid)

    plot_df_list[[m]] <- data.frame(
      N_Locais = x_grid,
      Value = y_grid,
      Model_Label = m
    )

    crossing_idx <- which(diff(y_grid >= target_value) != 0)

    if (y_grid[1] >= target_value) {
      exact_n <- min(x)
      status_msg <- "Already reached at the minimum number of locations."
    } else if (all(y_grid < target_value)) {
      exact_n <- NA_real_
      status_msg <- paste0("Not reached (max observed: ", round(max(y), 2), ")")
    } else {
      idx <- crossing_idx[1]
      x1 <- x_grid[idx]
      x2 <- x_grid[idx + 1]
      y1 <- y_grid[idx]
      y2 <- y_grid[idx + 1]
      exact_n <- if (abs(y2 - y1) < 1e-9) x1 else x1 + (target_value - y1) * (x2 - x1) / (y2 - y1)
      status_msg <- "Atingido."
    }

    opt_list[[m]] <- data.frame(
      Model_Label = m,
      Exact_N = exact_n,
      Sugerido_N = if (is.na(exact_n)) NA_integer_ else as.integer(ceiling(exact_n)),
      Status = status_msg,
      stringsAsFactors = FALSE
    )
  }

  opt_df <- do.call(rbind, opt_list)
  rownames(opt_df) <- NULL

  p <- NULL
  if (plot && length(plot_df_list) > 0) {
    plot_data <- do.call(rbind, plot_df_list)
    points_data <- summary_data %>%
      dplyr::filter(Model_Label %in% models) %>%
      dplyr::select(N_Locais, Value = dplyr::all_of(metric), Model_Label)

    p <- ggplot() +
      geom_line(data = plot_data, aes(x = N_Locais, y = Value, color = Model_Label), linewidth = 1.1) +
      geom_point(data = points_data, aes(x = N_Locais, y = Value, color = Model_Label), size = 3) +
      geom_hline(yintercept = target_value, linetype = "dashed", color = "grey40") +
      labs(
        title = paste("Location Optimization - Metric:", metric),
        subtitle = paste("Dashed line: Target of", target_value),
        x = "Number of Locations",
        y = metric,
        color = "Modelo"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom"
      )

    valid_crossings <- opt_df %>% dplyr::filter(!is.na(Exact_N))
    if (nrow(valid_crossings) > 0) {
      p <- p + geom_segment(
        data = valid_crossings,
        aes(x = Exact_N, xend = Exact_N, y = 0, yend = target_value, color = Model_Label),
        linetype = "dotted", linewidth = 0.8
      )
    }

    print(p)
  }

  return(list(opt_table = opt_df, plot = p))
}

#' Generate a consolidated sizing table with report MAE
#'
#' @param results             Result from run_scenario_grid() or summarize_results()
#' @param target_power        Statistical power target (default: 85)
#' @param target_acc          Global accuracy target (default: 80)
#' @param acc_metric          Accuracy metric: "pearson" (default) or "spearman"
#' @param target_cs           Coincidence-of-Selection CS10 target (default: 80)
#' @param target_top1_top2    Target for the Top 1 being within the predicted Top 2 (default: 80)
#' @param calibration_method  Output scale to read the targets from: "Pheno" (default), "BLUP", "REML" or "Cullis"
#' @return formatted data.frame
metSim_power_table <- function(results,
                               target_power = 85,
                               target_acc = 80,
                               acc_metric = "pearson",
                               target_cs = 80,
                               target_top1_top2 = 80,
                               calibration_method = "Pheno") {
  if ("Status" %in% names(results)) {
    summary_data <- summarize_results(results)
  } else {
    summary_data <- results
  }

  message("\n================================================================================")
  message("⚠️  [NOTE ON MAE & SHRINKAGE]")
  message("   The MAE (Mean Absolute Error) under the 'BLUP' calibration (raw BLUP) is artificially")
  message("   underestimated due to Ridge shrinkage. The model reduces the variance")
  message("   of the predictions to minimize the squared error, hiding the true physical error.")
  message("   To assess the real error magnitude on the field-data scale, always use")
  message("   the estimates under the 'Pheno' calibration.")
  message("================================================================================\n")

  pow_col <- "Power_D05_Glo"

  if (tolower(acc_metric) == "pearson") {
    acc_col <- "Acc_Glo_mean"
  } else {
    acc_col <- "Spearman_Glo_mean"
  }

  cs_col <- "CS_10_mean"
  top_col <- "P_Top1_in_Top2"
  mae_col <- paste0("MAE_Glo_", calibration_method, "_mean")

  opt_pow <- metSim_optimize(summary_data, metric = pow_col, target_value = target_power, plot = FALSE)$opt_table
  opt_acc <- metSim_optimize(summary_data, metric = acc_col, target_value = target_acc, plot = FALSE)$opt_table
  opt_cs <- metSim_optimize(summary_data, metric = cs_col, target_value = target_cs, plot = FALSE)$opt_table
  opt_top <- metSim_optimize(summary_data, metric = top_col, target_value = target_top1_top2, plot = FALSE)$opt_table

  models <- unique(summary_data$Model_Label)
  mae_at_opt <- numeric(length(models))
  names(mae_at_opt) <- models

  for (m in models) {
    sub_data <- summary_data %>%
      dplyr::filter(Model_Label == m) %>%
      arrange(N_Locais)
    opt_n <- opt_pow$Exact_N[opt_pow$Model_Label == m]

    if (is.na(opt_n) || nrow(sub_data) < 2) {
      mae_at_opt[m] <- NA_real_
    } else {
      x <- sub_data$N_Locais
      y <- sub_data[[mae_col]]
      f_mae <- tryCatch(
        {
          splinefun(x, y, method = "monoH.FC")
        },
        error = function(e) {
          splinefun(x, y)
        }
      )
      mae_at_opt[m] <- f_mae(opt_n)
    }
  }

  final_table <- data.frame(
    Modelo = models,
    Calibracao = calibration_method,
    Locais_Poder = opt_pow$Sugerido_N[match(models, opt_pow$Model_Label)],
    Locais_Acuracia = opt_acc$Sugerido_N[match(models, opt_acc$Model_Label)],
    Locais_CS10 = opt_cs$Sugerido_N[match(models, opt_cs$Model_Label)],
    Locais_Top1_in_Top2 = opt_top$Sugerido_N[match(models, opt_top$Model_Label)],
    MAE_no_Poder_Optimo = round(mae_at_opt, 4),
    stringsAsFactors = FALSE
  )

  colnames(final_table) <- c(
    "Modelo", "Calibration",
    paste0("Locations for ", target_power, "% Power (D05)"),
    paste0("Locations for ", target_acc, "% Accuracy (", toupper(acc_metric), ")"),
    paste0("Locations for ", target_cs, "% CS10"),
    paste0("Locations for ", target_top1_top2, "% Top1 in Top2"),
    "Global MAE at Optimal Power"
  )

  return(final_table)
}

#' Plot a multi-metric dashboard of MET network performance
#'
#' @param results             Result from run_scenario_grid() or summarize_results()
#' @param calibration_method  Output scale ("Pheno", "BLUP", "REML" or "Cullis")
#' @param acc_metric          Accuracy metric ("pearson" or "spearman")
#' @return ggplot object with 5 faceted sub-plots
plot_metSim_dashboard <- function(results, calibration_method = "Pheno", acc_metric = "pearson") {
  if ("Status" %in% names(results)) {
    summary_data <- summarize_results(results)
  } else {
    summary_data <- results
  }

  pow_col <- "Power_D05_Glo"
  if (tolower(acc_metric) == "pearson") {
    acc_col <- "Acc_Glo_mean"
  } else {
    acc_col <- "Spearman_Glo_mean"
  }
  cs_col <- "CS_10_mean"
  top_col <- "P_Top1_in_Top2"
  mae_col <- paste0("MAE_Glo_", calibration_method, "_mean")

  metrics_to_plot <- c(pow_col, acc_col, cs_col, top_col, mae_col)

  clean_labels <- c(
    "Power (D05, %)" = pow_col,
    "Accuracy (Global)" = acc_col,
    "Coincidence CS10 (%)" = cs_col,
    "Top1 in Top2 (%)" = top_col,
    "MAE Global" = mae_col
  )

  if ("Acc_Loc_Untested_mean" %in% names(summary_data) && any(!is.na(summary_data$Acc_Loc_Untested_mean))) {
    metrics_to_plot <- c(metrics_to_plot, "Acc_Loc_Untested_mean")
    clean_labels <- c(clean_labels, "Accuracy (Untested)" = "Acc_Loc_Untested_mean")
  }

  long_data <- summary_data %>%
    dplyr::select(N_Locais, Model_Label, dplyr::all_of(metrics_to_plot)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metrics_to_plot),
      names_to = "Metrica",
      values_to = "Valor"
    )
  lookup <- names(clean_labels)
  names(lookup) <- clean_labels

  long_data$Metrica_Clean <- lookup[long_data$Metrica]
  long_data$Metrica_Clean <- factor(long_data$Metrica_Clean, levels = names(clean_labels))

  p <- ggplot(long_data, aes(x = N_Locais, y = Valor, color = Model_Label)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    facet_wrap(~Metrica_Clean, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("MET Network Performance - Calibration: ", calibration_method),
      x = "Number of Locations in the Network",
      y = "Valor Estimado",
      color = "Modelo"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )

  return(p)
}

# ==============================================================================
# SECTION 12: SPARSE-TESTING-SPECIFIC PLOTS
# ==============================================================================

#' Plot Local Accuracy (Tested vs Untested) for Sparse Testing
#'
#' @param results Result from run_scenario_grid() or summarize_results()
#' @return ggplot object
plot_sparse_tested_vs_untested <- function(results) {
  if ("Status" %in% names(results)) {
    summary_data <- summarize_results(results)
  } else {
    summary_data <- results
  }

  if (!"Acc_Loc_Untested_mean" %in% names(summary_data)) {
    stop("The data.frame does not contain the metric Acc_Loc_Untested_mean. Did the simulation run with Sparse Testing (metSim v7.0+)?")
  }

  plot_data <- summary_data %>%
    dplyr::select(N_Locais, Model_Label, Acc_Loc_Tested_mean, Acc_Loc_Untested_mean)

  if ("Sparse_Frac" %in% names(summary_data)) {
    plot_data <- plot_data %>% dplyr::mutate(Sparse_Frac = summary_data$Sparse_Frac)
  } else {
    plot_data <- plot_data %>% dplyr::mutate(Sparse_Frac = 1.0)
  }

  long_data <- plot_data %>%
    tidyr::pivot_longer(
      cols = c(Acc_Loc_Tested_mean, Acc_Loc_Untested_mean),
      names_to = "Status_Teste",
      values_to = "Acuracia"
    ) %>%
    dplyr::mutate(
      Status_Teste = ifelse(Status_Teste == "Acc_Loc_Tested_mean", "Testados", "Nao-Testados")
    ) %>%
    dplyr::filter(!is.na(Acuracia))

  p <- ggplot(long_data, aes(x = N_Locais, y = Acuracia, color = Model_Label, linetype = Status_Teste)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    labs(
      title = "Local Accuracy: Tested vs. Untested (Sparse Testing)",
      subtitle = "FA models project untested cells; ID/DIAG omit them (line absent)",
      x = "Number of Locations in the Network",
      y = "Acuracia (Pearson)",
      color = "Modelo",
      linetype = "Grupo"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )

  if (length(unique(long_data$Sparse_Frac)) > 1) {
    p <- p + facet_wrap(~Sparse_Frac, labeller = label_both)
  }

  return(p)
}

# ==============================================================================
# END
# ==============================================================================
log_msg("metSim loaded successfully.")
