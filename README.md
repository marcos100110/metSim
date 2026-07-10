# metSim

R tools to simulate multi-environment breeding trials (MET) and compare models for
genotype-by-environment (GxE).

https://github.com/marcos100110/metSim/tree/main

## What it is

metSim makes simulated trial data where the true answer is known. It fits several
models to that data and checks how well each one finds the best genotypes.
Because the truth is set by the simulation, you can test model choices and trial
designs before running real trials.

It works in two steps: it builds a large pool of virtual environments (the TPE,
1000 by default) and then draws real trial locations from it. Models are fitted
with ASReml-R.

## What you can do with it

- Find how many locations and replicates you need for a target accuracy or power.
- Compare GxE models for low and high GxE.
- Compare new genotypes against standard checks, to help with variety pricing and
  positioning.
- Test sparse designs (each genotype in only some locations).
- See how the results change with more or fewer locations.
- Set up the simulation from your own real data.

## A few ideas (in plain words)

- **GxE**: how much the ranking of genotypes changes from one location to another.
- **TPE**: the big pool of virtual environments. The trial locations are drawn
  from it.
- **Environment / location**: used as synonyms throughout — each one is a single
  independent trial (a site drawn from the TPE). Some API names keep `Loc` /
  `location` (e.g. `n_locations`, `Acc_Loc`); they mean this same thing.
- **Checks**: reference genotypes with fixed, known true values. You use them as
  a benchmark to compare new genotypes against a standard and to measure error
  rates.
- **Output scales**: the fitted model reports each genotype's merit on different
  scales — BLUP, REML, Cullis and Pheno. They are rescalings of the same BLUP,
  so the ranking is identical; only the scale (and the size of MAE/RMSE) changes.
- **k_axes**: how complex the simulated GxE is. A higher value means more ways the
  rankings can differ between locations (the true structure has `k_axes + 1`
  factors).

## How it works

| Step | Function | What it does |
|---|---|---|
| 1. Inputs | `build_params` / `extract_params_from_asreml` | Set the inputs: correlation between locations, GxE complexity (`k_axes`), genetic and error CV, number of genotypes (`n_individuals`). Set them by hand, from a preset, or from a real fitted ASReml model. |
| 2. TPE | `generate_tpe` | Build the pool of virtual environments. |
| 3. Genetics | `simulate_genetics` | Draw the true genetic values for the genotypes. |
| 4. Trial | `realize_met` | Pick locations, add blocks and error, optional missing plots and sparse testing, and make the trial data. |
| 5. Analysis | `analyze_met` | Fit the models in ASReml-R and get the numbers. |
| 6. Grid | `run_scenario_grid` | Repeat many times over different network sizes. |
| 7. Output | `summarize_results`, `plot_*`, `metSim_power_table` | Summarize, plot and size the network. |

## Models

The fixed part is the location effect. The models differ in how they treat the
genotype across locations.

| Model | Random term | What it assumes |
|---|---|---|
| `MAIN` | `~ Hybrid` | same genotype value in every location |
| `ID` | `~ Env:Hybrid` | each location on its own (no sharing) |
| `DIAG` | `~ diag(Env):Hybrid` | each location on its own, different variance each |
| `CS` | `~ Hybrid + Env:Hybrid` | one shared value plus one interaction |
| `CORH` | `~ corh(Env):Hybrid` | shared correlation, different variance each |
| `FA(k)` | `~ fa(Env, k):Hybrid` | `k` factors that capture the GxE pattern |

FA models are fitted from 1 up to `max_k`.

## What it measures

| Column | What it means |
|---|---|
| `Acc_Loc`, `Acc_Glo` | how close the prediction is to the truth, per location and overall (a correlation) |
| `Spearman_Loc`, `Spearman_Glo` | same, but for rankings |
| `Acc_Loc_Tested`, `Acc_Loc_Untested` | accuracy where a genotype was tested and where it was not (sparse testing) |
| `MAE_Glo_*`, `RMSE_Glo_*` | error of the estimate vs the truth, one per output scale (BLUP/REML/Cullis/Pheno) |
| `GxE_Est_Cor` | estimated GxE, as a % |
| `Mean_GenCor` | average genetic correlation between locations |
| `CS_10`, `CS_20` | coincidence of selection: % of the true top 10% / 20% that the model also picks |
| `DG_10`, `DG_20` | realized genetic gain when selecting the model's top 10% / 20% |
| `Top1_Hit`, `Top2_Hit`, `Top3_Hit` | how often the true best genotype is in the model's top 1 / 2 / 3 |
| `Power_D15_Glo` … `Power_D02_Glo` | from `summarize_results`: % of simulations where a true difference (15 / 10 / 5 / 2 %) was called with the correct sign |
| `Status` | did the model fit converge |

## Requirements

- R (>= 4.0)
- **ASReml-R** (needs a paid license) — used to fit the models.
- R packages (installed automatically): `MASS`, `Matrix`, `dplyr`, `tidyr`,
  `ggplot2`, `sn`.

## Install

metSim depends on the commercial **ASReml-R**, so it is not on CRAN and
ASReml-R is not installed automatically. **Install ASReml-R first** (you need a
license) — metSim needs it present to load. Then install metSim from GitHub:

```r
# install.packages("remotes")
remotes::install_github("marcos100110/metSim")
library(metSim)
```

The other dependencies (MASS, Matrix, dplyr, tidyr, ggplot2, sn) install
automatically.

## Examples

All examples use 50 simulations (`n_reps_sim = 50`).

**1. Basic run**

```r
library(metSim)

params <- build_params(
  method   = "correlation_empirical",
  cor_mean = 0.6, k_axes = 3, psi_frac = 0,
  n_individuals = 50, mu_mean = 100, cv_env = 20,
  cvg_mean = 12, cvg_sd = 0,
  cve_mean = 12, cve_sd = 3
)

res <- run_scenario_grid(
  params,
  n_locations   = c(5, 10, 20),
  n_reps_sim    = 50,
  n_individuals = 50,
  models        = c("MAIN", "CS", "FA"),
  max_k         = 4,
  seed          = 42
)

summ <- summarize_results(res)
```

**2. Compare low and high GxE**

```r
# Low GxE: rankings stay about the same (high cor_mean)
low  <- build_params(method = "correlation_empirical", cor_mean = 0.85, k_axes = 2,
                     psi_frac = 0, n_individuals = 50, mu_mean = 100, cv_env = 20,
                     cvg_mean = 12, cvg_sd = 0, cve_mean = 12, cve_sd = 3)

# High GxE: rankings change a lot (low cor_mean)
high <- build_params(method = "correlation_empirical", cor_mean = 0.30, k_axes = 4,
                     psi_frac = 0, n_individuals = 50, mu_mean = 100, cv_env = 20,
                     cvg_mean = 12, cvg_sd = 0, cve_mean = 12, cve_sd = 3)

res_low  <- run_scenario_grid(low,  n_locations = c(10, 20), n_reps_sim = 50,
                              n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42)
res_high <- run_scenario_grid(high, n_locations = c(10, 20), n_reps_sim = 50,
                              n_individuals = 50, models = c("CS", "FA"), max_k = 5, seed = 42)
```

**3. Change the heritability**

```r
# More genetic CV and less error CV -> higher heritability -> better accuracy
h_low  <- build_params(method = "correlation_empirical", cor_mean = 0.55, k_axes = 3,
                       psi_frac = 0, n_individuals = 50, mu_mean = 100, cv_env = 20,
                       cvg_mean = 7,  cvg_sd = 0, cve_mean = 14, cve_sd = 3)

h_high <- build_params(method = "correlation_empirical", cor_mean = 0.55, k_axes = 3,
                       psi_frac = 0, n_individuals = 50, mu_mean = 100, cv_env = 20,
                       cvg_mean = 14, cvg_sd = 0, cve_mean = 7,  cve_sd = 3)

res_hlow  <- run_scenario_grid(h_low,  n_locations = c(10, 20), n_reps_sim = 50,
                               n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42)
res_hhigh <- run_scenario_grid(h_high, n_locations = c(10, 20), n_reps_sim = 50,
                               n_individuals = 50, models = c("CS", "FA"), max_k = 4, seed = 42)
```

**4. Sparse testing**

```r
res <- run_scenario_grid(
  params,
  n_locations   = c(10, 20, 30),
  n_reps_sim    = 50,
  n_individuals = 50,
  models        = c("MAIN", "CS", "FA"),
  max_k         = 4,
  sparse_fracs  = c(1.0, 0.5),
  seed          = 42
)
# Acc_Loc_Tested and Acc_Loc_Untested show how well each model predicts
# locations where a genotype was tested and where it was not.
```

**5. Save and read the results**

```r
res <- run_scenario_grid(params, n_locations = c(10, 20), n_reps_sim = 50,
                         n_individuals = 50, models = c("MAIN", "CS", "FA"),
                         max_k = 4, seed = 42)

saveRDS(res, "results.rds")       # save the raw output
summ <- summarize_results(res)    # average over the simulations
```

You can also start from a preset — `preset_low_gxe()`, `preset_moderate_gxe()`,
`preset_high_gxe()` — or set up the run from a real fitted ASReml model with
`extract_params_from_asreml()`.

## Functions

| Function | What it does |
|---|---|
| `build_params` | Set the simulation inputs by hand or from a preset. |
| `extract_params_from_asreml` | Set them from a real fitted ASReml model. |
| `preset_low_gxe` / `preset_moderate_gxe` / `preset_high_gxe` | Ready-made GxE scenarios. |
| `generate_tpe` | Build the pool of virtual environments. |
| `simulate_genetics` | Draw the true genetic values. |
| `realize_met` | Make the trial data. |
| `analyze_met` | Fit the models and get the numbers. |
| `run_scenario_grid` | Run many times over different network sizes. |
| `summarize_results` / `collect_scenarios` | Average and join results. |
| `metSim_optimize` / `metSim_power_table` | Size the network. |
| `plot_*` | Power curves, accuracy, model comparison, dashboard, sparse testing. |

## Notes

- You need ASReml-R to fit the models. The rest runs in base R plus the listed
  packages.
- Big runs (many locations, simulations and factors) can take a long time. Start
  small.
- The estimated GxE goes up with more factors, so compare GxE between models at
  the same number of factors.

## Status

**Version 1.0.0-beta** — beta, still in development. The core has been tested
across many scenarios and behaves as expected.

## License

metSim is free and open-source software, released under the **GPL-3** (GNU
General Public License, version 3 or later). Full text:
<https://www.gnu.org/licenses/gpl-3.0.html>.

Copyright © 2026 Marcos Filho.
