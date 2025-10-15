# r.ODE.RKHS

**r.ODE.RKHS** is an R package for non-parametric modeling of microbial dynamics from sparse longitudinal data. It implements the ODE-RKHS algorithm, which learns subject-specific dynamic features from microbial time series using a Reproducing Kernel Hilbert Space (RKHS) formulation of ordinary differential equations (ODEs). This method is particularly useful for reconstructing individual microbial trajectories and inferring dynamic signatures under irregular and sparse sampling.


## Table of Contents
1. [Overview](#Overview)
2. [Installation](#Installation)
3. [Quick Start (Simulated Data)](#Quick)

## Overview
What this package does

* Fits a non-parametric ODE in RKHS for each subject–taxon trajectory.  
* Handles irregular sampling, few time points, supports interpolation/extrapolation.
* Ships three simulated datasets and utilities to convert data into the required list-of-time-points format.

Why it matters
* In simulations and real cohorts, α-features separate distinct dynamic mechanisms and improve disease prediction compared with difference-based and spline-based approaches.

## Installation
```R
# From GitHub (recommend remotes)
install.packages("remotes")
remotes::install_github("ZhangXuan0429/r.ODE.RKHS")

# Local install (if you cloned the repo)
# setwd("/path/to/r.ODE.RKHS"); devtools::install()
```

## Quick Start (Simulated Data)
The package includes “wide” data frames: `simudata_linear_exp`, `simudata_predator_prey`, `simudata_periodic_random`.
Rows = samples; columns = t0..tK + Group.  
Convert to list of time-point matrices and fit.
```R
library(r.ODE.RKHS)

# 1) Load one dataset and keep a group (e.g., Linear)
data(simudata_linear_exp)
df_lin <- subset(simudata_linear_exp, Group == "Linear")

# 2) Build the list of time points (one matrix per column t*)
tp_cols   <- grep("^t\\d+$", names(df_lin), value = TRUE)
time_pts  <- as.numeric(sub("^t", "", tp_cols))
data_list <- lapply(tp_cols, function(tt) {
  m <- as.matrix(df_lin[, tt, drop = FALSE])
  rownames(m) <- rownames(df_lin)
  colnames(m) <- "value"   # one feature; pass this in `vertex`
  m
})

# 3) Fit ODE-RKHS (per subject–taxon)
fit_lin <- fit_oderkhs(
  data        = data_list,
  vertex      = "value",
  group_label = "Linear",
  sigma       = 1,
  lambda      = 1e-2,
  gamma_init  = 1,
  time_points = time_pts
)

# 4) Inspect outputs
head(fit_lin$alpha_df)   # α-features (per interval)
head(fit_lin$ode_result) # reconstructed trajectories
head(fit_lin$loss)       # losses and convergence flags
```

> Tip：For multi-taxon matrices at each time, keep all taxa as columns and pass a vector of column names to vertex.

