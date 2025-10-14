#' @title Fit ODE-RKHS trajectories and extract alpha-coefficients
#'
#' @description
#' Fits ODE-RKHS per subject–taxon on panel data organized as a list of time points.
#' For each taxon and subject, it alternates an Euler trajectory update and a
#' kernel-ridge update of the state-to-rate mapping coefficients (alpha), following
#' Eqs. (1)–(9). Returns alpha features, reconstructed trajectories, and loss summaries.
#'
#' @details
#' Alternating scheme: explicit Euler trajectory update (Eq. 4), finite-difference
#' targets (Eq. 5), kernel ridge update of coefficients (Eqs. 6–7), with geometric
#' schedule for the ODE penalty and termination by relative field change (Eq. 8).
#' The composite loss follows Eq. (9). The first alpha is obtained via kernel ridge
#' from finite differences of the observed trajectory; no manual alpha initialization
#' is required.
#'
#' @param data  A list of length \code{n_time}; each element is a numeric matrix or
#'   data frame of dimension \code{n_samples x n_taxa} for one time point, with row
#'   names = sample IDs and column names = taxa/features.
#' @param vertex Character vector of taxa/features (column names) to fit.
#' @param group_label Character scalar; a label written into outputs for bookkeeping
#'   (e.g., "Linear", "Exponential").
#' @param sigma Positive scalar; RBF kernel bandwidth.
#' @param lambda Positive scalar; base regularization strength.
#' @param gamma_init Positive scalar; initial ODE penalty weight.
#' @param time_points Numeric vector of length \code{n_time} with strictly increasing
#'   observation times.
#'
#' @return A list with three data frames:
#' \describe{
#'   \item{\code{alpha_df}}{Per-interval coefficients and derivative fits with columns:
#'         \code{TimePair}, \code{x_prev}, \code{x_current}, \code{dx_dt}, \code{Ypred},
#'         \code{Alpha}, \code{Sample}, \code{Species}, \code{Group}.}
#'   \item{\code{ode_result}}{Reconstructed trajectories with columns:
#'         \code{Time}, \code{Orig}, \code{Pred}, \code{Sample}, \code{Species}, \code{Group}.}
#'   \item{\code{loss}}{Per subject–taxon MSEs and total loss with columns:
#'         \code{mse_ypred}, \code{mse_ode}, \code{TotalLoss}, \code{Converged},
#'         \code{Sample}, \code{Species}, \code{Group}.}
#' }
#'
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' # Example using simulated data (see datasets in the package)
#' data(simudata_linear_exp)
#' df_lin <- subset(simudata_linear_exp, Group == "Linear")
#' tp_cols <- grep("^t\\d+$", names(df_lin), value = TRUE)
#' time_pts <- as.numeric(sub("^t", "", tp_cols))
#' data_list <- lapply(tp_cols, function(tt) {
#'   m <- as.matrix(df_lin[, tt, drop = FALSE])
#'   rownames(m) <- rownames(df_lin)
#'   colnames(m) <- "value"
#'   m
#' })
#' fit <- fit_oderkhs(
#'   data = data_list, vertex = "value", group_label = "Linear",
#'   sigma = 1, lambda = 1e-2, gamma_init = 1, time_points = time_pts
#' )
#' head(fit$alpha_df)
#'
fit_oderkhs <- function(data, vertex, group_label,
                        sigma, lambda, gamma_init, time_points) {

  alpha_list <- list()
  ode_results_list <- list()
  loss_list <- list()

  # Validate time settings
  if (length(time_points) != length(data)) {
    stop("Number of time points does not match the length of data list.")
  }
  n_time <- length(time_points)
  if (n_time < 2) stop("At least two time points are required.")
  h <- diff(time_points)

  # Alternating optimization parameters
  rho <- 1.5       # penalty growth factor
  max_iter <- 500  # maximum iterations
  tol <- 1e-3      # convergence threshold

  # Loop over all samples and taxa
  for (sample in rownames(data[[1]])) {
    for (species in vertex) {
      cat("Fitting sample:", sample, "species:", species, "group:", group_label, "\n")

      # ---- 1. Extract true trajectory ----
      z_true <- sapply(data, function(d) d[sample, species])

      # ---- 2. Initialization ----
      z_current <- z_true                 # initial trajectory
      f0 <- diff(z_true) / h              # initial derivative estimate
      f_prev <- f0
      gamma <- gamma_init
      converged <- FALSE
      final_total_loss <- NA

      # ---- 3. Alternating optimization ----
      for (iter in 1:max_iter) {

        # (3.1) Update trajectory given current f (Euler method)
        z_new <- numeric(n_time)
        z_new[1] <- z_current[1]
        for (i in 2:n_time) {
          z_new[i] <- z_new[i - 1] + h[i - 1] * f_prev[i - 1]
        }

        # (3.2) Compute finite differences and training pairs
        dx_dt <- diff(z_new) / h
        X_train <- cbind(
          z_prev = z_new[1:(n_time - 1)],
          z_current = z_new[2:n_time]
        )
        y_train <- dx_dt

        # (3.3) Kernel ridge regression update
        model <- try(krr(
          X_train = X_train,
          y_train = y_train,
          X_test = X_train,
          kernel_func = rbf_kernel,
          sigma = sigma,
          lambda = lambda,
          gamma = gamma
        ), silent = TRUE)

        if (inherits(model, "try-error")) next

        # (3.4) Compute total loss
        loss_data <- mean((z_new - z_true)^2)
        loss_ode  <- gamma * mean((z_new[-1] - z_new[-n_time] - h * f_prev)^2)
        K_train <- rbf_kernel(X_train, X_train, sigma)
        k <- nrow(X_train)
        h_mean <- mean(h)
        lambda_effective <- (lambda * k) / (gamma * h_mean^2)
        alpha0 <- solve(K_train + lambda_effective * diag(k), f0)
        reg_loss <- lambda * t(model$alpha - alpha0) %*% K_train %*% (model$alpha - alpha0)
        current_total_loss <- loss_data + loss_ode + reg_loss

        # (3.5) Update vector field
        f_new <- model$y_pred

        # (3.6) Check convergence
        delta_f <- sqrt(mean((f_new - f_prev)^2)) / sqrt(mean(f_prev^2))
        if (delta_f < tol) {
          converged <- TRUE
          final_total_loss <- current_total_loss
          cat("Sample", sample, "species", species,
              "converged at iteration", iter,
              "total loss:", round(final_total_loss, 4), "\n")
          break
        }

        # (3.7) Update for next iteration
        z_current <- z_new
        f_prev <- f_new
        gamma <- gamma * rho
        if (iter == max_iter) final_total_loss <- current_total_loss
      }

      # ---- 4. Store results ----
      time_pairs <- paste0("p", 1:(n_time - 1))
      alpha_df <- data.frame(
        TimePair = time_pairs,
        x_prev = X_train[, 1],
        x_current = X_train[, 2],
        dx_dt = y_train,
        Ypred = model$y_pred,
        Alpha = model$alpha,
        Sample = sample,
        Species = species,
        Group = group_label
      )
      alpha_list[[paste(sample, species)]] <- alpha_df

      ode_df <- data.frame(
        Time = time_points,
        Orig = z_true,
        Pred = z_current,
        Sample = sample,
        Species = species,
        Group = group_label
      )
      ode_results_list[[paste(sample, species)]] <- ode_df

      loss_list[[paste(sample, species)]] <- data.frame(
        mse_ypred = mean((f_new - y_train)^2),
        mse_ode = mean((z_new - z_true)^2),
        TotalLoss = final_total_loss,
        Converged = converged,
        Sample = sample,
        Species = species,
        Group = group_label
      )
    }
  }

  # ---- 5. Combine outputs ----
  return(list(
    alpha_df = dplyr::bind_rows(alpha_list),
    ode_result = dplyr::bind_rows(ode_results_list),
    loss = dplyr::bind_rows(loss_list)
  ))
}
