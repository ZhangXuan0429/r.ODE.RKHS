#' Simulated longitudinal dataset: linear vs exponential growth
#'
#' A wide data frame where rows are samples and columns are time points (t0, t1, ...),
#' plus a `Group` column indicating the dynamic model ("Linear" or "Exponential").
#'
#' @format A data frame with samples in rows and time points in columns.
#' @usage data(simudata_linear_exp)
#' @examples
#' data(simudata_linear_exp)
#' head(simudata_linear_exp)
#' @docType data
#' @keywords datasets
#' @name simudata_linear_exp
NULL

#' Simulated longitudinal dataset: predator–prey dynamics
#'
#' @format Same structure as \code{simudata_linear_exp}, but groups are
#'   "Predator" and "Prey".
#' @usage data(simudata_predator_prey)
#' @examples
#' data(simudata_predator_prey)
#' head(simudata_predator_prey)
#' @docType data
#' @keywords datasets
#' @name simudata_predator_prey
NULL

#' Simulated longitudinal dataset: periodic vs random fluctuations
#'
#' @format Same structure as \code{simudata_linear_exp}, but groups are
#'   "Periodic" and "Random".
#' @usage data(simudata_periodic_random)
#' @examples
#' data(simudata_periodic_random)
#' head(simudata_periodic_random)
#' @docType data
#' @keywords datasets
#' @name simudata_periodic_random
NULL
