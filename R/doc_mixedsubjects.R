#' mixedsubjects: Estimation and design tools for MSD-RCTs
#'
#' The mixedsubjects package provides estimators and design helpers for
#' mixed-subjects randomized controlled trials that combine a small labeled
#' (human) sample with inexpensive model predictions.
#'
#' @section Getting started:
#' Use the formula interface `Y ~ D | S` or `Y ~ D | S0 + S1` with
#' `msd_estimate()` to compute ATEs. If you only have predictions for assigned
#' arms, use the calibration family (DIM, GREG, PPI++/D–T). If you have
#' predictions for both arms, DiP estimators can further reduce variance.
#'
#' @section Key functions:
#' - `msd_estimate()`: estimate ATEs with analytic or bootstrap SEs.
#' - `msd_select()`: compare estimators by estimated SE.
#'
#' @docType package
#' @name mixedsubjects
NULL
