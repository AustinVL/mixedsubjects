#' Estimate ATEs in mixed-subjects design RCTs
#'
#' @param formula A formula of the form Y ~ D | S or Y ~ D | S1 - S0.
#' @param data A data.frame or `msd_data` containing outcome, treatment, and prediction columns.
#' @param estimator Estimator name.
#' @param crossfit Logical; whether to use cross-fitting when supported.
#' @param k Number of folds for cross-fitting (must be 2 for now).
#' @param folds Optional vector of fold assignments for labeled units.
#' @param se Standard error type.
#' @param n_boot Number of bootstrap replications (if bootstrap SE is requested).
#' @param bootstrap_parallel Logical; whether to run bootstrap in parallel.
#' @param bootstrap_seed Optional seed for bootstrap.
#' @param lambda_bounds Numeric vector of length 2 bounding tuning parameters.
#' @param tuning_fallback Strategy when tuning is unstable.
#' @param tuning_epsilon Threshold for near-zero denominators.
#' @param tuning_ridge Ridge factor used when `tuning_fallback = "ridge"`.
#' @param se_correction Whether to apply a delta-method correction for tuning.
#' @param fold_seed Optional seed for reproducible fold assignment.
#' @param alpha Significance level for confidence interval.
#' @param ... Reserved for future extensions.
#' @return An object of class `msd_estimate`.
#' @export
msd_estimate <- function(formula,
                         data,
                         estimator = c("dim", "greg", "ppi_pp", "dt", "dip", "dip_pp", "dt_dip"),
                         crossfit = TRUE,
                         k = 2,
                         folds = NULL,
                         se = c("analytic", "bootstrap"),
                         n_boot = 1000,
                         bootstrap_parallel = FALSE,
                         bootstrap_seed = NULL,
                         lambda_bounds = c(-5, 5),
                         tuning_fallback = c("clamp", "ridge", "zero"),
                         tuning_epsilon = 1e-8,
                         tuning_ridge = 1e-6,
                         se_correction = c("none", "delta"),
                         fold_seed = NULL,
                         alpha = 0.05,
                         ...) {
  estimator <- match.arg(estimator)
  se <- match.arg(se)
  tuning_fallback <- match.arg(tuning_fallback)
  se_correction <- match.arg(se_correction)
  if (inherits(data, "msd_data")) {
    if (missing(formula) || is.null(formula)) {
      formula <- data$formula
    }
    data <- data$data
  }
  parsed <- parse_msd_formula(formula, data)
  data_msd <- parsed$data
  needs_two_arm <- estimator %in% c("dip", "dip_pp", "dt_dip")
  if (!parsed$has_predictions && estimator != "dim") {
    stop("Predictions must be supplied for estimators other than dim.")
  }
  if (needs_two_arm && !parsed$has_two_arm) {
    stop("DiP estimators require two prediction columns (S1 - S0).")
  }

  if (estimator == "dim") {
    result <- dim_estimator(data_msd)
  } else if (estimator == "greg") {
    result <- greg_estimator(data_msd)
  } else if (estimator == "ppi_pp") {
    result <- tuned_estimator(
      data_msd,
      crossfit = crossfit,
      k = k,
      folds = folds,
      arm_specific = FALSE,
      fold_seed = fold_seed,
      lambda_bounds = lambda_bounds,
      tuning_fallback = tuning_fallback,
      tuning_epsilon = tuning_epsilon,
      tuning_ridge = tuning_ridge
    )
  } else if (estimator == "dt") {
    result <- tuned_estimator(
      data_msd,
      crossfit = crossfit,
      k = k,
      folds = folds,
      arm_specific = TRUE,
      fold_seed = fold_seed,
      lambda_bounds = lambda_bounds,
      tuning_fallback = tuning_fallback,
      tuning_epsilon = tuning_epsilon,
      tuning_ridge = tuning_ridge
    )
  } else if (estimator == "dip") {
    result <- dip_estimator(data_msd)
  } else if (estimator == "dip_pp") {
    result <- tuned_dip_estimator(
      data_msd,
      crossfit = crossfit,
      k = k,
      folds = folds,
      arm_specific = FALSE,
      fold_seed = fold_seed,
      lambda_bounds = lambda_bounds,
      tuning_fallback = tuning_fallback,
      tuning_epsilon = tuning_epsilon,
      tuning_ridge = tuning_ridge
    )
  } else {
    result <- tuned_dip_estimator(
      data_msd,
      crossfit = crossfit,
      k = k,
      folds = folds,
      arm_specific = TRUE,
      fold_seed = fold_seed,
      lambda_bounds = lambda_bounds,
      tuning_fallback = tuning_fallback,
      tuning_epsilon = tuning_epsilon,
      tuning_ridge = tuning_ridge
    )
  }

  se_info <- compute_se(
    result,
    data_msd,
    estimator,
    se,
    n_boot,
    bootstrap_parallel,
    bootstrap_seed,
    formula,
    se_correction
  )
  alpha <- as.numeric(alpha)
  crit <- stats::qnorm(1 - alpha / 2)
  conf_int <- result$estimate + c(-1, 1) * crit * se_info$std_error
  statistic <- result$estimate / se_info$std_error
  p_value <- 2 * (1 - stats::pnorm(abs(statistic)))

  output <- list(
    estimate = result$estimate,
    std.error = se_info$std_error,
    statistic = statistic,
    p.value = p_value,
    conf.int = conf_int,
    alpha = alpha,
    method = result$method,
    estimator = estimator,
    variance_decomp = se_info$variance_decomp,
    tuning = result$tuning,
    folds = result$folds,
    data_info = summarize_data(data_msd),
    call = match.call(),
    formula = formula
  )
  class(output) <- "msd_estimate"
  output
}

#' Select estimators based on analytic SE
#'
#' @param formula A formula of the form Y ~ D | S or Y ~ D | S1 - S0.
#' @param data A data.frame containing outcome, treatment, and prediction columns.
#' @param estimators Candidate estimators to compare.
#' @param se_type SE type (analytic only for now).
#' @param ... Additional arguments passed to `msd_estimate()`.
#' @return A data.frame sorted by standard error.
#' @export
msd_select <- function(formula,
                       data,
                       estimators = c("dim", "greg", "ppi_pp", "dt", "dip", "dip_pp", "dt_dip"),
                       se_type = c("analytic"),
                       ...) {
  se_type <- match.arg(se_type)
  results <- lapply(estimators, function(est) {
    fit <- msd_estimate(formula, data, estimator = est, se = se_type, ...)
    data.frame(
      estimator = est,
      estimate = fit$estimate,
      std.error = fit$std.error,
      method = fit$method,
      needs_two_arm_predictions = est %in% c("dip", "dip_pp", "dt_dip"),
      uses_crossfit = est %in% c("ppi_pp", "dt", "dip_pp", "dt_dip"),
      lambda_summary = summarize_lambda(fit$tuning$lambda_folds),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, results)
  out[order(out$std.error), , drop = FALSE]
}

#' @export
print.msd_estimate <- function(x, ...) {
  cat("Mixed-subjects ATE estimate\n")
  cat(sprintf("Estimator: %s\n", x$estimator))
  cat(sprintf("Estimate: %.4f\n", x$estimate))
  cat(sprintf("SE: %.4f\n", x$std.error))
  cat(sprintf("z: %.3f\n", x$statistic))
  cat(sprintf("p-value: %.4f\n", x$p.value))
  cat(sprintf("CI (%.0f%%): [%.4f, %.4f]\n", (1 - x$alpha) * 100, x$conf.int[1], x$conf.int[2]))
  if (!is.null(x$variance_decomp)) {
    cat("Variance decomposition:\n")
    for (name in names(x$variance_decomp)) {
      cat(sprintf("  %s: %.6f\n", name, x$variance_decomp[[name]]))
    }
  }
  invisible(x)
}

#' @export
summary.msd_estimate <- function(object, ...) {
  print(object)
  if (!is.null(object$tuning$lambda_folds)) {
    cat("\nTuning summary:\n")
    print(object$tuning$lambda_folds)
  }
  if (!is.null(object$folds)) {
    cat("\nFold assignment summary:\n")
    print(table(object$folds$fold, object$folds$arm))
  }
  invisible(object)
}

#' @export
tidy.msd_estimate <- function(x, ...) {
  data.frame(
    term = "ATE",
    estimate = x$estimate,
    std.error = x$std.error,
    statistic = x$statistic,
    p.value = x$p.value,
    conf.low = x$conf.int[1],
    conf.high = x$conf.int[2],
    method = x$method,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.msd_estimate <- function(x, ...) {
  info <- x$data_info
  data.frame(
    estimator = x$estimator,
    method = x$method,
    n_total = info$n_total,
    n_labeled = info$n_labeled,
    n_unlabeled = info$n_unlabeled,
    n1 = info$n1,
    n0 = info$n0,
    m1 = info$m1,
    m0 = info$m0,
    stringsAsFactors = FALSE
  )
}

#' @export
augment.msd_estimate <- function(x, data, ...) {
  parsed <- parse_msd_formula(x$formula, data)
  data_msd <- parsed$data
  fold <- rep(NA_integer_, nrow(data_msd))
  if (!is.null(x$folds)) {
    fold[x$folds$index] <- x$folds$fold
  }
  data.frame(
    data,
    fold = fold,
    residual = ifelse(data_msd$r == 1, data_msd$y - data_msd$s_assigned, NA_real_),
    stringsAsFactors = FALSE
  )
}

parse_msd_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("formula must be a valid formula.")
  }
  form_chr <- deparse(formula)
  split <- strsplit(form_chr, "\\|")[[1]]
  main_form <- stats::as.formula(trimws(split[1]))
  mf <- stats::model.frame(main_form, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  terms_obj <- stats::terms(main_form)
  term_labels <- attr(terms_obj, "term.labels")
  if (length(term_labels) < 1) {
    stop("Formula must include an outcome and treatment: Y ~ D | ...")
  }
  d_name <- term_labels[1]
  if (!d_name %in% colnames(mf)) {
    stop("Treatment variable not found in data: ", d_name)
  }
  d <- mf[[d_name]]
  if (is.logical(d)) {
    d <- as.integer(d)
  } else if (is.factor(d)) {
    if (nlevels(d) != 2) {
      stop("Treatment factor must have exactly two levels.")
    }
    d <- as.integer(d) - 1
  }
  if (!all(d %in% c(0, 1))) {
    stop("Treatment variable must be coded as 0/1 (or logical/two-level factor).")
  }
  if (length(split) == 1) {
    preds <- NULL
    pred_vars <- character()
    pred_expr <- NULL
  } else {
    pred_expr <- trimws(split[2])
    pred_form <- stats::as.formula(paste("~", pred_expr))
    preds <- stats::model.matrix(pred_form, data = data)
    preds <- preds[, colnames(preds) != "(Intercept)", drop = FALSE]
    pred_vars <- all.vars(pred_form)
  }

  if (is.null(preds)) {
    s0 <- s1 <- s_assigned <- rep(NA_real_, length(y))
    has_two_arm <- FALSE
    has_predictions <- FALSE
  } else if (length(pred_vars) >= 2 && grepl("-", pred_expr, fixed = TRUE)) {
    s1 <- data[[pred_vars[1]]]
    s0 <- data[[pred_vars[2]]]
    s_assigned <- ifelse(d == 1, s1, s0)
    has_two_arm <- TRUE
    has_predictions <- TRUE
  } else if (ncol(preds) == 1) {
    s_assigned <- preds[, 1]
    s0 <- ifelse(d == 0, s_assigned, NA_real_)
    s1 <- ifelse(d == 1, s_assigned, NA_real_)
    has_two_arm <- FALSE
    has_predictions <- TRUE
  } else if (ncol(preds) >= 2) {
    s0 <- preds[, 1]
    s1 <- preds[, 2]
    s_assigned <- ifelse(d == 1, s1, s0)
    has_two_arm <- TRUE
    has_predictions <- TRUE
  } else {
    stop("Prediction terms could not be parsed.")
  }

  r <- ifelse(is.na(y), 0, 1)
  data_msd <- data.frame(
    y = y,
    d = d,
    r = r,
    s0 = s0,
    s1 = s1,
    s_assigned = s_assigned,
    stringsAsFactors = FALSE
  )
  check_msd_inputs(data_msd, allow_missing_y = TRUE)
  list(data = data_msd, has_two_arm = has_two_arm, has_predictions = has_predictions)
}

check_msd_inputs <- function(data, allow_missing_y = FALSE) {
  required <- c("d", "r", "y", "s0", "s1", "s_assigned")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!all(data$d %in% c(0, 1))) {
    stop("Column 'd' must be coded as 0/1.")
  }
  if (!all(data$r %in% c(0, 1))) {
    stop("Column 'r' must be coded as 0/1.")
  }
  if (!allow_missing_y && any(is.na(data$y))) {
    stop("Outcome values must be provided for labeled units.")
  }
  invisible(TRUE)
}

dim_estimator <- function(data) {
  labeled <- data[data$r == 1, , drop = FALSE]
  estimate <- mean(labeled$y[labeled$d == 1]) - mean(labeled$y[labeled$d == 0])
  list(
    estimate = estimate,
    method = "DIM (labeled-only)",
    tuning = list(lambda = NULL, lambda_folds = NULL),
    folds = NULL
  )
}

greg_estimator <- function(data) {
  mu1 <- greg_arm_mean(data, arm = 1)
  mu0 <- greg_arm_mean(data, arm = 0)
  list(
    estimate = mu1 - mu0,
    method = "GREG (calibration)",
    tuning = list(lambda = 1, lambda_folds = NULL),
    folds = NULL
  )
}

greg_arm_mean <- function(data, arm) {
  labeled <- data[data$r == 1 & data$d == arm, , drop = FALSE]
  unlabeled <- data[data$r == 0 & data$d == arm, , drop = FALSE]
  y_mean <- mean(labeled$y)
  s_l <- mean(labeled$s_assigned)
  s_u <- mean(unlabeled$s_assigned)
  y_mean + (s_u - s_l)
}

tuned_estimator <- function(data,
                            crossfit,
                            k,
                            folds,
                            arm_specific,
                            fold_seed,
                            lambda_bounds,
                            tuning_fallback,
                            tuning_epsilon,
                            tuning_ridge) {
  if (!crossfit) {
    fold_info <- normalize_folds(data, rep(1L, sum(data$r == 1)), k = 1)
    k <- 1
  } else if (is.null(folds)) {
    fold_info <- normalize_folds(data, make_folds(data, k = k, seed = fold_seed), k = k)
  } else {
    fold_info <- normalize_folds(data, folds, k = k)
  }

  fold_vec <- fold_info$fold_vector
  fold_data <- fold_info$fold_data
  fold_levels <- sort(unique(fold_vec))
  labeled_idx <- fold_data$index
  mu1 <- mu0 <- numeric(length(fold_levels))
  lambda1 <- lambda0 <- numeric(length(fold_levels))

  for (i in seq_along(fold_levels)) {
    fold_id <- fold_levels[i]
    train_idx <- labeled_idx[fold_vec != fold_id]
    test_idx <- labeled_idx[fold_vec == fold_id]
    train_data <- data[train_idx, , drop = FALSE]
    test_data <- data[test_idx, , drop = FALSE]

    if (arm_specific) {
      lambda1[i] <- estimate_lambda_calibration(
        train_data,
        arm = 1,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
      lambda0[i] <- estimate_lambda_calibration(
        train_data,
        arm = 0,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
    } else {
      lambda <- estimate_lambda_pooled(
        train_data,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
      lambda1[i] <- lambda
      lambda0[i] <- lambda
    }

    mu1[i] <- arm_correction(test_data, data, arm = 1, lambda = lambda1[i])
    mu0[i] <- arm_correction(test_data, data, arm = 0, lambda = lambda0[i])
  }

  weights1 <- fold_weights(fold_vec[fold_data$arm == 1], fold_levels)
  weights0 <- fold_weights(fold_vec[fold_data$arm == 0], fold_levels)
  if (sum(weights1) == 0 || sum(weights0) == 0) {
    warning("Fold weights contain no labeled units in at least one arm; estimates may be unstable.", call. = FALSE)
  }

  mu1_hat <- sum(mu1 * as.numeric(weights1))
  mu0_hat <- sum(mu0 * as.numeric(weights0))

  tuning <- list(
    lambda = list(lambda1 = lambda1[length(lambda1)], lambda0 = lambda0[length(lambda0)]),
    lambda_folds = data.frame(
      fold = fold_levels,
      lambda1 = lambda1,
      lambda0 = lambda0,
      stringsAsFactors = FALSE
    )
  )

  list(
    estimate = mu1_hat - mu0_hat,
    method = if (arm_specific) "D-T (cross-fit, analytic SE)" else "PPI++ (cross-fit, analytic SE)",
    tuning = tuning,
    folds = fold_data
  )
}

dip_estimator <- function(data) {
  unlabeled <- data[data$r == 0, , drop = FALSE]
  labeled <- data[data$r == 1, , drop = FALSE]
  u_diff <- mean(unlabeled$s1 - unlabeled$s0)
  adj_t <- mean(labeled$y[labeled$d == 1] - labeled$s1[labeled$d == 1])
  adj_c <- mean(labeled$y[labeled$d == 0] - labeled$s0[labeled$d == 0])
  list(
    estimate = u_diff + adj_t - adj_c,
    method = "DiP (difference in predictions)",
    tuning = list(lambda = NULL, lambda_folds = NULL),
    folds = NULL
  )
}

tuned_dip_estimator <- function(data,
                                crossfit,
                                k,
                                folds,
                                arm_specific,
                                fold_seed,
                                lambda_bounds,
                                tuning_fallback,
                                tuning_epsilon,
                                tuning_ridge) {
  if (!crossfit) {
    fold_info <- normalize_folds(data, rep(1L, sum(data$r == 1)), k = 1)
    k <- 1
  } else if (is.null(folds)) {
    fold_info <- normalize_folds(data, make_folds(data, k = k, seed = fold_seed), k = k)
  } else {
    fold_info <- normalize_folds(data, folds, k = k)
  }

  fold_vec <- fold_info$fold_vector
  fold_data <- fold_info$fold_data
  fold_levels <- sort(unique(fold_vec))
  labeled_idx <- fold_data$index
  tau <- numeric(length(fold_levels))
  lambda1 <- lambda0 <- numeric(length(fold_levels))

  for (i in seq_along(fold_levels)) {
    fold_id <- fold_levels[i]
    train_idx <- labeled_idx[fold_vec != fold_id]
    test_idx <- labeled_idx[fold_vec == fold_id]
    train_data <- data[train_idx, , drop = FALSE]
    test_data <- data[test_idx, , drop = FALSE]

    if (arm_specific) {
      lambda1[i] <- estimate_lambda_dip(
        train_data,
        arm = 1,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
      lambda0[i] <- estimate_lambda_dip(
        train_data,
        arm = 0,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
    } else {
      lambda <- estimate_lambda_dip(
        train_data,
        arm = NULL,
        lambda_bounds = lambda_bounds,
        tuning_fallback = tuning_fallback,
        tuning_epsilon = tuning_epsilon,
        tuning_ridge = tuning_ridge
      )
      lambda1[i] <- lambda
      lambda0[i] <- lambda
    }

    tau[i] <- dip_fold_estimate(test_data, data, lambda1[i], lambda0[i])
  }

  weights <- fold_weights(fold_vec, fold_levels)
  if (sum(weights) == 0) {
    warning("Fold weights contain no labeled units; estimates may be unstable.", call. = FALSE)
  }
  estimate <- sum(tau * as.numeric(weights))

  tuning <- list(
    lambda = list(lambda1 = lambda1[length(lambda1)], lambda0 = lambda0[length(lambda0)]),
    lambda_folds = data.frame(
      fold = fold_levels,
      lambda1 = lambda1,
      lambda0 = lambda0,
      stringsAsFactors = FALSE
    )
  )

  list(
    estimate = estimate,
    method = if (arm_specific) "D-T DiP (cross-fit, analytic SE)" else "DiP++ (cross-fit, analytic SE)",
    tuning = tuning,
    folds = fold_data
  )
}

dip_fold_estimate <- function(labeled_fold, data, lambda1, lambda0) {
  unlabeled <- data[data$r == 0, , drop = FALSE]
  u_term <- mean(lambda1 * unlabeled$s1 - lambda0 * unlabeled$s0)
  adj_t <- mean(labeled_fold$y[labeled_fold$d == 1] - lambda1 * labeled_fold$s1[labeled_fold$d == 1])
  adj_c <- mean(labeled_fold$y[labeled_fold$d == 0] - lambda0 * labeled_fold$s0[labeled_fold$d == 0])
  u_term + adj_t - adj_c
}

arm_correction <- function(test_data, data, arm, lambda) {
  labeled <- test_data[test_data$d == arm, , drop = FALSE]
  unlabeled <- data[data$r == 0 & data$d == arm, , drop = FALSE]
  y_mean <- mean(labeled$y)
  s_l <- mean(labeled$s_assigned)
  s_u <- mean(unlabeled$s_assigned)
  y_mean + lambda * (s_u - s_l)
}

make_folds <- function(data, k = 2, seed = NULL) {
  if (k != 2) {
    stop("Only two-fold cross-fitting is supported.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  labeled_idx <- which(data$r == 1)
  fold <- rep(NA_integer_, length(labeled_idx))
  for (arm in c(0, 1)) {
    idx <- labeled_idx[data$d[labeled_idx] == arm]
    n_arm <- length(idx)
    fold_assign <- sample(rep(1:2, length.out = n_arm))
    fold[labeled_idx %in% idx] <- fold_assign
  }
  fold
}

estimate_lambda_calibration <- function(data,
                                        arm,
                                        lambda_bounds,
                                        tuning_fallback,
                                        tuning_epsilon,
                                        tuning_ridge) {
  labeled <- data[data$r == 1 & data$d == arm, , drop = FALSE]
  warn_small_n(labeled, arm = arm)
  s <- labeled$s_assigned
  y <- labeled$y
  m_d <- sum(data$r == 0 & data$d == arm)
  n_d <- nrow(labeled)
  var_s <- stats::var(s)
  cov_ys <- stats::cov(y, s)
  c_d <- var_s * (1 / m_d + 1 / n_d)
  d_d <- cov_ys / n_d
  lambda <- tune_ratio(d_d, c_d, var_s, lambda_bounds, tuning_fallback, tuning_epsilon, tuning_ridge)
  lambda
}

estimate_lambda_pooled <- function(data,
                                   lambda_bounds,
                                   tuning_fallback,
                                   tuning_epsilon,
                                   tuning_ridge) {
  labeled1 <- data[data$r == 1 & data$d == 1, , drop = FALSE]
  labeled0 <- data[data$r == 1 & data$d == 0, , drop = FALSE]
  if (nrow(labeled1) < 2 && nrow(labeled0) < 2) {
    return(0)
  }
  b1 <- lambda_component(data, arm = 1)$c_d
  b0 <- lambda_component(data, arm = 0)$c_d
  c1 <- lambda_component(data, arm = 1)$d_d
  c0 <- lambda_component(data, arm = 0)$d_d
  denom <- b1 + b0
  lambda <- tune_ratio(c1 + c0, denom, b1 + b0, lambda_bounds, tuning_fallback, tuning_epsilon, tuning_ridge)
  lambda
}

lambda_component <- function(data, arm) {
  labeled <- data[data$r == 1 & data$d == arm, , drop = FALSE]
  if (nrow(labeled) < 2) {
    return(list(c_d = 0, d_d = 0))
  }
  s <- labeled$s_assigned
  y <- labeled$y
  m_d <- sum(data$r == 0 & data$d == arm)
  n_d <- nrow(labeled)
  var_s <- stats::var(s)
  cov_ys <- stats::cov(y, s)
  c_d <- var_s * (1 / m_d + 1 / n_d)
  d_d <- cov_ys / n_d
  list(c_d = c_d, d_d = d_d)
}

estimate_lambda_dip <- function(data,
                                arm = NULL,
                                lambda_bounds,
                                tuning_fallback,
                                tuning_epsilon,
                                tuning_ridge) {
  labeled <- data[data$r == 1, , drop = FALSE]
  unlabeled <- data[data$r == 0, , drop = FALSE]
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  m <- nrow(unlabeled)
  if (m == 0 || n1 == 0 || n0 == 0) {
    return(0)
  }
  if (!is.null(arm)) {
    warn_small_n(labeled[labeled$d == arm, , drop = FALSE], arm = arm)
  }
  dip_var <- stats::var(unlabeled$s1 - unlabeled$s0) / m
  s1_var <- stats::var(labeled$s1[labeled$d == 1]) / n1
  s0_var <- stats::var(labeled$s0[labeled$d == 0]) / n0
  denom <- dip_var + s1_var + s0_var
  d1 <- stats::cov(labeled$y[labeled$d == 1], labeled$s1[labeled$d == 1]) / n1
  d0 <- stats::cov(labeled$y[labeled$d == 0], labeled$s0[labeled$d == 0]) / n0
  if (is.null(arm)) {
    numer <- d1 + d0
  } else if (arm == 1) {
    numer <- d1
  } else {
    numer <- d0
  }
  tune_ratio(numer, denom, denom, lambda_bounds, tuning_fallback, tuning_epsilon, tuning_ridge)
}

compute_se <- function(result,
                       data,
                       estimator,
                       se,
                       n_boot,
                       bootstrap_parallel,
                       bootstrap_seed,
                       formula,
                       se_correction) {
  if (se == "bootstrap") {
    std_error <- bootstrap_se(
      formula,
      data,
      estimator,
      n_boot,
      bootstrap_parallel,
      bootstrap_seed,
      result$folds
    )
    return(list(std_error = std_error, variance_decomp = NULL))
  }

  if (!is.null(result$folds) && estimator %in% c("ppi_pp", "dt", "dip_pp", "dt_dip")) {
    variance <- var_crossfit(data, estimator, result)
  } else {
    variance <- switch(
      estimator,
      dim = var_dim(data),
      greg = var_greg(data),
      dip = var_dip(data),
      ppi_pp = var_tuned(data, result$tuning),
      dt = var_tuned(data, result$tuning),
      dip_pp = var_dip_tuned(data, result$tuning),
      dt_dip = var_dip_tuned(data, result$tuning),
      var_dim(data)
    )
  }

  if (se_correction == "delta" && estimator %in% c("ppi_pp", "dt", "dip_pp", "dt_dip")) {
    variance$total <- variance$total + delta_lambda_variance(data, estimator, result, n_boot)
  }

  if (is.na(variance$total)) {
    warning("Analytic SE unavailable; falling back to bootstrap.", call. = FALSE)
    std_error <- bootstrap_se(
      formula,
      data,
      estimator,
      n_boot,
      bootstrap_parallel,
      bootstrap_seed,
      result$folds
    )
    return(list(std_error = std_error, variance_decomp = NULL))
  }

  list(
    std_error = sqrt(variance$total),
    variance_decomp = variance$decomp
  )
}

var_dim <- function(data) {
  labeled <- data[data$r == 1, , drop = FALSE]
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  var1 <- safe_div(stats::var(labeled$y[labeled$d == 1]), n1)
  var0 <- safe_div(stats::var(labeled$y[labeled$d == 0]), n0)
  total <- var1 + var0
  list(total = total, decomp = list(baseline = total))
}

var_greg <- function(data) {
  labeled <- data[data$r == 1, , drop = FALSE]
  m1 <- sum(data$r == 0 & data$d == 1)
  m0 <- sum(data$r == 0 & data$d == 0)
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  s1_var <- stats::var(labeled$s_assigned[labeled$d == 1])
  s0_var <- stats::var(labeled$s_assigned[labeled$d == 0])
  res1 <- labeled$y[labeled$d == 1] - labeled$s_assigned[labeled$d == 1]
  res0 <- labeled$y[labeled$d == 0] - labeled$s_assigned[labeled$d == 0]
  var1 <- safe_div(s1_var, m1) + safe_div(stats::var(res1), n1)
  var0 <- safe_div(s0_var, m0) + safe_div(stats::var(res0), n0)
  total <- var1 + var0
  list(total = total, decomp = list(baseline = total))
}

var_dip <- function(data) {
  labeled <- data[data$r == 1, , drop = FALSE]
  unlabeled <- data[data$r == 0, , drop = FALSE]
  m <- nrow(unlabeled)
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  dip_var <- safe_div(stats::var(unlabeled$s1 - unlabeled$s0), m)
  res1 <- labeled$y[labeled$d == 1] - labeled$s1[labeled$d == 1]
  res0 <- labeled$y[labeled$d == 0] - labeled$s0[labeled$d == 0]
  total <- dip_var + safe_div(stats::var(res1), n1) + safe_div(stats::var(res0), n0)
  list(total = total, decomp = list(unlabeled_variance = dip_var))
}

var_tuned <- function(data, tuning) {
  labeled <- data[data$r == 1, , drop = FALSE]
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  m1 <- sum(data$r == 0 & data$d == 1)
  m0 <- sum(data$r == 0 & data$d == 0)
  lambda1 <- tuning$lambda$lambda1
  lambda0 <- tuning$lambda$lambda0
  comps1 <- lambda_component(data, arm = 1)
  comps0 <- lambda_component(data, arm = 0)
  var1 <- safe_div(stats::var(labeled$y[labeled$d == 1]), n1) + lambda1^2 * comps1$c_d - 2 * lambda1 * comps1$d_d
  var0 <- safe_div(stats::var(labeled$y[labeled$d == 0]), n0) + lambda0^2 * comps0$c_d - 2 * lambda0 * comps0$d_d
  total <- var1 + var0
  list(total = total, decomp = list(prediction_rectifier = total - safe_div(stats::var(labeled$y[labeled$d == 1]), n1) - safe_div(stats::var(labeled$y[labeled$d == 0]), n0)))
}

var_dip_tuned <- function(data, tuning) {
  labeled <- data[data$r == 1, , drop = FALSE]
  unlabeled <- data[data$r == 0, , drop = FALSE]
  m <- nrow(unlabeled)
  n1 <- sum(labeled$d == 1)
  n0 <- sum(labeled$d == 0)
  lambda1 <- tuning$lambda$lambda1
  lambda0 <- tuning$lambda$lambda0
  dip_var <- safe_div(stats::var(lambda1 * unlabeled$s1 - lambda0 * unlabeled$s0), m)
  res1 <- labeled$y[labeled$d == 1] - lambda1 * labeled$s1[labeled$d == 1]
  res0 <- labeled$y[labeled$d == 0] - lambda0 * labeled$s0[labeled$d == 0]
  total <- dip_var + safe_div(stats::var(res1), n1) + safe_div(stats::var(res0), n0)
  list(total = total, decomp = list(unlabeled_variance = dip_var))
}

var_crossfit <- function(data, estimator, result) {
  folds <- result$folds
  fold_levels <- sort(unique(folds$fold))
  weights <- compute_fold_weights(folds, estimator)
  lambda1 <- result$tuning$lambda_folds$lambda1
  lambda0 <- result$tuning$lambda_folds$lambda0

  if (estimator %in% c("ppi_pp", "dt")) {
    var1 <- var_crossfit_arm(data, folds, fold_levels, lambda1, arm = 1)
    var0 <- var_crossfit_arm(data, folds, fold_levels, lambda0, arm = 0)
    total <- var1$total + var0$total
    decomp <- list(crossfit_covariance = var1$cov + var0$cov)
  } else {
    var <- var_crossfit_dip(data, folds, fold_levels, lambda1, lambda0, weights)
    total <- var$total
    decomp <- list(crossfit_covariance = var$cov)
  }

  list(total = total, decomp = decomp)
}

var_tuned_fold <- function(data, labeled_fold, lambda1, lambda0) {
  n1 <- sum(labeled_fold$d == 1)
  n0 <- sum(labeled_fold$d == 0)
  comps1 <- lambda_component_fold(labeled_fold, data, arm = 1)
  comps0 <- lambda_component_fold(labeled_fold, data, arm = 0)
  var1 <- safe_div(stats::var(labeled_fold$y[labeled_fold$d == 1]), n1) + lambda1^2 * comps1$c_d - 2 * lambda1 * comps1$d_d
  var0 <- safe_div(stats::var(labeled_fold$y[labeled_fold$d == 0]), n0) + lambda0^2 * comps0$c_d - 2 * lambda0 * comps0$d_d
  var1 + var0
}

var_dip_tuned_fold <- function(data, labeled_fold, lambda1, lambda0) {
  unlabeled <- data[data$r == 0, , drop = FALSE]
  m <- nrow(unlabeled)
  n1 <- sum(labeled_fold$d == 1)
  n0 <- sum(labeled_fold$d == 0)
  dip_var <- safe_div(stats::var(lambda1 * unlabeled$s1 - lambda0 * unlabeled$s0), m)
  res1 <- labeled_fold$y[labeled_fold$d == 1] - lambda1 * labeled_fold$s1[labeled_fold$d == 1]
  res0 <- labeled_fold$y[labeled_fold$d == 0] - lambda0 * labeled_fold$s0[labeled_fold$d == 0]
  dip_var + safe_div(stats::var(res1), n1) + safe_div(stats::var(res0), n0)
}

lambda_component_fold <- function(labeled_fold, data, arm) {
  labeled <- labeled_fold[labeled_fold$d == arm, , drop = FALSE]
  if (nrow(labeled) < 2) {
    return(list(c_d = 0, d_d = 0))
  }
  s <- labeled$s_assigned
  y <- labeled$y
  m_d <- sum(data$r == 0 & data$d == arm)
  n_d <- nrow(labeled)
  var_s <- stats::var(s)
  cov_ys <- stats::cov(y, s)
  c_d <- var_s * (1 / m_d + 1 / n_d)
  d_d <- cov_ys / n_d
  list(c_d = c_d, d_d = d_d)
}

var_crossfit_arm <- function(data, folds, fold_levels, lambda, arm) {
  weights <- compute_arm_weights(folds, arm)
  variances <- numeric(length(fold_levels))
  for (i in seq_along(fold_levels)) {
    fold_id <- fold_levels[i]
    labeled_fold <- data[folds$index[folds$fold == fold_id & folds$arm == arm], , drop = FALSE]
    comps <- lambda_component_fold(labeled_fold, data, arm = arm)
    n_d <- nrow(labeled_fold)
    var_y <- safe_div(stats::var(labeled_fold$y), n_d)
    variances[i] <- var_y + lambda[i]^2 * comps$c_d - 2 * lambda[i] * comps$d_d
  }
  cov_total <- 0
  if (length(fold_levels) > 1) {
    unlabeled <- data[data$r == 0 & data$d == arm, , drop = FALSE]
    m_d <- nrow(unlabeled)
    var_s <- safe_div(stats::var(unlabeled$s_assigned), m_d)
    for (i in seq_along(fold_levels)) {
      for (j in seq_along(fold_levels)) {
        if (i < j) {
          cov_total <- cov_total + 2 * weights[i] * weights[j] * lambda[i] * lambda[j] * var_s
        }
      }
    }
  }
  total <- sum((weights^2) * variances) + cov_total
  list(total = total, cov = cov_total)
}

var_crossfit_dip <- function(data, folds, fold_levels, lambda1, lambda0, weights) {
  variances <- numeric(length(fold_levels))
  for (i in seq_along(fold_levels)) {
    fold_id <- fold_levels[i]
    labeled_fold <- data[folds$index[folds$fold == fold_id], , drop = FALSE]
    variances[i] <- var_dip_tuned_fold(data, labeled_fold, lambda1[i], lambda0[i])
  }
  cov_total <- 0
  if (length(fold_levels) > 1) {
    unlabeled <- data[data$r == 0, , drop = FALSE]
    m <- nrow(unlabeled)
    s1 <- unlabeled$s1
    s0 <- unlabeled$s0
    for (i in seq_along(fold_levels)) {
      for (j in seq_along(fold_levels)) {
        if (i < j) {
          combo_i <- lambda1[i] * s1 - lambda0[i] * s0
          combo_j <- lambda1[j] * s1 - lambda0[j] * s0
          cov_ij <- safe_div(stats::cov(combo_i, combo_j), m)
          cov_total <- cov_total + 2 * weights[i] * weights[j] * cov_ij
        }
      }
    }
  }
  total <- sum((weights^2) * variances) + cov_total
  list(total = total, cov = cov_total)
}

compute_fold_weights <- function(folds, estimator) {
  fold_levels <- sort(unique(folds$fold))
  if (estimator %in% c("ppi_pp", "dt")) {
    weights1 <- compute_arm_weights(folds, arm = 1, fold_levels = fold_levels)
    weights0 <- compute_arm_weights(folds, arm = 0, fold_levels = fold_levels)
    weights <- (weights1 + weights0) / 2
  } else {
    weights <- fold_weights(folds$fold, fold_levels)
  }
  weights[as.character(fold_levels)]
}

compute_arm_weights <- function(folds, arm, fold_levels = NULL) {
  if (is.null(fold_levels)) {
    fold_levels <- sort(unique(folds$fold))
  }
  weights <- fold_weights(folds$fold[folds$arm == arm], fold_levels)
  weights[as.character(fold_levels)]
}

delta_lambda_variance <- function(data, estimator, result, n_boot) {
  folds <- result$folds
  if (is.null(folds)) {
    return(0)
  }
  fold_levels <- sort(unique(folds$fold))
  weights1 <- compute_arm_weights(folds, arm = 1, fold_levels = fold_levels)
  weights0 <- compute_arm_weights(folds, arm = 0, fold_levels = fold_levels)
  weights_all <- compute_fold_weights(folds, estimator)
  delta_var <- 0

  for (i in seq_along(fold_levels)) {
    fold_id <- fold_levels[i]
    train_idx <- folds$index[folds$fold != fold_id]
    test_idx <- folds$index[folds$fold == fold_id]
    train_data <- data[train_idx, , drop = FALSE]
    test_data <- data[test_idx, , drop = FALSE]
    if (estimator %in% c("ppi_pp", "dip_pp")) {
      lambda_var <- bootstrap_lambda(data, train_idx, estimator, n_boot = n_boot)
      if (estimator == "ppi_pp") {
        derivs <- lambda_derivative(test_data, data, estimator, arm_specific = TRUE)
        weight <- weights1[i] * derivs$d_tau_d_lambda1 + weights0[i] * derivs$d_tau_d_lambda0
        delta_var <- delta_var + weight^2 * lambda_var
      } else {
        derivs <- lambda_derivative(test_data, data, estimator)
        weight <- weights_all[i] * derivs$d_tau_d_lambda1
        delta_var <- delta_var + weight^2 * lambda_var
      }
    } else {
      lambda_var1 <- bootstrap_lambda(data, train_idx, estimator, arm = 1, n_boot = n_boot)
      lambda_var0 <- bootstrap_lambda(data, train_idx, estimator, arm = 0, n_boot = n_boot)
      derivs <- lambda_derivative(test_data, data, estimator, arm_specific = TRUE)
      delta_var <- delta_var + (weights1[i]^2) * derivs$d_tau_d_lambda1^2 * lambda_var1 +
        (weights0[i]^2) * derivs$d_tau_d_lambda0^2 * lambda_var0
    }
  }
  delta_var
}

lambda_derivative <- function(labeled_fold, data, estimator, arm_specific = FALSE) {
  unlabeled <- data[data$r == 0, , drop = FALSE]
  out <- list(d_tau_d_lambda1 = 0, d_tau_d_lambda0 = 0)
  if (estimator %in% c("ppi_pp", "dt")) {
    deriv1 <- mean(unlabeled$s1[unlabeled$d == 1]) - mean(labeled_fold$s_assigned[labeled_fold$d == 1])
    deriv0 <- mean(unlabeled$s0[unlabeled$d == 0]) - mean(labeled_fold$s_assigned[labeled_fold$d == 0])
    out$d_tau_d_lambda1 <- deriv1
    out$d_tau_d_lambda0 <- -deriv0
    return(out)
  }

  u_term1 <- mean(unlabeled$s1)
  u_term0 <- mean(unlabeled$s0)
  s1_mean <- if (any(labeled_fold$d == 1)) mean(labeled_fold$s1[labeled_fold$d == 1]) else 0
  s0_mean <- if (any(labeled_fold$d == 0)) mean(labeled_fold$s0[labeled_fold$d == 0]) else 0
  out$d_tau_d_lambda1 <- u_term1 - s1_mean
  out$d_tau_d_lambda0 <- -(u_term0 - s0_mean)
  out
}

bootstrap_lambda <- function(data, labeled_idx, estimator, arm = NULL, n_boot = 200) {
  labeled <- data[labeled_idx, , drop = FALSE]
  if (nrow(labeled) < 2) {
    return(0)
  }
  estimates <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    if (is.null(arm)) {
      resample <- sample(seq_len(nrow(labeled)), replace = TRUE)
      boot_labeled <- labeled[resample, , drop = FALSE]
    } else {
      arm_idx <- which(labeled$d == arm)
      if (length(arm_idx) == 0) {
        return(0)
      }
      resample <- sample(arm_idx, replace = TRUE)
      boot_labeled <- labeled[resample, , drop = FALSE]
    }
    boot_data <- data
    boot_data[labeled_idx, ] <- boot_labeled
    if (estimator == "ppi_pp") {
      estimates[b] <- estimate_lambda_pooled(
        boot_data,
        lambda_bounds = c(-Inf, Inf),
        tuning_fallback = "zero",
        tuning_epsilon = 1e-8,
        tuning_ridge = 1e-6
      )
    } else if (estimator == "dt") {
      estimates[b] <- estimate_lambda_calibration(
        boot_data,
        arm = arm,
        lambda_bounds = c(-Inf, Inf),
        tuning_fallback = "zero",
        tuning_epsilon = 1e-8,
        tuning_ridge = 1e-6
      )
    } else if (estimator == "dip_pp") {
      estimates[b] <- estimate_lambda_dip(
        boot_data,
        arm = NULL,
        lambda_bounds = c(-Inf, Inf),
        tuning_fallback = "zero",
        tuning_epsilon = 1e-8,
        tuning_ridge = 1e-6
      )
    } else {
      estimates[b] <- estimate_lambda_dip(
        boot_data,
        arm = arm,
        lambda_bounds = c(-Inf, Inf),
        tuning_fallback = "zero",
        tuning_epsilon = 1e-8,
        tuning_ridge = 1e-6
      )
    }
  }
  stats::var(estimates, na.rm = TRUE)
}

bootstrap_se <- function(formula, data, estimator, n_boot, bootstrap_parallel, bootstrap_seed, folds) {
  if (!is.null(bootstrap_seed)) {
    set.seed(bootstrap_seed)
  }
  labeled_idx <- which(data$r == 1)
  unlabeled_idx <- which(data$r == 0)

  bootstrap_once <- function(iter) {
    if (is.null(folds)) {
      resample_labeled <- sample(labeled_idx, replace = TRUE)
    } else {
      resample_labeled <- c()
      for (arm in c(0, 1)) {
        for (fold_id in unique(folds$fold)) {
          fold_idx <- folds$index[folds$arm == arm & folds$fold == fold_id]
          if (length(fold_idx) > 0) {
            resample_labeled <- c(resample_labeled, sample(fold_idx, replace = TRUE))
          }
        }
      }
    }
    resample_unlabeled <- c()
    for (arm in c(0, 1)) {
      arm_idx <- unlabeled_idx[data$d[unlabeled_idx] == arm]
      if (length(arm_idx) > 0) {
        resample_unlabeled <- c(resample_unlabeled, sample(arm_idx, replace = TRUE))
      }
    }
    boot_data <- data[c(resample_labeled, resample_unlabeled), , drop = FALSE]
    boot_folds <- NULL
    if (!is.null(folds)) {
      boot_folds <- folds$fold[match(resample_labeled, folds$index)]
    }
    boot_fit <- msd_estimate(
      formula,
      boot_data,
      estimator = estimator,
      se = "analytic",
      folds = boot_folds,
      crossfit = !is.null(boot_folds)
    )
    boot_fit$estimate
  }

  if (bootstrap_parallel) {
    if (.Platform$OS.type == "windows") {
      warning("bootstrap_parallel not supported on Windows; running serially.", call. = FALSE)
      boot_est <- vapply(seq_len(n_boot), bootstrap_once, numeric(1))
    } else {
      boot_est <- unlist(parallel::mclapply(seq_len(n_boot), bootstrap_once, mc.cores = parallel::detectCores(), mc.set.seed = TRUE))
    }
  } else {
    boot_est <- vapply(seq_len(n_boot), bootstrap_once, numeric(1))
  }
  stats::sd(boot_est)
}

summarize_data <- function(data) {
  n_total <- nrow(data)
  n_labeled <- sum(data$r == 1)
  n_unlabeled <- sum(data$r == 0)
  n1 <- sum(data$r == 1 & data$d == 1)
  n0 <- sum(data$r == 1 & data$d == 0)
  m1 <- sum(data$r == 0 & data$d == 1)
  m0 <- sum(data$r == 0 & data$d == 0)
  list(
    n_total = n_total,
    n_labeled = n_labeled,
    n_unlabeled = n_unlabeled,
    n1 = n1,
    n0 = n0,
    m1 = m1,
    m0 = m0
  )
}

summarize_lambda <- function(lambda_folds) {
  if (is.null(lambda_folds) || nrow(lambda_folds) == 0) {
    return(NA_character_)
  }
  paste0(
    "lambda1=[", paste(sprintf("%.3f", lambda_folds$lambda1), collapse = ", "), "], ",
    "lambda0=[", paste(sprintf("%.3f", lambda_folds$lambda0), collapse = ", "), "]"
  )
}

safe_div <- function(num, denom) {
  if (is.na(denom) || denom == 0) {
    return(NA_real_)
  }
  num / denom
}

`%||%` <- function(x, y) {
  if (!is.null(x)) {
    return(x)
  }
  y
}

fold_weights <- function(fold_assignments, fold_levels) {
  counts <- tabulate(match(fold_assignments, fold_levels), nbins = length(fold_levels))
  total <- sum(counts)
  if (total == 0) {
    return(setNames(rep(0, length(fold_levels)), fold_levels))
  }
  stats::setNames(counts / total, fold_levels)
}

normalize_folds <- function(data, folds, k) {
  labeled_idx <- which(data$r == 1)
  if (is.data.frame(folds)) {
    required <- c("index", "fold", "arm")
    missing <- setdiff(required, names(folds))
    if (length(missing) > 0) {
      stop("folds data.frame must contain columns: index, fold, arm.")
    }
    fold_data <- folds
    fold_data <- fold_data[order(fold_data$index), , drop = FALSE]
    fold_vec <- fold_data$fold
  } else {
    if (length(folds) != length(labeled_idx)) {
      stop("folds vector must have length equal to number of labeled units.")
    }
    fold_vec <- as.integer(folds)
    fold_data <- data.frame(index = labeled_idx, fold = fold_vec, arm = data$d[labeled_idx])
  }
  if (any(is.na(fold_vec))) {
    stop("folds must not contain NA.")
  }
  if (any(!fold_vec %in% seq_len(k))) {
    stop("folds must be coded within 1:k.")
  }
  list(fold_vector = fold_vec, fold_data = fold_data)
}

tune_ratio <- function(numer,
                       denom,
                       scale,
                       lambda_bounds,
                       tuning_fallback,
                       tuning_epsilon,
                       tuning_ridge) {
  if (is.na(denom) || abs(denom) < tuning_epsilon || denom <= 0) {
    if (tuning_fallback == "zero") {
      return(0)
    }
    if (tuning_fallback == "ridge") {
      denom <- denom + tuning_ridge * abs(scale)
    } else {
      return(clamp_lambda(0, lambda_bounds))
    }
  }
  lambda <- numer / denom
  clamp_lambda(lambda, lambda_bounds)
}

clamp_lambda <- function(lambda, lambda_bounds) {
  lower <- lambda_bounds[1]
  upper <- lambda_bounds[2]
  max(min(lambda, upper), lower)
}

warn_small_n <- function(labeled, arm) {
  if (nrow(labeled) < 10) {
    warning(sprintf("Small labeled sample in arm %s; tuning may be unstable.", arm), call. = FALSE)
  }
}
#' Create an MSD data object with auto-detected columns
#'
#' @param data A data.frame with outcome, treatment, and prediction columns.
#' @param outcome Name of outcome column (default: detect Y/y/outcome).
#' @param treatment Name of treatment column (default: detect D/d/treat/treatment).
#' @param s0 Name of control prediction column.
#' @param s1 Name of treated prediction column.
#' @param s_assigned Name of assigned-arm prediction column.
#' @return An object of class `msd_data`.
#' @export
msd_data <- function(data,
                     outcome = NULL,
                     treatment = NULL,
                     s0 = NULL,
                     s1 = NULL,
                     s_assigned = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  names_lower <- tolower(names(data))
  detect_name <- function(candidates) {
    idx <- match(candidates, names_lower)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) {
      return(NULL)
    }
    names(data)[idx[1]]
  }

  outcome <- outcome %||% detect_name(c("y", "outcome"))
  treatment <- treatment %||% detect_name(c("d", "treat", "treatment"))
  s0 <- s0 %||% detect_name(c("s0"))
  s1 <- s1 %||% detect_name(c("s1"))
  s_assigned <- s_assigned %||% detect_name(c("s"))

  if (is.null(outcome) || is.null(treatment)) {
    stop("Unable to detect outcome or treatment columns; please supply them explicitly.")
  }

  if (!is.null(s0) && !is.null(s1)) {
    formula <- stats::as.formula(sprintf("%s ~ %s | %s - %s", outcome, treatment, s1, s0))
  } else if (!is.null(s_assigned)) {
    formula <- stats::as.formula(sprintf("%s ~ %s | %s", outcome, treatment, s_assigned))
  } else {
    formula <- stats::as.formula(sprintf("%s ~ %s", outcome, treatment))
  }

  structure(
    list(
      data = data,
      outcome = outcome,
      treatment = treatment,
      s0 = s0,
      s1 = s1,
      s_assigned = s_assigned,
      formula = formula
    ),
    class = "msd_data"
  )
}

#' @export
print.msd_data <- function(x, ...) {
  cat("mixedsubjects data object\n")
  cat(sprintf("Outcome: %s\n", x$outcome))
  cat(sprintf("Treatment: %s\n", x$treatment))
  cat(sprintf("Prediction (assigned): %s\n", ifelse(is.null(x$s_assigned), "none", x$s_assigned)))
  cat(sprintf("Prediction (S0): %s\n", ifelse(is.null(x$s0), "none", x$s0)))
  cat(sprintf("Prediction (S1): %s\n", ifelse(is.null(x$s1), "none", x$s1)))
  cat(sprintf("Default formula: %s\n", deparse(x$formula)))
  info <- summarize_data(parse_msd_formula(x$formula, x$data)$data)
  cat(sprintf("N total: %d, labeled: %d, unlabeled: %d\n", info$n_total, info$n_labeled, info$n_unlabeled))
  invisible(x)
}

#' @export
summary.msd_data <- function(object, ...) {
  print(object)
  parsed <- parse_msd_formula(object$formula, object$data)
  data_msd <- parsed$data
  labeled <- data_msd[data_msd$r == 1, , drop = FALSE]
  unlabeled <- data_msd[data_msd$r == 0, , drop = FALSE]
  cat("\nSample sizes:\n")
  cat(sprintf("  Labeled: %d (D=1: %d, D=0: %d)\n", nrow(labeled), sum(labeled$d == 1), sum(labeled$d == 0)))
  cat(sprintf("  Unlabeled: %d (D=1: %d, D=0: %d)\n", nrow(unlabeled), sum(unlabeled$d == 1), sum(unlabeled$d == 0)))
  cat("\nOutcome summary (labeled):\n")
  print(summary(labeled$y))
  if (parsed$has_predictions) {
    cat("\nPrediction summaries (labeled):\n")
    if (!all(is.na(labeled$s_assigned))) {
      cat("Assigned prediction:\n")
      print(summary(labeled$s_assigned))
    }
    if (parsed$has_two_arm) {
      cat("S0 prediction:\n")
      print(summary(labeled$s0))
      cat("S1 prediction:\n")
      print(summary(labeled$s1))
    }
    cat("\nCorrelations (labeled):\n")
    corr_assigned <- stats::cor(labeled$y, labeled$s_assigned, use = "complete.obs")
    cat(sprintf("Corr(Y, S_assigned): %.3f\n", corr_assigned))
    if (parsed$has_two_arm) {
      corr_t <- stats::cor(labeled$y[labeled$d == 1], labeled$s1[labeled$d == 1], use = "complete.obs")
      corr_c <- stats::cor(labeled$y[labeled$d == 0], labeled$s0[labeled$d == 0], use = "complete.obs")
      cat(sprintf("Corr(Y, S1 | D=1): %.3f\n", corr_t))
      cat(sprintf("Corr(Y, S0 | D=0): %.3f\n", corr_c))
    }
  } else {
    cat("\nNo predictions detected; only DIM is available.\n")
  }
  invisible(object)
}
