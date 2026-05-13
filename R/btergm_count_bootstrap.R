



















# =============================================================================
# Bootstrap helpers for btergm_count()
# =============================================================================


# Combine several ERGMCntPrep objects into one stacked pseudo-likelihood object.
#
# Each selected time period contributes a set of conditional edge-value
# likelihood terms. Bootstrapping with replacement is handled by allowing the
# same period's ERGMCntPrep object to appear more than once in prep_list.
#
# @keywords internal
btergm_count_combine_preps <- function(prep_list) {
  if (length(prep_list) == 0L) {stop("'prep_list' must contain at least one ERGMCntPrep object.")}
  
  ok <- vapply(prep_list, inherits, logical(1), what = "ERGMCntPrep")
  if (!all(ok)) {stop("All elements of 'prep_list' must inherit from class 'ERGMCntPrep'.")}
  
  out <- list()
  
  out$rmr <- unlist(lapply(prep_list, `[[`, "rmr"), recursive = FALSE)
  out$cs <- unlist(lapply(prep_list, `[[`, "cs"), recursive = FALSE)
  out$ycwt <- unlist(lapply(prep_list, `[[`, "ycwt"), recursive = FALSE)
  out$iw <- unlist(lapply(prep_list, `[[`, "iw"), use.names = FALSE)
  
  # Optional diagnostics. These are not used by the optimizer, but they are
  # useful for checking whether offsets and sampling worked as intended.
  out$y <- unlist(lapply(prep_list, function(x) x$y %||% numeric(0)), use.names = FALSE)
  out$snd <- unlist(lapply(prep_list, function(x) x$snd %||% numeric(0)), use.names = FALSE)
  out$rec <- unlist(lapply(prep_list, function(x) x$rec %||% numeric(0)), use.names = FALSE)
  
  out$nev <- sum(vapply(prep_list, function(x) x$nev %||% length(x$iw), numeric(1)))
  out$nev_total <- sum(vapply(prep_list, function(x) x$nev_total %||% NA_real_, numeric(1)), na.rm = TRUE)
  
  out$count_offset_used <- any(vapply(
    prep_list,
    function(x) isTRUE(x$count_offset_used),
    logical(1)
  ))
  
  out$period_nev <- vapply(prep_list, function(x) x$nev %||% length(x$iw), numeric(1))
  out$period_nev_total <- vapply(prep_list, function(x) x$nev_total %||% NA_real_, numeric(1))
  
  class(out) <- "ERGMCntPrep"
  out
}


# Small infix helper for defaults.
#
# @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# Precompute one ERGMCntPrep object for each prepared time period.
#
# This is usually much faster than recomputing change scores during every
# bootstrap replicate.
#
# @keywords internal
btergm_count_precompute_period_preps <- function(form,
                                                 response,
                                                 count_offset,
                                                 time_steps,
                                                 reference,
                                                 max.count,
                                                 max.count.safety,
                                                 max.count.edgewise,
                                                 count.samples,
                                                 sample.size,
                                                 sample.method,
                                                 weight.type,
                                                 seed,
                                                 must.count,
                                                 WtSumAsSampSiz,
                                                 prep.cores,
                                                 env,
                                                 verbose = TRUE) {
  preps <- vector("list", time_steps)
  
  fit_env <- environment(form) # Make sure the formula has an environment.
  
  if (is.null(fit_env)) {
    fit_env <- env
    environment(form) <- fit_env
  }
  
  for (tt in seq_len(time_steps)) {
    if (verbose) {message("Preparing count-MPLE data for time step ", tt, " of ", time_steps, ".")}
    
    # Make the current time index visible wherever the formula may look.
    assign("i", tt, envir = env)
    assign("i", tt, envir = fit_env)
    
    preps[[tt]] <- ergmCntPrep_btergm(
      formula = form,
      response = response,
      reference = reference,
      max.count = max.count,
      max.count.safety = max.count.safety,
      max.count.edgewise = max.count.edgewise,
      count.samples = count.samples,
      sample.size = sample.size,
      sample.method = sample.method,
      weight.type = weight.type,
      cores = prep.cores, # prep.cores is used only here, before bootstrapping.
      seed = if (is.null(seed)) NULL else seed + tt,
      must.count = must.count,
      count_offset = count_offset[[tt]],
      WtSumAsSampSiz = WtSumAsSampSiz
    )
  }
  
  preps
}


# Fit one count-MPLE model to a selected set of time periods.
#
# This is the count-valued analogue of stacking MPLE contributions across
# time periods in binary btergm().
#
# @keywords internal
btergm_count_fit_indices <- function(indices,
                                     form,
                                     response,
                                     period_preps,
                                     reference,
                                     regularization,
                                     regularization.param,
                                     coef.init,
                                     optim.method,
                                     WtSumAsSampSiz,
                                     mple.cores = 1,
                                     estimate.cov = FALSE,
                                     env,
                                     verbose = FALSE,
                                     ...) {
  indices <- as.integer(indices)
  
  if (length(indices) == 0L) {stop("'indices' must contain at least one time step.")}
  
  # Set i to a valid time step so ergmCntMPLE_btergm() can parse the formula
  # and recover parameter names. The actual pseudo-likelihood data come from
  # the combined prep object below.
  
  ## OLD CODE (1):
  #assign("i", indices[[1]], envir = env)
  
  ## NEW CODE (2):
  fit_env <- environment(form)
  assign("i", indices[[1]], envir = fit_env)
  assign("i", indices[[1]], envir = env)
  
  environment(form) <- fit_env
  
  combined_prep <- btergm_count_combine_preps(period_preps[indices])
  
  fit <- ergmCntMPLE_btergm(
    formula = form,
    response = response,
    reference = reference,
    regularization = regularization,
    regularization.param = regularization.param,
    coef.init = coef.init,
    optim.method = optim.method,
    prep = combined_prep,
    WtSumAsSampSiz = WtSumAsSampSiz,
    cores = mple.cores,
    estimate.cov = estimate.cov,
    verbose = verbose,
    ...
  )
  
  fit
}


# Statistic function passed to boot::boot().
#
# @keywords internal
btergm_count_boot_stat <- function(data,
                                   indices,
                                   form,
                                   response,
                                   period_preps,
                                   reference,
                                   regularization,
                                   regularization.param,
                                   coef.init,
                                   optim.method,
                                   WtSumAsSampSiz,
                                   mple.cores = 1,
                                   estimate.cov = FALSE,
                                   env,
                                   verbose = FALSE,
                                   ...) {
  fit <- btergm_count_fit_indices(
    indices = indices,
    form = form,
    response = response,
    period_preps = period_preps,
    reference = reference,
    regularization = regularization,
    regularization.param = regularization.param,
    coef.init = coef.init,
    optim.method = optim.method,
    WtSumAsSampSiz = WtSumAsSampSiz,
    mple.cores = mple.cores,
    estimate.cov = estimate.cov,
    env = env,
    verbose = verbose,
    ...
  )
  
  ## OLD CODE (1):
  #stats::coef(fit)
  ## NEW CODE (2):
  btergm_count_extract_coef(fit)
}





btergm_count_extract_coef <- function(fit) {
  if (!is.null(fit$coef)) {
    return(fit$coef)
  }
  
  if (!is.null(fit$coefficients)) {
    return(fit$coefficients)
  }
  
  out <- try(stats::coef(fit), silent = TRUE)
  
  if (!inherits(out, "try-error") && !is.null(out)) {
    return(out)
  }
  
  stop(
    "Could not extract coefficients from the count-MPLE fit. ",
    "Check whether ergmCntMPLE_btergm() stores coefficients in fit$coef ",
    "or fit$coefficients."
  )
}

#' @export
coef.ergmCntMPLE <- function(object, ...) {
  if (!is.null(object$coef)) {
    object$coef
  } else {
    object$coefficients
  }
}