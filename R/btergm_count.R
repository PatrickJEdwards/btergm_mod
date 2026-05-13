#' Estimate a temporal ERGM for count-valued networks by bootstrapped count MPLE
#'
#' Estimate a temporal exponential-family random graph model for count-valued
#' networks using temporal bootstrapping and count-valued maximum
#' pseudolikelihood estimation.
#'
#' @param formula Formula for the temporal ERGM. The left-hand side should be
#'   a list of count-valued network objects or matrices in chronological order.
#' @param R Number of temporal bootstrap replications.
#' @param response Character string giving the name of the edge attribute that
#'   stores the count outcome. For example, \code{"count"}.
#' @param offset Logical. If \code{TRUE}, use structural-zero offset handling
#'   from \code{tergmprepare()}. For the first implementation, this should
#'   usually remain \code{FALSE}; user-supplied formula offsets will be handled
#'   separately.
#' @param returndata Logical. Return prepared data rather than estimating the
#'   model?
#' @param parallel Character. One of \code{"no"}, \code{"multicore"}, or
#'   \code{"snow"}.
#' @param ncpus Number of CPU cores for temporal bootstrap parallelization.
#' @param cl Optional cluster object for \code{parallel = "snow"}.
#' @param control.ergm Optional ergm control object. Reserved for compatibility
#'   with \code{btergm()}; not used directly by the count MPLE optimizer yet.
#' @param usefastglm Logical. Included for API compatibility with
#'   \code{btergm()}. Count MPLE does not use a binomial GLM, so this is ignored
#'   for now. The binary btergm() function eventually estimates a binomial 
#'   GLM-style pseudolikelihood, where fastglm makes sense. The Huang-Butts 
#'   count MPLE code instead builds a count pseudolikelihood and optimizes it 
#'   using trust() or optim(), so there is no immediate GLM step to replace with 
#'   fastglm.
#' @param verbose Logical. Print preprocessing and estimation messages?
#'
#' @param reference Reference measure for the count ERGM. One of
#'   \code{"poisson"}, \code{"geometric"}, \code{"uniform"}, or
#'   \code{"binomial"}.
#' @param max.count Maximum possible count value. Use \code{Inf} for unbounded
#'   count support.
#' @param max.count.safety Safety multiplier used when approximating the count
#'   support.
#' @param max.count.edgewise Logical. Use edgewise count support truncation?
#' @param count.samples Maximum number of support values to use per edge
#'   variable.
#' @param must.count Integer. Always include support values from 0 through this
#'   value.
#' @param sample.size Number of edge variables to sample for count MPLE. Use
#'   \code{Inf} to use all dyads.
#' @param sample.method Edge-variable sampling method. One of
#'   \code{"weighted"} or \code{"random"}.
#' @param weight.type Weighting scheme when \code{sample.method = "weighted"}.
#'   One of \code{"TNT"} or \code{"flatval"}.
#' @param regularization Regularization method. One of \code{"none"},
#'   \code{"L2"}, \code{"L1"}, \code{"pseudoHuber"}, or \code{"L1pH"}.
#' @param regularization.param Numeric vector of regularization parameters.
#' @param coef.init Optional starting values for the optimizer.
#' @param optim.method Optimizer. Use \code{"trust"} for trust-region
#'   optimization, or another method accepted by \code{optim()}.
#' @param WtSumAsSampSiz Logical. Whether inverse weights should sum to the
#'   sampled edge-variable count.
#' @param seed Optional random seed.
#' @param mple.cores Number of cores used inside each count-MPLE fit for
#'   edge-variable change-score calculations. Keep this at 1 when using outer
#'   bootstrap parallelization to avoid nested parallelism.
#' @param prep.cores Number of cores used during count-MPLE preparation,
#'   including count support construction, reference-measure ratios, and
#'   change-score computation. This happens before the temporal bootstrap.
#' @param prep.parallel Optional preparation-stage parallelization mode.
#'   Reserved for preparation workflows that parallelize across time periods.
#' @param prep.cl Optional cluster object for preparation-stage parallelization.
#' @param ... Additional arguments passed later to the count-MPLE optimizer.
#'
#' @return A fitted \code{btergm_count} object. In the current scaffold,
#'   estimation is not yet implemented.
#'
#' @importFrom ergm control.ergm
#' @export
btergm_count <- function(formula,
                         R = 500,
                         response = "tie_count",
                         offset = FALSE,
                         returndata = FALSE,
                         parallel = c("no", "multicore", "snow"),
                         ncpus = 1, # used for outer bootstrap
                         cl = NULL,
                         control.ergm = NULL,
                         usefastglm = FALSE,
                         verbose = TRUE,
                         reference = c("poisson", "geometric", "uniform", "binomial"), # Count-MPLE arguments from Huang and Butts code
                         max.count = Inf,
                         max.count.safety = 4,
                         max.count.edgewise = TRUE,
                         count.samples = Inf,
                         must.count = 5,
                         sample.size = Inf,
                         sample.method = c("weighted", "random"),
                         weight.type = c("TNT", "flatval"),
                         regularization = c("none", "L2", "L1", "pseudoHuber", "L1pH"),
                         regularization.param = c(1, 1e-4),
                         coef.init = NULL,
                         optim.method = "trust",
                         WtSumAsSampSiz = TRUE,
                         seed = NULL,
                         prep.cores = 1, # used before bootstrap for ergmCntPrep_btergm()
                         prep.parallel = c("auto", "fork", "psock", "serial"),
                         prep.cl = NULL,
                         mple.cores = 1, # used inside MPLE calls, kept low for memory safety
                         ...) {
  
  # ------------------------------------------------------------
  # 1. Match and validate arguments
  # ------------------------------------------------------------
  
  parallel <- match.arg(parallel)
  prep.parallel <- match.arg(prep.parallel)
  reference <- match.arg(reference)
  sample.method <- match.arg(sample.method)
  weight.type <- match.arg(weight.type)
  regularization <- match.arg(regularization)
  
  
  if (missing(formula) || !inherits(formula, "formula")) {stop("'formula' must be provided as a formula object.")}
  
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 0) {stop("'R' must be a non-negative numeric value of length 1.")}
  R <- as.integer(R)
  
  if (!is.character(response) || length(response) != 1L || is.na(response)) {stop("'response' must be a character string naming the count edge attribute.")}
  
  if (!is.logical(offset) || length(offset) != 1L || is.na(offset)) {stop("'offset' must be TRUE or FALSE.")}
  
  if (!is.logical(returndata) || length(returndata) != 1L || is.na(returndata)) {stop("'returndata' must be TRUE or FALSE.")}
  
  if (!is.numeric(ncpus) || length(ncpus) != 1L || is.na(ncpus) || ncpus < 1) {stop("'ncpus' must be a positive numeric value of length 1.")}
  ncpus <- as.integer(ncpus)
  
  if (!is.numeric(prep.cores) || length(prep.cores) != 1L || is.na(prep.cores) || prep.cores < 1) {stop("'prep.cores' must be a positive numeric value of length 1.")}
  prep.cores <- as.integer(prep.cores)
  
  if (!is.numeric(mple.cores) || length(mple.cores) != 1L || is.na(mple.cores) || mple.cores < 1) {stop("'mple.cores' must be a positive numeric value of length 1.")}
  mple.cores <- as.integer(mple.cores)
  
  if (is.null(control.ergm)) {control.ergm <- ergm::control.ergm()}
  
  if (isTRUE(usefastglm)) {warning("'usefastglm = TRUE' is ignored for btergm_count() at this stage. ", "The Huang-Butts count MPLE routine uses direct optimization, not a ", "binomial GLM.")}
  
  
  
  # ------------------------------------------------------------
  # 2. Prepare temporal data using the same outer structure as btergm()
  # ------------------------------------------------------------
  
  #l <- tergmprepare(formula = formula, offset = offset, verbose = verbose)
  l <- btergm_count_prepare(
    formula = formula,
    response = response,
    offset = offset,
    blockdiag = FALSE,
    verbose = verbose,
    temporal_covariates = "count"
  )
  
  # Make prepared networks and covariates visible to the prepared formula.
  # This mirrors the local assignment pattern in btergm().
  for (i in seq_along(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]], envir = environment())
  }
  
  assign("offsmat", l$offsmat, envir = environment())
  
  form <- stats::as.formula(l$form, env = environment())
  
  
  
  
  # ------------------------------------------------------------
  # 3. Basic temporal checks and messages
  # ------------------------------------------------------------
  
  if (l$time.steps == 1L) {warning(paste("The confidence intervals and standard errors are meaningful only", "with temporal replication. Only one time step was provided."))}
  
  if (verbose && !returndata) {
    if (parallel == "no") {parallel.msg <- "on a single computing core"
    } else if (parallel == "multicore") {parallel.msg <- paste("using multicore forking on", ncpus, "cores")
    } else if (parallel == "snow") {parallel.msg <- paste("using parallel processing on", ncpus, "cores")}
    
    if (offset) {offset.msg <- "with offset matrix and "
    } else {offset.msg <- "with "}
    
    message("\nStarting count pseudolikelihood estimation ", offset.msg, R, " temporal bootstrap replications ", parallel.msg, "...")
  } else if (verbose && returndata) {message("\nReturning prepared count-TERGM data structure.")}
  
  
  # ------------------------------------------------------------
  # 4. Store count-MPLE controls for the next implementation step
  # ------------------------------------------------------------
  
  count_control <- list(
    response = response,
    reference = reference,
    max.count = max.count,
    max.count.safety = max.count.safety,
    max.count.edgewise = max.count.edgewise,
    count.samples = count.samples,
    must.count = must.count,
    sample.size = sample.size,
    sample.method = sample.method,
    weight.type = weight.type,
    regularization = regularization,
    regularization.param = regularization.param,
    coef.init = coef.init,
    optim.method = optim.method,
    WtSumAsSampSiz = WtSumAsSampSiz,
    seed = seed,
    prep.cores = prep.cores,
    prep.parallel = prep.parallel,
    mple.cores = mple.cores
  )
  
  if (returndata) {
    return(
      structure(
        list(
          call = match.call(),
          formula = formula,
          formula2 = form,
          prepared = l,
          count_control = count_control,
          parallel = parallel,
          ncpus = ncpus,
          cl = cl,
          control.ergm = control.ergm,
          dots = list(...)
        ),
        class = "btergm_count_prepared"
      )
    )
  }
  
  
  
  # ------------------------------------------------------------
  # 5. Precompute count-MPLE components by time period
  # ------------------------------------------------------------
  
  # Important:
  # prep.cores is used ONLY for ergmCntPrep_btergm(), which computes
  # count support values, reference-measure ratios, and change scores.
  # This happens before the temporal bootstrap begins.
  #
  # mple.cores is reserved for count-MPLE calls inside the full-sample
  # fit and bootstrap workers. In the current design, those calls receive
  # precomputed prep objects, so mple.cores should usually remain 1.
  
  if (verbose) {message("\nPrecomputing count-MPLE components by time period using ", prep.cores, " core(s) for ergmCntPrep_btergm().")}
  
  assign("count_offset", l$count_offset, envir = environment())
  
  period_preps <- btergm_count_precompute_period_preps(
    form = form,
    response = response,
    count_offset = l$count_offset,
    time_steps = l$time.steps,
    reference = reference,
    max.count = max.count,
    max.count.safety = max.count.safety,
    max.count.edgewise = max.count.edgewise,
    count.samples = count.samples,
    sample.size = sample.size,
    sample.method = sample.method,
    weight.type = weight.type,
    seed = seed,
    must.count = must.count,
    WtSumAsSampSiz = WtSumAsSampSiz,
    prep.cores = prep.cores,
    prep.parallel = prep.parallel,
    prep.cl = prep.cl,
    env = environment(),
    verbose = verbose
  )
  
  
  
  
  
  # ------------------------------------------------------------
  # 6. Estimate original model using all prepared time periods
  # ------------------------------------------------------------
  
  if (verbose) {message("\nEstimating count-MPLE model on the original time sequence.")}
  
  fit0 <- btergm_count_fit_indices(
    indices = seq_len(l$time.steps),
    form = form,
    response = response,
    period_preps = period_preps,
    reference = reference,
    regularization = regularization,
    regularization.param = regularization.param,
    coef.init = coef.init,
    optim.method = optim.method,
    WtSumAsSampSiz = WtSumAsSampSiz,
    mple.cores = mple.cores,                      # Full-sample MPLE call. Usually keep this at 1 for memory safety.
    estimate.cov = TRUE,                          # Full-sample fit can compute covariance for diagnostics.
    env = environment(),
    verbose = verbose,
    ...
  )
  ## OLD CODE (1):
  #theta0 <- stats::coef(fit0)
  ## NEW CODE (2):
  theta0 <- btergm_count_extract_coef(fit0)
  
  
  clean_coef_names <- btergm_count_clean_coef_names(
    parnam = names(theta0),
    prepared = l
  )
  
  names(theta0) <- clean_coef_names
  
  # Keep the internal full-sample count-MPLE fit aligned with the outer object.
  fit0$coef <- theta0
  fit0$coefficients <- theta0
  
  if (!is.null(fit0$cov)) {
    dimnames(fit0$cov) <- list(clean_coef_names, clean_coef_names)
  }
  
  if (!is.null(fit0$pll.hessian)) {
    dimnames(fit0$pll.hessian) <- list(clean_coef_names, clean_coef_names)
  }
  
  
  if (verbose) {message("\nOriginal count-MPLE estimate completed.")}
  
  # ------------------------------------------------------------
  # 7. Temporal bootstrap using boot::boot()
  # ------------------------------------------------------------
  
  if (R > 0L) {
    if (verbose) {message("\nStarting temporal bootstrap for count-MPLE estimates.")}
    
    boot_cl <- cl
    own_boot_cluster <- FALSE
    
    if (parallel == "snow") {
      if (is.null(boot_cl)) {
        if (verbose) {
          message("\nStarting snow cluster with ", ncpus, " worker(s) for temporal bootstrap.")
        }
        
        boot_cl <- parallel::makeCluster(ncpus, type = "PSOCK")
        own_boot_cluster <- TRUE
        
        on.exit({
          if (own_boot_cluster && !is.null(boot_cl)) {
            parallel::stopCluster(boot_cl)
          }
        }, add = TRUE)
      }
      
      # Load required packages on workers.
      parallel::clusterEvalQ(
        boot_cl,
        {
          pkgs <- c(
            "ergm",
            "network",
            "sna",
            "statnet.common",
            "sampling",
            "trust",
            "MASS",
            "btergm"
          )
          
          ok <- vapply(
            pkgs,
            requireNamespace,
            quietly = TRUE,
            FUN.VALUE = logical(1)
          )
          
          if (!all(ok)) {
            stop(
              "The following packages could not be loaded on the snow workers: ",
              paste(pkgs[!ok], collapse = ", ")
            )
          }
          
          NULL
        }
      )
      
      # Export internal helper functions needed by the bootstrap statistic.
      helper_env <- environment(btergm_count_boot_stat)
      
      parallel::clusterExport(
        boot_cl,
        varlist = c(
          "btergm_count_boot_stat",
          "btergm_count_fit_indices",
          "btergm_count_combine_preps",
          "btergm_count_extract_coef",
          "ergmCntMPLE_btergm",
          "ergmCntNLPL",
          "ergmCntNLPLDeriv",
          "%||%"
        ),
        envir = helper_env
      )
    }
    
    boot_obj <- boot::boot(
      data = seq_len(l$time.steps),
      statistic = btergm_count_boot_stat,
      R = R,
      sim = "ordinary",
      stype = "i",
      parallel = parallel,
      ncpus = ncpus,
      cl = boot_cl,
      form = form,
      response = response,
      period_preps = period_preps,
      reference = reference,
      regularization = regularization,
      regularization.param = regularization.param,
      coef.init = theta0,
      optim.method = optim.method,
      WtSumAsSampSiz = WtSumAsSampSiz,
      mple.cores = mple.cores,
      estimate.cov = FALSE,
      env = environment(),
      verbose = FALSE,
      ...
    )
    
    # Make sure boot$t0 equals the full-sample estimate we just computed.
    # boot::boot() computes t0 internally too, but this avoids tiny differences
    # caused by different starting values.
    boot_obj$t0 <- theta0
    
  } else {
    boot_obj <- structure(
      list(
        t0 = theta0,
        t = matrix(numeric(0), nrow = 0L, ncol = length(theta0)),
        R = 0L,
        data = seq_len(l$time.steps),
        seed = .Random.seed,
        statistic = btergm_count_boot_stat,
        sim = "ordinary",
        call = match.call()
      ),
      class = "boot"
    )
    
    colnames(boot_obj$t) <- names(theta0)
  }
  
  if (!is.null(boot_obj$t) && ncol(boot_obj$t) == length(theta0)) {
    colnames(boot_obj$t) <- names(theta0)
  }
  
  # ------------------------------------------------------------
  # 8. Return fitted btergm_count object
  # ------------------------------------------------------------
  
  out <- structure(
    list(
      call = match.call(),
      coef = theta0,
      coefficients = theta0,
      boot = boot_obj,
      R = R,
      nobs = sum(vapply(period_preps, function(x) length(x$iw), numeric(1))),
      time.steps = l$time.steps,
      formula = formula,
      formula2 = form,
      response = response,
      fit = fit0,
      period_preps = period_preps,
      count_control = count_control,
      data = l,
      auto.adjust = l$auto.adjust,
      offset = offset,
      directed = l$directed,
      bipartite = l$bipartite,
      nvertices = l$nvertices
    ),
    class = "btergm_count"
  )
  
  out
}
  
  