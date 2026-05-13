# =============================================================================
# Count-valued ERGM MPLE helpers for btergm_count()
# =============================================================================
#
# This file adapts the count-valued ERGM MPLE implementation from Huang and
# Butts for use inside btergm_count().
#
# Main modifications relative to the original replication script:
#
#   1. The public-facing function is btergm_count(), not ergmCntMPLE_btergm().
#   2. ergmCntMPLE_btergm() and ergmCntPrep_btergm() are internal helpers.
#   3. User-supplied offset(edgecov(...)) terms are removed upstream by
#      btergm_count_prepare().
#   4. Dyads marked by count_offset matrices are excluded before count-MPLE
#      edge-variable sampling.
#   5. The C++ objective functions are retained as the optimization backend.
#
# Convention for count_offset:
#
#   count_offset[i, j] == 0     means dyad i -> j is active.
#   count_offset[i, j] != 0     means dyad i -> j is inactive.
#   is.na(count_offset[i, j])   means dyad i -> j is inactive.
#
# This distinction matters because inactive dyads should not enter the
# conditional pseudo-likelihood at all. They are not merely dyads with a low
# expected count.
#
# =============================================================================
# Original Huang-Butts count-MPLE notes, lightly adapted
# =============================================================================
#
# Maximum Pseudo-likelihood Estimation for General Count ERGMs
#
# The ergm package does not directly support MPLE for count ERGMs. It instead
# relies on contrastive divergence for valued models. CD can be slow or
# unreliable in some settings, while full MLE can be too expensive for large
# graphs with high-valued counts. In that regime, MPLE may be useful.
#
# For a general count ERGM, the conditional edge value probability is:
#
#   Pr(Y_i = y | Y_i^c)
#     =
#   h(y, Y_i^c) exp(theta * t(y, Y_i^c))
#   ------------------------------------------------------------
#   sum_{k = 0}^Inf h(k, Y_i^c) exp(theta * t(k, Y_i^c))
#
# Equivalently:
#
#   Pr(Y_i = y | Y_i^c)
#     =
#   1 /
#   sum_{k = 0}^Inf
#     [h(k, Y_i^c) / h(y, Y_i^c)]
#     exp(theta * [t(k, Y_i^c) - t(y, Y_i^c)])
#
# The pseudo-likelihood is the product of these conditional edge-value
# probabilities across edge variables.
#
# Useful reference-measure ratios:
#
#   Uniform/Geometric:
#     h(k, Y_i^c) / h(y, Y_i^c) = 1
#
#   Binomial:
#     h(k, Y_i^c) / h(y, Y_i^c)
#       = k!(m - k)! / [y!(m - y)!]
#
#   Poisson:
#     h(k, Y_i^c) / h(y, Y_i^c) = y! / k!
#
# Implementation idea:
#
#   1. Precompute reference-measure ratios.
#   2. Precompute change scores for each sampled edge variable and each
#      candidate count value in that edge variable's local support.
#   3. Optimize the resulting count pseudo-likelihood using trust-region
#      optimization or optim().
#
# Memory warning:
#
#   The precomputed object can be very large. Memory use scales roughly with:
#
#     number of edge variables
#       x number of model statistics
#       x number of candidate count values per edge variable
#
#   The count support can be controlled with max.count, max.count.safety,
#   max.count.edgewise, must.count, and count.samples.
#
# Edge-variable subsampling:
#
#   The count pseudo-likelihood can be approximated by sampling edge variables.
#   Weighted sampling can oversample rare values, especially nonzero or high
#   count dyads. Random sampling samples active dyads uniformly.
#
# Important modification for btergm_count():
#
#   In the original standalone count-MPLE code, all dyads returned by
#   network.dyadcount() and as.sociomatrix(nw, response) are potential edge
#   variables. In btergm_count(), this is not correct when the user supplies
#   offset(edgecov(receiver_ineligible_mat_list)) or similar terms. Those dyads
#   must be removed from the sampling frame before any support construction or
#   change-score calculation occurs.
#
#
# Notes in this section are adapted from the Huang and Butts count-MPLE
# replication script, with modifications for btergm_count().
#
# =============================================================================




# Count-valued ERGM MPLE helpers for btergm_count
#
# Internal helper functions for fitting count-valued ERGMs by maximum
# pseudolikelihood inside btergm_count(). These functions adapt the count-MPLE
# implementation from Huang and Butts for use with temporally prepared
# btergm-style network lists.
#
# The main package-specific modification is that dyads marked as structurally
# inactive by count_offset matrices are excluded before edge-variable sampling,
# count-support construction, and change-score precomputation.
#
# These functions are not intended to be called directly by end users.
#
# Reference:
# Huang, Peng, and Carter T. Butts. 2024. "Parameter Estimation Procedures for
# Exponential-Family Random Graph Models on Count-Valued Networks: A Comparative
# Simulation Study." Social Networks 76: 51-67.








#' Fit a count-valued ERGM by maximum pseudolikelihood
#'
#' Internal count-valued ERGM MPLE fitting function used by
#' \code{btergm_count()}. This function wraps the Huang-Butts count-MPLE
#' optimizer and adds support for dyad-level exclusion through a
#' \code{count_offset} matrix.
#'
#' @details
#' This function estimates a count-valued ERGM by maximizing the product of
#' conditional edge-value probabilities. For each sampled edge variable, the
#' method approximates the local normalizing factor by evaluating candidate
#' count values over a finite support. The preparation step precomputes the
#' reference-measure ratios and change scores needed by the optimizer.
#'
#' The computational burden can be high because the preparation object scales
#' with the number of sampled edge variables, the number of model statistics,
#' and the number of candidate count values evaluated for each edge variable.
#' The arguments \code{max.count}, \code{max.count.safety},
#' \code{max.count.edgewise}, \code{must.count}, and \code{count.samples}
#' control the count support approximation.
#'
#' The argument \code{sample.size} controls edge-variable subsampling. If
#' \code{sample.size} is smaller than the number of active dyads, the
#' pseudo-likelihood is approximated using a sample of edge variables. The
#' \code{sample.method} and \code{weight.type} arguments control how this
#' sample is drawn.
#'
#' For \code{btergm_count()}, \code{count_offset} is the crucial addition. Dyads
#' marked as inactive are excluded before edge-variable sampling and before
#' change-score precomputation. This is how user-supplied terms such as
#' \code{offset(edgecov(receiver_ineligible_mat_list))} are translated into
#' structural dyad exclusions for count MPLE.
#'
#' @param formula An ERGM formula for a single count-valued network.
#' @param response Character string giving the edge attribute that stores the
#'   count response.
#' @param reference Reference measure for the count ERGM. One of
#'   \code{"uniform"}, \code{"binomial"}, \code{"geometric"}, or
#'   \code{"poisson"}.
#' @param max.count Maximum possible count value. Use \code{Inf} for unbounded
#'   count support.
#' @param max.count.safety Safety multiplier used when approximating count
#'   support.
#' @param max.count.edgewise Logical. Should support truncation be done
#'   edgewise?
#' @param count.samples Maximum number of support values to evaluate per edge
#'   variable.
#' @param regularization Regularization method. One of \code{"none"},
#'   \code{"L2"}, \code{"L1"}, \code{"pseudoHuber"}, or \code{"L1pH"}.
#' @param regularization.param Numeric vector of regularization parameters.
#' @param sample.size Number of edge variables to sample. Use \code{Inf} to use
#'   all active dyads.
#' @param sample.method Edge-variable sampling method. One of
#'   \code{"weighted"} or \code{"random"}.
#' @param weight.type Weighting scheme for weighted sampling. One of
#'   \code{"TNT"} or \code{"flatval"}.
#' @param coef.init Optional initial coefficient vector.
#' @param cores Number of cores used for count-MPLE change-score computation.
#' @param optim.method Optimizer. Use \code{"trust"} for trust-region
#'   optimization, or another method accepted by \code{optim()}.
#' @param seed Optional random seed for edge-variable sampling.
#' @param must.count Integer. Always include support values from 0 through this
#'   value when building the local count support.
#' @param prep Optional precomputed \code{ERGMCntPrep} object.
#' @param count_offset Optional matrix marking inactive dyads. Entries equal to
#'   0 are treated as active. Nonzero or \code{NA} entries are excluded from the
#'   pseudo-likelihood.
#' @param WtSumAsSampSiz Logical. Should inverse inclusion weights sum to the
#'   sampled active edge-variable count?
#' @param ... Additional arguments passed to \code{trust::trust()} or
#'   \code{stats::optim()}.
#'
#' @return An object of class \code{"ergmCntMPLE"}.
#'
#' @importFrom ergm ergm_model nparam param_names
#' @importFrom stats as.formula optim optimHess rnorm pnorm
#' @importFrom MASS ginv
#' @keywords internal
ergmCntMPLE_btergm <- function(formula,
                               response,
                               reference = c("uniform", "binomial", "geometric", "poisson"),
                               max.count = Inf,
                               max.count.safety = 4,
                               max.count.edgewise = TRUE,
                               count.samples = Inf,
                               regularization = c("none", "L2", "L1", "pseudoHuber", "L1pH"),
                               regularization.param = c(1, 1e-4),
                               sample.size = Inf,
                               sample.method = c("weighted", "random"),
                               weight.type = c("TNT", "flatval"),
                               coef.init = NULL,
                               cores = 1,
                               optim.method = "trust",
                               seed = NULL,
                               must.count = 5,
                               prep = NULL,
                               count_offset = NULL,
                               WtSumAsSampSiz = TRUE,
                               estimate.cov = TRUE,
                               verbose = TRUE,
                               ...) {
  
  if (verbose) message("Starting count-valued ERGM MPLE.")
  if (verbose) message("Extracting network and preprocessing valued response.")
  
  nw0 <- ergm.getnetwork(formula);response2 <- response #here we use response2 because ergm_preprocess_response consumes the response string, but we need it later for ergmCntPrep (submitted as Issue#463 on github/ergm)
  ergm_preprocess_response_fun <- getFromNamespace("ergm_preprocess_response", "ergm")
  ergm_preprocess_response_fun(nw = nw0, response = response2)
  #generate new formula that use preprocessed networks
  
  ## OLD CODE (1):
  #mod<-ergm_model(formula=as.formula(paste0("nw0~",as.character(formula)[3])))
  
  ## NEW CODE (2): Current ergm versions expect the formula as the first argument to ergm_model(), not as a named argument called formula.
  #rhs <- paste(deparse(formula[[3]]), collapse = "")
  #form0 <- stats::as.formula(paste("nw0 ~", rhs), env = environment())
  #mod <- ergm::ergm_model(form0)
  
  ## NEW CODE (3): creates a new formula environment that contains nw0 but can still find memory, delrecip, and i through the parent environment.
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  # Use the formula's environment, because this is where btergm_count()
  # stores period-indexed temporal covariates such as memory, delrecip, and i.
  model_env <- new.env(parent = environment(formula))
  assign("nw0", nw0, envir = model_env)
  form0 <- stats::as.formula(paste("nw0 ~", rhs), env = model_env)
  
  mod <- ergm::ergm_model(form0)
  
  np<-nparam(mod)
  parnam<- param_names(mod)
  
  if (verbose) {message("Model parsed: ", np, " parameter(s).")}
  
  rm(nw0);gc()
  if(is.null(prep)){
    if (verbose) message("Preparing count pseudo-likelihood components.")
    prep <- ergmCntPrep_btergm(
      formula = formula,
      response = response,
      reference = reference, 
      max.count = max.count, 
      max.count.safety = max.count.safety, 
      max.count.edgewise = max.count.edgewise, 
      count.samples = count.samples, 
      sample.size = sample.size,
      sample.method = sample.method,
      cores = cores,
      weight.type = weight.type,
      seed = seed,
      must.count = must.count,
      count_offset = count_offset,
      WtSumAsSampSiz = WtSumAsSampSiz
    )
    
    if (verbose) message("Finished preparing count pseudo-likelihood components.")
  } else {
    if (verbose) message("Using precomputed count pseudo-likelihood components.")
  }
  rfun<-switch(match.arg(regularization),
               none=0,
               L1=1,
               L2=2,
               pseudoHuber=3,
               L1pH=3
  )
  if (optim.method == "trust") {
    if (!requireNamespace("trust", quietly = TRUE)) {stop("Package 'trust' is required for optim.method = 'trust'.", call. = FALSE)}
    usetrust <- TRUE
    ps <- rep(1, np)
  } else {
    usetrust <- FALSE
  }
  if(match.arg(regularization)=="L1pH"){
    regularization.param[2]<-0.1
  }
  
  if (verbose) {message("Starting optimization using method: ", optim.method, ".")}
  
  if(usetrust)
    fit<-trust::trust(objfun=ergmCntNLPLDeriv, parinit=coef.init, rinit=1, rmax=100, parscale=ps, obj=prep, rtype=rfun, rparam=regularization.param, ...)
  else
    fit<-stats::optim(coef.init, fn=ergmCntNLPL, obj=prep, rtype=rfun, rparam=regularization.param, method=optim.method, ...)
  while((match.arg(regularization)=="L1pH")&&(regularization.param[2]>1e-7)){
    regularization.param[2]<-regularization.param[2]/2
    
    if (verbose) {message("Reducing pseudo-Huber radius to ", regularization.param[2], ".") }
    
    #cat("\tReducing pseudo-Huber radius to",regularization.param[2],"\n")
    if(usetrust)
      fit<-trust::trust(objfun=ergmCntNLPLDeriv, parinit=fit$argument, rinit=1, rmax=100, parscale=ps, obj=prep, rtype=rfun, rparam=regularization.param, ...)
    else
      fit<-stats::optim(fit$par, fn=ergmCntNLPL, obj=prep, rtype=rfun, rparam=regularization.param, method=optim.method, ...)
  }
  
  if (verbose) {message("Optimization finished.")}
  
  if (usetrust){
    fit$par<-fit$argument
  }
  
  fit$coef<-fit$par
  names(fit$coef)<-parnam
  
  # Also store coefficients using the conventional name expected by coef.default().
  fit$coefficients <- fit$coef
  names(fit$coefficients) <- parnam
  
  fit$pseudo.deviance<-2*fit$value
  
  if (isTRUE(estimate.cov)) {
    
    if (verbose) {message("Computing Hessian and covariance matrix.")}
    
    if (usetrust) {
      fit$pll.hessian <- -fit$hessian
    } else {
      fit$pll.hessian <- -stats::optimHess(
        fit$par,
        fn = ergmCntNLPL,
        obj = prep,
        rtype = rfun,
        rparam = regularization.param
      )
    }
    
    fit$cov <- tryCatch(
      solve(-fit$pll.hessian),
      error = function(e) {
        warning("Could not invert the count-MPLE Hessian. ", "Using MASS::ginv() instead. Original error: ", conditionMessage(e))
        
        if (!requireNamespace("MASS", quietly = TRUE)) {
          matrix(
            NA_real_,
            nrow = length(fit$coef),
            ncol = length(fit$coef),
            dimnames = list(names(fit$coef), names(fit$coef))
          )
        } else {
          MASS::ginv(-fit$pll.hessian)
        }
      }
    )
    
    if (is.null(dimnames(fit$cov))) {
      dimnames(fit$cov) <- list(names(fit$coef), names(fit$coef))
    }
    
  } else {
    
    # During temporal bootstrap, we only need coefficients.
    # Skipping this avoids failures from singular replicate-specific Hessians.
    fit$pll.hessian <- NULL
    fit$cov <- NULL
    
    if (verbose) {message("Skipping Hessian/covariance calculation.")}
  }
  
  
  fit$formula<-formula
  fit$reference<-reference
  fit$response<-response
  fit$regularization<-match.arg(regularization)
  fit$regularization.param<-regularization.param
  fit$sample.method<-match.arg(sample.method)
  if(match.arg(sample.method)=="weighted")
    fit$weight.type<-match.arg(weight.type)
  class(fit)<-"ergmCntMPLE"
  
  if (verbose) {message("Finished count-valued ERGM MPLE.")}
  
  fit
}












#' Precompute count-valued ERGM pseudo-likelihood components
#'
#' Internal preparation function for count-valued ERGM MPLE. This function
#' computes reference-measure ratios, count support values, change scores, and
#' inverse inclusion weights for active edge variables. Compared with the
#' original Huang-Butts preparation function, this version allows inactive dyads
#' to be excluded through a \code{count_offset} matrix before edge-variable
#' sampling.
#'
#' @param formula An ERGM formula for a single count-valued network.
#' @param nw Optional network object. Retained for compatibility with the
#'   original Huang-Butts implementation.
#' @param response Character string naming the count edge attribute.
#' @param reference Reference measure for the count ERGM.
#' @param max.count Maximum possible count value.
#' @param max.count.safety Safety multiplier for approximating count support.
#' @param max.count.edgewise Logical. Should support truncation be edgewise?
#' @param count.samples Maximum number of support values to evaluate per edge.
#' @param sample.size Number of active edge variables to sample.
#' @param sample.method Edge-variable sampling method.
#' @param weight.type Weighting scheme for weighted sampling.
#' @param cores Number of cores used for change-score computation.
#' @param seed Optional random seed.
#' @param must.count Integer. Always include support values from 0 through this
#'   value.
#' @param count_offset Optional matrix marking inactive dyads. Entries equal to
#'   0 are active. Nonzero or \code{NA} entries are excluded.
#' @param WtSumAsSampSiz Logical. Should inverse inclusion weights sum to the
#'   sampled active edge-variable count?
#'
#' @return An object of class \code{"ERGMCntPrep"} containing precomputed
#'   pseudo-likelihood components.
#'
#' @importFrom ergm ergm.getnetwork ergm_model nparam param_names
#' @importFrom network network.dyadcount is.directed
#' @importFrom network as.sociomatrix 
#' @importFrom sna gvectorize
#' @importFrom sampling inclusionprobabilities UPpoisson
#' @importFrom parallel mclapply parLapply parLapplyLB makeCluster stopCluster clusterCall
#' @keywords internal
ergmCntPrep_btergm <- function(formula,
                               nw,
                               response,
                               reference = c("uniform", "binomial", "geometric", "poisson"),
                               max.count = Inf,
                               max.count.safety = 4,
                               max.count.edgewise = TRUE,
                               count.samples = Inf,
                               sample.size = Inf,
                               sample.method = c("weighted", "random"),
                               weight.type = c("TNT", "flatval"),
                               cores = 1,
                               prep.parallel = c("auto", "fork", "psock", "serial"),
                               prep.cl = NULL,
                               seed = NULL,
                               must.count = 1,
                               count_offset = NULL,
                               WtSumAsSampSiz = TRUE,
                               verbose = TRUE) {
  
  if (verbose) message("Preparing count-MPLE object.")
  if (verbose) message("Extracting network and response matrix.")
  
  prep.parallel <- match.arg(prep.parallel)
  
  if (prep.parallel == "auto") {prep.parallel <- if (.Platform$OS.type == "windows") "psock" else "fork"}
  
  if (prep.parallel == "fork" && .Platform$OS.type == "windows") {stop("prep.parallel = 'fork' is not available on Windows.")}
  
  #Get the network, and set things up
  nw<-ergm.getnetwork(formula)
  isdir<-is.directed(nw)  
  n<-network.size(nw)
  
  ## NEW CODE (1): # Make a self-contained formula for worker processes. PSOCK workers do not inherit objects like `g1` from the main R session, so we bind the extracted network to a local formula environment.
  #rhs <- paste(deparse(formula[[3]]), collapse = "")
  #formula_env <- new.env(parent = environment(formula))
  #formula_env$nw_for_godfather <- nw
  #formula_gf <- stats::as.formula(paste("nw_for_godfather ~", rhs), env = formula_env)
  
  ## NEW CODE (2): 
  rhs_expr <- formula[[3]]
  rhs <- paste(deparse(rhs_expr), collapse = "")
  
  formula_env <- new.env(parent = environment(formula))
  formula_env$nw_for_godfather <- nw
  
  # Copy non-function objects used in the RHS formula into the formula environment.
  # This is important for Windows PSOCK workers, which do not inherit objects from
  # the main R session.
  rhs_vars <- all.vars(rhs_expr)
  copied_vars <- character(0)
  
  for (nm in rhs_vars) {
    if (exists(nm, envir = environment(formula), inherits = TRUE)) {
      obj <- get(nm, envir = environment(formula), inherits = TRUE)
      
      # Do not copy functions or ERGM term names. We only need matrix/list/vector
      # objects that appear inside terms like edgecov(memory).
      if (!is.function(obj)) {
        assign(nm, obj, envir = formula_env)
        copied_vars <- c(copied_vars, nm)
      }
    }
  }
  
  if (verbose && length(copied_vars) > 0L) {message("Copied formula-side object(s) into worker formula environment: ", paste(copied_vars, collapse = ", "), ".")}
  
  formula_gf <- stats::as.formula(paste("nw_for_godfather ~", rhs), env = formula_env)
  
  
  
  ## OLD CODE (1):
  #nev<-network.dyadcount(nw)              #Number of edge variables
  #ss<-min(sample.size,nev)                #Sample size
  #lrmr<-switch(match.arg(reference),      #Function to calculate log ref mes rat
  #             uniform = function(k,y) {rep(0,length(k))},
  #             binomial = function(k,y) {lchoose(max.count,k)-lchoose(max.count,y)},
  #             geometric = function(k,y) {rep(0,length(k))},
  #             poisson = function(k,y) {lfactorial(y)-lfactorial(k)}
  #)
  ##Set up edge variable sampling - use inverse freq weighting to smooth a bit
  #ally<-as.sociomatrix(nw,response)       #All y values
  #alli<-row(ally)                         #Senders
  #allj<-col(ally)                         #Receivers
  #ally<-gvectorize(ally,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  #alli<-gvectorize(alli,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  #allj<-gvectorize(allj,mode=ifelse(isdir,"digraph","graph"),censor.as.na=FALSE)
  
  ## NEW CODE (2): replaces the old edge-variable setup block.
  ##               The original code computes nev and ss before knowing which
  ##               dyads are active. The modified version computes nev and ss
  ##               after removing inactive dyads.
  count_offset_supplied <- !is.null(count_offset)
  
  nev_total <- network.dyadcount(nw)  # Number of edge variables before applying count_offset.
  
  lrmr <- switch(
    match.arg(reference),
    uniform = function(k, y) { rep(0, length(k)) },
    binomial = function(k, y) { lchoose(max.count, k) - lchoose(max.count, y) },
    geometric = function(k, y) { rep(0, length(k)) },
    poisson = function(k, y) { lfactorial(y) - lfactorial(k) }
  )
  
  # Set up edge-variable matrices.
  ally_mat <- as.sociomatrix(nw, response) # All y values
  alli_mat <- row(ally_mat)                # Senders
  allj_mat <- col(ally_mat)                # Receivers
  
  # If supplied, count_offset must mark inactive dyads.
  # Convention:
  #   0 means active/eligible
  #   nonzero or NA means inactive/excluded
  if (!is.null(count_offset)) {
    count_offset <- as.matrix(count_offset)
    
    if (!identical(dim(count_offset), dim(ally_mat))) {
      if (!is.null(rownames(count_offset)) && !is.null(colnames(count_offset)) &&
          !is.null(rownames(ally_mat)) && !is.null(colnames(ally_mat))) {
        
        count_offset <- count_offset[
          rownames(ally_mat),
          colnames(ally_mat),
          drop = FALSE
        ]
        
      } else {
        stop(
          "'count_offset' is not conformable with the network response matrix. ",
          "It must have the same dimensions as as.sociomatrix(nw, response), ",
          "or it must have row and column names that can be used for alignment."
        )
      }
    }
  } else {
    count_offset <- matrix(
      0,
      nrow = nrow(ally_mat),
      ncol = ncol(ally_mat),
      dimnames = dimnames(ally_mat)
    )
  }
  
  # Vectorize outcome, sender indices, receiver indices, and offset indicators
  # using exactly the same mode, so all vectors refer to the same edge variables.
  vec_mode <- ifelse(isdir, "digraph", "graph")
  
  ally <- gvectorize(ally_mat, mode = vec_mode, censor.as.na = FALSE)
  alli <- gvectorize(alli_mat, mode = vec_mode, censor.as.na = FALSE)
  allj <- gvectorize(allj_mat, mode = vec_mode, censor.as.na = FALSE)
  offv <- gvectorize(count_offset, mode = vec_mode, censor.as.na = FALSE)
  
  # Keep only active dyads.
  active_idx <- which(!is.na(offv) & offv == 0)
  
  if (length(active_idx) == 0L) {stop("All dyads are excluded by 'count_offset'. ","Check the offset matrix: 0 should mean active, and nonzero/NA should mean inactive.")}
  
  ally <- ally[active_idx]
  alli <- alli[active_idx]
  allj <- allj[active_idx]
  
  # Number of edge variables after removing inactive dyads.
  nev <- length(ally) # total number of active edge variables
  
  # Sample size must be computed after excluding inactive dyads.
  ss <- min(sample.size, nev) # number of sampled active edge variables

  if (verbose) {message("Dyad filtering complete: ",nev, " active dyads out of ", nev_total, " total dyads.")}
  
  if (verbose) {message("Sampling edge variables: target sample size = ",ss, " active dyads.")}
  
  ## OLD CODE (1):
  #if(ss==nev){                            #If ss==nev, no sampling
  #  samp<-1:nev                              #Everyone's in the sample
  #  iw<-rep(1,nev)                           #Inclusion weight is 1
  #}else{                                  #Else, sample EVs
  #  if(match.arg(sample.method)=="random"){
  #    set.seed(seed)                         #Set seed for sampling
  #    samp<-sample(1:nev,ss)                 #Choose at random
  #    iw <- rep(1,ss)                        #Inclusion weight, before regularization
  #  }else{
  #    if(match.arg(weight.type)=="TNT"){          #"Tie/No-Tie" style weighting
  #      taby<-table(ally>0)  #c(FALSE, TRUE)
  #      wght<-(1/(2*taby))[1+(ally>0)]
  #      wght<-inclusionprobabilities(wght,ss)
  #      set.seed(seed)                             #Set seed for sampling
  #      samp<-which(UPpoisson(wght)>0)             #Draw the sample
  #      iw<- 1/(wght[samp])                        #Inclusion weight, before regularization
  #    }else if(match.arg(weight.type)=="flatval"){ #"Flat" value distribution weighting
  #      if(ss/nev>0.15)
  #        warning("Target sample size is ",round(ss/nev*100),"% of the total EV count.  Weighted sampling may be unreliable here - you may want to consider random sampling.")
  #      taby<-table(ally)
  #      wght<-(1/length(taby)/taby)[match(ally,names(taby))] #Ideal weights
  #      wght<-inclusionprobabilities(wght,ss)
  #      set.seed(seed)                           #Set seed for sampling
  #      samp<-which(UPpoisson(wght)>0)           #Draw the sample
  #      iw<- 1/(wght[samp])                        #Inclusion weight, before regularization
  #    }else{
  #      stop("Unknown weighting method ",weight.type," in ergmCntPrep.  Cannot go on like this.\n")
  #    }
  #  }
  #}
  
  ## NEW CODE (2): This does two things.
  ##               First, it samples only from the active dyad vector. 
  ##               Second, it makes the TNT weighting more robust if all active 
  ##               dyads are zero or all active dyads are nonzero. The input data 
  ##               probably have both, but this safeguard prevents a frustrating 
  ##               failure later.
  if (ss == nev) {
    # No sampling: all active dyads are included.
    samp <- seq_len(nev)
    iw <- rep(1, nev)
  } else {
    if (match.arg(sample.method) == "random") {
      
      set.seed(seed)
      samp <- sample(seq_len(nev), ss)
      iw <- rep(1, ss)
      
    } else {
      
      if (match.arg(weight.type) == "TNT") {
        
        # Tie/no-tie style weighting among active dyads only.
        grp <- ally > 0
        taby <- table(grp)
        
        wght <- numeric(length(grp))
        
        if ("FALSE" %in% names(taby)) {
          wght[!grp] <- 1 / (length(taby) * as.numeric(taby[["FALSE"]]))
        }
        
        if ("TRUE" %in% names(taby)) {
          wght[grp] <- 1 / (length(taby) * as.numeric(taby[["TRUE"]]))
        }
        
        # Fallback if all active dyads are in a single group.
        if (any(wght <= 0) || any(!is.finite(wght))) {
          wght <- rep(1 / length(grp), length(grp))
        }
        
        wght <- sampling::inclusionprobabilities(wght, ss)
        
        set.seed(seed)
        samp <- which(sampling::UPpoisson(wght) > 0)
        iw <- 1 / wght[samp]
        
      } else if (match.arg(weight.type) == "flatval") {
        
        if (ss / nev > 0.15) {
          warning(
            "Target sample size is ",
            round(ss / nev * 100),
            "% of the active edge-variable count. Weighted sampling may be unreliable here; ",
            "you may want to consider random sampling."
          )
        }
        
        taby <- table(ally)
        
        wght <- 1 / (length(taby) * as.numeric(taby[as.character(ally)]))
        
        if (any(wght <= 0) || any(!is.finite(wght))) {
          wght <- rep(1 / length(ally), length(ally))
        }
        
        wght <- sampling::inclusionprobabilities(wght, ss)
        
        set.seed(seed)
        samp <- which(sampling::UPpoisson(wght) > 0)
        iw <- 1 / wght[samp]
        
      } else {
        stop(
          "Unknown weighting method ",
          weight.type,
          " in ergmCntPrep_btergm()."
        )
      }
    }
  }
  
  ss<-length(samp)                        #Should not have changed, but can...
  
  if (verbose) {message("Edge-variable sampling complete: ", ss, " dyads selected.")}
  
  #Inclusion weight: regularize it to sample size or total dyads.
  if (WtSumAsSampSiz) {
    iw <- iw * ss / sum(iw)
  } else {
    iw <- iw * nev / sum(iw)
  }
  #There may be some ignoreable numerical difference
  if(WtSumAsSampSiz && !isTRUE(all.equal(sum(iw), ss))) warning("WtSumAsSampSiz = TRUE, but sum of inclusion weights is unequal to sample size; diff = ",sum(iw) - ss)
  if(!WtSumAsSampSiz && !isTRUE(all.equal(sum(iw), nev))) warning("WtSumAsSampSiz = FALSE, but sum of inclusion weights is unequal to active edge-variable count; diff = ",sum(iw) - nev)
  
  y<-ally[samp]                           # Observed sampled y values
  snd<-alli[samp]                         # Sampled senders
  rec<-allj[samp]                         # Sampled receivers
  
  if (verbose) message("Constructing local count supports.")
  
  #Walk through edge variables and calculate exciting things
  if(max.count.edgewise){        #Use per-edge values
    yub<-pmin(pmax(must.count+1,ceiling(y+max.count.safety*4*sqrt(y))),max.count)  #Upper bounds
    ylb<-pmax(must.count+1,floor(y-max.count.safety*4*sqrt(y)))    #Lower bounds (0:must.count always included)
    yrng<-vector(mode="list",length=ss)        #y values for each edge
    ycwt<-vector(mode="list",length=ss)        #Weights for approximation
    for(i in 1:ss){
      if(count.samples>=yub[i]-ylb[i]){
        yrng[[i]]<-c(0:must.count,ylb[i]:yub[i])
        ycwt[[i]]<-rep(1,length(yrng[[i]]))
      }else{
        yrng[[i]]<-c(0:must.count,round(seq(from=ylb[i],to=yub[i],length=count.samples)))
        ycwt[[i]]<-c(rep(1, must.count+1),1,diff(round(seq(from=ylb[i],to=yub[i],length=count.samples))))
      }
    } 
  }else{                         #Use uniform edge value ranges
    maxy<-max(ally)                         #Maximum y value
    yrng<-vector(mode="list",length=ss)     #y values for each edge
    ycwt<-vector(mode="list",length=ss)     #Weights for approximation
    if(max.count<Inf){                      #Range of y values for normalization
      if(count.samples>=max.count+1){
        for(i in 1:ss){
          yrng[[i]]<-0:max.count
          ycwt[[i]]<-rep(1,1+max.count)
        }
      }else{
        for(i in 1:ss){
          yrng[[i]]<-c(0,round(seq(from=1,to=max.count,length=count.samples)))
          ycwt[[i]]<-c(1,1,diff(yrng[[i]]))
        }
      }
    }else{
      if(count.samples>=max.count.safety*maxy+1){
        for(i in 1:ss){
          yrng[[i]]<-0:(max.count.safety*maxy)
          ycwt[[i]]<-rep(1,1+max.count.safety*maxy)
        }
      }else{
        for(i in 1:ss){
          yrng[[i]]<-round(seq(from=0,to=ceiling(max.count.safety*maxy),length=count.samples))
          ycwt[[i]]<-c(diff(yrng[[i]]),1)
        }
      }
    }
  }
  
  if (verbose) message("Finished constructing local count supports.")
  
  
  ## NEW CODE (1): create one Windows PSOCK cluster once and reuse it for both `rmr` and `cs`
  #win_cl <- NULL
  #if (.Platform$OS.type == "windows" && cores > 1L) {
  #  
  #  if (verbose) {message("Starting Windows PSOCK cluster with ", cores, " cores.")}
  #  
  #  win_cl <- parallel::makeCluster(cores, type = "PSOCK")
  #  on.exit(parallel::stopCluster(win_cl), add = TRUE)
  #  
  #  parallel::clusterCall(win_cl, function() {
  #    library(ergm)
  #    library(network)
  #    library(sna)
  #    library(statnet.common)
  #    NULL
  #  })
  #}
  
  ## NEW CODE (2):
  prep_cl <- prep.cl
  own_prep_cluster <- FALSE
  if (cores > 1L && prep.parallel == "psock") {
    if (is.null(prep_cl)) {
      if (verbose) {message("Starting PSOCK cluster with ", cores, " cores for count-MPLE preparation.")}
      prep_cl <- parallel::makeCluster(cores, type = "PSOCK")
      own_prep_cluster <- TRUE
      on.exit({if (own_prep_cluster && !is.null(prep_cl)) {parallel::stopCluster(prep_cl)}}, add = TRUE)
    }
    parallel::clusterCall(prep_cl, function() {
      pkgs <- c(
        "ergm",
        "network",
        "sna",
        "statnet.common",
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
          "The following packages could not be loaded on the preparation workers: ",
          paste(pkgs[!ok], collapse = ", ")
        )
      }
      
      NULL
    })
  } else if (cores > 1L && prep.parallel == "fork") {
    if (verbose) {message("Using forked parallelism with ", cores, " cores for count-MPLE preparation.")}
  } else {
    if (verbose) {message("Using serial computation for count-MPLE preparation.")}
  }
  
  
  if (verbose) message("Computing reference-measure ratios.")
  
  ## OLD CODE (1): #Ref meas ratio (one entry per EV)
  #rmr<-mclapply(1:ss,function(i){lrmr(yrng[[i]],y[i])},mc.cores=cores) 
  ## NEW CODE (2):
  #rmr <- safe_mclapply(seq_len(ss), function(i) {lrmr(yrng[[i]], y[i])}, cores = cores, cl = win_cl)
  ## NEW CODE (3):
  rmr_fun <- function(i) {lrmr(yrng[[i]], y[i])}
  # For PSOCK, avoid serializing the full local function environment if possible.
  if (prep.parallel == "psock" && cores > 1L) {environment(rmr_fun) <- .GlobalEnv}
  rmr <- safe_mclapply(
    seq_len(ss),
    rmr_fun,
    cores = cores,
    cl = prep_cl,
    parallel.type = prep.parallel,
    packages = character(),
    export = c("lrmr", "yrng", "y"),
    exportenv = environment()
  )
  
  if (verbose) message("Finished computing reference-measure ratios.")
  
  
  
  if (verbose) message("Computing change scores. This is usually the slowest step.")
  
  ## OLD CODE (1): #Changescore list (one entry per EV)
  #cs<-mclapply(1:ss,function(i){ergm.godfather(formula=formula, changes=lapply(yrng[[i]],function(z){matrix(c(snd[i],rec[i],z),nrow=1)}), response=response, changes.only=TRUE)},mc.cores=cores)
  ## NEW CODE (2):
  #cs <- safe_mclapply(seq_len(ss), function(i) {ergm.godfather(formula = formula, changes = lapply(yrng[[i]], function(z) {matrix(c(snd[i], rec[i], z), nrow = 1)}), response = response, changes.only = TRUE)}, cores = cores)
  ## NEW CODE (3): 
  #cs <- safe_mclapply(seq_len(ss), function(i) {ergm::ergm.godfather(object = formula, changes = lapply(yrng[[i]], function(z) {matrix(c(snd[i], rec[i], z), nrow = 1)}), response = response, changes.only = TRUE)}, cores = cores, cl = win_cl, packages = c("ergm", "network", "sna", "statnet.common"))
  ## NEW CODE (4): 
  #cs <- safe_mclapply(seq_len(ss), function(i) {ergm::ergm.godfather(object = formula_gf, changes = lapply(yrng[[i]], function(z) {matrix(c(snd[i], rec[i], z), nrow = 1)}), response = response, changes.only = TRUE)}, cores = cores, cl = win_cl, packages = c("ergm", "network", "sna", "statnet.common"))
  ## NEW CODE (5):
  cs_fun <- function(i) {ergm::ergm.godfather(object = formula_gf, changes = lapply(yrng[[i]], function(z) {matrix(c(snd[i], rec[i], z), nrow = 1)}), response = response, changes.only = TRUE)}
  # For PSOCK, avoid capturing the full local parent environment.
  if (prep.parallel == "psock" && cores > 1L) {environment(cs_fun) <- .GlobalEnv}
  cs <- safe_mclapply(
    seq_len(ss),
    cs_fun,
    cores = cores,
    cl = prep_cl,
    parallel.type = prep.parallel,
    packages = c("ergm", "network", "sna", "statnet.common", "btergm"),
    export = c("formula_gf", "yrng", "snd", "rec", "response"),
    exportenv = environment()
  )
  
  if (verbose) message("Finished computing change scores.")
  
  #Return a list with all the goodies
  
  ## OLD CODE (1): 
  #out<-list(snd=snd,rec=rec,y=y,rmr=rmr,cs=cs,ycwt=ycwt,iw=iw)
  ## NEW CODE (2): Modify prepared object so it stores basic active-dyad information
  out <- list(
    rmr = rmr,
    cs = cs,
    iw = iw,
    ycwt = ycwt,
    y = y,
    snd = snd,
    rec = rec,
    nev = nev,
    nev_total = nev_total,
    active_idx = active_idx,
    count_offset_used = count_offset_supplied
  )
  
  class(out)<-"ERGMCntPrep"
  
  if (verbose) message("Count-MPLE preparation complete.")
  
  out
}

























#' Safe parallel lapply helper
#'
#' Internal helper that uses forked parallelism on Unix-like systems and PSOCK
#' clusters on Windows.
#'
#' @param X List or vector to iterate over.
#' @param FUN Function to apply.
#' @param cores Number of cores.
#' @param ... Additional arguments passed to \code{FUN}.
#' @param cl Optional cluster object.
#' @param packages Character vector of packages to load on PSOCK workers.
#' @param load_balance Logical. Use load-balanced parallel apply?
#'
#' @return A list.
#'
#' @keywords internal
safe_mclapply <- function(X,
                          FUN,
                          cores = 1L,
                          ...,
                          cl = NULL,
                          packages = character(),
                          load_balance = TRUE,
                          parallel.type = c("auto", "fork", "psock", "serial"),
                          export = character(),
                          exportenv = parent.frame()) {
  
  cores <- as.integer(cores)
  parallel.type <- match.arg(parallel.type)
  
  # Windows: use PSOCK cluster; Unix/macOS/Linux: use forked parallelism
  if (parallel.type == "auto") {parallel.type <- if (.Platform$OS.type == "windows") "psock" else "fork"}
  
  # Serial fallback:
  if (cores <= 1L || parallel.type == "serial") {return(lapply(X, FUN, ...))}
  
  
  if (parallel.type == "fork") {
    if (.Platform$OS.type == "windows") {stop("Forked parallelism is not available on Windows.")}
    return(parallel::mclapply(X, FUN, ..., mc.cores = cores))
  }
  
  if (parallel.type == "psock") {
    own_cluster <- is.null(cl)
    
    if (own_cluster) {
      cl <- parallel::makeCluster(cores, type = "PSOCK")
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
    # Load needed packages on Windows workers
    if (length(packages) > 0L) {parallel::clusterCall(cl, function(pkgs) {invisible(lapply(pkgs, library, character.only = TRUE))}, packages)}
    if (length(export) > 0L) {parallel::clusterExport(cl, varlist = export, envir = exportenv)}
    
    # Load-balanced version is usually better because some edge variables have much larger count-support sets than others.
    if (load_balance) {
      return(parallel::parLapplyLB(cl, X, FUN, ...))
    } else {
      return(parallel::parLapply(cl, X, FUN, ...))
    }
  }
  
  stop("Unknown parallel.type: ", parallel.type)
}























# Count-MPLE objective functions.
#
# These C++ functions are called internally by ergmCntMPLE_btergm().
# They operate on an already prepared ERGMCntPrep object, so they do not
# need to know about temporal structure, bootstrap resampling, or offset
# matrices. Dyad exclusion happens upstream in ergmCntPrep_btergm().

#Negative log pseudo-likelihood for ERGM count models, using precomputed quantities.
Rcpp::cppFunction('double ergmCntNLPL(NumericVector coef, List obj, int rtype, NumericVector rparam){
  double lpl=0.0,ils,pd,reg=0.0;
  int i,j,k;

  //Make sure that we were passed a legitimate input, lest we crash
  if(!obj.inherits("ERGMCntPrep"))
    stop("Must be called with an ERGMCntPrep object.");

  //Go ahead and coerce the importance weights
  NumericVector iw = as<NumericVector>(obj["iw"]);

  //Gonna add it up
  for(i=0;i<iw.size();i++){
    //Set up the variables we will need
    ils=0.0;
    NumericVector rmr = as<NumericVector>(as<List>(obj["rmr"])[i]);
    NumericMatrix cs = as<NumericMatrix>(as<List>(obj["cs"])[i]);
    NumericVector ycwt = as<NumericVector>(as<List>(obj["ycwt"])[i]);
    //Find the inverse log sum
    for(j=0;j<rmr.size();j++){
      pd=0.0;
      for(k=0;k<coef.size();k++)
        pd+=cs(j,k)*coef[k];
      if(j==0)
        ils=rmr[j]+pd+log(ycwt[j]);
      else
        ils=R::logspace_add(ils,rmr[j]+pd+log(ycwt[j]));
    }
    //Add to the total (multiplying by inverse inclusion weight)
    lpl-=iw[i]*ils;
  }
  
  //Check process
  //Rcout << "Coef:";
  //for(i=0;i<coef.size();i++)
  //  Rcout << " " << coef[i];
  //Rcout << " LPL: " << lpl << "\\n";

  //Compute the regularization penalty, if any
  if(rtype==1){                              //L1 regularization
    for(i=0;i<coef.length();i++)
      reg+=fabs(coef[i]);
    reg*=rparam[0];
  }else if(rtype==2){                        //L2 regularization
    for(i=0;i<coef.length();i++)
      reg+=coef[i]*coef[i];
    reg*=rparam[0];
  }else if(rtype==3){                        //pseudo-Huber regularization
    for(i=0;i<coef.length();i++)
      reg+=rparam[1]*(sqrt(1.0+coef[i]*coef[i]/(rparam[1]*rparam[1]))-1.0);
    reg*=rparam[0];
  }

  //Return the negative log pseudo-likelihood (plus any regularization penalty)
  return -lpl+reg;
}')


#Negative log pseudo-likelihood and derivatives for ERGM count models, using precomputed 
#quantities.  Note that this output is formatted for the trust package.
Rcpp::cppFunction('Rcpp::List ergmCntNLPLDeriv(NumericVector coef, List obj, int rtype, NumericVector rparam){
  double lpl=0.0,ils,pd,reg=0.0;
  int i,j,k,l,p=coef.size();
  NumericVector gr(p),igr(p);                //Gradient of the lpl (penalized)
  NumericMatrix hess(p,p),ihess(p,p);        //Hessian of the lpl (penalized)

  //Make sure that we were passed a legitimate input, lest we crash
  if(!obj.inherits("ERGMCntPrep"))
    stop("Must be called with an ERGMCntPrep object.");

  //Go ahead and coerce the importance weights
  NumericVector iw = as<NumericVector>(obj["iw"]);

  //Gonna add it up
  for(i=0;i<iw.size();i++){
    //Set up the variables we will need
    ils=0.0;
    NumericVector rmr = as<NumericVector>(as<List>(obj["rmr"])[i]);
    NumericMatrix cs = as<NumericMatrix>(as<List>(obj["cs"])[i]);
    NumericVector ycwt = as<NumericVector>(as<List>(obj["ycwt"])[i]);
    for(j=0;j<p;j++){
      igr[j]=0.0;
      for(k=0;k<p;k++){
        ihess(j,k)=0.0;
      }
    }
    //Find the inverse log sum
    for(j=0;j<rmr.size();j++){
      pd=0.0;
      for(k=0;k<p;k++)
        pd+=cs(j,k)*coef[k];
      if(j==0)
        ils=rmr[j]+pd+log(ycwt[j]);
      else
        ils=R::logspace_add(ils,rmr[j]+pd+log(ycwt[j]));
    }
    //Have to make a second pass to get derivative elements (this is slower than one
    //pass, but we need to precompute ils to avoid overflow issues)
    for(j=0;j<rmr.size();j++){
      pd=0.0;
      for(k=0;k<p;k++)
        pd+=cs(j,k)*coef[k];
      for(k=0;k<p;k++){
        igr[k]-=cs(j,k)*exp(rmr[j]+pd-ils)*ycwt[j];   //Local gradient contribution
        for(l=0;l<p;l++){                             //Local Hessian contribution
          ihess(k,l)+=cs(j,k)*cs(j,l)*exp(rmr[j]+pd-ils)*ycwt[j];
        }
      }
    }
    //Add to the total (multiplying by inverse inclusion weight)
    lpl-=iw[i]*ils;                                                     //LPL
    for(j=0;j<p;j++){
      gr[j]-=iw[i]*igr[j];                                              //Gradient NLPL
      for(k=0;k<p;k++){
        hess(j,k)-=iw[i]*(igr[j]*igr[k]-ihess(j,k));                    //Hessian NLPL
      }
    }
  }

  //Compute the regularization penalty, if any
  if(rtype==1){                              //L1 regularization
    for(i=0;i<p;i++){
      reg+=fabs(coef[i]);
      gr[i]+=(coef[i] > 0.0 ? 1.0 : (coef[i] < 0.0 ? -1.0 : 0.0))*rparam[0];
    }
    reg*=rparam[0];
  }else if(rtype==2){                        //L2 regularization
    for(i=0;i<p;i++){
      reg+=coef[i]*coef[i];
      gr[i]+=2.0*coef[i]*rparam[0];
      hess(i,i)+=2.0*rparam[0];
    }
    reg*=rparam[0];
  }else if(rtype==3){                        //pseudo-Huber regularization
    for(i=0;i<p;i++){
      reg+=rparam[1]*(sqrt(1.0+coef[i]*coef[i]/(rparam[1]*rparam[1]))-1.0);
      gr[i]+=rparam[0]*coef[i]/(rparam[1]*sqrt(1.0+coef[i]*coef[i]/(rparam[1]*rparam[1])));
      hess(i,i)+=rparam[0]*rparam[1]/((rparam[1]*rparam[1]+coef[i]*coef[i])*sqrt(1.0+coef[i]*coef[i]/(rparam[i]*rparam[i])));
    }
    reg*=rparam[0];
  }

  //Return the penalized negative log pseudo-likelihood, gradient, and Hessian
  List out = List::create(Named("value") = -lpl+reg, Named("gradient") = gr, Named("hessian") = hess);
  return out;
}')