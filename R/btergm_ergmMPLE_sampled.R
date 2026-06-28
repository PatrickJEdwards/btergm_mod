#' Sample informative dyads for ERGM MPLE
#'
#' Internal helper used by \code{btergm_ergmMPLE_sampled()} to reduce the
#' informative dyad set before calling \code{ergm}'s MPLE change-statistic
#' machinery. The function supports simple random dyad sampling, approximate
#' tie/no-tie balanced sampling, and a mode that keeps all observed ties while
#' sampling non-ties.
#'
#' This function is intended for internal use. The main user-facing wrapper is
#' \code{btergm_ergmMPLE_sampled()}.
#'
#' @param nw A \pkg{network} object.
#' @param fd An \code{ergm} informative-dyad mask, usually an \code{rlebdm}
#'   object created by \code{ergm}'s internal \code{as.rlebdm()}.
#' @param sample.size Target number of informative dyads to retain. Use
#'   \code{Inf} for no sampling.
#' @param sample.method Dyad sampling method. \code{"none"} preserves ordinary
#'   \code{ergmMPLE()} behavior. \code{"random"} samples informative dyads
#'   uniformly. \code{"TNT"} approximately balances ties and non-ties.
#'   \code{"all_ties"} keeps all observed ties and samples non-ties.
#' @param seed Optional random seed.
#' @param WtSumAsSampSiz Logical. If \code{TRUE}, rescale inverse-probability
#'   weights so that their sum equals the sampled dyad count.
#' @param as.rlebdm_fun Function used to convert a logical matrix back to an
#'   \code{rlebdm} object.
#'
#' @return A list containing the sampled informative-dyad mask, sampling weights,
#'   and sampling diagnostics.
#'
#' @keywords internal
sample_mple_fd <- function(
    nw,
    fd,
    sample.size = Inf,
    sample.method = c("none", "random", "TNT", "all_ties"),
    seed = NULL,
    WtSumAsSampSiz = TRUE,
    as.rlebdm_fun = getFromNamespace("as.rlebdm", "ergm")
  ) {
  
  sample.method <- match.arg(sample.method)
  
  # No sampling.
  if (is.infinite(sample.size) || sample.method == "none") {
    return(list(
      fd = fd,
      ipw = NULL,
      sampled = FALSE
    ))
  }
  
  # Convert informative-dyad mask to a logical matrix.
  # This assumes fd has an as.matrix() method. If not, use the relevant
  # statnet/ergm internal converter for rlebdm objects.
  fd_mat <- as.matrix(fd)
  fd_mat <- !is.na(fd_mat) & fd_mat != 0
  
  dy <- which(fd_mat, arr.ind = TRUE)
  colnames(dy) <- c("tail", "head")
  
  nev <- nrow(dy)
  ss <- min(sample.size, nev)
  
  if (ss == nev) {
    ipw <- rep(1, nev)
    names(ipw) <- paste(dy[, "tail"], dy[, "head"], sep = "->")
    return(list(
      fd = fd,
      ipw = ipw,
      sampled = FALSE
    ))
  }
  
  # Get observed binary tie values for informative dyads.
  ymat <- as.matrix(nw, matrix.type = "adjacency")
  
  # dy contains full ergm vertex indices.
  # For bipartite networks, fd usually indexes heads as b1 + 1, ..., b1 + b2,
  # while as.matrix.network.adjacency(nw) may return a rectangular b1 x b2 matrix.
  dy_y <- dy
  
  n_total <- network::network.size(nw)
  is_bip <- network::is.bipartite(nw)
  
  if (isTRUE(is_bip)) {
    b1 <- network::get.network.attribute(nw, "bipartite")
    b2 <- n_total - b1
    
    # If ymat is rectangular b1 x b2, convert full head vertex IDs to
    # rectangular column indices.
    if (nrow(ymat) == b1 && ncol(ymat) == b2) {
      dy_y[, "head"] <- dy_y[, "head"] - b1
    }
  }
  
  # Safety check before indexing.
  if (
    any(dy_y[, "tail"] < 1L) ||
    any(dy_y[, "tail"] > nrow(ymat)) ||
    any(dy_y[, "head"] < 1L) ||
    any(dy_y[, "head"] > ncol(ymat))
  ) {
    stop(
      "Dyad indices do not align with the response matrix. ",
      "fd matrix dimension is ", paste(dim(as.matrix(fd)), collapse = " x "), "; ",
      "response matrix dimension is ", paste(dim(ymat), collapse = " x "), "."
    )
  }
  
  y <- ymat[cbind(dy_y[, "tail"], dy_y[, "head"])]
  tie <- y != 0
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (sample.method == "random") {
    
    samp <- sample(seq_len(nev), ss, replace = FALSE)
    
    pi <- rep(ss / nev, nev)
    iw <- 1 / pi[samp]
    
  } else if (sample.method == "TNT") {
    
    # Approximate tie/no-tie balanced sample.
    tie_idx <- which(tie)
    notie_idx <- which(!tie)
    
    n_tie <- length(tie_idx)
    n_notie <- length(notie_idx)
    
    # Target half ties, half non-ties, subject to availability.
    ss_tie <- min(n_tie, ceiling(ss / 2))
    ss_notie <- min(n_notie, ss - ss_tie)
    
    # If one group is too small, fill remaining sample from the other group.
    remaining <- ss - ss_tie - ss_notie
    if (remaining > 0 && n_tie > ss_tie) {
      add <- min(remaining, n_tie - ss_tie)
      ss_tie <- ss_tie + add
      remaining <- remaining - add
    }
    if (remaining > 0 && n_notie > ss_notie) {
      add <- min(remaining, n_notie - ss_notie)
      ss_notie <- ss_notie + add
      remaining <- remaining - add
    }
    
    samp_tie <- if (ss_tie > 0) sample(tie_idx, ss_tie, replace = FALSE) else integer(0)
    samp_notie <- if (ss_notie > 0) sample(notie_idx, ss_notie, replace = FALSE) else integer(0)
    samp <- c(samp_tie, samp_notie)
    
    pi <- numeric(nev)
    if (n_tie > 0) pi[tie_idx] <- ss_tie / n_tie
    if (n_notie > 0) pi[notie_idx] <- ss_notie / n_notie
    
    iw <- 1 / pi[samp]
    
  } else if (sample.method == "all_ties") {
    
    tie_idx <- which(tie)
    notie_idx <- which(!tie)
    
    n_tie <- length(tie_idx)
    n_notie <- length(notie_idx)
    
    # "all_ties" should literally keep all observed ties.
    # If sample.size is smaller than the number of ties, stop rather than
    # silently sampling ties away.
    if (n_tie > sample.size) {
      stop(
        "sample.method = 'all_ties' requires sample.size >= number of observed ties. ",
        "There are ", n_tie, " observed ties, but sample.size = ", sample.size, ". ",
        "Increase sample.size or use sample.method = 'TNT'."
      )
    }
    
    # Keep all ties and sample as many non-ties as needed, up to sample.size.
    ss_notie <- min(n_notie, sample.size - n_tie)
    
    samp_notie <- if (ss_notie > 0) {
      sample(notie_idx, ss_notie, replace = FALSE)
    } else {
      integer(0)
    }
    
    samp <- c(tie_idx, samp_notie)
    
    pi <- numeric(nev)
    pi[tie_idx] <- 1
    
    if (n_notie > 0 && ss_notie > 0) {
      pi[notie_idx] <- ss_notie / n_notie
    }
    
    iw <- 1 / pi[samp]
    
  } else {
    stop("Unknown sample.method.")
  }
  
  # Optional rescaling. Multiplying all weights by a constant should not change
  # coefficient estimates in the GLM, but it can help numerical behavior.
  if (WtSumAsSampSiz) {
    iw <- iw * length(iw) / sum(iw)
  } else {
    iw <- iw * nev / sum(iw)
  }
  
  dy_samp <- dy[samp, , drop = FALSE]
  
  fd_samp_mat <- matrix(
    FALSE,
    nrow = nrow(fd_mat),
    ncol = ncol(fd_mat),
    dimnames = dimnames(fd_mat)
  )
  
  fd_samp_mat[cbind(dy_samp[, "tail"], dy_samp[, "head"])] <- TRUE
  
  # This is the line most likely to need adjustment depending on namespace.
  # If you are modifying ergm::ergmMPLE() inside the ergm namespace, as.rlebdm()
  # should be visible. If you are doing this inside btergm, you may need
  # getFromNamespace().
  fd_samp <- as.rlebdm_fun(fd_samp_mat)
  
  ipw <- iw
  names(ipw) <- paste(dy_samp[, "tail"], dy_samp[, "head"], sep = "->")
  
  list(
    fd = fd_samp,
    ipw = ipw,
    sampled = TRUE,
    sampled_dyads = dy_samp,
    nev_total = nev,
    nev_sampled = length(samp)
  )
}







#' Evaluate ergm_conlist inside the ergm namespace
#'
#' Internal helper used because \code{ergm_conlist()} is an internal
#' \pkg{ergm} generic whose methods may not dispatch correctly from another
#' package namespace.
#'
#' @param x Constraint term list.
#' @param nw A \pkg{network} object.
#' @param term.options ERGM term options.
#'
#' @return An \code{ergm_conlist} object.
#'
#' @keywords internal
.btergm_make_ergm_conlist <- function(x, nw, term.options) {
  eval_env <- new.env(parent = asNamespace("ergm"))
  
  assign("x", x, envir = eval_env)
  assign("nw", nw, envir = eval_env)
  assign("term.options", term.options, envir = eval_env)
  
  eval(
    quote(ergm_conlist(x, nw, term.options = term.options)),
    envir = eval_env
  )
}


#' Retrieve a formula updater compatible with the installed statnet version
#'
#' @return A function for updating ERGM formulas.
#'
#' @keywords internal
.btergm_nonsimp_update_formula <- function() {
  tryCatch(
    getFromNamespace("nonsimp_update.formula", "ergm"),
    error = function(e1) {
      tryCatch(
        getFromNamespace("nonsimp_update.formula", "statnet.common"),
        error = function(e2) {
          stats::update.formula
        }
      )
    }
  )
}


#' Retrieve arr_from_coo from the available statnet namespace
#'
#' @return The internal \code{arr_from_coo()} helper.
#'
#' @keywords internal
.btergm_arr_from_coo <- function() {
  tryCatch(
    getFromNamespace("arr_from_coo", "ergm"),
    error = function(e) {
      getFromNamespace("arr_from_coo", "statnet.common")
    }
  )
}












#' ERGM MPLE with optional informative-dyad sampling
#'
#' This function is a wrapper around \code{\link[ergm]{ergmMPLE}} that preserves
#' ordinary \code{ergmMPLE()} behavior when no dyad sampling is requested, but
#' optionally samples the informative dyad set before \code{ergm}'s MPLE design
#' matrix is constructed. This is useful for large bipartite TERGM applications
#' where constructing MPLE rows for every informative dyad can be memory
#' intensive.
#'
#' When \code{sample.method = "none"} or \code{sample.size = Inf}, the function
#' simply calls \code{ergm::ergmMPLE()} and returns its result. When sampling is
#' active, the function constructs the informative dyad mask, samples dyads using
#' \code{sample_mple_fd()}, calls \code{ergm}'s internal MPLE design machinery on
#' the sampled dyads, and applies inverse-probability weights to the resulting
#' MPLE weights.
#'
#' The sampled implementation currently supports \code{output = "matrix"},
#' \code{output = "dyadlist"}, and \code{output = "array"}. It does not support
#' \code{output = "fit"} when sampling is active.
#'
#' @param formula An ERGM formula.
#' @param constraints ERGM constraints formula. Defaults to \code{~.}.
#' @param obs.constraints Observational constraints formula. Defaults to
#'   \code{~-observed}.
#' @param output Output type. See \code{\link[ergm]{ergmMPLE}}. With sampling,
#'   \code{"matrix"}, \code{"dyadlist"}, and \code{"array"} are supported.
#' @param expand.bipartite Logical. Whether to expand bipartite output as in
#'   \code{ergmMPLE()}.
#' @param control Control settings from \code{\link[ergm]{control.ergm}}.
#' @param verbose Logical. Print diagnostic output from ERGM internals.
#' @param ... Additional arguments passed to \code{ergm} internals.
#' @param basis Network basis. Defaults to \code{ergm::ergm.getnetwork(formula)}.
#' @param sample.size Target number of informative dyads to retain. Use
#'   \code{Inf} for no sampling.
#' @param sample.method Dyad sampling method. One of \code{"none"},
#'   \code{"random"}, \code{"TNT"}, or \code{"all_ties"}.
#' @param WtSumAsSampSiz Logical. If \code{TRUE}, rescale inverse-probability
#'   weights so that their sum equals the sampled dyad count.
#' @param seed Optional random seed for dyad sampling.
#'
#' @return A list with \code{response}, \code{predictor}, and \code{weights},
#'   with an \code{etamap} attribute, matching the structure used by
#'   \code{ergmMPLE()} for matrix-like outputs.
#'
#' @keywords internal
#'
#' @seealso \code{\link[ergm]{ergmMPLE}}
#' @export
btergm_ergmMPLE_sampled <- function(
    formula,
    constraints = ~.,
    obs.constraints = ~-observed,
    output = c("matrix", "array", "dyadlist", "fit"),
    expand.bipartite = FALSE,
    control = ergm::control.ergm(),
    verbose = FALSE,
    ...,
    basis = ergm::ergm.getnetwork(formula),
    
    sample.size = Inf,
    sample.method = c("none", "random", "TNT", "all_ties"),
    WtSumAsSampSiz = TRUE,
    seed = NULL
) {
  
  output <- match.arg(output)
  sample.method <- match.arg(sample.method)

  # Preserve original ergm behavior when no dyad sampling is requested.
  if (sample.method == "none" || is.infinite(sample.size)) {
    return(ergm::ergmMPLE(
      formula = formula,
      constraints = constraints,
      obs.constraints = obs.constraints,
      output = output,
      expand.bipartite = expand.bipartite,
      control = control,
      verbose = verbose,
      ...,
      basis = basis
    ))
  }
  
  # Sampled MPLE does not currently support output = "fit".
  if (output == "fit") {
    stop(
      "btergm_ergmMPLE_sampled() does not support output = 'fit' ",
      "when sample.method != 'none'. Use output = 'matrix' or 'dyadlist'."
    )
  }
  
  # Pull ergm internals explicitly. This is safer inside your btergm package
  # than relying on unqualified internal ergm functions being visible.
  #check_dots_used <- getFromNamespace("check_dots_used", "statnet.common")
  #unused_dots_warning <- getFromNamespace("unused_dots_warning", "statnet.common")
  #check.control.class <- getFromNamespace("check.control.class", "ergm")
  #handle.control.toplevel <- getFromNamespace("handle.control.toplevel", "ergm")
  nonsimp_update.formula <- .btergm_nonsimp_update_formula()
  ergm_model <- getFromNamespace("ergm_model", "ergm")
  handle_auto_constraints <- getFromNamespace(".handle.auto.constraints", "ergm")
  
  prune.ergm_conlist <- getFromNamespace("prune.ergm_conlist", "ergm")
  as.rlebdm <- getFromNamespace("as.rlebdm", "ergm")
  ergm.pl <- getFromNamespace("ergm.pl", "ergm")
  ergm_state <- getFromNamespace("ergm_state", "ergm")
  
  #make_ergm_conlist <- function(x, nw, term.options) {
  #  eval_env <- new.env(parent = asNamespace("ergm"))
  #  assign("x", x, envir = eval_env)
  #  assign("nw", nw, envir = eval_env)
  #  assign("term.options", term.options, envir = eval_env)
  #  
  #  eval(
  #    quote(ergm_conlist(x, nw, term.options = term.options)),
  #    envir = eval_env
  #  )
  #}
  
  if (output == "array") {
    arr_from_coo <- .btergm_arr_from_coo()
  }
  
  #check_dots_used(error = unused_dots_warning)
  #check.control.class("ergm", "ergmMPLE")
  #handle.control.toplevel("ergm", ...)
  
  nw <- basis
  
  # If sampled, force dyad indices internally so we can align IPW weights.
  # If the user requested matrix output, we will drop tail/head before returning.
  need_indices <- output %in% c("array", "dyadlist") || sample.method != "none"
  
  if (need_indices) {
    formula <- nonsimp_update.formula(formula, . ~ indices + .)
  }
  
  model <- ergm_model(
    formula,
    nw,
    ...,
    term.options = control$term.options
  )
  
  tmp <- handle_auto_constraints(nw, constraints, obs.constraints)
  nw <- tmp$nw
  conterms <- tmp$conterms
  conterms.obs <- tmp$conterms.obs

  if ("constraints" %in% names(control$MCMC.prop.args)) {
    conlist <- prune.ergm_conlist(control$MCMC.prop.args$constraints)
    class(conlist) <- "ergm_conlist"
  } else {
    conlist <- .btergm_make_ergm_conlist(
      conterms,
      nw,
      term.options = control$term.options
    )
  }
  
  if ("constraints" %in% names(control$obs.MCMC.prop.args)) {
    conlist.obs <- prune.ergm_conlist(control$obs.MCMC.prop.args$constraints)
    class(conlist.obs) <- "ergm_conlist"
  } else {
    conlist.obs <- .btergm_make_ergm_conlist(
      conterms.obs,
      nw,
      term.options = control$term.options
    )
  }
  
  fd <- as.rlebdm(conlist, conlist.obs, which = "informative")

  sample_info <- sample_mple_fd(
    nw = nw,
    fd = fd,
    sample.size = sample.size,
    sample.method = sample.method,
    seed = seed,
    WtSumAsSampSiz = WtSumAsSampSiz,
    as.rlebdm_fun = as.rlebdm
  )
  
  pl <- ergm.pl(
    ergm_state(nw, model = model),
    sample_info$fd,
    verbose = verbose,
    control = control,
    ignore.offset = TRUE,
    ...
  )
  
  # Apply inverse-probability weights once.
  if (isTRUE(sample_info$sampled)) {
    if (!all(c("tail", "head") %in% colnames(pl$xmat.full))) {
      stop(
        "Sampled MPLE requires 'tail' and 'head' columns in pl$xmat.full, ",
        "but they were not found. Check that indices were added to the formula."
      )
    }
    key <- paste(
      pl$xmat.full[, "tail"],
      pl$xmat.full[, "head"],
      sep = "->"
    )
    ipw <- sample_info$ipw[key]
    if (anyNA(ipw)) {
      stop(
        "Could not align inverse-probability weights with ergm.pl() rows. ",
        "This usually means the dyad keys returned by sample_mple_fd() do not ",
        "match the tail/head columns returned by ergm.pl()."
      )
    }
    pl$wend <- pl$wend * unname(ipw)
  }
  
  out <- switch(
    output,
    
    matrix = {
      predictor <- pl$xmat.full
      
      # Drop dyad ID columns so the output matches ordinary ergm::ergmMPLE().
      predictor <- predictor[
        ,
        !colnames(predictor) %in% c("tail", "head"),
        drop = FALSE
      ]
      
      list(
        response = pl$zy,
        predictor = predictor,
        weights = pl$wend
      )
    },
    
    dyadlist = {
      o <- order(pl$xmat.full[, "tail"], pl$xmat.full[, "head"])
      
      list(
        response = pl$zy[o],
        predictor = pl$xmat.full[o, , drop = FALSE],
        weights = pl$wend[o]
      )
    },
    
    array = {
      bip <- if (!expand.bipartite) ergm::b1.size(nw) else 0
      
      vertex_names <- network::get.vertex.attribute(nw, "vertex.names")
      
      vn <- if (is.null(vertex_names) || all(is.na(vertex_names))) {
        seq_len(network::network.size(nw))
      } else {
        vertex_names
      }
      
      t_names <- if (bip) vn[seq_len(bip)] else vn
      h_names <- if (bip) vn[-seq_len(bip)] else vn
      
      dyads <- pl$xmat.full[, 1:2, drop = FALSE]
      storage.mode(dyads) <- "integer"
      
      terms <- colnames(pl$xmat.full)[-(1:2)]
      dn <- list(tail = t_names, head = h_names)
      
      if (bip) {
        dyads[, 2] <- dyads[, 2] - bip
      }
      
      xa <- arr_from_coo(
        pl$xmat.full[, -(1:2), drop = FALSE],
        cbind(
          dyads[, 1],
          dyads[, 2],
          rep(seq_along(terms), each = nrow(dyads))
        ),
        dimnames = c(dn, list(term = terms))
      )
      
      ym <- arr_from_coo(pl$zy, dyads, dimnames = dn)
      wm <- arr_from_coo(pl$wend, dyads, x0 = 0, dimnames = dn)
      
      list(
        response = ym,
        predictor = xa,
        weights = wm
      )
    }
  )
  
  eta <- model$etamap
  
  # For matrix output, btergm_ergmMPLE_sampled() drops the internally added
  # tail/head columns before returning pl$predictor. The etamap must be
  # realigned so eta$offsettheta has the same length/order as out$predictor.
  if (output == "matrix" && !is.null(eta$offsettheta)) {
    
    old_offsettheta <- eta$offsettheta
    pred_names <- colnames(out$predictor)
    
    if (!is.null(names(old_offsettheta)) &&
        any(names(old_offsettheta) %in% pred_names)) {
      
      # Preferred path: map offset flags by column name.
      offset_names <- names(old_offsettheta)[old_offsettheta]
      eta$offsettheta <- pred_names %in% offset_names
      
    } else if (length(old_offsettheta) == length(pred_names)) {
      
      # Already aligned.
      eta$offsettheta <- old_offsettheta
      
    } else if (!is.null(colnames(pl$xmat.full)) &&
               length(old_offsettheta) == ncol(pl$xmat.full)) {
      
      # Common sampled-MPLE case: old etamap matches pl$xmat.full,
      # but matrix output dropped tail/head.
      keep <- !colnames(pl$xmat.full) %in% c("tail", "head")
      eta$offsettheta <- old_offsettheta[keep]
      
    } else {
      
      stop(
        "Internal sampled MPLE error: etamap$offsettheta has length ",
        length(old_offsettheta),
        ", but the returned predictor matrix has ",
        length(pred_names),
        " columns. This usually means etamap was not realigned after ",
        "dropping the internally added tail/head columns."
      )
    }
  }
  
  structure(out, etamap = eta)
}
