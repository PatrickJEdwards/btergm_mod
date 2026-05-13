# Register S3 list class for S4 method dispatch.
methods::setOldClass("btergm_count")
methods::setOldClass("summary.btergm_count")

#' Methods for btergm_count objects
#'
#' Basic methods for extracting coefficients, printing fitted objects,
#' computing observation counts, summarizing fitted models, and constructing
#' bootstrap confidence intervals for count-valued bootstrapped TERGM objects.
#'
#' @name btergm_count_methods
#' @aliases coef.btergm_count
#' @aliases coef.ergmCntMPLE
#' @aliases print.btergm_count
#' @aliases nobs.btergm_count
#' @aliases confint.btergm_count
#' @aliases summary.btergm_count
#' @aliases print.summary.btergm_count
#' @aliases coef,btergm_count-method
#' @aliases confint,btergm_count-method
#' @aliases nobs,btergm_count-method
#' @aliases summary,btergm_count-method
#' @keywords internal
NULL




btergm_count_coef <- function(object) {
  if (!is.null(object$coef)) {
    return(object$coef)
  }
  
  if (!is.null(object$coefficients)) {
    return(object$coefficients)
  }
  
  stop("Could not find coefficients in this btergm_count object.")
}

#' @rdname btergm_count_methods
#' @method coef btergm_count
#' @export
coef.btergm_count <- function(object, ...) {
  btergm_count_coef(object)
}
methods::setMethod(f = "coef", signature = "btergm_count", definition = function(object, ...) {
    btergm_count_coef(object)
  })

#' @rdname btergm_count_methods
#' @method print btergm_count
#' @export
print.btergm_count <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCount-valued bootstrapped TERGM by MPLE\n")
  cat("\nFormula:   ")
  print(x$formula)
  
  cat("\nResponse variable:", x$response, "\n")
  cat("Reference measure:", x$count_control$reference, "\n")
  cat("Time steps:", x$time.steps, "\n")
  cat("Bootstrap replications:", x$R, "\n")
  cat("Sampled active edge-variable contributions:", x$nobs, "\n")
  
  cat("\nMPLE Coefficients:\n")
  print.default(
    format(x$coef, digits = digits),
    print.gap = 2,
    quote = FALSE
  )
  
  invisible(x)
}

#' @rdname btergm_count_methods
#' @method nobs btergm_count
#' @export
nobs.btergm_count <- function(object, ...) {
  c(
    "Number of time steps" = object$time.steps,
    "Sampled active edge-variable contributions" = object$nobs,
    "Bootstrap replications" = object$R
  )
}
methods::setMethod(f = "nobs", signature = "btergm_count", definition = function(object, ...) {
    nobs.btergm_count(object, ...)
  })

#' @rdname btergm_count_methods
#' @method confint btergm_count
#' @export
confint.btergm_count <- function(object, parm, level = 0.95, type = c("perc", "basic", "norm"), ...) {
  type <- match.arg(type)
  
  cf <- btergm_count_coef(object)
  pnames <- names(cf)
  
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  
  if (object$R == 0L || is.null(object$boot$t) || nrow(object$boot$t) == 0L) {
    stop("No bootstrap replications are available. Refit with R > 0.")
  }
  
  bt <- object$boot$t
  
  if (is.null(colnames(bt))) {
    colnames(bt) <- pnames
  }
  
  keep <- stats::complete.cases(bt)
  bt0 <- bt
  bt <- bt[keep, , drop = FALSE]
  
  if (nrow(bt) < nrow(bt0)) {
    warning(
      nrow(bt0) - nrow(bt),
      " bootstrap replication(s) dropped because of missing or failed estimates."
    )
  }
  
  alpha <- 1 - level
  
  ci <- matrix(NA_real_, nrow = length(pnames), ncol = 4)
  rownames(ci) <- pnames
  colnames(ci) <- c(
    "Estimate",
    "Boot mean",
    paste0(100 * alpha / 2, "%"),
    paste0(100 * (1 - alpha / 2), "%")
  )
  
  ci[, "Estimate"] <- cf
  ci[, "Boot mean"] <- colMeans(bt, na.rm = TRUE)
  
  for (jj in seq_along(pnames)) {
    b <- bt[, jj]
    
    if (type == "perc") {
      ci[jj, 3:4] <- stats::quantile(
        b,
        probs = c(alpha / 2, 1 - alpha / 2),
        na.rm = TRUE,
        names = FALSE
      )
    } else if (type == "basic") {
      q <- stats::quantile(
        b,
        probs = c(1 - alpha / 2, alpha / 2),
        na.rm = TRUE,
        names = FALSE
      )
      ci[jj, 3:4] <- 2 * cf[jj] - q
    } else if (type == "norm") {
      se <- stats::sd(b, na.rm = TRUE)
      z <- stats::qnorm(1 - alpha / 2)
      ci[jj, 3:4] <- c(cf[jj] - z * se, cf[jj] + z * se)
    }
  }
  
  ci[parm, , drop = FALSE]
}
methods::setMethod(f = "confint", signature = "btergm_count", definition = function(object, parm, level = 0.95, type = c("perc", "basic", "norm"),...) {
    confint.btergm_count(
      object = object,
      parm = parm,
      level = level,
      type = type,
      ...
    )
  })

#' @rdname btergm_count_methods
#' @method summary btergm_count
#' @export
summary.btergm_count <- function(object, level = 0.95, type = c("perc", "basic", "norm"), include.mple.se = TRUE, ...) {
  type <- match.arg(type)
  
  out <- list()
  out$call <- object$call
  out$formula <- object$formula
  out$response <- object$response
  out$reference <- object$count_control$reference
  out$regularization <- object$count_control$regularization
  out$regularization.param <- object$count_control$regularization.param
  out$time.steps <- object$time.steps
  out$R <- object$R
  out$nobs <- object$nobs
  out$level <- level
  out$type <- type
  out$coef <- btergm_count_coef(object)
  
  out$ci <- confint(object, level = level, type = type, ...)
  
  bt <- object$boot$t
  if (!is.null(bt) && nrow(bt) > 0L) {
    cf <- btergm_count_coef(object)
    colnames(bt) <- names(cf)
    out$boot.se <- apply(bt, 2, stats::sd, na.rm = TRUE)
  } else {
    out$boot.se <- rep(NA_real_, length(object$coef))
    names(out$boot.se) <- names(object$coef)
  }
  
  if (isTRUE(include.mple.se) &&
      !is.null(object$fit) &&
      !is.null(object$fit$cov)) {
    out$mple.se <- sqrt(diag(object$fit$cov))
  } else {
    out$mple.se <- NULL
  }
  
  class(out) <- "summary.btergm_count"
  out
}
methods::setMethod(f = "summary", signature = "btergm_count", definition = function(object, level = 0.95, type = c("perc", "basic", "norm"), include.mple.se = TRUE, ...) {
    summary.btergm_count(
      object = object,
      level = level,
      type = type,
      include.mple.se = include.mple.se,
      ...
    )
  })

#' @rdname btergm_count_methods
#' @method print summary.btergm_count
#' @export
print.summary.btergm_count <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n==============================\n")
  cat("Summary of count TERGM fit\n")
  cat("==============================\n\n")
  
  cat("Formula:   ")
  print(x$formula)
  cat("\n")
  
  cat("Response variable:", x$response, "\n")
  cat("Reference measure:", x$reference, "\n")
  cat("Time steps:", x$time.steps, "\n")
  cat("Bootstrap replications:", x$R, "\n")
  cat("Sampled active edge-variable contributions:", x$nobs, "\n")
  
  cat("\nRegularization:", x$regularization)
  if (!is.null(x$regularization) && x$regularization != "none") {
    cat(
      " (regularization parameter ",
      paste(x$regularization.param, collapse = ", "),
      ")\n",
      sep = ""
    )
  } else {
    cat("\n")
  }
  
  tab <- x$ci
  
  if (!is.null(x$boot.se)) {
    tab <- cbind(
      tab[, "Estimate", drop = FALSE],
      "Boot Std.Err" = x$boot.se[rownames(tab)],
      tab[, setdiff(colnames(tab), "Estimate"), drop = FALSE]
    )
  }
  
  if (!is.null(x$mple.se)) {
    common <- intersect(rownames(tab), names(x$mple.se))
    mple_col <- rep(NA_real_, nrow(tab))
    names(mple_col) <- rownames(tab)
    mple_col[common] <- x$mple.se[common]
    
    tab <- cbind(
      tab[, "Estimate", drop = FALSE],
      "MPLE Std.Err" = mple_col,
      tab[, setdiff(colnames(tab), "Estimate"), drop = FALSE]
    )
  }
  
  cat(
    "\nMPLE estimates with ",
    100 * x$level,
    "% temporal bootstrap confidence intervals:\n",
    sep = ""
  )
  
  printCoefmat(
    tab,
    digits = digits,
    signif.stars = FALSE,
    has.Pvalue = FALSE,
    P.values = FALSE,
    na.print = "NA",
    ...
  )
  
  if (!is.null(x$mple.se)) {
    cat(
      "\nNote: MPLE Std.Err is the Hessian-based standard error from the ",
      "full-sample count-MPLE fit. Boot Std.Err and confidence intervals ",
      "come from the temporal bootstrap.\n",
      sep = ""
    )
  }
  
  invisible(x)
}


btergm_count_clean_coef_names <- function(parnam, prepared) {
  if (is.null(parnam)) {
    return(parnam)
  }
  
  rhs_terms <- prepared$rhs.terms
  temporal_labels <- prepared$count_temporal_labels %||% list()
  
  clean <- parnam
  
  rhs_labels <- vapply(
    rhs_terms,
    btergm_count_rhs_label,
    character(1),
    temporal_labels = temporal_labels
  )
  
  # Best case: one parameter per RHS term.
  # This is true for your current model:
  # sum, nonzero, nodeofactor(...), nodeifactor(...), edgecov(memory), edgecov(delrecip)
  if (length(rhs_labels) == length(clean)) {
    replace <- !is.na(rhs_labels) & nzchar(rhs_labels)
    clean[replace] <- rhs_labels[replace]
    return(make.unique(clean, sep = "_"))
  }
  
  # Fallback: replace duplicated generic edgecov names in order.
  edge_labels <- rhs_labels[!is.na(rhs_labels) & grepl("^(memory|delayed_reciprocity|edgecov)", rhs_labels)]
  
  edge_idx <- which(
    grepl("^edgecov", clean) &
      (duplicated(clean) | duplicated(clean, fromLast = TRUE) | clean %in% c("edgecov.sum.NULL", "edgecov.NULL"))
  )
  
  if (length(edge_idx) == length(edge_labels)) {
    clean[edge_idx] <- edge_labels
  }
  
  make.unique(clean, sep = "_")
}


btergm_count_rhs_label <- function(term, temporal_labels = list()) {
  x <- gsub("\\s+", "", term)
  
  if (grepl("^edgecov\\(memory\\[\\[i\\]\\]\\)", x)) {
    return(temporal_labels$memory %||% "memory")
  }
  
  if (grepl("^edgecov\\(delrecip\\[\\[i\\]\\]\\)", x)) {
    return(temporal_labels$delrecip %||% "delayed_reciprocity")
  }
  
  # Optional: label other edgecov terms more cleanly if ergm gives a NULL name.
  if (grepl("^edgecov\\(", x)) {
    inside <- sub("^edgecov\\((.*)\\)$", "\\1", x)
    inside <- sub("\\[\\[i\\]\\]", "", inside)
    inside <- sub("\\[i\\]", "", inside)
    return(paste0("edgecov.", inside))
  }
  
  NA_character_
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}