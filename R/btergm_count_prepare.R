#' Prepare data for bootstrapped count-valued TERGM MPLE
#'
#' Internal helper for btergm_count().
#'
#' This function reuses tergmprepare() for temporal formula handling,
#' covariate indexing, and dimension adjustment, then adds count-specific
#' handling for user-supplied offset terms and temporal dependence covariates.
#'
#' @param formula Original temporal ERGM formula.
#' @param response Character string naming the count edge attribute.
#' @param offset Logical. Passed to tergmprepare().
#' @param blockdiag Logical. Currently should be FALSE for btergm_count().
#' @param verbose Logical. Print preparation messages?
#' @param temporal_covariates Character. "count" keeps lagged counts for
#'   autoregression and delayed reciprocity; "binary" uses lagged nonzero
#'   indicators.
#'
#' @return A prepared list like tergmprepare(), with additional elements:
#'   count_offset, count_offset_terms, and response.
#'
#' @keywords internal
#' 
#' @importFrom network is.network list.edge.attributes set.edge.attribute set.vertex.attribute as.edgelist network



btergm_count_prepare <- function(formula,
                                 response = "count",
                                 offset = FALSE,
                                 blockdiag = FALSE,
                                 verbose = TRUE,
                                 temporal_covariates = c("count", "binary")) {
  
  temporal_covariates <- match.arg(temporal_covariates)
  
  if (blockdiag) {
    stop("btergm_count_prepare() currently supports blockdiag = FALSE only.")
  }
  
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.")
  }
  
  if (!is.character(response) || length(response) != 1L || is.na(response)) {
    stop("'response' must be a single character string.")
  }
  
  # 1. Reuse existing btergm preparation machinery.
  #l <- tergmprepare(formula = formula, offset = offset, blockdiag = blockdiag, verbose = verbose)
  l <- tergmprepare_count(
    formula = formula,
    response = response,
    offset = offset,
    blockdiag = blockdiag,
    verbose = verbose
  )
  
  l$response <- response
  
  # 2. Build count-specific offset matrices from all user-supplied offset(...)
  #    terms, then remove those terms from the formula used by ergmCntMPLE().
  l <- btergm_count_extract_offsets(
    prepared = l,
    formula = formula,
    offset = offset,
    verbose = verbose
  )
  
  # 3. Rebuild temporal covariates in count-aware form where needed.
  l <- btergm_count_rebuild_temporal_covariates(
    prepared = l,
    formula = formula,
    response = response,
    temporal_covariates = temporal_covariates,
    verbose = verbose
  )
  
  # 4. Make sure the dependent networks are usable by ergmCntMPLE().
  #    ergmCntMPLE() expects a network object with a response edge attribute.
  l$networks <- lapply(
    l$networks,
    btergm_count_ensure_network_response,
    response = response,
    directed = l$directed,
    bipartite = l$bipartite
  )
  
  l
}





## Helper for extracting count matrices:
## ergmCntPrep() reads the outcome using as.sociomatrix(nw, response), not the ordinary binary adjacency matrix
btergm_count_response_matrix <- function(x, response = "count") {
  if (network::is.network(x)) {
    edge_attrs <- network::list.edge.attributes(x)
    
    if (!response %in% edge_attrs) {
      stop(
        "The network does not contain edge attribute '",
        response,
        "'. Available edge attributes are: ",
        paste(edge_attrs, collapse = ", "),
        "."
      )
    }
    
    out <- try(
      as.matrix(
        x,
        matrix.type = "adjacency",
        attrname = response
      ),
      silent = TRUE
    )
    
    if (inherits(out, "try-error")) {
      stop(
        "Could not extract count response edge attribute '",
        response,
        "' using as.matrix(..., attrname = response)."
      )
    }
  } else {
    out <- as.matrix(x)
  }
  
  storage.mode(out) <- "numeric"
  out[is.na(out)] <- 0
  out
}


## Helper for converting matrices back to count networks:
## Lets btergm_count_prepare() return networks that ergmCntMPLE() can use.
btergm_count_ensure_network_response <- function(x, response = "count", directed = TRUE, bipartite = FALSE) {
  if (network::is.network(x)) {
    edge_attrs <- network::list.edge.attributes(x)
    
    if (response %in% edge_attrs) {
      return(x)
    }
    
    # If the network lacks the response attribute, fall through and rebuild it
    # from its matrix representation.
    mat <- as.matrix(x)
  } else {
    mat <- as.matrix(x)
  }
  
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  nw <- network::network(
    mat != 0,
    directed = directed,
    bipartite = if (isTRUE(bipartite)) nrow(mat) else FALSE,
    loops = FALSE,
    matrix.type = "adjacency"
  )
  
  if (!is.null(rownames(mat))) {
    network::set.vertex.attribute(nw, "vertex.names", rownames(mat))
  }
  
  el <- network::as.edgelist(nw)
  
  if (nrow(el) > 0) {
    vals <- mat[cbind(el[, 1], el[, 2])]
    network::set.edge.attribute(nw, response, vals)
  }
  
  nw
}



## Helpers for user-supplied offset terms:
## In binary btergm(), offset columns are detected after ergmMPLE() creates the MPLE predictor matrix. 
## The code then removes offset columns and excludes any dyad marked by an offset term.
## For count MPLE, we need to do this before calling ergmCntMPLE() or ergmCntPrep().
btergm_count_extract_offsets <- function(prepared, formula, offset = FALSE, verbose = TRUE) {
  l <- prepared
  
  offset_idx <- grepl("^\\s*offset\\s*\\(", l$rhs.terms)
  
  l$count_offset_terms <- l$rhs.terms[offset_idx]
  
  # Start with all dyads active.
  count_offset <- vector("list", l$time.steps)
  
  for (tt in seq_len(l$time.steps)) {
    base_mat <- matrix(
      0,
      nrow = nrow(as.matrix(l$networks[[tt]])),
      ncol = ncol(as.matrix(l$networks[[tt]]))
    )
    
    rownames(base_mat) <- rownames(as.matrix(l$networks[[tt]]))
    colnames(base_mat) <- colnames(as.matrix(l$networks[[tt]]))
    
    count_offset[[tt]] <- base_mat
  }
  
  # Include automatic structural-zero offsets from tergmprepare() if requested.
  # In count MPLE, these should not remain in the formula. They should become
  # dyad exclusions.
  if (isTRUE(offset)) {
    for (tt in seq_len(l$time.steps)) {
      count_offset[[tt]][as.matrix(l$offsmat[[tt]]) != 0] <- 1
    }
  }
  
  # Evaluate all user-supplied offset(edgecov(...)) or offset(dyadcov(...))
  # terms after tergmprepare() has indexed them by [[i]].
  if (length(l$count_offset_terms) > 0L) {
    eval_env <- new.env(parent = environment(formula))
    
    for (nm in l$covnames) {
      assign(nm, l[[nm]], envir = eval_env)
    }
    
    assign("offsmat", l$offsmat, envir = eval_env)
    
    for (tt in seq_len(l$time.steps)) {
      assign("i", tt, envir = eval_env)
      
      for (term in l$count_offset_terms) {
        cov_expr <- btergm_count_extract_offset_covariate_expression(term)
        
        mat <- eval(parse(text = cov_expr), envir = eval_env)
        mat <- as.matrix(mat)
        
        mat <- btergm_count_align_matrix(
          mat = mat,
          target = count_offset[[tt]],
          term = term,
          time = tt
        )
        
        count_offset[[tt]][mat != 0 | is.na(mat)] <- 1
      }
    }
  }
  
  l$count_offset <- count_offset
  
  # Remove offset(...) terms from the formula that will be passed to ergmCntMPLE().
  if (any(offset_idx)) {
    l$rhs.terms <- l$rhs.terms[!offset_idx]
  }
  
  if (length(l$rhs.terms) == 0L) {
    stop(
      "After removing offset(...) terms, the model has no RHS terms. ",
      "Add at least one count ERGM term, such as sum or nonzero."
    )
  }
  
  ## OLD CODE (1):
  #rhs <- paste(l$rhs.terms, collapse = " + ")
  #lhs <- deparse(l$form[[2]])
  #l$form <- stats::as.formula(paste(lhs, "~", rhs), env = environment(l$form))
  
  ## NEW CODE (2):
  rhs <- paste(l$rhs.terms, collapse = " + ")
  # tergmprepare() returns l$form as a character string, not a formula object.
  # The prepared LHS is always networks[[i]] when blockdiag = FALSE.
  l$form <- paste("networks[[i]] ~", rhs)
  
  if (verbose && length(l$count_offset_terms) > 0L) {
    message(
      "\nRemoved ",
      length(l$count_offset_terms),
      " user-supplied offset term(s) from the count ERGM formula ",
      "and stored them as count_offset matrices."
    )
  }
  
  l
}

## Parsing helper:
btergm_count_extract_offset_covariate_expression <- function(term) {
  # Example input:
  # offset(edgecov(receiver_ineligible_mat_list[[i]]))
  # offset(dyadcov(receiver_ineligible_mat_list[[i]]))
  
  x <- trimws(term)
  
  x <- sub("^offset\\s*\\((.*)\\)\\s*$", "\\1", x, perl = TRUE)
  x <- trimws(x)
  
  x <- sub("^edgecov\\s*\\((.*)\\)\\s*$", "\\1", x, perl = TRUE)
  x <- sub("^dyadcov\\s*\\((.*)\\)\\s*$", "\\1", x, perl = TRUE)
  
  # Keep only the first argument to edgecov()/dyadcov().
  # This handles edgecov(x, attr = "...") without letting the attr argument
  # become part of the object expression.
  btergm_count_first_top_level_arg(x)
}

## First-argument helper:
btergm_count_first_top_level_arg <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  
  depth_paren <- 0L
  depth_square <- 0L
  
  for (ii in seq_along(chars)) {
    ch <- chars[[ii]]
    
    if (ch == "(") depth_paren <- depth_paren + 1L
    if (ch == ")") depth_paren <- depth_paren - 1L
    if (ch == "[") depth_square <- depth_square + 1L
    if (ch == "]") depth_square <- depth_square - 1L
    
    if (ch == "," && depth_paren == 0L && depth_square == 0L) {
      return(trimws(substr(x, 1L, ii - 1L)))
    }
  }
  
  trimws(x)
}

## Matrix-alignment helper:
btergm_count_align_matrix <- function(mat, target, term = NULL, time = NULL) {
  mat <- as.matrix(mat)
  target <- as.matrix(target)
  
  if (identical(dim(mat), dim(target))) {
    if (!is.null(rownames(target)) && !is.null(rownames(mat))) {
      mat <- mat[rownames(target), , drop = FALSE]
    }
    
    if (!is.null(colnames(target)) && !is.null(colnames(mat))) {
      mat <- mat[, colnames(target), drop = FALSE]
    }
    
    return(mat)
  }
  
  if (!is.null(rownames(mat)) && !is.null(colnames(mat)) &&
      !is.null(rownames(target)) && !is.null(colnames(target))) {
    return(mat[rownames(target), colnames(target), drop = FALSE])
  }
  
  stop(
    "Offset matrix is not conformable",
    if (!is.null(term)) paste0(" for term: ", term) else "",
    if (!is.null(time)) paste0(" at time step ", time) else "",
    ". Provide row and column names so it can be aligned."
  )
}



## TEMPORAL COVARIATES:
#### `memory(type = "autoregression", lag = 1)`:
####    This is the cleanest count-valued temporal dependence term. It uses the 
####    previous count as a dyadic covariate. A positive coefficient means dyads 
####    with larger prior counts tend to have larger current counts.
#### `delrecip(lag = 1)`: 
####    This is also count-compatible. With temporal_covariates = "count", it 
####    uses the previous reciprocal count. A positive coefficient means prior 
####    support from j to i predicts current support from i to j.
#### `memory(type = "stability")`
#### `memory(type = "innovation")` 
#### `memory(type = "loss")`
####    These are basically binary presence/absence concepts. In the code above, 
####    I deliberately make them depend on whether the prior count is nonzero. I 
####    would not use raw counts for these unless you write a new theoretical definition.
#### `timecov()`
####    Usually this is fine, because it becomes an edgecov() time trend. But if 
####    you use timecov() in a way that transforms the dependent network itself, 
####    we should make it response-aware later. For the first implementation, I 
####    would avoid dependent-network timecov() and instead use explicit time 
####    covariate matrices when needed.

## Count-aware temporal covariate rebuilding:
btergm_count_rebuild_temporal_covariates <- function(
    prepared, formula, response = "count", temporal_covariates = c("count", "binary"), verbose = TRUE) {
  
  temporal_covariates <- match.arg(temporal_covariates)
  
  l <- prepared
  
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  
  rhs_terms <- btergm_count_split_rhs(rhs)
  
  memory_terms <- rhs_terms[grepl("^\\s*memory\\s*\\(", rhs_terms)]
  delrecip_terms <- rhs_terms[grepl("^\\s*delrecip\\s*\\(", rhs_terms)]
  
  if (length(memory_terms) == 0L && length(delrecip_terms) == 0L) {
    return(l)
  }
  
  original_networks <- eval(formula[[2]], envir = environment(formula))
  
  if (!("list" %in% class(original_networks)) &&
      !("network.list" %in% class(original_networks))) {
    original_networks <- list(original_networks)
  }
  
  original_count_mats <- lapply(
    original_networks,
    btergm_count_response_matrix,
    response = response
  )
  
  max_lag <- length(original_count_mats) - l$time.steps
  
  if (max_lag < 0L) {
    stop("Could not align temporal covariates with prepared networks.")
  }
  
  # Existing tergmprepare() stores only one object called memory.
  # We follow that convention here.
  if (length(memory_terms) > 0L && "memory" %in% l$covnames) {
    mem_term <- memory_terms[[1]]
    
    mem_type <- btergm_count_extract_character_arg(
      term = mem_term,
      arg = "type",
      default = "stability"
    )
    
    mem_lag <- btergm_count_extract_integer_arg(
      term = mem_term,
      arg = "lag",
      default = 1L
    )
    
    memory <- vector("list", l$time.steps)
    
    for (tt in seq_len(l$time.steps)) {
      source_time <- max_lag + tt - mem_lag
      
      if (source_time < 1L) {
        stop("Could not align memory() term at prepared time step ", tt, ".")
      }
      
      lag_mat <- original_count_mats[[source_time]]
      
      lag_mat <- btergm_count_align_matrix(
        mat = lag_mat,
        target = as.matrix(l$networks[[tt]]),
        term = mem_term,
        time = tt
      )
      
      memory[[tt]] <- btergm_count_make_memory_matrix(
        lag_mat = lag_mat,
        type = mem_type,
        temporal_covariates = temporal_covariates
      )
    }
    
    l[["memory"]] <- memory
    
    if (verbose) {
      message(
        "\nRebuilt memory() as a count-aware temporal covariate using type = '",
        mem_type,
        "' and lag = ",
        mem_lag,
        "."
      )
    }
    
    if (is.null(l$count_temporal_labels)) {
      l$count_temporal_labels <- list()
    }
    
    l$count_temporal_labels$memory <- paste0(
      "memory_", mem_type, "_lag", mem_lag
    )
  }
  
  if (length(delrecip_terms) > 0L && "delrecip" %in% l$covnames) {
    del_term <- delrecip_terms[[1]]
    
    del_lag <- btergm_count_extract_integer_arg(
      term = del_term,
      arg = "lag",
      default = 1L
    )
    
    mutuality <- btergm_count_extract_logical_arg(
      term = del_term,
      arg = "mutuality",
      default = FALSE
    )
    
    delrecip <- vector("list", l$time.steps)
    
    for (tt in seq_len(l$time.steps)) {
      source_time <- max_lag + tt - del_lag
      
      if (source_time < 1L) {
        stop("Could not align delrecip() term at prepared time step ", tt, ".")
      }
      
      lag_mat <- original_count_mats[[source_time]]
      
      lag_mat <- btergm_count_align_matrix(
        mat = lag_mat,
        target = as.matrix(l$networks[[tt]]),
        term = del_term,
        time = tt
      )
      
      lag_mat <- t(lag_mat)
      
      if (isTRUE(mutuality)) {
        delrecip[[tt]] <- ifelse(lag_mat > 0, 1, -1)
      } else if (temporal_covariates == "count") {
        delrecip[[tt]] <- lag_mat
      } else {
        delrecip[[tt]] <- 1 * (lag_mat > 0)
      }
    }
    
    l[["delrecip"]] <- delrecip
    
    if (verbose) {
      message(
        "\nRebuilt delrecip() as a count-aware temporal covariate using lag = ",
        del_lag,
        "."
      )
    }
    
    if (is.null(l$count_temporal_labels)) {
      l$count_temporal_labels <- list()
    }
    
    l$count_temporal_labels$delrecip <- if (isTRUE(mutuality)) {
      paste0("delayed_reciprocity_mutuality_lag", del_lag)
    } else {
      paste0("delayed_reciprocity_lag", del_lag)
    }
    
  }
  
  l
}


## Temporal parsing helpers:
btergm_count_make_memory_matrix <- function(lag_mat, type = "stability", temporal_covariates = c("count", "binary")) {
  temporal_covariates <- match.arg(temporal_covariates)
  
  lag_mat <- as.matrix(lag_mat)
  storage.mode(lag_mat) <- "numeric"
  lag_mat[is.na(lag_mat)] <- 0
  
  if (type == "autoregression") {
    if (temporal_covariates == "count") {
      return(lag_mat)
    } else {
      return(1 * (lag_mat > 0))
    }
  }
  
  if (type == "stability") {
    return(ifelse(lag_mat > 0, 1, -1))
  }
  
  if (type == "innovation") {
    return(ifelse(lag_mat == 0, 1, 0))
  }
  
  if (type == "loss") {
    return(ifelse(lag_mat > 0, -1, 0))
  }
  
  stop("'type' argument in memory() not recognized.")
}

btergm_count_extract_integer_arg <- function(term, arg, default) {
  pattern <- paste0(arg, "\\s*=\\s*([0-9]+)")
  
  if (!grepl(pattern, term, perl = TRUE)) {
    return(as.integer(default))
  }
  
  as.integer(sub(paste0(".*", pattern, ".*"), "\\1", term, perl = TRUE))
}

btergm_count_extract_logical_arg <- function(term, arg, default) {
  pattern <- paste0(arg, "\\s*=\\s*(TRUE|FALSE|T|F)")
  
  if (!grepl(pattern, term, perl = TRUE)) {
    return(default)
  }
  
  val <- sub(paste0(".*", pattern, ".*"), "\\1", term, perl = TRUE)
  as.logical(val)
}

btergm_count_extract_character_arg <- function(term, arg, default) {
  pattern <- paste0(arg, "\\s*=\\s*['\"]([^'\"]+)['\"]")
  
  if (!grepl(pattern, term, perl = TRUE)) {
    return(default)
  }
  
  sub(paste0(".*", pattern, ".*"), "\\1", term, perl = TRUE)
}

btergm_count_split_rhs <- function(rhs) {
  rhsterms <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
  
  if (length(rhsterms) > 1L) {
    for (ii in length(rhsterms):2L) {
      left <- gregexpr("\\(", rhsterms[ii])[[1]]
      left <- left[left != -1L]
      left <- length(left)
      
      right <- gregexpr("\\)", rhsterms[ii])[[1]]
      right <- right[right != -1L]
      right <- length(right)
      
      if (left != right) {
        rhsterms[ii - 1L] <- paste(rhsterms[ii - 1L], rhsterms[ii], sep = " + ")
        rhsterms <- rhsterms[-ii]
      }
    }
  }
  
  trimws(rhsterms)
}