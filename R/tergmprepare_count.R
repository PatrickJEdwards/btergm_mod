tergmprepare_count_dim_matrix <- function(obj, name, n_periods = NULL) {
  if (is.null(obj)) {stop("Object '", name, "' is NULL in tergmprepare_count().")}
  
  if (network::is.network(obj) || is.matrix(obj)) {
    obj <- list(obj)
  }
  
  if (!is.list(obj)) {
    obj <- list(obj)
  }
  
  if (!is.null(n_periods) && length(obj) < n_periods) {stop("Object '", name, "' has length ", length(obj), ", but tergmprepare_count() expects ", n_periods, " time steps.")}
  
  dims <- lapply(seq_along(obj), function(tt) {
    z <- obj[[tt]]
    
    if (is.null(z)) {
      stop(
        "Object '", name, "' contains NULL at time index ", tt,
        ". This usually means a temporal covariate list was not trimmed or ",
        "constructed correctly after memory(), delrecip(), or an offset term."
      )
    }
    
    dz <- dim(as.matrix(z))
    
    if (is.null(dz) || length(dz) != 2L) {
      stop(
        "Object '", name, "' at time index ", tt,
        " could not be converted to a two-dimensional matrix."
      )
    }
    
    dz
  })
  
  out <- do.call(cbind, dims)
  rownames(out) <- c("row", "col")
  out
}




#' Prepare data structure for count-valued TERGM estimation
#'
#' Prepare a count-valued temporal ERGM data structure for use with
#' \code{btergm_count()}, including composition change, temporal covariates,
#' and structural-zero matrices.
#'
#' This is a count-valued variant of \code{tergmprepare()}. It follows the same
#' general preparation logic as the binary TERGM preparation routine, but it
#' treats the dependent networks as count-valued networks. When the dependent
#' networks are \code{network} objects, the count outcome is extracted from the
#' edge attribute named by \code{response}. The internal preparation then uses
#' count-valued matrices rather than binary adjacency matrices.
#'
#' The function adjusts the dimensions of networks and covariates within each
#' time step so that they are conformable for estimation. Depending on
#' \code{offset}, nodes that are absent from some objects within a time step are
#' either retained with structural-zero offset matrices or removed from the
#' corresponding time step. This behavior matches the original
#' \code{tergmprepare()} logic.
#'
#' The function also constructs temporal dyadic covariates from count-valued
#' dependent networks. In particular, \code{memory(type = "autoregression")}
#' uses the lagged count value, and \code{delrecip()} uses the lagged reciprocal
#' count value when \code{mutuality = FALSE}. Other \code{memory()} types are
#' interpreted as prior tie presence or absence, with nonzero counts treated as
#' prior tie presence.
#'
#' This function is intended for internal use by \code{btergm_count()} and is
#' not intended to be called directly by end users.
#'
#' @param formula The original temporal ERGM formula provided by the user.
#'   The left-hand side should be a list of count-valued networks or matrices
#'   in chronological order.
#' @param response Character string giving the name of the edge attribute that
#'   stores the count-valued outcome when the left-hand side consists of
#'   \code{network} objects. For example, \code{"tie_count"}.
#' @param offset Logical. Indicates whether absent nodes should be added where
#'   they are missing and represented through structural-zero offset matrices
#'   \code{offset = TRUE}, or removed where they are not shared across all
#'   relevant objects within a time step \code{offset = FALSE}.
#' @param blockdiag Logical. Should time steps be arranged in a block-diagonal
#'   structure? This is retained for compatibility with the original preparation
#'   function, but count-valued bootstrapped MPLE should usually use
#'   \code{blockdiag = FALSE}.
#' @param verbose Logical. Print details about dimension adjustment and temporal
#'   covariate construction?
#'
#' @return A list with the following slots:
#' \describe{
#'   \item{lhs.original}{A character object containing the original name of the
#'     object on the left-hand side of the formula.}
#'   \item{networks}{A list of adjusted \code{network} objects whose edge
#'     attribute named by \code{response} stores the count-valued outcome.}
#'   \item{response}{The name of the count-valued edge attribute.}
#'   \item{num.vertices}{The maximum number of nodes across time points after
#'     preparation.}
#'   \item{directed}{Logical. Are the dependent networks directed?}
#'   \item{bipartite}{Logical. Are the dependent networks bipartite?}
#'   \item{form}{The prepared formula, with \code{networks[[i]]} on the
#'     left-hand side and time-indexed covariates on the right-hand side.}
#'   \item{time.steps}{The number of time steps after accounting for temporal
#'     lag terms.}
#'   \item{rhs.terms}{A character vector containing the prepared right-hand-side
#'     model terms.}
#'   \item{covnames}{A character vector containing the names of the dependent
#'     network list and covariate objects stored in the returned list.}
#'   \item{auto.adjust}{Logical. Did the function adjust network or covariate
#'     dimensions?}
#'   \item{nvertices}{A matrix containing the number of rows and columns at each
#'     prepared time step.}
#'   \item{offsmat}{A list of structural-zero matrices. If \code{offset = FALSE},
#'     these matrices contain only zeros after node removal. If
#'     \code{offset = TRUE}, they contain ones for structurally unavailable
#'     dyads.}
#'   \item{original_networks_count}{The original dependent network list before
#'     count-valued matrix conversion. Used internally to restore vertex
#'     attributes.}
#'   \item{count_time_index}{The original time indices represented by the
#'     prepared network list after lagged temporal terms have been handled.}
#' }
#'
#' @author Patrick J. Edwards, adapted from Philip Leifeld's
#'   \code{tergmprepare()} implementation in \pkg{btergm}.
#'
#' @importFrom network is.network as.sociomatrix network.vertex.names
#'   list.vertex.attributes list.edge.attributes get.vertex.attribute
#'   set.vertex.attribute set.edge.attribute as.edgelist is.directed
#'   is.bipartite get.network.attribute network
#' @keywords internal
#' @noRd
tergmprepare_count <- function(formula, response = "count", offset = TRUE, blockdiag = FALSE, verbose = TRUE) {
  
  # extract response networks and make sure they are saved in a list
  l <- list()
  l$lhs.original <- deparse(formula[[2]])  # for reporting purposes later on
  l$networks <- eval(parse(text = deparse(formula[[2]])), envir = environment(formula))
  if ("list" %in% class(l$networks) || "network.list" %in% class(l$networks)) {
    # do nothing
  } else {
    l$networks <- list(l$networks)
  }
  
  # Count-model metadata. We keep this because the internal dimension
  # adjustment will use count-valued matrices, but the final objects passed
  # to ergmCntMPLE() should be network objects with vertex attributes and a
  # count-valued edge attribute.
  l$response <- response
  l$original_networks_count <- l$networks
  
  l$count_vertex_attributes <- vector("list", length(l$networks))
  l$count_vertex_names <- vector("list", length(l$networks))
  
  for (ii in seq_along(l$networks)) {
    if (network::is.network(l$networks[[ii]])) {
      l$count_vertex_names[[ii]] <- network::network.vertex.names(l$networks[[ii]])
      
      vertex_attr_names <- network::list.vertex.attributes(l$networks[[ii]])
      vertex_attr_names <- setdiff(vertex_attr_names, "na")
      
      l$count_vertex_attributes[[ii]] <- lapply(vertex_attr_names, function(attr) {
        network::get.vertex.attribute(l$networks[[ii]], attr)
      })
      
      names(l$count_vertex_attributes[[ii]]) <- vertex_attr_names
    } else {
      l$count_vertex_names[[ii]] <- rownames(as.matrix(l$networks[[ii]]))
      l$count_vertex_attributes[[ii]] <- list()
    }
  }
  
  # extract additional information
  l$num.vertices <- max(sapply(l$networks, function(x)
    network::get.network.attribute(network::network(x), "n")))  # number nodes
  if (is.network(l$networks[[1]])) {
    l$directed <- network::is.directed(l$networks[[1]])  # directed?
    l$bipartite <- network::is.bipartite(l$networks[[1]])  # bipartite?
  } else {
    if (is.mat.directed(as.matrix(l$networks[[1]]))) {
      l$directed <- TRUE
    } else {
      l$directed <- FALSE
    }
    if (is.mat.onemode(as.matrix(l$networks[[1]]))) {
      l$bipartite <- FALSE
    } else {
      l$bipartite <- TRUE
    }
  }
  
  # Convert dependent networks to count-valued matrices for internal
  # preparation. This is the key count-specific change. From this point until
  # the final reconstruction step, l$networks is a list of count matrices.
  l$networks <- lapply(seq_along(l$networks), function(ii) {
    x <- l$networks[[ii]]
    
    if (network::is.network(x)) {
      edge_attrs <- network::list.edge.attributes(x)
      
      if (!response %in% edge_attrs) {
        stop(
          "The dependent network at time step ",
          ii,
          " does not contain edge attribute '",
          response,
          "'. Available edge attributes are: ",
          paste(edge_attrs, collapse = ", "),
          "."
        )
      }
      
      mat <- network::as.sociomatrix(x, attrname = response)
      mat <- as.matrix(mat)
      storage.mode(mat) <- "numeric"
      mat[is.na(mat)] <- 0
      
      vn <- network::network.vertex.names(x)
      
      if (!is.null(vn) && length(vn) == nrow(mat)) {
        rownames(mat) <- vn
        colnames(mat) <- vn
      }
      
      diag(mat) <- 0
      
      mat
    } else {
      mat <- as.matrix(x)
      storage.mode(mat) <- "numeric"
      mat[is.na(mat)] <- 0
      diag(mat) <- 0
      mat
    }
  })
  
  
  
  # convert list elements to matrices if unknown data type
  for (i in 1:length(l$networks)) {
    if (!is.network(l$networks[[i]]) && !is.matrix(l$networks[[i]]) && !"list" %in% class(l$networks[[i]])) {
      tryCatch(
        {
          l$networks[[i]] <- as.matrix(l$networks[[i]])
        },
        error = function(cond) {
          stop(paste("Object", i, "could not be converted to a matrix."))
        }
      )
    }
  }
  
  # adjust and disassemble formula
  l$form <- update.formula(formula, networks[[i]] ~ .)
  l$time.steps <- length(l$networks)
  tilde <- deparse(l$form[[1]])
  lhs <- deparse(l$form[[2]])
  rhs <- paste(deparse(l$form[[3]]), collapse = "")  # for long formulae
  rhs <- gsub("\\s+", " ", rhs)
  
  # parse rhs of formula and add indices to edgecov and dyadcov terms
  rhsterms <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
  if (length(rhsterms) > 1) {  # 'transform' in 'timecov' may contain '+'
    for (i in length(rhsterms):2) {
      leftbracketmatches <- gregexpr("\\(", rhsterms[i])[[1]]
      leftbracketmatches <- leftbracketmatches[leftbracketmatches != -1]
      leftbracketmatches <- length(leftbracketmatches)
      rightbracketmatches <- gregexpr("\\)", rhsterms[i])[[1]]
      rightbracketmatches <- rightbracketmatches[rightbracketmatches != -1]
      rightbracketmatches <- length(rightbracketmatches)
      if (leftbracketmatches != rightbracketmatches) {
        rhsterms[i - 1] <- paste(rhsterms[i - 1], rhsterms[i], sep = " + ")
        rhsterms <- rhsterms[-i]
      }
    }
  }
  l$rhs.terms <- rhsterms
  rhs.operators <- rep("+", length(l$rhs.terms) - 1)
  
  # preprocess dyadcov and edgecov terms, memory terms, and timecov terms
  covnames <- character()
  for (k in 1:length(l$rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", l$rhs.terms[k])) {  # edgecov or dyadcov
      
      # split up into components
      if (grepl(",\\s*?((attr)|\\\")", l$rhs.terms[k])) { # with attrib arg.
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)"
      } else { # without attribute argument
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)"
      }
      x1 <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)  # before the covariate
      x2 <- sub(s, "\\5", l$rhs.terms[k], perl = TRUE)  # name of the cov.
      if (grepl("\\[.*\\]", x2)) {
        stop(paste0("Covariate names are not allowed to have indices: ", x2,
                    ". Please prepare a list object before estimation."))
      }
      if (grepl("^\"", x2)) next  # ignore built-in matrices b/c conformable
      x3 <- sub(s, "\\6", l$rhs.terms[k], perl = TRUE)  # after the covariate
      x.current <- eval(parse(text = x2), envir = environment(formula))
      type <- class(x.current)[1]
      l$covnames <- c(l$covnames, x2)
      l[[x2]] <- x.current
      if (grepl("\\[i\\]+$", x2)) {
        stop(paste0("Error in the following model term: ", l$rhs.terms[k],
                    ". The index 'i' is used internally by btergm. Please use a ",
                    "different index, for example 'j'."))
      }
      
      # add brackets if necessary, convert to list, and reassemble rhs term
      if (grepl("[^\\]]\\]$", x2)) {
        # time-varying covariate with given indices (e.g., formula[1:5])
        l$rhs.terms[k] <- paste0(x1, x2, x3)
        if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix",
                        "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
          x.current <-list(x.current)
          l[[x2]] <- x.current
        }
        if (length(x.current) != l$time.steps) {
          stop(paste(x2, "has", length(x.current), "elements, but there are",
                     l$time.steps, "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
      } else if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix",
                             "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
        # time-independent covariate
        if (!type %in% c("matrix", "network")) {
          x.current <- as.matrix(x.current)
        }
        l[[x2]] <- list()
        for (i in 1:l$time.steps) {
          l[[x2]][[i]] <- x.current
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
        l$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "list" || type == "network.list") {
        # time-varying covariate
        if (length(x.current) != l$time.steps) {
          stop(paste(x2, "has", length(get(x2)), "elements, but there are", l$time.steps, "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
        l$rhs.terms[k] <- paste0(x1, x2, x3)
      } else {  # something else --> try to convert to matrix list
        tryCatch(
          {
            l[[x2]] <- list(rep(as.matrix(x.current)), l$time.steps)
          },
          error = function(cond) {stop(paste0("Object '", x2, "' could not be converted to a matrix."))}
        )
      }
    } else if (grepl("memory", l$rhs.terms[k])) {  # memory terms
      
      # extract type argument
      s <- "(?:memory\\((?:.*type\\s*=\\s*)?(?:\"|'))(\\w+)(?:(\"|').*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        type <- "stability"
      } else {
        type <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
      }
      
      # extract lag argument
      s <- "(?:memory\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        lag <- 1
      } else {
        lag <- as.integer(sub(s, "\\1", l$rhs.terms[k], perl = TRUE))
      }
      if (lag > length(l$networks) - 1) {
        stop("The 'lag' argument in the 'memory' term is too large.")
      }
      
      # process dependent list of networks
      mem <- l$networks[-(length(l$networks):(length(l$networks) - lag + 1))]
      mem <- lapply(mem, as.matrix)
      memory <- list()
      for (i in 1:length(mem)) {
        mem[[i]] <- as.matrix(mem[[i]])
        storage.mode(mem[[i]]) <- "numeric"
        mem[[i]][is.na(mem[[i]])] <- 0
        
        if (type == "autoregression") {
          # Count-valued lagged tie strength.
          memory[[i]] <- mem[[i]]
        } else if (type == "stability") {
          # Binary prior presence/absence stability.
          memory[[i]] <- ifelse(mem[[i]] > 0, 1, -1)
        } else if (type == "innovation") {
          # Current tie formation where there was no prior tie.
          memory[[i]] <- ifelse(mem[[i]] == 0, 1, 0)
        } else if (type == "loss") {
          # Current absence/loss where there was a prior tie.
          memory[[i]] <- ifelse(mem[[i]] > 0, -1, 0)
        } else {
          stop("'type' argument in the 'memory' term not recognized.")
        }
      }
      rm(mem)
      
      # re-introduce as edgecov and name of model term including brackets
      l[["memory"]] <- memory
      if (blockdiag == TRUE) {
        l$rhs.terms[k] <- "edgecov(memory)"
      } else {
        l$rhs.terms[k] <- "edgecov(memory[[i]])"
      }
      l$covnames <- c(l$covnames, "memory")
    } else if (grepl("delrecip", l$rhs.terms[k])) {  # delayed reciprocity
      
      # extract mutuality argument
      s <- "(?:delrecip\\((?:.*mutuality\\s*=\\s*)?)((TRUE)|(FALSE)|T|F)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        mutuality <- FALSE
      } else {
        mutuality <- as.logical(sub(s, "\\1", l$rhs.terms[k], perl = TRUE))
      }
      
      # extract lag argument
      s <- "(?:delrecip\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"  # get lag
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        lag <- 1
      } else {
        lag <- as.integer(sub(s, "\\1", l$rhs.terms[k], perl = TRUE))
      }
      if (lag > length(l$networks) - 1) {
        stop("The 'lag' argument in the 'delrecip' term is too large.")
      }
      
      # process dependent list of networks
      dlr <- l$networks[-(length(l$networks):(length(l$networks) - lag + 1))]
      dlr <- lapply(dlr, function(x) t(as.matrix(x)))
      delrecip <- list()
      for (i in 1:length(dlr)) {
        dlr[[i]] <- as.matrix(dlr[[i]])
        storage.mode(dlr[[i]]) <- "numeric"
        dlr[[i]][is.na(dlr[[i]])] <- 0
        
        if (mutuality == TRUE) {
          delrecip[[i]] <- ifelse(dlr[[i]] > 0, 1, -1)
        } else {
          # Count-valued prior reciprocal tie strength.
          delrecip[[i]] <- dlr[[i]]
        }
      }
      rm(dlr)
      
      # re-introduce as edgecov and name of model term including brackets
      l[["delrecip"]] <- delrecip
      if (blockdiag == TRUE) {
        l$rhs.terms[k] <- "edgecov(delrecip)"
      } else {
        l$rhs.terms[k] <- "edgecov(delrecip[[i]])"
      }
      l$covnames <- c(l$covnames, "delrecip")
    } else if (grepl("timecov", l$rhs.terms[k])) {  # time covariate
      
      # extract x argument
      s <- "(?:timecov\\((?:.*x\\s*=\\s*)?)(\\w+)(?:.*\\))"
      if (sub(s, "\\1", l$rhs.terms[k], perl = TRUE) %in% c("minimum",
                                                            "maximum", "transform", "min", "max", "trans")) {
        s <- "(?:timecov\\(?:.*x\\s*=\\s*)(\\w+)(?:.*\\))"
      }
      
      # ensure there are no duplicate model term names
      countprevtc <- 1
      if (k > 1) {
        for (i in (k - 1):1) {
          if (grepl("timecov", l$rhs.terms[i])) {
            countprevtc <- countprevtc + 1
          }
        }
      }
      if (countprevtc > 0) {
        countprevtc <- as.character(countprevtc)
      } else {
        countprevtc <- ""
      }
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        x <- NULL
        label <- paste0("timecov", countprevtc)
      } else {
        x <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
        label <- paste0("timecov", countprevtc, ".", x)
      }
      
      # extract minimum argument
      s <- "(?:timecov\\(.*minimum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        minimum <- 1
      } else {
        minimum <- as.integer(sub(s, "\\1", l$rhs.terms[k], perl = TRUE))
      }
      
      # extract maximum argument
      s <- "(?:timecov\\(.*maximum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        maximum <- l$time.steps
      } else {
        maximum <- as.integer(sub(s, "\\1", l$rhs.terms[k], perl = TRUE))
      }
      
      # extract transform argument
      s <- "(?:timecov\\(.*transform\\s*=\\s*)(.+?)(?:(?:,|\\)$)]*.*)"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        transform <- function(t) t
      } else {
        transform <- eval(parse(text = sub(s, "\\1", l$rhs.terms[k],
                                           perl = TRUE)))
      }
      
      # process dependent list of networks
      if (is.null(x)) {
        covariate <- l[["networks"]]
        onlytime <- TRUE
      } else {
        onlytime <- FALSE
        covariate <- get(x)
      }
      tc <- timecov(covariate = covariate, minimum = minimum,
                    maximum = maximum, transform = transform, onlytime = onlytime)
      
      # re-introduce as edgecov and name of model term including brackets
      l[[label]] <- tc
      labelsuffix <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
      labelsuffix <-
        if (blockdiag == TRUE) {
          l$rhs.terms[k] <- paste0("edgecov(", label, ")")
        } else {
          l$rhs.terms[k] <- paste0("edgecov(", label, "[[i]])")
        }
      l$covnames <- c(l$covnames, label)
      
      # reporting
      if (verbose == TRUE) {
        timecovreporting <- matrix(sapply(tc, function(x) mean(x[1, 2])),
                                   nrow = 1)
        colnames(timecovreporting) <- paste0("t=", 1:length(l$networks))
        rownames(timecovreporting) <- ""
        message("Mean transformed timecov values:")
        print(timecovreporting)
      }
    }
  }
  l$covnames <- c("networks", l$covnames)
  
  missing_covs <- l$covnames[vapply(l$covnames, function(cn) is.null(l[[cn]]), logical(1))]
  
  if (length(missing_covs) > 0L) {
    stop(
      "The following covariates are listed in l$covnames but are NULL: ",
      paste(missing_covs, collapse = ", ")
    )
  }
  
  lengths <- sapply(l$covnames, function(cn) length(l[[cn]]))
  
  if (any(lengths == 0L)) {
    stop(
      "The following covariates have length 0: ",
      paste(names(lengths)[lengths == 0L], collapse = ", ")
    )
  }
  
  mn <- min(lengths)
  t.end <- max(lengths)
  t.start <- t.end - mn + 1
  
  # Trim all longer lists to the last mn elements, matching the original
  # btergm convention for lagged temporal covariates.
  for (ii in seq_along(l$covnames)) {
    cn <- l$covnames[[ii]]
    ll <- l[[cn]]
    
    if (length(ll) > mn) {
      l[[cn]] <- utils::tail(ll, mn)
    }
    
    if (length(l[[cn]]) != mn) {
      stop(
        "Covariate '", cn, "' has length ", length(l[[cn]]),
        " after trimming, but expected length ", mn, "."
      )
    }
    
    if (is.list(l[[cn]]) && !network::is.network(l[[cn]])) {
      null_idx <- which(vapply(l[[cn]], is.null, logical(1)))
      
      if (length(null_idx) > 0L) {
        stop(
          "Covariate '", cn, "' contains NULL elements after trimming at indices: ",
          paste(null_idx, collapse = ", ")
        )
      }
    }
  }
  
  l$time.steps <- mn
  
  # determine and report initial dimensions of networks and covariates
  if (verbose == TRUE) {
    if (length(l$covnames) > 1) {
      dimensions <- lapply(
        l$covnames,
        function(nm) tergmprepare_count_dim_matrix(
          obj = l[[nm]],
          name = nm,
          n_periods = l$time.steps
        )
      )
      rownames(dimensions[[1]]) <- paste(l$lhs.original, c("(row)", "(col)"))
      for (i in 2:length(dimensions)) {
        rownames(dimensions[[i]]) <- c(paste(l$covnames[i], "(row)"),
                                       paste(l$covnames[i], "(col)"))
      }
      dimensions <- do.call(rbind, dimensions)
      colnames(dimensions) <- paste0("t=", t.start:t.end) #1:length(l$networks))
      message("\nInitial dimensions of the network and covariates:")
      print(dimensions)
    } else {
      message("\nNo covariates provided.")
    }
  }
  
  # determine whether covariate dimensions need to be automatically adjusted
  l$auto.adjust <- FALSE
  if (length(l$covnames) > 1) {
    # check number of rows and columns
    dim_list <- lapply(
      l$covnames,
      function(nm) tergmprepare_count_dim_matrix(
        obj = l[[nm]],
        name = nm,
        n_periods = l$time.steps
      )
    )
    
    nr <- do.call(rbind, lapply(dim_list, function(x) x["row", , drop = FALSE]))
    nc <- do.call(rbind, lapply(dim_list, function(x) x["col", , drop = FALSE]))
    
    rownames(nr) <- l$covnames
    rownames(nc) <- l$covnames
    for (i in 1:ncol(nr)) {
      if (length(unique(nr[, i])) > 1) {
        l$auto.adjust <- TRUE
      }
    }
    for (i in 1:ncol(nc)) {
      if (length(unique(nc[, i])) > 1) {
        l$auto.adjust <- TRUE
      }
    }
    if (verbose == TRUE && l$auto.adjust == TRUE) {
      message(paste("\nDimensions differ across networks within time steps."))
    }
    # check if labels are present
    if (l$auto.adjust == TRUE) {
      for (i in 1:length(l$covnames)) {
        for (t in 1:l$time.steps) {
          if (is.null(rownames(as.matrix(l[[l$covnames[i]]][[t]]))) ||
              is.null(colnames(as.matrix(l[[l$covnames[i]]][[t]])))) {
            stop(paste0("The dimensions of the covariates differ, but ",
                        "covariate '", l$covnames[i],
                        " does not have node labels at t = ", t,
                        ". Automatic adjustment of dimensions is therefore not ",
                        "possible."))
          }
        }
      }
    }
    # check if there are different labels despite identical dimensions
    if (l$auto.adjust == FALSE) {
      for (t in 1:l$time.steps) {
        rlabels.i <- list()
        clabels.i <- list()
        for (i in 1:length(l$covnames)) {
          rlabels.i[[i]] <- rownames(as.matrix(l[[l$covnames[i]]][[t]]))
          clabels.i[[i]] <- colnames(as.matrix(l[[l$covnames[i]]][[t]]))
        }
        rlabels.i <- do.call(rbind, rlabels.i)
        clabels.i <- do.call(rbind, clabels.i)
        flag <- FALSE
        if (!is.null(rlabels.i)) {
          for (j in 1:ncol(rlabels.i)) {
            if (length(unique(rlabels.i[, j])) > 1) {
              l$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
        if (!is.null(clabels.i)) {
          for (j in 1:ncol(clabels.i)) {
            if (length(unique(clabels.i[, j])) > 1) {
              l$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
      }
      if (verbose == TRUE && flag == TRUE) {
        message(paste("\nSame dimensions but different labels across",
                      "networks within time steps."))
      }
    }
  }
  if (verbose == TRUE && l$auto.adjust == TRUE) {
    message("Trying to auto-adjust the dimensions of the networks. ",
            "If this fails, provide conformable matrices or network objects.")
  } else if (verbose == TRUE) {
    message("\nAll networks are conformable.")
  }
  
  # do mutual adjustment of networks and covariates at each time step
  structzero.df <- data.frame(label = character(), time = integer(),
                              object = character(), where = character())
  if (length(l$covnames) > 0 && l$auto.adjust == TRUE) {
    for (i in 1:l$time.steps) {
      for (j in 1:length(l$covnames)) {
        for (k in 1:length(l$covnames)) {
          if (j != k) {
            nw.j <- l[[l$covnames[j]]][[i]]
            rn.j <- rownames(as.matrix(nw.j))
            cn.j <- colnames(as.matrix(nw.j))
            nr.j <- nrow(as.matrix(nw.j))
            nc.j <- ncol(as.matrix(nw.j))
            nw.k <- l[[l$covnames[k]]][[i]]
            rn.k <- rownames(as.matrix(nw.k))
            cn.k <- colnames(as.matrix(nw.k))
            nr.k <- nrow(as.matrix(nw.k))
            nc.k <- ncol(as.matrix(nw.k))
            if (is.null(rn.j) || is.null(cn.j)) {
              stop(paste0("Missing row or column labels in object '",
                          l$covnames[j], "'. Provide row and column ",
                          "labels for all networks and covariates."))
            } else if (is.null(rn.k) || is.null(cn.k)) {
              stop(paste0("Missing row or column labels in object '",
                          l$covnames[k], "'. Provide row and column ",
                          "labels for all networks and covariates."))
            } else {
              if (is.null(rn.j) && !is.null(rn.k) && nr.j == nr.k) {
                if (is.network(nw.j)) {
                  network::set.vertex.attribute(nw.j, "vertex.names", rn.k)
                } else {
                  rownames(nw.j) <- rn.k
                }
              } else if (is.null(rn.k) && !is.null(rn.j) && nr.j == nr.k) {
                if (is.network(nw.k)) {
                  network::set.vertex.attribute(nw.k, "vertex.names", rn.j)
                } else {
                  rownames(nw.k) <- rn.j
                }
              } else if ((is.null(rn.k) || is.null(rn.j)) && nr.j != nr.k) {
                stop(paste0("Object '", l$covnames[j],
                            "' is incompatible with object '", l$covnames[k],
                            "' at t = ", i, "."))
              }
              # adjust j to k
              nw.j.labels <- adjust(nw.j, nw.k, remove = FALSE,
                                    value = 1, returnlabels = TRUE)
              nw.j <- adjust(nw.j, nw.k, remove = FALSE, value = 1)
              l[[l$covnames[j]]][[i]] <- nw.j
              ro <- nw.j.labels$added.row
              co <- nw.j.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)),
                                 object = rep(l$covnames[j], length(ro)),
                                 where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)),
                                 object = rep(l$covnames[j], length(co)),
                                 where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
              # adjust k back to j
              nw.k.labels <- adjust(nw.k, nw.j, remove = FALSE,
                                    value = 1, returnlabels = TRUE)
              nw.k <- adjust(nw.k, nw.j, remove = FALSE, value = 1)
              l[[l$covnames[k]]][[i]] <- nw.k
              ro <- nw.k.labels$added.row
              co <- nw.k.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)),
                                 object = rep(l$covnames[j], length(ro)),
                                 where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)),
                                 object = rep(l$covnames[j], length(co)),
                                 where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
            }
          }
        }
      }
    }
  }
  
  # check whether all dimensions are cross-sectionally conformable now
  nr.net <- sapply(l$networks, function(x) nrow(as.matrix(x)))
  for (i in 1:length(l$covnames)) {
    nr <- sapply(l[[l$covnames[i]]], function(x) {
      nrow(as.matrix(x))
    })
    for (j in 1:l$time.steps) {
      if (nr[j] != nr.net[j]) {
        stop(paste0("Covariate object '", l$covnames[i],
                    "' does not have the same number of rows as the dependent ",
                    "network at time step ", j, "."))
      }
    }
  }
  nc.net <- sapply(l$networks, function(x) ncol(as.matrix(x)))
  for (i in 1:length(l$covnames)) {
    nc <- sapply(l[[l$covnames[i]]], function(x) {
      ncol(as.matrix(x))
    })
    for (j in 1:l$time.steps) {
      if (nc[j] != nc.net[j]) {
        stop(paste0("Covariate object '", l$covnames[i],
                    "' does not have the same number of columns as the dependent ",
                    "network at time step ", j, "."))
      }
    }
  }
  
  # reporting
  if (verbose == TRUE) {
    if (l$auto.adjust == TRUE) {
      sz.row <- unique(structzero.df[structzero.df$where == "row", -3])
      szrownum <- numeric(length(l$networks))
      for (i in 1:length(l$networks)) {
        szrownum[i] <- nrow(sz.row[sz.row$time == i, ])
      }
      sz.col <- unique(structzero.df[structzero.df$where == "col", -3])
      szcolnum <- numeric(length(l$networks))
      for (i in 1:length(l$networks)) {
        szcolnum[i] <- nrow(sz.col[sz.col$time == i, ])
      }
      totrow <- sapply(l$networks, function(x) nrow(as.matrix(x)))
      totcol <- sapply(l$networks, function(x) ncol(as.matrix(x)))
      if (offset == TRUE) {
        dimensions <- rbind(totrow, totcol, szrownum, szcolnum,
                            totrow - szrownum, totcol - szcolnum)
        rownames(dimensions) <- c("total number of rows",
                                  "total number of columns", "row-wise structural zeros",
                                  "column-wise structural zeros", "remaining rows",
                                  "remaining columns")
      } else {
        dimensions <- rbind(szrownum, szcolnum, totrow - szrownum,
                            totcol - szcolnum)
        rownames(dimensions) <- c("maximum deleted nodes (row)",
                                  "maximum deleted nodes (col)", "remaining rows",
                                  "remaining columns")
      }
      colnames(dimensions) <- paste0("t=", t.start:t.end)
      if (nrow(structzero.df) > 0) {
        if (offset == TRUE) {
          message("\nNodes affected completely by structural zeros:")
        } else {
          message("\nAbsent nodes:")
        }
        szcopy <- structzero.df
        szcopy$time <- szcopy$time - 1 + t.start  # correct lagged starting time
        print(unique(szcopy))
      } else {
        message("\nAll nodes are retained.")
      }
      
      message("\nNumber of nodes per time step after adjustment:")
      print(dimensions)
    }
  }
  
  l$nvertices <- sapply(l$networks, function(x) c(nrow(as.matrix(x)),
                                                  ncol(as.matrix(x))))
  rownames(l$nvertices) <- c("row", "col")
  colnames(l$nvertices) <- paste0("t=", t.start:t.end)
  
  # create list of offset matrices (required both for offset and node removal)
  l$offsmat <- list()
  for (i in 1:l$time.steps) {
    mat <- matrix(0, nrow = nrow(as.matrix(l$networks[[i]])),
                  ncol = ncol(as.matrix(l$networks[[i]])))
    rownames(mat) <- rownames(as.matrix(l$networks[[i]]))
    colnames(mat) <- colnames(as.matrix(l$networks[[i]]))
    l$offsmat[[i]] <- mat
  }
  if (nrow(structzero.df) > 0) {
    for (i in 1:nrow(structzero.df)) {
      if (structzero.df$where[i] == "row") {
        index <- which(rownames(l$offsmat[[structzero.df$time[i]]]) ==
                         structzero.df$label[i])
        l$offsmat[[structzero.df$time[i]]][index, ] <- 1
      } else {
        index <- which(colnames(l$offsmat[[structzero.df$time[i]]]) ==
                         structzero.df$label[i])
        l$offsmat[[structzero.df$time[i]]][, index] <- 1
      }
    }
  }
  
  # offset preparation or node removal for MPLE
  if (offset == TRUE) {
    # add offset to formula and reassemble formula
    l$rhs.terms[length(l$rhs.terms) + 1] <- "offset(edgecov(offsmat[[i]]))"
    rhs.operators[length(rhs.operators) + 1] <- "+"
  } else {
    # delete nodes with structural zeros
    if (l$auto.adjust == TRUE) {
      if (l$bipartite) {
        l$offsmat <- lapply(l$offsmat, function(x) {
          x_ind <- which(apply(x, 1, function(y) all(y == 1)))
          y_ind <- which(apply(x, 2, function(y) all(y == 1)))
          if (length(x_ind) > 0 && length(y_ind) > 0) {
            x[-x_ind, -y_ind]
          } else if (length(x_ind) > 0) {
            x[-x_ind, ]
          } else if (length(y_ind) > 0) {
            x[, -y_ind]
          } else {
            x
          }
        })
      } else {
        l$offsmat <- suppressMessages(handleMissings(l$offsmat, na = 1, method = "remove"))
      }
      for (j in 1:length(l$covnames)) {
        l[[l$covnames[j]]] <- adjust(l[[l$covnames[j]]], l$offsmat)
      }
    }
  }
  
  # determine and report initial dimensions of networks and covariates
  if (verbose == TRUE && length(l$covnames) > 1) {
    dimensions <- lapply(lapply(l$covnames, function(x) l[[x]]),
                         function(y) sapply(y, function(z) dim(as.matrix(z))))
    rownames(dimensions[[1]]) <- paste(l$lhs.original, c("(row)", "(col)"))
    for (i in 2:length(dimensions)) {
      rownames(dimensions[[i]]) <- c(paste(l$covnames[i], "(row)"),
                                     paste(l$covnames[i], "(col)"))
    }
    dimensions <- do.call(rbind, dimensions)
    colnames(dimensions) <- paste0("t=", t.start:t.end) #1:length(l$networks))
    message("\nDimensions of the network and covariates after adjustment:")
    print(dimensions)
  }
  
  # Rebuild the adjusted count matrices as network objects with the response
  # edge attribute and vertex attributes. This keeps ergm count terms and
  # node-level terms like nodeofactor() working after preparation.
  l$count_time_index <- t.start:t.end
  
  l$networks <- lapply(seq_along(l$networks), function(ii) {
    original_time <- l$count_time_index[ii]
    
    btergm_count_matrix_to_network(
      mat = as.matrix(l$networks[[ii]]),
      original_network = l$original_networks_count[[original_time]],
      response = response,
      directed = l$directed,
      bipartite = l$bipartite
    )
  })
  
  # assemble formula
  rhs <- l$rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], l$rhs.terms[i + 1])
    }
  }
  f <- paste(lhs, tilde, rhs)
  l$form <- as.formula(f, env = environment())  # l$form <- as.formula(f)
  
  # for mtergm estimation using MCMC: create block-diagonal matrices
  if (blockdiag == TRUE) {
    if (l$bipartite == TRUE) {
      stop(paste("MCMC estimation is currently only supported for one-mode",
                 "networks. Use the btergm function instead."))
    }
    # also save formula without time indices for ergm estimation
    l$form <- update.formula(l$form, networks ~ .)
    l$form <- paste(deparse(l$form), collapse = "")
    l$form <- paste(l$form, "+ offset(edgecov(offsmat))")
    l$form <- as.formula(l$form, env = environment())
    # make covariates block-diagonal
    if (length(l$covnames) > 1) {
      for (j in 2:length(l$covnames)) {
        l[[l$covnames[j]]] <- as.matrix(Matrix::bdiag(lapply(
          l[[l$covnames[j]]], as.matrix)))
      }
    }
    # create block-diagonal offset matrix and merge with existing offsmat term
    l$offsmat <- as.matrix(Matrix::bdiag(l$offsmat))  # make block-diagonal
    bdoffset <- lapply(l$networks, as.matrix)
    for (i in 1:length(bdoffset)) {
      bdoffset[[i]][, ] <- 1
    }
    bdoffset <- as.matrix((Matrix::bdiag(bdoffset) - 1) * -1)  # off-diagonal
    l$offsmat <- l$offsmat + bdoffset
    rm(bdoffset)
    l$offsmat[l$offsmat > 0] <- 1
    # make dependent network block-diagonal
    if (is.network(l$networks[[1]])) {  # network
      attrnames <- network::list.vertex.attributes(l$networks[[1]])
      attributes <- list()
      for (i in 1:length(l$networks)) {
        attrib <- list()
        for (j in 1:length(attrnames)) {
          attrib[[j]] <- network::get.vertex.attribute(l$networks[[i]],
                                                       attrnames[j])
        }
        attributes[[i]] <- attrib
        l$networks[[i]] <- as.matrix(l$networks[[i]])
      }
      l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                     directed = l$directed, bipartite = l$bipartite)
      for (i in 1:length(attrnames)) {  # collate attributes and merge back in
        attrib <- unlist(lapply(attributes, function(x) x[[i]]))
        network::set.vertex.attribute(l$networks, attrnames[i], attrib)
      }
    } else {  # matrix
      l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                     directed = l$directed, bipartite = l$bipartite)
    }
    if (verbose == TRUE) {
      cat("\n")  # to get a blank line before the MCMC MLE output starts
    }
  }
  # convert formula back to character object because environment invalid anyway
  form3 <- paste(deparse(l$form[[3]]), collapse = "")  # for long formulae
  form3 <- gsub("\\s+", " ", form3)
  l$form <- paste(deparse(l$form[[2]]), deparse(l$form[[1]]), form3)
  return(l)  # return the environment with all the data
}











#' Convert an adjusted count matrix to a count-valued network object
#'
#' Internal helper used by \code{tergmprepare_count()} to rebuild adjusted
#' count-valued matrices as \code{network} objects after composition adjustment.
#' The resulting network stores nonzero entries as edges and stores the original
#' count values in the edge attribute named by \code{response}. Vertex
#' attributes are restored from the original network by matching vertex names.
#'
#' @param mat A count-valued adjacency matrix after dimension adjustment.
#' @param original_network The original \code{network} object for the same time
#'   point, used to restore vertex attributes.
#' @param response Character string naming the count-valued edge attribute to
#'   create.
#' @param directed Logical. Should the output network be directed?
#' @param bipartite Logical. Should the output network be bipartite?
#'
#' @return A \code{network} object with count values stored as an edge attribute.
#'
#' @importFrom network network set.vertex.attribute list.vertex.attributes
#'   get.vertex.attribute network.vertex.names as.edgelist set.edge.attribute
#'   is.network
#' @keywords internal
#' @noRd
btergm_count_matrix_to_network <- function(mat,
                                           original_network = NULL,
                                           response = "count",
                                           directed = TRUE,
                                           bipartite = FALSE) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  diag(mat) <- 0
  
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
  
  # Restore vertex attributes by matching vertex.names.
  if (!is.null(original_network) && network::is.network(original_network)) {
    old_names <- network::network.vertex.names(original_network)
    new_names <- network::network.vertex.names(nw)
    
    vertex_attr_names <- network::list.vertex.attributes(original_network)
    vertex_attr_names <- setdiff(vertex_attr_names, c("na", "vertex.names"))
    
    for (attr in vertex_attr_names) {
      old_values <- network::get.vertex.attribute(original_network, attr)
      names(old_values) <- old_names
      
      new_values <- old_values[new_names]
      
      network::set.vertex.attribute(nw, attr, unname(new_values))
    }
  }
  
  # Set count response edge attribute.
  el <- network::as.edgelist(nw)
  
  if (nrow(el) > 0L) {
    vals <- mat[cbind(el[, 1], el[, 2])]
    vals[is.na(vals)] <- 0
    network::set.edge.attribute(nw, response, vals)
  }
  
  nw
}