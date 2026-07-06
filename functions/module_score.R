#' Calculate bulk RNA-seq module scores
#'
#' Calculates an AddModuleScore-style score for one or more gene sets in a bulk
#' RNA-seq expression matrix. Genes are binned by average expression, control
#' genes are sampled from matching bins, and each sample is scored as the module
#' mean minus the matched-control mean.
#'
#' @param x Numeric expression matrix with genes in rows and samples in columns.
#' @param genes_to_score Character vector of genes or a named list of character
#'   vectors. Each vector is one module.
#' @param universe Optional character vector of genes eligible as controls.
#'   Defaults to all genes in `x`.
#' @param nbin Number of average-expression bins used for control matching.
#' @param ctrl Maximum number of control genes sampled per module gene.
#' @param seed Random seed. Use `NULL` to leave the random state unmanaged.
#' @param min_genes Minimum number of module genes that must be present in `x`.
#' @param warn_missing If `TRUE`, warn when requested genes are absent from `x`.
#' @param verbose If `TRUE`, print a module audit summary and return score,
#'   gene-set, and control details.
#'
#' @return If `verbose = FALSE`, a data frame with samples in rows and modules
#'   in columns. If `verbose = TRUE`, a list with `scores`, `gene_set_summary`,
#'   and `details`.
calc_bulk_module_score <- function(x,
                                   genes_to_score,
                                   universe = NULL,
                                   nbin = 24,
                                   ctrl = 100,
                                   seed = 1,
                                   min_genes = 1,
                                   warn_missing = TRUE,
                                   verbose = FALSE) {
  x <- .validate_expression_matrix(x)
  genes_to_score <- .normalize_genes_to_score(genes_to_score)
  .check_score_args(nbin = nbin, ctrl = ctrl, seed = seed, min_genes = min_genes)
  if (!isTRUE(verbose) && !identical(verbose, FALSE)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  nbin <- as.integer(nbin)
  ctrl <- as.integer(ctrl)
  min_genes <- as.integer(min_genes)

  all_genes <- rownames(x)
  if (length(all_genes) < 2L) {
    stop("x must contain at least two genes.", call. = FALSE)
  }

  if (is.null(universe)) {
    universe <- all_genes
  } else {
    universe <- unique(as.character(stats::na.omit(universe)))
    missing_universe <- setdiff(universe, all_genes)
    if (length(missing_universe) > 0L) {
      stop("universe contains genes not present in x.", call. = FALSE)
    }
  }

  module_genes <- unique(unlist(genes_to_score, use.names = FALSE))
  control_universe <- setdiff(universe, module_genes)
  if (length(control_universe) == 0L) {
    stop("control universe is empty after removing module genes.", call. = FALSE)
  }

  gene_bins <- .get_expression_bins(x, nbin = nbin)
  old_seed <- .save_random_seed()
  on.exit(.restore_random_seed(old_seed), add = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  scores <- vector("list", length(genes_to_score))
  names(scores) <- names(genes_to_score)
  details <- if (verbose) scores else NULL
  gene_set_summary <- data.frame(
    module = names(genes_to_score),
    n_requested = lengths(genes_to_score),
    n_present = integer(length(genes_to_score)),
    n_missing = integer(length(genes_to_score)),
    n_controls = integer(length(genes_to_score)),
    stringsAsFactors = FALSE
  )

  for (module_name in names(genes_to_score)) {
    requested <- genes_to_score[[module_name]]
    present <- intersect(requested, all_genes)
    missing <- setdiff(requested, present)

    if (length(present) < min_genes) {
      stop("module has fewer than min_genes present in x: ", module_name, call. = FALSE)
    }
    if (warn_missing && length(missing) > 0L) {
      warning(
        "module has missing genes: ", module_name,
        " (", length(missing), " missing)",
        call. = FALSE
      )
    }

    score_calc <- .score_one_module(
      x = x,
      genes = present,
      universe = control_universe,
      gene_bins = gene_bins,
      ctrl = ctrl
    )
    scores[[module_name]] <- score_calc$scores

    gene_set_summary[gene_set_summary$module == module_name, "n_present"] <- length(present)
    gene_set_summary[gene_set_summary$module == module_name, "n_missing"] <- length(missing)
    gene_set_summary[gene_set_summary$module == module_name, "n_controls"] <- length(score_calc$pooled_controls)

    if (verbose) {
      details[[module_name]] <- list(
        requested_genes = requested,
        present_genes = present,
        missing_genes = missing,
        controls_by_gene = score_calc$controls_by_gene,
        pooled_controls = score_calc$pooled_controls
      )
    }
  }

  score_df <- as.data.frame(scores, check.names = FALSE)
  rownames(score_df) <- colnames(x)
  attr(score_df, "gene_set_summary") <- gene_set_summary
  attr(score_df, "params") <- list(
    nbin = length(unique(gene_bins)),
    ctrl = ctrl,
    seed = seed,
    universe_size = length(control_universe)
  )

  if (!verbose) {
    return(score_df)
  }

  .print_gene_set_summary(gene_set_summary)

  list(scores = score_df, gene_set_summary = gene_set_summary, details = details)
}

.print_gene_set_summary <- function(gene_set_summary) {
  cat("\nBulk module score audit\n")
  print(gene_set_summary, row.names = FALSE)
  invisible(gene_set_summary)
}

.validate_expression_matrix <- function(x, name = "x") {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (!is.numeric(x) || is.null(rownames(x)) || is.null(colnames(x))) {
    stop(name, " must be a numeric matrix with gene rownames and sample colnames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) {
    stop(name, " has duplicated gene rownames.", call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) {
    stop(name, " has duplicated sample colnames.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " contains non-finite values.", call. = FALSE)
  }

  x
}

.normalize_genes_to_score <- function(genes_to_score) {
  if (is.character(genes_to_score)) {
    genes_to_score <- list(bulk_module_score = genes_to_score)
  }
  if (!is.list(genes_to_score) || length(genes_to_score) == 0L) {
    stop("genes_to_score must be a non-empty character vector or named list of character vectors.", call. = FALSE)
  }

  genes_to_score <- lapply(genes_to_score, function(genes) {
    genes <- as.character(stats::na.omit(genes))
    unique(genes[nzchar(genes)])
  })
  genes_to_score <- genes_to_score[lengths(genes_to_score) > 0L]
  if (length(genes_to_score) == 0L) {
    stop("genes_to_score is empty after removing missing or blank genes.", call. = FALSE)
  }

  gene_set_names <- names(genes_to_score)
  if (is.null(gene_set_names)) {
    gene_set_names <- rep("", length(genes_to_score))
  }
  empty_names <- is.na(gene_set_names) | !nzchar(gene_set_names)
  gene_set_names[empty_names] <- paste0("module_", which(empty_names))
  names(genes_to_score) <- make.unique(gene_set_names, sep = "_")

  genes_to_score
}

.check_score_args <- function(nbin, ctrl, seed, min_genes) {
  if (!.is_single_number(nbin) || nbin < 2L) {
    stop("nbin must be a single number >= 2.", call. = FALSE)
  }
  if (!.is_single_number(ctrl) || ctrl < 1L) {
    stop("ctrl must be a single number >= 1.", call. = FALSE)
  }
  if (!.is_single_number(min_genes) || min_genes < 1L) {
    stop("min_genes must be a single number >= 1.", call. = FALSE)
  }
  if (!is.null(seed) && !.is_single_number(seed)) {
    stop("seed must be NULL or a single finite number.", call. = FALSE)
  }
}

.check_perm_args <- function(n_perm, ctrl, seed, module_name) {
  if (!.is_single_number(n_perm) || n_perm < 1L) {
    stop("n_perm must be a single positive number.", call. = FALSE)
  }
  if (!.is_single_number(ctrl) || ctrl < 1L) {
    stop("ctrl must be a single positive number.", call. = FALSE)
  }
  if (!is.null(seed) && !.is_single_number(seed)) {
    stop("seed must be NULL or a single finite number.", call. = FALSE)
  }
  if (!is.character(module_name) || length(module_name) != 1L || !nzchar(module_name)) {
    stop("module_name must be a single non-empty string.", call. = FALSE)
  }
}

.is_single_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

# Snapshot and restore the caller's RNG state so that setting `seed` for
# reproducible control sampling does not disturb random draws elsewhere in the
# user's session.
.save_random_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
}

.restore_random_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
}

.get_expression_bins <- function(x, nbin = 24) {
  x <- .validate_expression_matrix(x)
  if (!.is_single_number(nbin) || nbin < 2L) {
    stop("nbin must be a single number >= 2.", call. = FALSE)
  }

  all_genes <- rownames(x)
  nbin <- min(as.integer(nbin), length(all_genes))
  gene_means <- rowMeans(x)
  breaks <- stats::quantile(gene_means, probs = seq(0, 1, length.out = nbin + 1L), na.rm = TRUE)

  # Prefer quantile breaks so bins hold roughly equal gene counts. When many
  # genes share a value (ties collapse the breaks), fall back to rank-based
  # binning so every bin is still populated.
  if (length(unique(breaks)) == length(breaks)) {
    gene_bins <- cut(gene_means, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  } else {
    gene_ranks <- rank(gene_means, ties.method = "average", na.last = "keep")
    gene_bins <- ceiling(gene_ranks / max(gene_ranks, na.rm = TRUE) * nbin)
    gene_bins <- pmin(pmax(gene_bins, 1L), nbin)
  }
  names(gene_bins) <- all_genes
  gene_bins
}

.score_one_module <- function(x, genes, universe, gene_bins, ctrl) {
  present <- intersect(genes, rownames(x))
  if (length(present) == 0L) {
    stop("genes has no genes present in x.", call. = FALSE)
  }

  control_universe <- setdiff(universe, present)
  if (length(control_universe) == 0L) {
    stop("control universe is empty after removing module genes.", call. = FALSE)
  }

  # For each module gene, draw up to `ctrl` controls from its own expression
  # bin. If that bin has no eligible genes (e.g. a restricted universe), fall
  # back to the whole control universe for that gene.
  controls_by_gene <- lapply(present, function(gene) {
    bin_genes <- control_universe[gene_bins[control_universe] == gene_bins[[gene]]]
    if (length(bin_genes) == 0L) {
      bin_genes <- control_universe
    }
    sample(bin_genes, min(ctrl, length(bin_genes)), replace = FALSE)
  })
  names(controls_by_gene) <- present

  # Pool controls across module genes and de-duplicate so a gene sampled for
  # several module genes is counted once in the control mean.
  pooled_controls <- unique(unlist(controls_by_gene, use.names = FALSE))
  list(
    scores = colMeans(x[present, , drop = FALSE]) -
      colMeans(x[pooled_controls, , drop = FALSE]),
    controls_by_gene = controls_by_gene,
    pooled_controls = pooled_controls
  )
}
