#' Permutation test for a bulk RNA-seq module score
#'
#' Tests one observed gene set against random gene sets of the same size. The
#' default null keeps genes expression-matched by drawing each null gene from
#' the same expression bin as the corresponding observed gene.
#'
#' @param x Numeric expression matrix with genes in rows and samples in columns.
#' @param gene_list Character vector of genes in the observed module.
#' @param metadata Optional sample metadata.
#' @param sample_col Column in `metadata` containing sample IDs.
#' @param group_col Optional metadata column used to summarize scores by group.
#' @param module_name Name used in output tables and plot titles.
#' @param n_perm Number of null gene sets to score.
#' @param null_method `"matched_bins"` samples null genes from matching
#'   expression bins. `"random"` samples from the full candidate universe.
#' @param random_genome If `TRUE`, use genome-wide random null gene sets from
#'   `universe`. This is a convenience switch equivalent to
#'   `null_method = "random"`.
#' @param universe Optional character vector of genes eligible for null sets and
#'   controls. Defaults to all genes in `x`.
#' @param nbin Number of expression bins.
#' @param ctrl Maximum number of control genes sampled per module gene.
#' @param seed Random seed. Use `NULL` to leave the random state unmanaged.
#' @param summary Summary used for grouped scores: `"mean"` or `"median"`.
#' @param trajectory_stat Group-level trajectory statistic: `"last_minus_first"`,
#'   `"slope"`, or `"none"`.
#' @param alternative Empirical p-value direction: `"two.sided"`, `"greater"`,
#'   or `"less"`.
#' @param make_plot If `TRUE`, include a trajectory plot in the result.
#' @param make_histogram If `TRUE`, include a null distribution histogram when a
#'   trajectory statistic is calculated.
#'
#' @return A list containing observed scores, permutation scores, summaries,
#'   empirical p-values, optional plots, and parameters.
perm_bulk_module_score <- function(x,
                                   gene_list,
                                   metadata = NULL,
                                   sample_col = "sample_id",
                                   group_col = NULL,
                                   module_name = "bulk_module_score",
                                   n_perm = 1000,
                                   null_method = c("matched_bins", "random"),
                                   random_genome = FALSE,
                                   universe = NULL,
                                   nbin = 24,
                                   ctrl = 100,
                                   seed = 1,
                                   summary = c("mean", "median"),
                                   trajectory_stat = c("last_minus_first", "slope", "none"),
                                   alternative = c("two.sided", "greater", "less"),
                                   make_plot = TRUE,
                                   make_histogram = TRUE) {
  x <- .validate_expression_matrix(x)
  null_method <- match.arg(null_method)
  if (!is.logical(random_genome) || length(random_genome) != 1L || is.na(random_genome)) {
    stop("random_genome must be TRUE or FALSE.", call. = FALSE)
  }
  if (random_genome) {
    null_method <- "random"
  }
  summary <- match.arg(summary)
  trajectory_stat <- match.arg(trajectory_stat)
  alternative <- match.arg(alternative)
  .check_perm_args(n_perm = n_perm, ctrl = ctrl, seed = seed, module_name = module_name)

  n_perm <- as.integer(n_perm)
  ctrl <- as.integer(ctrl)
  summary_fun <- switch(summary,
    mean = mean,
    median = stats::median
  )

  gene_list <- .normalize_genes_to_score(list(module = gene_list))$module
  observed_genes <- intersect(gene_list, rownames(x))
  missing_genes <- setdiff(gene_list, observed_genes)
  if (length(observed_genes) == 0L) {
    stop("gene_list has no genes present in x.", call. = FALSE)
  }

  if (is.null(universe)) {
    universe <- rownames(x)
  } else {
    universe <- unique(as.character(stats::na.omit(universe)))
    missing_universe <- setdiff(universe, rownames(x))
    if (length(missing_universe) > 0L) {
      stop("universe contains genes not present in x.", call. = FALSE)
    }
  }

  candidate_genes <- setdiff(universe, observed_genes)
  if (length(candidate_genes) < length(observed_genes)) {
    stop("universe has too few non-module genes for permutation.", call. = FALSE)
  }

  metadata_df <- .prepare_sample_metadata(
    sample_ids = colnames(x),
    metadata = metadata,
    sample_col = sample_col
  )

  old_seed <- .save_random_seed()
  on.exit(.restore_random_seed(old_seed), add = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  gene_bins <- .get_expression_bins(x, nbin = nbin)
  observed_scores <- .score_one_module(
    x = x,
    genes = observed_genes,
    universe = universe,
    gene_bins = gene_bins,
    ctrl = ctrl
  )$scores
  observed_score_df <- data.frame(
    permutation = "observed",
    sample_id = colnames(x),
    module_score = observed_scores,
    stringsAsFactors = FALSE
  )

  permutation_gene_sets <- vector("list", n_perm)
  permutation_scores <- vector("list", n_perm)
  for (perm_idx in seq_len(n_perm)) {
    null_genes <- if (null_method == "matched_bins") {
      .sample_matched_null_genes(observed_genes, candidate_genes, gene_bins)
    } else {
      sample(candidate_genes, length(observed_genes), replace = FALSE)
    }
    permutation_gene_sets[[perm_idx]] <- null_genes

    null_scores <- .score_one_module(
      x = x,
      genes = null_genes,
      universe = universe,
      gene_bins = gene_bins,
      ctrl = ctrl
    )$scores
    permutation_scores[[perm_idx]] <- data.frame(
      permutation = paste0("perm_", perm_idx),
      sample_id = colnames(x),
      module_score = null_scores,
      stringsAsFactors = FALSE
    )
  }

  permutation_score_df <- do.call(rbind, permutation_scores)
  rownames(permutation_score_df) <- NULL

  observed_summary <- .summarize_perm_scores(
    observed_score_df,
    metadata_df = metadata_df,
    group_col = group_col,
    summary_fun = summary_fun
  )
  permutation_summary <- .summarize_perm_scores(
    permutation_score_df,
    metadata_df = metadata_df,
    group_col = group_col,
    summary_fun = summary_fun
  )
  null_mean_summary <- stats::aggregate(
    score ~ .,
    data = permutation_summary[, setdiff(colnames(permutation_summary), "permutation"), drop = FALSE],
    FUN = mean
  )
  null_mean_summary$permutation <- "permutation_mean"
  null_median_summary <- stats::aggregate(
    score ~ .,
    data = permutation_summary[, setdiff(colnames(permutation_summary), "permutation"), drop = FALSE],
    FUN = stats::median
  )
  null_median_summary$permutation <- "permutation_median"

  sample_p_values <- .calc_sample_p_values(
    observed_score_df = observed_score_df,
    permutation_score_df = permutation_score_df,
    alternative = alternative
  )
  group_p_values <- .calc_group_p_values(
    observed_summary = observed_summary,
    permutation_summary = permutation_summary,
    group_col = group_col,
    alternative = alternative
  )
  trajectory <- .calc_trajectory_result(
    observed_summary = observed_summary,
    permutation_summary = permutation_summary,
    group_col = group_col,
    statistic = trajectory_stat,
    alternative = alternative
  )

  result <- list(
    observed_scores = observed_score_df,
    permutation_scores = permutation_score_df,
    observed_summary = observed_summary,
    permutation_summary = permutation_summary,
    null_mean_summary = null_mean_summary,
    null_median_summary = null_median_summary,
    sample_p_values = sample_p_values,
    group_p_values = group_p_values,
    trajectory = trajectory$summary,
    trajectory_null = trajectory$null,
    permutation_gene_sets = permutation_gene_sets,
    missing_genes = missing_genes,
    plot_data = .build_perm_plot_data(observed_summary, permutation_summary, null_median_summary),
    params = list(
      module_name = module_name,
      n_present_genes = length(observed_genes),
      n_requested_genes = length(gene_list),
      n_perm = n_perm,
      null_method = null_method,
      random_genome = random_genome,
      nbin = nbin,
      ctrl = ctrl,
      seed = seed,
      summary = summary,
      group_col = group_col,
      plot_x = if (is.null(group_col)) "sample_index" else "group",
      alternative = alternative
    )
  )

  result <- c(
    result,
    list(
      plot = if (make_plot) .plot_perm_trajectory(result) else NULL,
      histogram = if (make_histogram && !is.null(result$trajectory)) {
        .plot_perm_distribution(result)
      } else {
        NULL
      }
    )
  )

  class(result) <- c("bulk_module_score_perm", class(result))
  result
}

.prepare_sample_metadata <- function(sample_ids, metadata, sample_col) {
  if (is.null(metadata)) {
    return(data.frame(
      sample_id = sample_ids,
      sample_index = seq_along(sample_ids),
      stringsAsFactors = FALSE
    ))
  }

  metadata <- as.data.frame(metadata)
  if (!sample_col %in% colnames(metadata)) {
    stop("metadata is missing sample_col: ", sample_col, call. = FALSE)
  }
  if (anyNA(metadata[[sample_col]]) || anyDuplicated(metadata[[sample_col]])) {
    stop("metadata sample identifiers must be non-missing and unique.", call. = FALSE)
  }

  sample_match <- match(sample_ids, as.character(metadata[[sample_col]]))
  if (anyNA(sample_match)) {
    stop("metadata is missing samples present in x.", call. = FALSE)
  }

  metadata <- metadata[sample_match, , drop = FALSE]
  metadata$sample_id <- sample_ids
  metadata$sample_index <- seq_along(sample_ids)
  metadata
}

.sample_matched_null_genes <- function(observed_genes, candidate_genes, gene_bins) {
  selected <- character(0)
  all_bins <- sort(unique(gene_bins[candidate_genes]))

  # Draw one null gene per observed gene from the observed gene's expression
  # bin, without reuse, so the null set is size- and expression-matched.
  for (gene in observed_genes) {
    target_bin <- gene_bins[[gene]]
    available <- candidate_genes[
      gene_bins[candidate_genes] == target_bin & !candidate_genes %in% selected
    ]

    # If the matching bin is exhausted, widen to the nearest bins by distance
    # until an unused candidate is found.
    if (length(available) == 0L) {
      for (bin in all_bins[order(abs(all_bins - target_bin))]) {
        available <- candidate_genes[
          gene_bins[candidate_genes] == bin & !candidate_genes %in% selected
        ]
        if (length(available) > 0L) {
          break
        }
      }
    }

    if (length(available) == 0L) {
      stop("not enough candidate genes to draw a matched null set.", call. = FALSE)
    }

    selected <- c(selected, sample(available, 1L))
  }

  selected
}

.summarize_perm_scores <- function(score_df, metadata_df, group_col, summary_fun) {
  if (is.null(group_col)) {
    out <- data.frame(
      permutation = score_df$permutation,
      sample_id = score_df$sample_id,
      sample_index = metadata_df$sample_index[match(score_df$sample_id, metadata_df$sample_id)],
      score = score_df$module_score,
      stringsAsFactors = FALSE
    )
    return(out[order(out$permutation, out$sample_index), , drop = FALSE])
  }
  if (!group_col %in% colnames(metadata_df)) {
    stop("metadata is missing group_col: ", group_col, call. = FALSE)
  }

  score_df <- merge(
    score_df,
    metadata_df[, c("sample_id", group_col), drop = FALSE],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  score_df <- score_df[!is.na(score_df[[group_col]]), , drop = FALSE]

  group_values <- score_df[[group_col]]
  group_order <- if (is.factor(group_values)) {
    levels(group_values)
  } else {
    sort(unique(group_values))
  }

  pieces <- split(score_df, list(score_df$permutation, group_values), drop = TRUE)
  out <- do.call(rbind, lapply(pieces, function(piece) {
    data.frame(
      permutation = piece$permutation[[1]],
      group = piece[[group_col]][[1]],
      score = summary_fun(piece$module_score),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out$group <- factor(out$group, levels = group_order)
  out[order(out$permutation, out$group), , drop = FALSE]
}

.calc_sample_p_values <- function(observed_score_df, permutation_score_df, alternative) {
  do.call(rbind, lapply(unique(observed_score_df$sample_id), function(sample_id) {
    observed <- observed_score_df$module_score[observed_score_df$sample_id == sample_id][[1L]]
    null <- permutation_score_df$module_score[permutation_score_df$sample_id == sample_id]
    data.frame(
      sample_id = sample_id,
      observed_score = observed,
      null_mean = mean(null),
      null_sd = stats::sd(null),
      empirical_p = .empirical_p_value(observed, null, alternative = alternative),
      stringsAsFactors = FALSE
    )
  }))
}

.calc_group_p_values <- function(observed_summary, permutation_summary, group_col, alternative) {
  if (is.null(group_col)) {
    return(NULL)
  }

  do.call(rbind, lapply(levels(observed_summary$group), function(group_value) {
    observed <- observed_summary$score[observed_summary$group == group_value][[1L]]
    null <- permutation_summary$score[permutation_summary$group == group_value]
    data.frame(
      group = group_value,
      observed_score = observed,
      null_mean = mean(null),
      null_sd = stats::sd(null),
      empirical_p = .empirical_p_value(observed, null, alternative = alternative),
      stringsAsFactors = FALSE
    )
  }))
}

.calc_trajectory_result <- function(observed_summary,
                                    permutation_summary,
                                    group_col,
                                    statistic,
                                    alternative) {
  if (is.null(group_col) || statistic == "none") {
    return(list(summary = NULL, null = NULL))
  }

  observed_stat <- .calc_trajectory_stat(observed_summary, statistic = statistic)
  null_stats <- vapply(
    split(permutation_summary, permutation_summary$permutation),
    .calc_trajectory_stat,
    numeric(1),
    statistic = statistic
  )

  list(
    summary = data.frame(
      statistic = statistic,
      observed = observed_stat,
      null_mean = mean(null_stats),
      null_sd = stats::sd(null_stats),
      empirical_p = .empirical_p_value(observed_stat, null_stats, alternative = alternative),
      stringsAsFactors = FALSE
    ),
    null = data.frame(
      permutation = names(null_stats),
      statistic = statistic,
      statistic_value = as.numeric(null_stats),
      stringsAsFactors = FALSE
    )
  )
}

.calc_trajectory_stat <- function(summary_df, statistic = "last_minus_first") {
  statistic <- match.arg(statistic, c("last_minus_first", "slope"))
  summary_df <- summary_df[order(summary_df$group), , drop = FALSE]

  if (statistic == "last_minus_first") {
    return(utils::tail(summary_df$score, 1L) - summary_df$score[[1L]])
  }

  x <- seq_len(nrow(summary_df))
  stats::coef(stats::lm(summary_df$score ~ x))[[2L]]
}

.empirical_p_value <- function(observed, null, alternative = "two.sided") {
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) == 0L) {
    return(NA_real_)
  }

  # Add-one (Phipson & Smyth 2010) empirical p-value: count null statistics at
  # least as extreme as the observed, then add 1 to numerator and denominator so
  # the p-value is never exactly 0.
  if (alternative == "greater") {
    return((sum(null >= observed) + 1) / (length(null) + 1))
  }
  if (alternative == "less") {
    return((sum(null <= observed) + 1) / (length(null) + 1))
  }

  # Two-sided: extremeness is distance from the null center.
  null_center <- mean(null)
  (sum(abs(null - null_center) >= abs(observed - null_center)) + 1) / (length(null) + 1)
}

.build_perm_plot_data <- function(observed_summary, permutation_summary, null_median_summary) {
  observed_for_plot <- observed_summary
  observed_for_plot$type <- "observed"
  permutation_for_plot <- permutation_summary
  permutation_for_plot$type <- "permutation"
  null_median_for_plot <- null_median_summary
  null_median_for_plot$type <- "permutation_median"

  plot_cols <- Reduce(intersect, list(
    colnames(permutation_for_plot),
    colnames(null_median_for_plot),
    colnames(observed_for_plot)
  ))
  rbind(
    permutation_for_plot[, plot_cols, drop = FALSE],
    null_median_for_plot[, plot_cols, drop = FALSE],
    observed_for_plot[, plot_cols, drop = FALSE]
  )
}

.plot_perm_trajectory <- function(perm_result,
                                  observed_color = "#D7191C",
                                  null_color = "grey70",
                                  null_alpha = 0.35,
                                  null_line_width = 0.18,
                                  median_line_width = 0.8,
                                  observed_line_width = 0.9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to create the permutation plot.", call. = FALSE)
  }

  plot_df <- perm_result$plot_data
  if (is.null(plot_df) || nrow(plot_df) == 0L) {
    stop("perm_result does not contain plot_data.", call. = FALSE)
  }

  x_col <- perm_result$params$plot_x
  x_label <- if (identical(x_col, "group")) "Timepoint (in days)" else "Sample"
  perm_df <- plot_df[plot_df$type == "permutation", , drop = FALSE]
  median_df <- plot_df[plot_df$type == "permutation_median", , drop = FALSE]
  obs_df <- plot_df[plot_df$type == "observed", , drop = FALSE]
  x_values <- plot_df[[x_col]]
  if (is.factor(x_values)) {
    x_breaks <- seq_along(levels(x_values))
    x_labels <- levels(x_values)
    plot_df$.plot_x_value <- as.numeric(x_values)
  } else if (is.character(x_values)) {
    x_factor <- factor(x_values, levels = unique(x_values))
    x_breaks <- seq_along(levels(x_factor))
    x_labels <- levels(x_factor)
    plot_df$.plot_x_value <- as.numeric(x_factor)
  } else {
    x_breaks <- sort(unique(x_values))
    x_labels <- x_breaks
    plot_df$.plot_x_value <- x_values
  }
  perm_df <- plot_df[plot_df$type == "permutation", , drop = FALSE]
  median_df <- plot_df[plot_df$type == "permutation_median", , drop = FALSE]
  obs_df <- plot_df[plot_df$type == "observed", , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_smooth(
      data = perm_df,
      ggplot2::aes(x = .data$.plot_x_value, y = .data$score, group = .data$permutation),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.55,
      color = null_color,
      alpha = null_alpha,
      linewidth = null_line_width
    ) +
    ggplot2::geom_smooth(
      data = median_df,
      ggplot2::aes(x = .data$.plot_x_value, y = .data$score, group = 1),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.55,
      color = "black",
      linewidth = median_line_width
    ) +
    ggplot2::geom_smooth(
      data = obs_df,
      ggplot2::aes(x = .data$.plot_x_value, y = .data$score, group = 1),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.55,
      color = observed_color,
      linewidth = observed_line_width
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = base::range(x_breaks),
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::labs(
      title = paste0(perm_result$params$module_name, " module-score permutation test"),
      x = x_label,
      y = "Module Score"
    ) +
    ggplot2::theme_classic(base_size = 9)
}

.plot_perm_distribution <- function(perm_result,
                                    fill = "grey75",
                                    observed_color = "#D7191C") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to create the permutation distribution plot.", call. = FALSE)
  }
  if (is.null(perm_result$trajectory) || is.null(perm_result$trajectory_null)) {
    stop(
      "perm_result must contain trajectory and trajectory_null. ",
      "Run with trajectory_stat other than 'none'.",
      call. = FALSE
    )
  }

  observed <- perm_result$trajectory$observed[[1L]]
  p_value <- perm_result$trajectory$empirical_p[[1L]]
  null_values <- perm_result$trajectory_null$statistic_value
  x_range <- base::range(null_values, observed, na.rm = TRUE)
  histogram_counts <- graphics::hist(null_values, breaks = "FD", plot = FALSE)$counts
  p_label_x <- x_range[[1L]] + base::diff(x_range) * 0.03
  p_label_y <- base::max(histogram_counts, na.rm = TRUE) * 0.84

  ggplot2::ggplot(perm_result$trajectory_null, ggplot2::aes(x = .data$statistic_value)) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = fill,
      color = "grey30",
      alpha = 0.85,
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(
      xintercept = observed,
      color = observed_color,
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    ggplot2::annotate(
      "text",
      x = p_label_x,
      y = p_label_y,
      label = paste0("empirical p = ", signif(p_value, 3)),
      hjust = 0,
      vjust = 1,
      family = "Nimbus Sans",
      size = 6 / ggplot2::.pt
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::labs(
      title = paste0(perm_result$params$module_name, " permutation null distribution"),
      subtitle = paste0(
        "Gene set size = ", perm_result$params$n_present_genes,
        "; permutations = ", perm_result$params$n_perm
      ),
      x = "Null statistic value",
      y = "Count"
    ) +
    ggplot2::theme_classic(base_size = 9)
}
