#!/usr/bin/env Rscript
# Created:
# 2026-05-28
#
# Inputs:
# - functions/ and load_functions.R
# - data/GSE122380_metadata.csv
# - data/GSE122380_vst.csv
#
# Outputs:
# - results/GSE122380_modulescore_public.csv
# - results/GSE122380_modulescore_custom.csv
# - results/GSE122380_modulescore_custom_marker_permutation.csv
# - report/assets/figures/reference_pca_and_sample_correlation.png
# - report/assets/figures/reference_pca_and_sample_correlation.svg
# - report/assets/figures/public_gene_set_module_score_run.png
# - report/assets/figures/public_gene_set_module_score_run.svg
# - report/assets/figures/custom_module_score_run.png
# - report/assets/figures/custom_module_score_run.svg
# - report/assets/figures/custom_marker_score_heatmap.png
# - report/assets/figures/custom_marker_score_heatmap.svg
# - report/assets/figures/custom_marker_permutation_run.png
# - report/assets/figures/custom_marker_permutation_run.svg
# - report/assets/figures/public_score_relationship.png
# - report/assets/figures/public_score_relationship.svg
# - report/assets/figures/expression_bin_visualization.png
# - report/assets/figures/expression_bin_visualization.svg
# - report/assets/figures/module_score_simulation.png
# - report/assets/figures/module_score_simulation.svg
#
# Purpose:
# Run the public tutorial analysis for bulk RNA-seq module scoring.
#
# Notes:
# Uses a processed GSE122380 cardiomyocyte differentiation matrix. The workflow
# scores public GO/MSigDB modules, scores custom developmental marker modules,
# and creates the tutorial figures used by the bookdown report.

# 1.1 load scoring functions and tutorial dependencies -----------------

required_packages <- base::c(
  'dplyr',
  'ggplot2',
  'GSVA',
  'magrittr',
  'msigdbr',
  'patchwork',
  'readr',
  'base64enc',
  'ggrastr',
  'ggtext',
  'scales',
  'singscore',
  'svglite',
  'tibble',
  'tidyr',
  'viridisLite')

missing_packages <- required_packages[
  !base::vapply(required_packages, base::requireNamespace, logical(1), quietly = TRUE)]
if (base::length(missing_packages) > 0L) {
  base::stop(
    'Install required packages before running the workflow: ',
    base::paste(missing_packages, collapse = ', '),
    call. = FALSE)
}

base::source('load_functions.R', local = FALSE)
base::suppressPackageStartupMessages(base::library(magrittr))

# 1.2 define tutorial figure helper functions -----------------

make_module_specs <- function(module_order, labels, source_ids, gene_sets, colors, gradients) {
  module_specs <- tibble::tibble(
    module_id = module_order,
    module_label = base::unname(labels[module_order]),
    source_id = base::unname(source_ids[module_order]),
    n_genes = base::vapply(gene_sets[module_order], base::length, integer(1)),
    fill = base::unname(colors[module_order]),
    panel_label = letters[base::seq_along(module_order)],
    gradient = gradients[module_order])
  module_specs %>%
    dplyr::mutate(
      title = dplyr::if_else(
        .data$source_id == 'custom marker set',
        base::paste0(.data$module_label, ' (n=', .data$n_genes, ')'),
        base::paste0(.data$source_id, ' ', base::tolower(.data$module_label), ' (n=', .data$n_genes, ')')))
}

format_module_score_labels <- function(x) {
  formatted <- base::ifelse(
    base::abs(x) < 1e-10,
    '0',
    base::ifelse(
      base::abs(x - base::round(x)) < 1e-10,
      base::as.character(base::round(x)),
      base::formatC(x, format = 'f', digits = 2)))
  base::sub('\\.?0+$', '', formatted)
}

gradient_end_mid_breaks <- function(x) {
  x <- x[base::is.finite(x)]
  if (base::length(x) == 0L) {
    return(NULL)
  }
  x_range <- gradient_end_limits(x)
  step <- 0.05
  middle <- base::round(base::mean(x_range) / step) * step
  base::unique(base::round(base::c(x_range[[1]], middle, x_range[[2]]), digits = 2))
}

gradient_end_limits <- function(x) {
  x <- x[base::is.finite(x)]
  if (base::length(x) == 0L) {
    return(NULL)
  }
  step <- 0.05
  x_range <- base::range(x)
  base::round(
    base::c(
      base::floor(x_range[[1]] / step) * step,
      base::ceiling(x_range[[2]] / step) * step),
    digits = 2)
}

score_to_long <- function(score_df, module_specs) {
  score_df %>%
    dplyr::select(dplyr::all_of(base::c('sample_id', 'day', module_specs$module_id))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(module_specs$module_id),
      names_to = 'module_id',
      values_to = 'module_score') %>%
    dplyr::left_join(module_specs[, base::c('module_id', 'module_label', 'title')], by = 'module_id') %>%
    dplyr::mutate(
      module_id = base::factor(.data$module_id, levels = module_specs$module_id),
      title = base::factor(.data$title, levels = module_specs$title))
}

plot_module_boxplot <- function(score_data, module_id, module_specs) {
  module_spec <- module_specs[module_specs$module_id == module_id, , drop = FALSE]
  plot_data <- score_data %>%
    dplyr::filter(.data$module_id == !!module_id)

  plot_bulk_module_score_boxplot(
    score_df = plot_data,
    score_col = 'module_score',
    x_col = 'day',
    x_order = day_order,
    fill = module_spec$fill[[1]],
    title = module_spec$title[[1]],
    tag = module_spec$panel_label[[1]],
    x_label = 'Timepoint (in days)',
    y_label = 'Bulk module score',
    base_size = fs(7),
    base_family = figure_family) +
    ggplot2::labs(
      title = module_spec$title[[1]],
      tag = module_spec$panel_label[[1]],
      x = 'Timepoint (in days)',
      y = 'Bulk module score') +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'plain', size = fs(7), lineheight = 1.05),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.01, 0.99),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.title.y = ggplot2::element_text(size = fs(7), margin = ggplot2::margin(r = 0)),
      axis.text = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(8, 6, 14, 6))
}

plot_pca_score <- function(pca_data, day_path, module_id, module_specs, x_label, y_label) {
  module_spec <- module_specs[module_specs$module_id == module_id, , drop = FALSE]
  pca_x_range <- base::range(pca_data$PC1, na.rm = TRUE)
  pca_y_range <- base::range(pca_data$PC2, na.rm = TRUE)
  pca_x_span <- base::diff(pca_x_range)
  pca_y_span <- base::diff(pca_y_range)
  dev_arrow <- tibble::tibble(
    x = pca_x_range[[1]] + 0.06 * pca_x_span,
    xend = pca_x_range[[2]] - 0.06 * pca_x_span,
    y = pca_y_range[[1]] + 0.06 * pca_y_span,
    label_x = base::mean(pca_x_range),
    label_y = pca_y_range[[1]] + 0.135 * pca_y_span)

  ggplot2::ggplot(pca_data, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
    ggplot2::geom_segment(
      data = dev_arrow,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      inherit.aes = FALSE,
      color = 'black',
      linewidth = gs(0.28),
      arrow = grid::arrow(length = grid::unit(0.065, 'in'), type = 'closed')) +
    ggplot2::geom_text(
      data = dev_arrow,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = 'Developmental time'),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      color = 'black') +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[[module_id]]),
      size = gs(1.05),
      alpha = 0.96) +
    ggplot2::scale_color_gradientn(
      colors = module_spec$gradient[[1]],
      name = 'Bulk module score',
      limits = gradient_end_limits(pca_data[[module_id]]),
      breaks = gradient_end_mid_breaks,
      labels = format_module_score_labels,
      guide = ggplot2::guide_colorbar(
        direction = 'horizontal',
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(0.95, 'in'),
        barheight = grid::unit(0.06, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
        theme = ggplot2::theme(
          legend.title = ggplot2::element_text(face = 'bold'),
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')))) +
    ggplot2::labs(
      title = module_spec$title[[1]],
      x = x_label,
      y = y_label) +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'plain', size = fs(7), lineheight = 1.05),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      legend.position = base::c(0.97, 0.95),
      legend.justification = base::c(1, 1),
      legend.direction = 'horizontal',
      legend.title = ggplot2::element_text(size = fs(6), face = 'bold'),
      legend.title.align = 0.5,
      legend.text = ggplot2::element_text(size = fs(5.5)),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.78), color = NA),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(8, 6, 14, 6))
}

plot_reference_pca <- function(pca_data, day_path, x_label, y_label, pca_gene_fraction, pca_gene_count) {
  reference_day_colors <- stats::setNames(reference_day_palette, 1:15)
  phase_colors <- base::c(
    'Pluripotent' = reference_day_colors[['1']],
    'Mesoderm' = reference_day_colors[['4']],
    'Immature cardiomyocyte' = reference_day_colors[['8']],
    'Mature cardiomyocyte' = reference_day_colors[['15']])
  phase_label_fills <- phase_colors
  phase_label_fills[['Mature cardiomyocyte']] <- '#B8860B'

  pca_data <- pca_data %>%
    dplyr::mutate(
      phase = dplyr::case_when(
        .data$day <= 2 ~ 'Pluripotent',
        .data$day <= 5 ~ 'Mesoderm',
        .data$day <= 9 ~ 'Immature cardiomyocyte',
        TRUE ~ 'Mature cardiomyocyte'),
      phase = base::factor(.data$phase, levels = base::names(phase_colors)))

  ellipse_df <- dplyr::bind_rows(base::lapply(base::names(phase_colors), function(phase_name) {
    phase_data <- pca_data[pca_data$phase == phase_name, base::c('PC1', 'PC2'), drop = FALSE]
    if (base::nrow(phase_data) < 3L) {
      return(NULL)
    }

    center <- base::colMeans(phase_data)
    covariance <- stats::cov(phase_data)
    ellipse_angle <- base::seq(0, 2 * base::pi, length.out = 160)
    ellipse_circle <- base::rbind(base::cos(ellipse_angle), base::sin(ellipse_angle))
    ellipse_coords <- base::t(center + base::sqrt(stats::qchisq(0.80, df = 2)) *
      base::t(base::chol(covariance)) %*% ellipse_circle)

    tibble::tibble(
      phase = phase_name,
      PC1 = ellipse_coords[, 1],
      PC2 = ellipse_coords[, 2])
  }))

  ellipse_layers <- base::lapply(base::names(phase_colors), function(phase_name) {
    ggplot2::geom_path(
      data = ellipse_df[ellipse_df$phase == phase_name, , drop = FALSE],
      ggplot2::aes(x = .data$PC1, y = .data$PC2),
      inherit.aes = FALSE,
      color = phase_colors[[phase_name]],
      linewidth = gs(0.50),
      alpha = 0.9)
  })

  day_label_offsets <- tibble::tibble(
    day = 1:15,
    nudge_x = base::c(-2, -4, -2, 0, 0, 0, -2, -2, -8, 5, -10, 8, -9, 10, 0),
    nudge_y = base::c(-6, 5, -5, 5, -5, 6, -7, 6, -5, 5, -9, 0, -2, -5, -12))

  day_path <- day_path %>%
    dplyr::left_join(day_label_offsets, by = 'day') %>%
    dplyr::mutate(
      day_label = base::paste0('D', .data$day),
      label_x = .data$PC1 + .data$nudge_x,
      label_y = .data$PC2 + .data$nudge_y)

  phase_labels <- tibble::tibble(
    day = base::c(1, 4, 8, 14),
    phase = base::c('Pluripotent', 'Mesoderm', 'Immature cardiomyocyte', 'Mature cardiomyocyte'),
    nudge_x = base::c(-4, -22, 8, 2),
    nudge_y = base::c(-12, 4, 19, -16))
  phase_labels <- phase_labels %>%
    dplyr::left_join(day_path[, base::c('day', 'PC1', 'PC2')], by = 'day') %>%
    dplyr::mutate(
      label_x = .data$PC1 + .data$nudge_x,
      label_y = .data$PC2 + .data$nudge_y)

  plot_bounds <- dplyr::bind_rows(
    pca_data %>% dplyr::transmute(PC1 = .data$PC1, PC2 = .data$PC2),
    ellipse_df %>% dplyr::transmute(PC1 = .data$PC1, PC2 = .data$PC2),
    day_path %>% dplyr::transmute(PC1 = .data$label_x, PC2 = .data$label_y),
    phase_labels %>% dplyr::transmute(PC1 = .data$label_x, PC2 = .data$label_y))
  pca_xlim <- base::range(plot_bounds$PC1, na.rm = TRUE) + base::c(-8, 8)
  pca_ylim <- base::range(plot_bounds$PC2, na.rm = TRUE) + base::c(-8, 8)

  ggplot2::ggplot(pca_data, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
    ellipse_layers +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$day),
      size = gs(0.95),
      alpha = 0.92) +
    ggplot2::geom_text(
      data = day_path,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$day_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      color = 'black') +
    ggplot2::geom_label(
      data = phase_labels,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$phase, fill = .data$phase),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      color = 'white',
      linewidth = 0,
      label.r = grid::unit(0.08, 'lines'),
      label.padding = grid::unit(0.18, 'lines')) +
    ggplot2::scale_color_gradientn(
      colors = reference_day_palette,
      name = 'Differentiation day',
      breaks = base::c(1, 5, 10, 15),
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0,
        barwidth = grid::unit(0.78, 'in'),
        barheight = grid::unit(0.06, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')))) +
    ggplot2::scale_fill_manual(
      values = phase_label_fills,
      guide = 'none') +
    ggplot2::labs(
      title = 'PCA: GSE122380',
      subtitle = base::paste0('Top 10% variable genes (n=', scales::comma(pca_gene_count), ')'),
      tag = 'a',
      x = x_label,
      y = y_label) +
    ggplot2::coord_cartesian(
      xlim = pca_xlim,
      ylim = pca_ylim,
      clip = 'off') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'plain', size = fs(7)),
      plot.subtitle = ggplot2::element_text(face = 'plain', size = fs(6)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      legend.position = base::c(0.98, 0.98),
      legend.justification = base::c(1, 1),
      legend.direction = 'horizontal',
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
      legend.margin = ggplot2::margin(1, 1, 1, 1),
      legend.title = ggplot2::element_text(size = fs(6), face = 'bold'),
      legend.text = ggplot2::element_text(size = fs(5.5)),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 6, 10))
}

make_timepoint_correlation_matrix <- function(expression_mat, metadata) {
  day_order <- base::sort(base::unique(metadata$day))
  day_mean_mat <- base::sapply(day_order, function(day_value) {
    sample_ids <- metadata$sample_id[metadata$day == day_value]
    base::rowMeans(expression_mat[, sample_ids, drop = FALSE], na.rm = TRUE)
  })
  base::colnames(day_mean_mat) <- base::paste0('D', day_order)

  stats::cor(day_mean_mat, method = 'pearson', use = 'pairwise.complete.obs')
}

get_day_colors <- function(days) {
  day_palette <- annotation_day_palette
  day_index <- base::round(scales::rescale(days, to = base::c(1, 256), from = base::range(days)))
  stats::setNames(day_palette[day_index], base::paste0('D', days))
}

plot_timepoint_correlation_heatmap <- function(correlation_mat) {
  ordered_labels <- base::paste0('D', base::seq_len(base::ncol(correlation_mat)))
  y_labels <- base::rev(ordered_labels)
  n_labels <- base::length(ordered_labels)
  gap_after <- base::ceiling(n_labels / 2)
  gap_size <- 0.08
  add_heatmap_gap <- function(index) {
    index + base::ifelse(index > gap_after, gap_size, 0)
  }
  x_positions <- add_heatmap_gap(base::seq_along(ordered_labels))
  y_positions <- add_heatmap_gap(base::seq_along(y_labels))
  max_x_position <- base::max(x_positions)
  max_y_position <- base::max(y_positions)

  heatmap_data <- base::as.data.frame(base::as.table(correlation_mat[ordered_labels, y_labels]))
  base::names(heatmap_data) <- base::c('x_label', 'y_label', 'correlation')
  heatmap_data <- heatmap_data %>%
    dplyr::mutate(
      x_index = add_heatmap_gap(base::match(.data$x_label, ordered_labels)),
      y_index = add_heatmap_gap(base::match(.data$y_label, y_labels)))

  ordered_days <- base::as.integer(base::sub('^D', '', ordered_labels))
  y_days <- base::as.integer(base::sub('^D', '', y_labels))
  day_colors <- get_day_colors(base::sort(base::unique(ordered_days)))
  top_annotation <- tibble::tibble(
    x_index = x_positions,
    y_index = max_y_position + 0.82,
    fill = base::unname(day_colors[ordered_labels]))
  left_annotation <- tibble::tibble(
    x_index = 0.15,
    y_index = y_positions,
    fill = base::unname(day_colors[base::paste0('D', y_days)]))
  top_labels <- tibble::tibble(
    x_index = x_positions,
    y_index = max_y_position + 1.37,
    label = ordered_labels)

  ggplot2::ggplot() +
    ggplot2::annotate(
      'segment',
      x = 1,
      xend = max_x_position,
      y = max_y_position + 2.45,
      yend = max_y_position + 2.45,
      linewidth = gs(0.20),
      arrow = grid::arrow(length = grid::unit(0.055, 'in'), type = 'closed')) +
    ggplot2::annotate(
      'text',
      x = (max_x_position + 1) / 2,
      y = max_y_position + 2.84,
      label = 'Developmental time',
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5)) +
    ggplot2::geom_tile(
      data = heatmap_data,
      ggplot2::aes(x = .data$x_index, y = .data$y_index, fill = .data$correlation),
      width = 1.01,
      height = 1.01,
      color = NA,
      linewidth = 0) +
    ggplot2::geom_tile(
      data = top_annotation,
      ggplot2::aes(x = .data$x_index, y = .data$y_index),
      fill = top_annotation$fill,
      width = 1,
      height = 0.55,
      color = NA) +
    ggplot2::geom_tile(
      data = left_annotation,
      ggplot2::aes(x = .data$x_index, y = .data$y_index),
      fill = left_annotation$fill,
      width = 0.55,
      height = 1,
      color = NA) +
    ggplot2::geom_text(
      data = top_labels,
      ggplot2::aes(x = .data$x_index, y = .data$y_index, label = .data$label),
      inherit.aes = FALSE,
      family = figure_family,
      size = gfs(5.5),
      angle = 45,
      hjust = 0,
      vjust = 0.5,
      color = 'grey25') +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      position = 'top',
      limits = base::c(-0.16, max_x_position + 0.51),
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_positions,
      labels = y_labels,
      limits = base::c(0.49, max_y_position + 3.09),
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_fill_gradientn(
      colors = correlation_palette,
      limits = base::c(0, 1),
      name = 'Pearson r',
      breaks = base::c(0, 0.5, 1),
      labels = function(x) base::ifelse(
        x %% 1 == 0,
        base::formatC(x, format = 'f', digits = 0),
        base::formatC(x, format = 'f', digits = 1)),
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(1.25, 'in'),
        barheight = grid::unit(0.08, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')))) +
    ggplot2::labs(
      title = NULL,
      tag = 'b',
      x = NULL,
      y = NULL) +
    ggplot2::coord_fixed(clip = 'off') +
    ggplot2::theme_minimal(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.text.x = ggplot2::element_text(size = fs(6), angle = 45, hjust = 0, vjust = 0.5, margin = ggplot2::margin(b = -2)),
      axis.text.y = ggplot2::element_text(size = fs(6)),
      legend.position = 'bottom',
      legend.justification = 'center',
      legend.title = ggplot2::element_text(size = fs(6), face = 'bold'),
      legend.title.align = 0.5,
      legend.text = ggplot2::element_text(size = fs(5.5)),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-24, 0, 0, 0),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 0, 10))
}

plot_module_pairs <- function(score_data, pca_data, day_path, module_specs, x_label, y_label) {
  paired_plots <- base::unlist(base::lapply(module_specs$module_id, function(module_id) {
    list(
      plot_module_boxplot(score_data = score_data, module_id = module_id, module_specs = module_specs),
      plot_pca_score(
        pca_data = pca_data,
        day_path = day_path,
        module_id = module_id,
        module_specs = module_specs,
        x_label = x_label,
        y_label = y_label))
  }), recursive = FALSE)

  patchwork::wrap_plots(paired_plots, ncol = 2)
}

plot_score_heatmap <- function(score_data, module_specs) {
  row_gap_size <- 0
  col_gap_size <- 0.16
  collapsed_score_data <- score_data %>%
    dplyr::group_by(.data$day, .data$module_id) %>%
    dplyr::summarise(module_score = base::mean(.data$module_score, na.rm = TRUE), .groups = 'drop')

  day_order <- base::sort(base::unique(collapsed_score_data$day))
  score_matrix <- base::matrix(
    NA_real_,
    nrow = base::nrow(module_specs),
    ncol = base::length(day_order),
    dimnames = list(module_specs$module_id, base::paste0('D', day_order)))
  for (score_idx in base::seq_len(base::nrow(collapsed_score_data))) {
    score_matrix[
      collapsed_score_data$module_id[[score_idx]],
      base::paste0('D', collapsed_score_data$day[[score_idx]])] <- collapsed_score_data$module_score[[score_idx]]
  }
  if (base::any(!base::is.finite(score_matrix))) {
    stop('custom score heatmap clustering matrix contains non-finite values')
  }
  column_cluster <- stats::hclust(stats::dist(base::t(score_matrix)), method = 'complete')
  clustered_labels <- column_cluster$labels[column_cluster$order]
  clustered_days <- base::as.integer(base::sub('^D', '', clustered_labels))
  column_clusters <- stats::cutree(column_cluster, k = 3)

  sample_order <- score_data %>%
    dplyr::distinct(.data$day) %>%
    dplyr::mutate(day_label = base::paste0('D', .data$day)) %>%
    dplyr::filter(.data$day_label %in% clustered_labels) %>%
    dplyr::arrange(base::match(.data$day, clustered_days)) %>%
    dplyr::mutate(
      column_cluster = base::unname(column_clusters[.data$day_label]),
      cluster_boundary = .data$column_cluster != dplyr::lag(.data$column_cluster, default = dplyr::first(.data$column_cluster)),
      cumulative_gap = base::cumsum(.data$cluster_boundary) * col_gap_size,
      sample_index_raw = dplyr::row_number(),
      sample_index = .data$sample_index_raw + .data$cumulative_gap)

  row_info <- module_specs %>%
    dplyr::mutate(
      y_index = base::rev(base::seq_along(.data$module_id)) * (1 + row_gap_size),
      title = base::as.character(.data$title))

  heatmap_data <- collapsed_score_data %>%
    dplyr::left_join(sample_order, by = 'day') %>%
    dplyr::left_join(row_info[, base::c('module_id', 'y_index')], by = 'module_id') %>%
    dplyr::mutate(module_score_winsor = base::pmax(base::pmin(.data$module_score, 4), -4))

  day_colors <- get_day_colors(base::sort(base::unique(sample_order$day)))
  top_annotation <- sample_order %>%
    dplyr::mutate(
      y_index = base::max(row_info$y_index) + 0.48,
      fill = base::unname(day_colors[base::paste0('D', .data$day)]))

  get_dendrogram_segments <- function(cluster, y_base, y_height) {
    dendrogram <- stats::as.dendrogram(cluster)
    leaf_x <- stats::setNames(sample_order$sample_index, sample_order$day_label)
    max_height <- base::max(cluster$height, na.rm = TRUE)
    if (!base::is.finite(max_height) || max_height == 0) {
      max_height <- 1
    }

    collect_segments <- function(node) {
      node_height <- base::attr(node, 'height')
      if (base::isTRUE(base::attr(node, 'leaf'))) {
        node_label <- base::attr(node, 'label')
        return(list(
          x = base::unname(leaf_x[[node_label]]),
          y = node_height,
          segments = tibble::tibble()))
      }

      children <- base::lapply(node, collect_segments)
      child_x <- base::vapply(children, `[[`, numeric(1), 'x')
      child_y <- base::vapply(children, `[[`, numeric(1), 'y')
      node_x <- base::mean(child_x)
      child_segments <- dplyr::bind_rows(base::lapply(children, `[[`, 'segments'))
      vertical_segments <- tibble::tibble(
        x = child_x,
        xend = child_x,
        y = child_y,
        yend = node_height)
      horizontal_segment <- tibble::tibble(
        x = base::min(child_x),
        xend = base::max(child_x),
        y = node_height,
        yend = node_height)

      list(
        x = node_x,
        y = node_height,
        segments = dplyr::bind_rows(child_segments, vertical_segments, horizontal_segment))
    }

    collect_segments(dendrogram)$segments %>%
      dplyr::mutate(
        y = y_base + (.data$y / max_height) * y_height,
        yend = y_base + (.data$yend / max_height) * y_height)
  }

  dendrogram_segments <- get_dendrogram_segments(
    cluster = column_cluster,
    y_base = base::max(row_info$y_index) + 0.72,
    y_height = 0.82)

  max_sample_x <- base::max(sample_order$sample_index)
  annotation_x <- max_sample_x + 0.72
  label_x <- annotation_x + 0.28
  legend_values <- base::seq(-4, 4, length.out = 180)
  legend_x <- base::seq(4.35, 6.85, length.out = base::length(legend_values))
  legend_data <- tibble::tibble(
    x = legend_x,
    y = -0.66,
    module_score_winsor = legend_values)
  timepoint_legend_values <- base::seq(
    base::min(sample_order$day),
    base::max(sample_order$day),
    length.out = 180)
  timepoint_legend_x <- base::seq(9.45, 12.95, length.out = base::length(timepoint_legend_values))
  timepoint_legend_data <- tibble::tibble(
    x = timepoint_legend_x,
    y = -0.66,
    day = timepoint_legend_values,
    fill = viridisLite::viridis(180))

  ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$sample_index, y = .data$y_index, fill = .data$module_score_winsor)) +
    ggplot2::geom_segment(
      data = dendrogram_segments,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      inherit.aes = FALSE,
      color = 'grey20',
      linewidth = gs(0.28),
      lineend = 'square') +
    ggplot2::geom_tile(
      data = top_annotation,
      ggplot2::aes(x = .data$sample_index, y = .data$y_index),
      inherit.aes = FALSE,
      fill = top_annotation$fill,
      color = NA,
      width = 1.01,
      height = 0.58) +
    ggplot2::geom_tile(
      data = legend_data,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$module_score_winsor),
      inherit.aes = FALSE,
      width = 0.04,
      height = 0.18,
      color = NA) +
    ggplot2::geom_tile(
      data = timepoint_legend_data,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = timepoint_legend_data$fill,
      width = 0.04,
      height = 0.18,
      color = NA) +
    ggplot2::annotate(
      'text',
      x = 4.10,
      y = -0.66,
      label = 'Bulk module score',
      family = figure_family,
      fontface = 'bold',
      hjust = 1,
      size = gfs(6),
      color = 'black') +
    ggplot2::annotate(
      'text',
      x = base::c(4.35, 5.60, 6.85),
      y = -0.93,
      label = base::c('-4', '0', '4'),
      family = figure_family,
      size = gfs(5.5),
      color = 'black') +
    ggplot2::annotate(
      'text',
      x = 9.20,
      y = -0.66,
      label = 'Day',
      family = figure_family,
      fontface = 'bold',
      hjust = 1,
      size = gfs(6),
      color = 'black') +
    ggplot2::annotate(
      'text',
      x = base::c(9.45, 11.20, 12.95),
      y = -0.93,
      label = base::c('1', '8', '15'),
      family = figure_family,
      size = gfs(5.5),
      color = 'black') +
    ggplot2::geom_tile(
      data = row_info,
      ggplot2::aes(x = annotation_x, y = .data$y_index),
      inherit.aes = FALSE,
      fill = row_info$fill,
      color = NA,
      width = 0.36,
      height = 1) +
    ggplot2::geom_tile(color = NA, width = 1.01, height = 1.01) +
    ggplot2::annotate(
      'text',
      x = 0.28,
      y = base::max(row_info$y_index) + 1.70,
      label = 'a',
      family = panel_tag_family,
      fontface = 'bold',
      size = gfs(8),
      color = 'black') +
    ggplot2::geom_text(
      data = row_info,
      ggplot2::aes(x = label_x, y = .data$y_index, label = .data$title, color = .data$module_id),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      hjust = 0,
      size = gfs(6)) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      position = 'top',
      limits = base::c(-0.55, max_sample_x + 5.8),
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      expand = ggplot2::expansion(mult = 0, add = base::c(0.18, 0.72))) +
    ggplot2::scale_fill_gradientn(
      colors = correlation_palette,
      limits = base::c(-4, 4),
      breaks = base::c(-4, 0, 4),
      labels = format_module_score_labels,
      guide = 'none') +
    ggplot2::scale_color_manual(
      values = stats::setNames(module_specs$fill, module_specs$module_id),
      guide = 'none') +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1, clip = 'off') +
    ggplot2::theme_minimal(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x.top = ggplot2::element_text(size = fs(7), face = 'bold', margin = ggplot2::margin(b = 0)),
      axis.text.x = ggplot2::element_text(size = fs(6), margin = ggplot2::margin(b = -10)),
      axis.text.y = ggplot2::element_text(size = fs(6), face = 'bold', color = row_info$fill),
      legend.position = 'none',
      plot.margin = ggplot2::margin(14, 12, 2, 8))
}

plot_score_arrow_relationship <- function(data,
                                          x_col,
                                          y_col,
                                          title,
                                          x_label,
                                          y_label,
                                          x_short,
                                          y_short,
                                          legend_title = 'Differentiation day',
                                          show_legend = TRUE,
                                          panel_aspect_ratio = NULL,
                                          stats_label_y = 1,
                                          title_size = fs(7),
                                          axis_title_size = fs(7),
                                          axis_text_size = fs(6),
                                          legend_title_size = fs(6),
                                          legend_text_size = fs(5.5),
                                          quadrant_label_size = gfs(5.5),
                                          stats_label_size = gfs(6),
                                          plot_base_size = fs(7),
                                          tag = NULL) {
  complete_data <- data %>%
    dplyr::filter(base::is.finite(.data[[x_col]]), base::is.finite(.data[[y_col]]))
  r_value <- stats::cor(complete_data[[x_col]], complete_data[[y_col]])
  r_squared <- r_value^2
  fit <- stats::lm(stats::reformulate(x_col, response = y_col), data = complete_data)
  day_medians <- complete_data %>%
    dplyr::group_by(.data$day) %>%
    dplyr::summarise(x = stats::median(.data[[x_col]], na.rm = TRUE), .groups = 'drop') %>%
    dplyr::arrange(.data$day)
  arrow_start_x <- day_medians$x[[1L]]
  arrow_end_x <- day_medians$x[[base::nrow(day_medians)]]
  arrow_data <- tibble::tibble(
    x = arrow_start_x,
    xend = arrow_end_x,
    y = stats::predict(fit, newdata = stats::setNames(data.frame(arrow_start_x), x_col)),
    yend = stats::predict(fit, newdata = stats::setNames(data.frame(arrow_end_x), x_col)))
  if (base::isTRUE(base::all.equal(arrow_data$x, arrow_data$xend))) {
    x_quantiles <- stats::quantile(complete_data[[x_col]], probs = base::c(0.10, 0.90), na.rm = TRUE)
    arrow_data$x <- x_quantiles[[1L]]
    arrow_data$xend <- x_quantiles[[2L]]
    arrow_data$y <- stats::predict(fit, newdata = stats::setNames(data.frame(arrow_data$x), x_col))
    arrow_data$yend <- stats::predict(fit, newdata = stats::setNames(data.frame(arrow_data$xend), x_col))
  }
  x_limits <- base::range(base::c(complete_data[[x_col]], arrow_data$x, arrow_data$xend, 0), na.rm = TRUE)
  y_limits <- base::range(base::c(complete_data[[y_col]], arrow_data$y, arrow_data$yend, 0), na.rm = TRUE)
  x_limit <- base::max(base::abs(x_limits), na.rm = TRUE) * 1.08
  y_limit <- base::max(base::abs(y_limits), na.rm = TRUE) * 1.08
  if (!base::is.finite(x_limit) || x_limit == 0) {
    x_limit <- 1
  }
  if (!base::is.finite(y_limit) || y_limit == 0) {
    y_limit <- 1
  }
  x_limits <- base::c(-x_limit, x_limit)
  y_limits <- base::c(-y_limit, y_limit)
  x_span <- base::diff(x_limits)
  y_span <- base::diff(y_limits)

  quadrant_labels <- tibble::tibble(
    x = base::c(
      x_limits[[2]] - 0.135 * x_span,
      x_limits[[1]] + 0.135 * x_span,
      x_limits[[1]] + 0.135 * x_span,
      x_limits[[2]] - 0.135 * x_span),
    y = base::c(
      y_limits[[2]] - 0.14 * y_span,
      y_limits[[2]] - 0.14 * y_span,
      y_limits[[1]] + 0.14 * y_span,
      y_limits[[1]] + 0.14 * y_span),
    label = base::c(
      base::paste('High', x_short, '\nHigh', y_short),
      base::paste('Low', x_short, '\nHigh', y_short),
      base::paste('Low', x_short, '\nLow', y_short),
      base::paste('High', x_short, '\nLow', y_short)))

  ggplot2::ggplot(complete_data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data$day)) +
    ggplot2::geom_hline(yintercept = 0, color = 'black', linewidth = gs(0.24), linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = 0, color = 'black', linewidth = gs(0.24), linetype = 'dashed') +
    ggplot2::geom_point(size = gs(0.85), alpha = 1) +
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = gs(0.34)) +
    ggplot2::geom_label(
      data = quadrant_labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      family = figure_family,
      size = quadrant_label_size,
      lineheight = 0.9,
      color = 'black',
      fill = scales::alpha('white', 0.82),
      linewidth = gs(0.22),
      label.r = grid::unit(0.04, 'lines'),
      label.padding = grid::unit(0.12, 'lines')) +
    ggplot2::annotate(
      'text',
      x = x_limits[[2]] - 0.04 * x_span,
      y = y_limits[[2]] - 0.04 * y_span,
      label = base::sprintf('italic(r) == %.2f*","~~italic(R)^2 == %.2f', r_value, r_squared),
      parse = TRUE,
      hjust = 1,
      vjust = 0.5,
      size = stats_label_size,
      family = figure_family,
      color = 'black') +
    ggplot2::scale_color_viridis_c(
      option = 'D',
      name = legend_title,
      breaks = base::c(1, 8, 15),
      guide = ggplot2::guide_colorbar(
        direction = 'horizontal',
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(1.25, 'in'),
        barheight = grid::unit(0.08, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')))) +
    ggplot2::labs(
      title = title,
      tag = tag,
      x = x_label,
      y = y_label) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = y_limits) +
    ggplot2::theme_classic(base_size = plot_base_size, base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'plain', size = title_size),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.position = base::ifelse(show_legend, 'bottom', 'none'),
      legend.title = ggplot2::element_text(size = legend_title_size, face = 'bold'),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-4, 0, 0, 0),
      aspect.ratio = panel_aspect_ratio,
      plot.margin = ggplot2::margin(6, 7, 6, 7))
}

plot_public_score_correlation <- function(score_df) {
  p_mesoderm_pluripotency <- plot_score_arrow_relationship(
    data = score_df,
    x_col = 'embryonic_stem_cell_signature',
    y_col = 'paraxial_mesoderm',
    title = 'Pluripotency signature versus mesoderm',
    x_label = 'Pluripotency signature bulk module score',
    y_label = 'Mesoderm bulk module score',
    x_short = 'pluripotency',
    y_short = 'mesoderm',
    show_legend = FALSE,
    panel_aspect_ratio = 1,
    title_size = fs(7),
    axis_title_size = fs(7),
    axis_text_size = fs(6),
    legend_title_size = fs(6),
    legend_text_size = fs(5.5),
    quadrant_label_size = gfs(5.5),
    stats_label_size = gfs(6),
    plot_base_size = fs(7),
    tag = 'a')

  p_cardiac_mesoderm <- plot_score_arrow_relationship(
    data = score_df,
    x_col = 'paraxial_mesoderm',
    y_col = 'kegg_cardiac_muscle_contraction',
    title = 'Mesoderm versus cardiac muscle contraction',
    x_label = 'Mesoderm bulk module score',
    y_label = 'Cardiac muscle contraction bulk module score',
    x_short = 'mesoderm',
    y_short = 'cardiac',
    panel_aspect_ratio = 1,
    title_size = fs(7),
    axis_title_size = fs(7),
    axis_text_size = fs(6),
    legend_title_size = fs(6),
    legend_text_size = fs(5.5),
    quadrant_label_size = gfs(5.5),
    stats_label_size = gfs(6),
    plot_base_size = fs(7),
    tag = 'b')

  p_mesoderm_pluripotency + p_cardiac_mesoderm +
    patchwork::plot_layout(ncol = 2, guides = 'collect') &
    ggplot2::theme(
      legend.position = 'bottom',
      legend.justification = 'center',
      legend.box.margin = ggplot2::margin(-2, 0, 0, 0))
}

make_expression_bins <- function(expression_mat, nbin = 24) {
  gene_means <- base::rowMeans(expression_mat, na.rm = TRUE)
  breaks <- stats::quantile(gene_means, probs = base::seq(0, 1, length.out = nbin + 1L), na.rm = TRUE)
  if (base::length(base::unique(breaks)) == base::length(breaks)) {
    gene_bins <- base::cut(gene_means, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  } else {
    gene_ranks <- base::rank(gene_means, ties.method = 'average', na.last = 'keep')
    gene_bins <- base::ceiling(gene_ranks / base::max(gene_ranks, na.rm = TRUE) * nbin)
    gene_bins <- base::pmin(base::pmax(gene_bins, 1L), nbin)
  }
  tibble::tibble(
    gene = base::names(gene_means),
    mean_expression = base::as.numeric(gene_means),
    mean_cpm = base::pmax(2^base::as.numeric(gene_means), 0),
    bin = base::as.integer(gene_bins))
}

plot_expression_bin_visualization <- function(cpm_mat, nbin = 24, selected_bin = 12) {
  bin_data <- make_expression_bins(cpm_mat, nbin = nbin) %>%
    dplyr::filter(base::is.finite(.data$mean_expression))
  expressed_gene_count <- base::nrow(bin_data)
  average_genes_per_bin <- base::mean(base::table(bin_data$bin))
  stats_label <- base::sprintf(
    'Expressed genes: %s\nExpression bins: %s\nAvg genes/bin: %.0f',
    scales::comma(expressed_gene_count),
    nbin,
    average_genes_per_bin)
  selected_genes <- bin_data %>%
    dplyr::filter(.data$bin == selected_bin) %>%
    dplyr::arrange(.data$mean_expression) %>%
    dplyr::mutate(gene_index = dplyr::row_number())
  selected_gene_count <- base::nrow(selected_genes)
  total_gene_count <- 39239L
  expression_quantiles <- stats::quantile(
    bin_data$mean_expression,
    probs = base::c(0.25, 0.75),
    na.rm = TRUE)
  selected_stats_label <- base::sprintf(
    paste0(
      'Total genes: %s<br>',
      'Bottom 25%%: %.2f log<sub>2</sub>CPM<br>',
      'Top 75%%: %.2f log<sub>2</sub>CPM'),
    scales::comma(total_gene_count),
    expression_quantiles[[1L]],
    expression_quantiles[[2L]])
  selected_quantiles <- stats::quantile(
    selected_genes$mean_cpm,
    probs = base::c(0.25, 0.50, 0.75),
    na.rm = TRUE)
  selected_quantile_lines <- tibble::tibble(
    yintercept = base::as.numeric(selected_quantiles),
    label = base::c('25%', '50%', '75%'))
  bin_summary <- bin_data %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      median_expression = stats::median(.data$mean_expression, na.rm = TRUE),
      lower = base::min(.data$mean_expression, na.rm = TRUE),
      upper = base::max(.data$mean_expression, na.rm = TRUE),
      .groups = 'drop')
  bin_ranges <- bin_data %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      xmin = base::min(.data$mean_expression, na.rm = TRUE),
      xmax = base::max(.data$mean_expression, na.rm = TRUE),
      .groups = 'drop')
  histogram_breaks <- base::seq(
    base::min(bin_data$mean_expression, na.rm = TRUE),
    base::max(bin_data$mean_expression, na.rm = TRUE),
    length.out = 61L)
  histogram_counts <- graphics::hist(
    bin_data$mean_expression,
    breaks = histogram_breaks,
    plot = FALSE)
  histogram_data <- tibble::tibble(
    xmin = utils::head(histogram_counts$breaks, -1L),
    xmax = utils::tail(histogram_counts$breaks, -1L),
    xmid = histogram_counts$mids,
    count = histogram_counts$counts)
  histogram_data$bin <- base::vapply(histogram_data$xmid, function(xmid) {
    bin_match <- bin_ranges$bin[xmid >= bin_ranges$xmin & xmid <= bin_ranges$xmax]
    if (base::length(bin_match) == 0L) {
      bin_match <- bin_ranges$bin[base::which.min(base::abs(xmid - (bin_ranges$xmin + bin_ranges$xmax) / 2))]
    }
    bin_match[[1L]]
  }, integer(1))
  soften_colors <- function(colors, amount = 0.24, alpha = 0.88) {
    rgb_mat <- grDevices::col2rgb(colors) / 255
    rgb_mat <- rgb_mat + (1 - rgb_mat) * amount
    scales::alpha(grDevices::rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ]), alpha)
  }
  expression_palette <- soften_colors(correlation_palette, amount = 0.16, alpha = 1)
  histogram_palette <- soften_colors(correlation_palette, amount = 0.18, alpha = 0.94)
  histogram_x_limits <- base::range(bin_data$mean_expression, na.rm = TRUE)
  histogram_y_max <- base::max(histogram_data$count, na.rm = TRUE)
  legend_x_min <- histogram_x_limits[[1L]] + base::diff(histogram_x_limits) * 0.055
  legend_x_max <- legend_x_min + base::diff(histogram_x_limits) * 0.22
  legend_y_max <- histogram_y_max * 0.93
  legend_y_min <- histogram_y_max * 0.885
  legend_bar_data <- tibble::tibble(
    index = base::seq_len(80L),
    xmin = base::seq(legend_x_min, legend_x_max, length.out = 81L)[-81L],
    xmax = base::seq(legend_x_min, legend_x_max, length.out = 81L)[-1L],
    ymin = legend_y_min,
    ymax = legend_y_max,
    bin_value = base::seq(0, nbin, length.out = 80L))
  legend_label_data <- tibble::tibble(
    x = base::c(legend_x_min, (legend_x_min + legend_x_max) / 2, legend_x_max),
    y = legend_y_min - histogram_y_max * 0.055,
    label = base::c('0', '12', '24'))

  p_bins <- ggplot2::ggplot(bin_summary, ggplot2::aes(x = .data$bin, y = .data$median_expression)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$median_expression),
      shape = 16,
      size = gs(2.0),
      stroke = 0) +
    ggplot2::geom_point(
      data = bin_summary[bin_summary$bin == selected_bin, , drop = FALSE],
      ggplot2::aes(color = .data$median_expression),
      shape = 16,
      size = gs(2.8),
      stroke = 0) +
    ggplot2::annotate(
      'segment',
      x = selected_bin + 5.6,
      xend = selected_bin + 0.55,
      y = bin_summary$median_expression[bin_summary$bin == selected_bin] - 1.05,
      yend = bin_summary$median_expression[bin_summary$bin == selected_bin] - 0.08,
      linewidth = gs(0.24),
      color = 'black') +
    ggplot2::annotate(
      'text',
      x = selected_bin + 5.8,
      y = bin_summary$median_expression[bin_summary$bin == selected_bin] - 1.08,
      label = 'Bin 12',
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      hjust = 0,
      color = 'black') +
    ggplot2::annotate(
      'text',
      x = 1,
      y = base::max(bin_summary$upper, na.rm = TRUE),
      label = stats_label,
      family = figure_family,
      size = gfs(5.5),
      hjust = 0,
      vjust = 1,
      lineheight = 0.95,
      color = 'black') +
    ggplot2::scale_x_continuous(breaks = base::c(1, 6, 12, 18, 24)) +
    ggplot2::scale_color_gradientn(
      colors = expression_palette,
      guide = 'none') +
    ggplot2::labs(
      title = 'Expression matching into 24 bins',
      tag = 'b',
      x = 'Bin expression number',
      y = bquote('Mean log'[2] * '(CPM)')) +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = fs(7)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(8, 8, 8, 8))

  selected_color_limits <- base::range(selected_genes$mean_cpm, na.rm = TRUE)
  selected_color_breaks <- base::round(
    base::c(
      selected_color_limits[[1]],
      base::mean(selected_color_limits),
      selected_color_limits[[2]]),
    digits = 0)
  selected_color_breaks <- base::unique(selected_color_breaks)
  selected_color_labels <- base::formatC(selected_color_breaks, format = 'f', digits = 0)
  selected_y_limits <- base::range(selected_genes$mean_cpm, selected_quantile_lines$yintercept, na.rm = TRUE)
  selected_y_padding <- base::diff(selected_y_limits) * 0.18
  if (!base::is.finite(selected_y_padding) || selected_y_padding <= 0) {
    selected_y_padding <- 0.5
  }

  p_selected <- ggplot2::ggplot(selected_genes, ggplot2::aes(x = .data$gene_index, y = .data$mean_cpm)) +
    ggplot2::geom_hline(
      data = selected_quantile_lines,
      ggplot2::aes(yintercept = .data$yintercept),
      linewidth = gs(0.24),
      linetype = 'dashed',
      color = 'black',
      show.legend = FALSE) +
    ggplot2::geom_text(
      data = selected_quantile_lines,
      ggplot2::aes(
        x = selected_gene_count * 0.985,
        y = .data$yintercept,
        label = .data$label),
      inherit.aes = FALSE,
      parse = FALSE,
      family = figure_family,
      size = gfs(5),
      hjust = 1,
      vjust = -0.35,
      color = 'black') +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$mean_cpm),
      size = gs(0.62),
      alpha = 0.9) +
    ggtext::geom_richtext(
      data = tibble::tibble(
        x = selected_gene_count * 0.985,
        y = selected_y_limits[[2L]] + selected_y_padding * 0.92,
        label = selected_stats_label),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      family = figure_family,
      size = gfs(5),
      hjust = 1,
      vjust = 1,
      lineheight = 0.95,
      fill = scales::alpha('white', 0.82),
      label.color = NA,
      label.padding = grid::unit(base::c(0.08, 0.08, 0.08, 0.08), 'lines'),
      color = 'black') +
    ggplot2::scale_x_continuous(
      breaks = base::c(0, selected_gene_count),
      expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::scale_y_continuous(
      limits = base::c(
        selected_y_limits[[1L]] - selected_y_padding * 0.18,
        selected_y_limits[[2L]] + selected_y_padding),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_color_gradientn(
      colors = expression_palette,
      name = 'CPM',
      limits = selected_color_limits,
      breaks = selected_color_breaks,
      labels = selected_color_labels,
      guide = ggplot2::guide_colorbar(
        direction = 'horizontal',
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(0.70, 'in'),
        barheight = grid::unit(0.055, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.16),
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')))) +
    ggplot2::labs(
      title = base::paste0('Bin ', selected_bin, ': Genes ranked by CPM'),
      tag = 'c',
      x = base::paste0('Rank of genes within Bin ', selected_bin),
      y = 'CPM') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = fs(7)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      legend.position = base::c(0.02, 0.98),
      legend.justification = base::c(0, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
      legend.title = ggplot2::element_text(size = fs(6), face = 'bold', hjust = 0.5),
      legend.text = ggplot2::element_text(size = fs(5.5)),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(8, 8, 8, 8))

  histogram_overlap <- base::diff(base::range(histogram_breaks)) / base::length(histogram_breaks) * 0.035
  p_histogram <- ggplot2::ggplot(histogram_data) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin - histogram_overlap,
        xmax = .data$xmax + histogram_overlap,
        ymin = 0,
        ymax = .data$count,
        fill = .data$bin),
      color = NA,
      linewidth = 0) +
    ggrastr::rasterise(
      ggplot2::geom_rect(
        data = legend_bar_data,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data$bin_value),
        inherit.aes = FALSE,
        color = NA),
      dpi = figure_dpi) +
    ggplot2::annotate(
      'rect',
      xmin = legend_x_min,
      xmax = legend_x_max,
      ymin = legend_y_min,
      ymax = legend_y_max,
      fill = NA,
      color = 'black',
      linewidth = gs(0.18)) +
    ggplot2::annotate(
      'text',
      x = (legend_x_min + legend_x_max) / 2,
      y = legend_y_max + histogram_y_max * 0.055,
      label = 'Expression bin',
      family = figure_family,
      fontface = 'bold',
      size = gfs(6),
      hjust = 0.5,
      color = 'black') +
    ggplot2::geom_text(
      data = legend_label_data,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      family = figure_family,
      size = gfs(5.5),
      color = 'black') +
    ggplot2::scale_x_reverse(
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = base::c(0, 0.05))) +
    ggplot2::scale_fill_gradientn(
      colors = base::rev(histogram_palette),
      limits = base::c(0, nbin),
      guide = 'none') +
    ggplot2::labs(
      title = 'Distribution of gene expression, colored by expression bin',
      tag = 'a',
      x = bquote('Mean log'[2] * '(CPM)'),
      y = 'Genes') +
    ggplot2::coord_cartesian(
      xlim = base::rev(base::range(bin_data$mean_expression, na.rm = TRUE)),
      clip = 'off') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = fs(7)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(2, 74, 8, 74))

  p_histogram
}

simulate_module_score_example <- function() {
  base::set.seed(229)
  n_genes <- 6000L
  n_module <- 60L
  n_samples <- 48L
  sample_group <- base::rep(base::c('Condition A', 'Condition B'), each = n_samples / 2)
  gene_ids <- base::paste0('gene_', base::seq_len(n_genes))
  gene_baseline <- base::pmin(base::pmax(stats::rnorm(n_genes, mean = 5.1, sd = 1.25), 0.25), 10.5)
  background_loading <- 0.10 + 0.92 * stats::plogis((gene_baseline - 5.35) * 1.45)
  module_pool <- base::which(
    gene_baseline >= stats::quantile(gene_baseline, 0.72) &
      gene_baseline <= stats::quantile(gene_baseline, 0.90))
  module_genes <- gene_ids[base::sample(module_pool, n_module, replace = FALSE)]
  sample_background <- base::ifelse(
    sample_group == 'Condition B',
    stats::rnorm(n_samples, mean = 0.95, sd = 0.11),
    stats::rnorm(n_samples, mean = 0, sd = 0.11))
  sample_offset <- stats::rnorm(n_samples, mean = 0, sd = 0.08)
  expression_mat <- base::matrix(
    stats::rnorm(n_genes * n_samples, sd = 0.22),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, base::paste0('sample_', base::seq_len(n_samples))))
  expression_mat <- expression_mat +
    gene_baseline +
    base::rep(sample_offset, each = n_genes) +
    background_loading %*% base::t(sample_background)

  score_result <- calc_bulk_module_score(
    x = expression_mat,
    genes_to_score = list(module = module_genes),
    nbin = 24,
    ctrl = 80,
    seed = 7,
    warn_missing = FALSE,
    verbose = TRUE)
  module_score <- score_result$scores[['module']]
  pooled_controls <- score_result$details$module$pooled_controls
  gsva_param <- GSVA::gsvaParam(
    exprData = expression_mat,
    geneSets = list(module = module_genes),
    kcdf = 'Gaussian',
    minSize = 5,
    maxSize = Inf,
    verbose = FALSE)
  gsva_score <- base::as.numeric(GSVA::gsva(gsva_param, verbose = FALSE)['module', ])
  ranked_expression <- singscore::rankGenes(expression_mat)
  singscore_score <- base::suppressWarnings(singscore::simpleScore(
    ranked_expression,
    upSet = module_genes,
    centerScore = TRUE,
    knownDirection = TRUE)$TotalScore)
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expression_mat,
    geneSets = list(module = module_genes),
    minSize = 5,
    maxSize = Inf,
    normalize = TRUE,
    verbose = FALSE)
  ssgsea_score <- base::as.numeric(GSVA::gsva(ssgsea_param, verbose = FALSE)['module', ])
  plage_param <- GSVA::plageParam(
    exprData = expression_mat,
    geneSets = list(module = module_genes),
    minSize = 5,
    maxSize = Inf,
    verbose = FALSE)
  raw_mean_score <- base::colMeans(expression_mat[module_genes, , drop = FALSE])
  orient_to_reference <- function(score, reference) {
    score_correlation <- stats::cor(score, reference, use = 'complete.obs')
    if (base::is.finite(score_correlation) && score_correlation < 0) {
      return(-score)
    }
    score
  }
  plage_score <- orient_to_reference(
    base::as.numeric(GSVA::gsva(plage_param, verbose = FALSE)['module', ]),
    raw_mean_score)

  score_data <- tibble::tibble(
    sample_id = base::colnames(expression_mat),
    condition = sample_group,
    raw_mean = raw_mean_score,
    gsva = gsva_score,
    singscore = singscore_score,
    ssgsea = ssgsea_score,
    plage = plage_score,
    matched_score = module_score) %>%
    tidyr::pivot_longer(
      cols = base::c('raw_mean', 'gsva', 'singscore', 'ssgsea', 'plage', 'matched_score'),
      names_to = 'method',
      values_to = 'score') %>%
    dplyr::mutate(
      method = base::factor(
        .data$method,
        levels = base::c('raw_mean', 'gsva', 'singscore', 'ssgsea', 'plage', 'matched_score'),
        labels = base::c(
          'Raw module mean',
          'GSVA score',
          'singscore',
          'ssGSEA',
          'PLAGE',
          'Expression-matched score')))

  gene_shift_data <- tibble::tibble(
    gene = gene_ids,
    mean_expression = gene_baseline,
    condition_delta = base::rowMeans(expression_mat[, sample_group == 'Condition B', drop = FALSE]) -
      base::rowMeans(expression_mat[, sample_group == 'Condition A', drop = FALSE]),
    category = 'Other genes') %>%
    dplyr::mutate(
      category = dplyr::case_when(
        .data$gene %in% pooled_controls ~ 'Matched controls',
        .data$gene %in% module_genes ~ 'Module genes',
        TRUE ~ .data$category),
      category = base::factor(.data$category, levels = base::c('Other genes', 'Matched controls', 'Module genes')))

  effect_data <- score_data %>%
    dplyr::group_by(.data$method, .data$condition) %>%
    dplyr::summarise(mean_score = base::mean(.data$score), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = 'condition', values_from = 'mean_score') %>%
    dplyr::mutate(delta = .data[['Condition B']] - .data[['Condition A']])

  list(
    score_data = score_data,
    gene_shift_data = gene_shift_data,
    effect_data = effect_data)
}

plot_module_score_simulation <- function() {
  sim <- simulate_module_score_example()
  base::set.seed(230)
  background_sample <- sim$gene_shift_data %>%
    dplyr::filter(.data$category == 'Other genes') %>%
    dplyr::slice_sample(n = 850)
  matched_control_genes <- sim$gene_shift_data %>%
    dplyr::filter(.data$category == 'Matched controls')
  module_genes <- sim$gene_shift_data %>%
    dplyr::filter(.data$category == 'Module genes')
  category_colors <- base::c(
    'Other genes' = '#787878',
    'Matched controls' = '#2F6DB3',
    'Module genes' = '#D95F24')
  condition_colors <- base::c(
    'Condition A' = '#1F5A83',
    'Condition B' = '#B54A1F')
  score_method_labels <- base::c(
    'Raw module mean' = 'Raw mean',
    'GSVA score' = 'GSVA',
    'singscore' = 'singscore',
    'ssGSEA' = 'ssGSEA',
    'PLAGE' = 'PLAGE',
    'Expression-matched score' = 'Module score')
  score_data <- sim$score_data %>%
    dplyr::mutate(
      method_label = base::factor(
        base::unname(score_method_labels[.data$method]),
        levels = base::unname(score_method_labels)),
      condition_label = base::factor(
        .data$condition,
        levels = base::names(condition_colors),
        labels = base::c('A', 'B')))
  zero_line_data <- tibble::tibble(
    method_label = base::factor(
      base::c('GSVA', 'PLAGE', 'Module score'),
      levels = base::unname(score_method_labels)),
    yintercept = 0)
  effect_plot_data <- sim$effect_data %>%
    dplyr::arrange(.data$method) %>%
    dplyr::mutate(
      method_index = dplyr::row_number(),
      method_label = base::unname(score_method_labels[base::as.character(.data$method)]))

  p_gene_shift <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = background_sample,
      ggplot2::aes(x = .data$mean_expression, y = .data$condition_delta),
      color = category_colors[['Other genes']],
      alpha = 1,
      size = gs(0.74)) +
    ggplot2::geom_point(
      data = matched_control_genes,
      ggplot2::aes(x = .data$mean_expression, y = .data$condition_delta, color = .data$category),
      alpha = 0.82,
      size = gs(0.74)) +
    ggplot2::geom_point(
      data = module_genes,
      ggplot2::aes(x = .data$mean_expression, y = .data$condition_delta, color = .data$category),
      alpha = 0.88,
      size = gs(0.80)) +
    ggplot2::geom_smooth(
      data = sim$gene_shift_data,
      ggplot2::aes(x = .data$mean_expression, y = .data$condition_delta),
      method = 'loess',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = gs(0.32),
      span = 0.62) +
    ggplot2::scale_color_manual(
      values = category_colors[base::c('Matched controls', 'Module genes')],
      name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(size = gs(1.35), alpha = 1),
        keywidth = grid::unit(0.18, 'in'))) +
    ggplot2::labs(
      title = 'Module genes follow an expression-matched background shift',
      tag = 'a',
      x = bquote('Baseline mean log'[2] * '(CPM)'),
      y = 'Condition B - A') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = fs(7)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      legend.position = base::c(0.02, 0.98),
      legend.justification = base::c(0, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
      legend.text = ggplot2::element_text(size = fs(6.5)),
      legend.key.width = grid::unit(0.18, 'in'),
      legend.spacing.x = grid::unit(0.02, 'in'),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(8, 7, 8, 8))

  p_scores <- ggplot2::ggplot(score_data, ggplot2::aes(x = .data$condition_label, y = .data$score, fill = .data$condition)) +
    ggplot2::geom_hline(
      data = zero_line_data,
      ggplot2::aes(yintercept = .data$yintercept),
      inherit.aes = FALSE,
      color = 'black',
      linewidth = gs(0.22),
      linetype = 'dashed') +
    ggplot2::geom_boxplot(
      width = 0.38,
      outlier.shape = NA,
      linewidth = gs(0.30),
      alpha = 0.70,
      staplewidth = 0.3) +
    ggplot2::facet_wrap(~method_label, scales = 'free_y', nrow = 1) +
    ggplot2::scale_fill_manual(values = condition_colors, guide = 'none') +
    ggplot2::scale_color_manual(values = condition_colors, guide = 'none') +
    ggplot2::labs(
      title = NULL,
      tag = 'b',
      x = NULL,
      y = 'Score') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = fs(7), face = 'bold'),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      plot.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      panel.spacing.x = grid::unit(0.18, 'lines'),
      plot.margin = ggplot2::margin(8, 6, 8, 8))

  p_effect <- ggplot2::ggplot(effect_plot_data, ggplot2::aes(x = .data$method_index, y = .data$delta, fill = .data$method)) +
    ggplot2::geom_hline(yintercept = 0, color = 'black', linewidth = gs(0.22)) +
    ggplot2::geom_col(width = 0.46, color = NA, linewidth = 0, alpha = 0.82) +
    ggplot2::scale_x_continuous(
      breaks = effect_plot_data$method_index,
      labels = effect_plot_data$method_label,
      limits = base::c(0.58, base::max(effect_plot_data$method_index) + 0.42),
      expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_fill_manual(
      values = base::c(
        'Raw module mean' = '#D95F24',
        'GSVA score' = '#7C4D9E',
        'singscore' = '#4B8E6A',
        'ssGSEA' = '#20A5A6',
        'PLAGE' = '#C49A00',
        'Expression-matched score' = '#2F6DB3'),
      guide = 'none') +
    ggplot2::labs(
      title = 'Simulated condition effect per single-sample scoring method',
      tag = 'c',
      x = NULL,
      y = 'Mean difference (B - A)') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = fs(7)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text.x = ggplot2::element_text(size = fs(6), angle = 28, hjust = 1),
      axis.text.y = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(8, 8, 8, 8))

  p_gene_shift / p_scores / p_effect +
    patchwork::plot_layout(heights = base::c(1.04, 0.95, 0.76))
}

style_permutation_plot <- function(plot, tag = NULL) {
  plot +
    ggplot2::labs(tag = tag) +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      plot.title = ggplot2::element_text(face = 'plain', size = fs(7)),
      plot.subtitle = ggplot2::element_text(face = 'plain', size = fs(6)),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(8, 6, 14, 6))
}

make_permutation_summary_table <- function(permutation_results, module_labels) {
  dplyr::bind_rows(base::lapply(base::names(permutation_results), function(module_id) {
    result <- permutation_results[[module_id]]
    module_label <- module_labels[[module_id]]

    dplyr::bind_rows(
      result$observed_summary %>%
        dplyr::mutate(record_type = 'observed_summary'),
      result$null_median_summary %>%
        dplyr::mutate(record_type = 'null_median_summary'),
      result$group_p_values %>%
        dplyr::mutate(record_type = 'group_p_values'),
      result$trajectory %>%
        dplyr::mutate(record_type = 'trajectory')) %>%
      dplyr::mutate(
        module_id = module_id,
        module = module_label,
        .before = 1)
  }))
}

save_figure <- function(path, plot, width, height) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = figure_dpi,
    bg = 'white',
    limitsize = FALSE)
}

inline_svg_font <- function(path, svg_font_url, font_file, font_name) {
  svg_text <- base::readLines(path, warn = FALSE)
  encoded_font <- base64enc::base64encode(font_file)
  font_data_uri <- base::paste0('data:font/otf;base64,', encoded_font)
  svg_text <- base::gsub(
    pattern = base::paste0('url\\("', svg_font_url, '"\\) format\\("opentype"\\)'),
    replacement = base::paste0('url("', font_data_uri, '") format("opentype")'),
    x = svg_text,
    fixed = FALSE)
  if (!base::any(base::grepl(font_data_uri, svg_text, fixed = TRUE))) {
    base::stop('failed to inline ', font_name, ' in ', path, call. = FALSE)
  }
  base::writeLines(svg_text, path, useBytes = TRUE)
}

save_svg <- function(path, plot, width, height) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = svglite::svglite,
    web_fonts = list(
      svglite::font_face(
        'Nimbus Sans',
        otf = '../fonts/NimbusSans-Regular.otf',
        weight = 400),
      svglite::font_face(
        'Nimbus Sans',
        otf = '../fonts/NimbusSans-Bold.otf',
        weight = 700),
      svglite::font_face(
        'Nimbus Sans',
        otf = '../fonts/NimbusSans-Italic.otf',
        style = 'italic',
        weight = 400),
      svglite::font_face(
        'Nimbus Sans',
        otf = '../fonts/NimbusSans-BoldItalic.otf',
        style = 'italic',
        weight = 700)),
    fix_text_size = FALSE,
    bg = 'white',
    limitsize = FALSE)

  inline_svg_font(
    path,
    '../fonts/NimbusSans-Regular.otf',
    'report/assets/fonts/NimbusSans-Regular.otf',
    'Nimbus Sans Regular')
  inline_svg_font(
    path,
    '../fonts/NimbusSans-Bold.otf',
    'report/assets/fonts/NimbusSans-Bold.otf',
    'Nimbus Sans Bold')
  inline_svg_font(
    path,
    '../fonts/NimbusSans-Italic.otf',
    'report/assets/fonts/NimbusSans-Italic.otf',
    'Nimbus Sans Italic')
  inline_svg_font(
    path,
    '../fonts/NimbusSans-BoldItalic.otf',
    'report/assets/fonts/NimbusSans-BoldItalic.otf',
    'Nimbus Sans Bold Italic')
}

save_report_figure <- function(path_stub, plot, width, height) {
  width <- fd(width)
  height <- fd(height)

  save_figure(
    path = base::paste0(path_stub, '.png'),
    plot = plot,
    width = width,
    height = height)
  save_svg(
    path = base::paste0(path_stub, '.svg'),
    plot = plot,
    width = width,
    height = height)
}

# 1.3 load canonical inputs -----------------

metadata <- readr::read_csv(
  'data/GSE122380_metadata.csv',
  show_col_types = FALSE) %>%
  dplyr::mutate(
    sample_id = base::as.character(.data$sample_id),
    day = base::as.integer(.data$day),
    cell_line = base::as.character(.data$cell_line))

vst_tbl <- readr::read_csv(
  'data/GSE122380_vst.csv',
  show_col_types = FALSE)

vst_mat <- vst_tbl %>%
  tibble::column_to_rownames('gene_symbol') %>%
  base::as.matrix()

cpm_tbl <- readr::read_csv(
  'data/GSE122380_cpm.csv',
  show_col_types = FALSE)

cpm_mat <- cpm_tbl %>%
  tibble::column_to_rownames('gene_symbol') %>%
  base::as.matrix()

metadata <- metadata %>%
  dplyr::filter(.data$sample_id %in% base::colnames(vst_mat)) %>%
  dplyr::arrange(base::match(.data$sample_id, base::colnames(vst_mat)))
vst_mat <- vst_mat[, metadata$sample_id, drop = FALSE]
cpm_mat <- cpm_mat[, metadata$sample_id, drop = FALSE]
base::stopifnot(base::identical(base::colnames(vst_mat), metadata$sample_id))
base::stopifnot(base::identical(base::colnames(cpm_mat), metadata$sample_id))

# 1.4 define parameters -----------------

score_ctrl = 100
score_nbin = 24
score_seed = 1
figure_dpi = 600
figure_scale = 1
figure_font_scale = 1
figure_geom_scale = 1
figure_family = 'Nimbus Sans'
panel_tag_family = 'Nimbus Sans'
fs <- function(size) size * figure_font_scale
fd <- function(size) size * figure_scale
gs <- function(size) size * figure_geom_scale
gfs <- function(size) fs(size) / ggplot2::.pt
custom_permutation_n = 500
pca_gene_fraction = 0.10
day_order = 1:15
reference_day_palette <- viridisLite::viridis(15)
reference_day_palette[[15]] <- '#D8B11E'
annotation_day_palette <- grDevices::colorRampPalette(reference_day_palette)(256)
correlation_palette <- base::c(
  '#093F60', '#176086', '#2C83AA', '#56A5B8', '#82B6BB',
  '#AECFC0', '#D5E3BB', '#F6E699', '#FAD171', '#F5B14A',
  '#EA832A', '#D95F24', '#C43C22', '#A92325', '#831026')

public_module_order <- base::c(
  'embryonic_stem_cell_signature',
  'paraxial_mesoderm',
  'kegg_cardiac_muscle_contraction')

custom_module_order <- base::c(
  'progenitor',
  'mesoderm',
  'cardiogenic',
  'cardiomyocyte')

custom_marker_genes <- list(
  progenitor = base::c('ESRG', 'NANOG', 'SOX2', 'TDGF1', 'POU5F1'),
  mesoderm = base::c('EOMES', 'MESP1', 'MESP2', 'MIXL1', 'TBXT', 'FOXF1', 'MSX1', 'MSX2', 'TWIST1', 'SNAI1', 'HAND1'),
  cardiogenic = base::c('TBX5', 'GATA4', 'NKX2-5', 'MEF2C'),
  cardiomyocyte = base::c('ACTN2', 'MYH6', 'MYH7', 'MYL7', 'PLN', 'RYR2', 'TNNC1', 'TNNI3', 'TNNT2', 'TTN'))

public_module_colors <- base::c(
  embryonic_stem_cell_signature = '#20A5A6',
  paraxial_mesoderm = '#D95F02',
  kegg_cardiac_muscle_contraction = '#2F6DB3')

custom_colors <- base::c(
  progenitor = '#20A5A6',
  mesoderm = '#D95F02',
  cardiogenic = '#D9B300',
  cardiomyocyte = '#2F6DB3')

score_gradients <- list(
  embryonic_stem_cell_signature = viridisLite::mako(256),
  paraxial_mesoderm = viridisLite::inferno(256),
  kegg_cardiac_muscle_contraction = viridisLite::cividis(256),
  progenitor = viridisLite::mako(256),
  mesoderm = viridisLite::inferno(256),
  cardiogenic = viridisLite::viridis(256),
  cardiomyocyte = viridisLite::cividis(256))

# 1.5 define local helper functions -----------------

get_msigdb_gene_set <- function(gene_set_name,
                                matrix_genes,
                                collection = 'C2',
                                subcollection = NULL) {
  msigdb_sets <- msigdbr::msigdbr(
    db_species = 'HS',
    species = 'Homo sapiens',
    collection = collection,
    subcollection = subcollection)

  genes <- msigdb_sets %>%
    dplyr::filter(.data$gs_name == .env$gene_set_name, !base::is.na(.data$gene_symbol)) %>%
    dplyr::distinct(.data$gene_symbol) %>%
    dplyr::pull(.data$gene_symbol) %>%
    base::intersect(matrix_genes)

  list(
    genes = genes,
    map = msigdb_sets %>% dplyr::filter(.data$gs_name == .env$gene_set_name))
}

# 1.6 create directories -----------------

base::dir.create('results', recursive = TRUE, showWarnings = FALSE)
base::dir.create('report/assets/figures', recursive = TRUE, showWarnings = FALSE)

# 2.0 build public gene sets -----------------

public_gene_set_sources <- list(
  embryonic_stem_cell_signature = get_msigdb_gene_set(
    'BENPORATH_ES_2',
    base::rownames(vst_mat),
    collection = 'C2',
    subcollection = 'CGP'),
  paraxial_mesoderm = get_msigdb_gene_set(
    'GOBP_PARAXIAL_MESODERM_FORMATION',
    base::rownames(vst_mat),
    collection = 'C5',
    subcollection = 'GO:BP'),
  kegg_cardiac_muscle_contraction = get_msigdb_gene_set(
    'KEGG_CARDIAC_MUSCLE_CONTRACTION',
    base::rownames(vst_mat),
    collection = 'C2',
    subcollection = 'CP:KEGG_LEGACY'))

public_genes_to_score <- base::lapply(public_gene_set_sources, `[[`, 'genes')
base::stopifnot(base::all(base::vapply(public_genes_to_score, base::length, integer(1)) > 0L))

public_module_labels <- base::c(
  embryonic_stem_cell_signature = 'Embryonic stem cell signature',
  paraxial_mesoderm = 'Paraxial mesoderm',
  kegg_cardiac_muscle_contraction = 'Cardiac muscle contraction')

public_source_ids <- base::c(
  embryonic_stem_cell_signature = 'Ben-Porath ES2',
  paraxial_mesoderm = 'GO:0048341',
  kegg_cardiac_muscle_contraction = 'KEGG')

public_module_specs <- make_module_specs(
  module_order = public_module_order,
  labels = public_module_labels,
  source_ids = public_source_ids,
  gene_sets = public_genes_to_score,
  colors = public_module_colors,
  gradients = score_gradients)

# 3.0 score public and custom modules -----------------

public_score_result <- calc_bulk_module_score(
  x = vst_mat,
  genes_to_score = public_genes_to_score[public_module_order],
  ctrl = score_ctrl,
  nbin = score_nbin,
  seed = score_seed,
  warn_missing = FALSE,
  verbose = TRUE)

public_score_df <- public_score_result$scores %>%
  tibble::rownames_to_column('sample_id') %>%
  dplyr::left_join(metadata[, base::c('sample_id', 'day', 'cell_line')], by = 'sample_id')

custom_genes_to_score <- base::lapply(custom_marker_genes, base::intersect, base::rownames(vst_mat))
custom_genes_to_score <- custom_genes_to_score[custom_module_order]
base::stopifnot(base::all(base::vapply(custom_genes_to_score, base::length, integer(1)) > 0L))

custom_labels <- base::c(
  progenitor = 'Progenitor',
  mesoderm = 'Mesoderm',
  cardiogenic = 'Cardiogenic',
  cardiomyocyte = 'Cardiomyocyte')

custom_module_specs <- make_module_specs(
  module_order = custom_module_order,
  labels = custom_labels,
  source_ids = stats::setNames(base::rep('custom marker set', base::length(custom_module_order)), custom_module_order),
  gene_sets = custom_genes_to_score,
  colors = custom_colors,
  gradients = score_gradients)

custom_score_result <- calc_bulk_module_score(
  x = vst_mat,
  genes_to_score = custom_genes_to_score,
  ctrl = score_ctrl,
  nbin = score_nbin,
  seed = score_seed,
  warn_missing = FALSE,
  verbose = TRUE)

custom_score_df <- custom_score_result$scores %>%
  tibble::rownames_to_column('sample_id') %>%
  dplyr::left_join(metadata[, base::c('sample_id', 'day', 'cell_line')], by = 'sample_id')

public_score_long <- score_to_long(score_df = public_score_df, module_specs = public_module_specs)
custom_score_long <- score_to_long(score_df = custom_score_df, module_specs = custom_module_specs)

# 4.0 build top-variable-gene pca coordinates and score plots -----------------

pca_gene_vars <- base::apply(vst_mat, 1, stats::var, na.rm = TRUE)
pca_gene_count <- base::ceiling(base::length(pca_gene_vars) * pca_gene_fraction)
pca_genes <- base::names(base::sort(pca_gene_vars, decreasing = TRUE))[base::seq_len(pca_gene_count)]
pca_input_mat <- vst_mat[pca_genes, , drop = FALSE]

pca <- stats::prcomp(base::t(pca_input_mat), center = TRUE, scale. = FALSE)
pca_df <- tibble::as_tibble(pca$x[, 1:2], rownames = 'sample_id') %>%
  dplyr::left_join(public_score_df, by = 'sample_id') %>%
  dplyr::left_join(
    custom_score_df %>% dplyr::select(dplyr::all_of(base::c('sample_id', custom_module_order))),
    by = 'sample_id')

pc1_flip <- base::ifelse(
  stats::median(pca_df$PC1[pca_df$day == base::max(pca_df$day)], na.rm = TRUE) <
    stats::median(pca_df$PC1[pca_df$day == base::min(pca_df$day)], na.rm = TRUE),
  -1,
  1)
pca_df <- pca_df %>%
  dplyr::mutate(PC1 = .data$PC1 * pc1_flip)

day_path <- pca_df %>%
  dplyr::group_by(.data$day) %>%
  dplyr::summarise(
    PC1 = stats::median(.data$PC1, na.rm = TRUE),
    PC2 = stats::median(.data$PC2, na.rm = TRUE),
    .groups = 'drop') %>%
  dplyr::arrange(.data$day)
var_explained <- (pca$sdev^2) / base::sum(pca$sdev^2)
x_label <- base::sprintf('PC1 (%.1f%%)', var_explained[[1]] * 100)
y_label <- base::sprintf('PC2 (%.1f%%)', var_explained[[2]] * 100)

p_reference_pca <- plot_reference_pca(
  pca_data = pca_df,
  day_path = day_path,
  x_label = x_label,
  y_label = y_label,
  pca_gene_fraction = pca_gene_fraction,
  pca_gene_count = pca_gene_count)

timepoint_correlation_mat <- make_timepoint_correlation_matrix(
  expression_mat = pca_input_mat,
  metadata = metadata)
p_reference_timepoint_correlation <- plot_timepoint_correlation_heatmap(
  correlation_mat = timepoint_correlation_mat)
p_reference_overview <- p_reference_pca + p_reference_timepoint_correlation +
  patchwork::plot_layout(ncol = 2, widths = base::c(1, 1.08)) +
  patchwork::plot_annotation(
    caption = base::paste0(
      '<b>a.</b> PCA of bulk RNA-seq samples using the top 10% most variable VST genes.',
      '\n',
      '<b>b.</b> Pearson correlation matrix for day-collapsed VST profiles using the same genes; row and column strips encode differentiation day.')) &
  ggplot2::theme(
    plot.caption = ggtext::element_markdown(
      family = figure_family,
      size = fs(7),
      color = 'grey20',
      hjust = 0,
      lineheight = 1.18,
      margin = ggplot2::margin(t = 8)))

p_public_summary <- plot_module_pairs(
  score_data = public_score_long,
  pca_data = pca_df,
  day_path = day_path,
  module_specs = public_module_specs,
  x_label = x_label,
  y_label = y_label)

p_custom_summary <- plot_module_pairs(
  score_data = custom_score_long,
  pca_data = pca_df,
  day_path = day_path,
  module_specs = custom_module_specs,
  x_label = x_label,
  y_label = y_label)

p_custom_heatmap <- plot_score_heatmap(
  score_data = custom_score_long,
  module_specs = custom_module_specs)

p_public_score_correlation <- plot_public_score_correlation(score_df = public_score_df)
p_expression_bin_visualization <- plot_expression_bin_visualization(
  cpm_mat = cpm_mat,
  nbin = score_nbin,
  selected_bin = 12)
p_module_score_simulation <- plot_module_score_simulation()

# 5.0 run custom marker permutation tests -----------------

custom_permutation_results <- base::lapply(custom_module_order, function(module_id) {
  perm_bulk_module_score(
    x = vst_mat,
    gene_list = custom_genes_to_score[[module_id]],
    metadata = metadata,
    sample_col = 'sample_id',
    group_col = 'day',
    module_name = custom_labels[[module_id]],
    n_perm = custom_permutation_n,
    null_method = 'matched_bins',
    random_genome = FALSE,
    ctrl = score_ctrl,
    nbin = score_nbin,
    seed = 100 + base::match(module_id, custom_module_order),
    summary = 'median',
    trajectory_stat = 'last_minus_first',
    alternative = 'two.sided',
    make_plot = TRUE,
    make_histogram = TRUE)
}) %>%
  stats::setNames(custom_module_order)

custom_permutation_plots <- base::unlist(base::lapply(base::seq_along(custom_module_order), function(module_idx) {
  module_id <- custom_module_order[[module_idx]]
  result <- custom_permutation_results[[module_id]]
  list(
    style_permutation_plot(result$plot, tag = letters[[module_idx]]),
    style_permutation_plot(result$histogram))
}), recursive = FALSE)

p_custom_permutation <- patchwork::wrap_plots(custom_permutation_plots, ncol = 2) &
  ggplot2::theme(plot.margin = ggplot2::margin(10, 8, 26, 8))

custom_permutation_csv <- make_permutation_summary_table(
  permutation_results = custom_permutation_results,
  module_labels = custom_labels)

# 6.0 save outputs -----------------

readr::write_csv(
  public_score_df[, base::c('sample_id', 'day', 'cell_line', public_module_order)],
  'results/GSE122380_modulescore_public.csv')
readr::write_csv(
  custom_score_df[, base::c('sample_id', 'day', 'cell_line', custom_module_order)],
  'results/GSE122380_modulescore_custom.csv')
readr::write_csv(
  custom_permutation_csv,
  'results/GSE122380_modulescore_custom_marker_permutation.csv')

save_report_figure(
  path_stub = 'report/assets/figures/reference_pca_and_sample_correlation',
  plot = p_reference_overview,
  width = 7.2,
  height = 3.81)
save_report_figure(
  path_stub = 'report/assets/figures/public_gene_set_module_score_run',
  plot = p_public_summary,
  width = 7.2,
  height = 8.52)
save_report_figure(
  path_stub = 'report/assets/figures/custom_module_score_run',
  plot = p_custom_summary,
  width = 7.2,
  height = 11.52)
save_report_figure(
  path_stub = 'report/assets/figures/custom_marker_score_heatmap',
  plot = p_custom_heatmap,
  width = 6.48,
  height = 1.83)
save_report_figure(
  path_stub = 'report/assets/figures/public_score_relationship',
  plot = p_public_score_correlation,
  width = 7.2,
  height = 3.65)
save_report_figure(
  path_stub = 'report/assets/figures/expression_bin_visualization',
  plot = p_expression_bin_visualization,
  width = 7.2,
  height = 2.72)
save_report_figure(
  path_stub = 'report/assets/figures/module_score_simulation',
  plot = p_module_score_simulation,
  width = 7.2,
  height = 5.75)
save_report_figure(
  path_stub = 'report/assets/figures/custom_marker_permutation_run',
  plot = p_custom_permutation,
  width = 7.2,
  height = 9.12)

output_paths <- base::c(
  'results/GSE122380_modulescore_public.csv',
  'results/GSE122380_modulescore_custom.csv',
  'results/GSE122380_modulescore_custom_marker_permutation.csv',
  'report/assets/figures/reference_pca_and_sample_correlation.png',
  'report/assets/figures/reference_pca_and_sample_correlation.svg',
  'report/assets/figures/public_gene_set_module_score_run.png',
  'report/assets/figures/public_gene_set_module_score_run.svg',
  'report/assets/figures/custom_module_score_run.png',
  'report/assets/figures/custom_module_score_run.svg',
  'report/assets/figures/custom_marker_score_heatmap.png',
  'report/assets/figures/custom_marker_score_heatmap.svg',
  'report/assets/figures/custom_marker_permutation_run.png',
  'report/assets/figures/custom_marker_permutation_run.svg',
  'report/assets/figures/public_score_relationship.png',
  'report/assets/figures/public_score_relationship.svg',
  'report/assets/figures/expression_bin_visualization.png',
  'report/assets/figures/expression_bin_visualization.svg',
  'report/assets/figures/module_score_simulation.png',
  'report/assets/figures/module_score_simulation.svg')

output_paths

# 7.0 remove script scratch objects -----------------

base::rm(
  list = base::c(
    'annotation_day_palette',
    'correlation_palette',
    'public_module_colors',
    'public_genes_to_score',
    'public_module_labels',
    'public_module_order',
    'public_module_specs',
    'public_score_df',
    'public_score_long',
    'public_score_result',
    'public_source_ids',
    'custom_colors',
    'custom_genes_to_score',
    'custom_marker_genes',
    'custom_labels',
    'custom_module_order',
    'custom_module_specs',
    'custom_permutation_csv',
    'custom_permutation_n',
    'custom_permutation_plots',
    'custom_permutation_results',
    'custom_score_df',
    'custom_score_long',
    'custom_score_result',
    'day_order',
    'day_path',
    'figure_dpi',
    'figure_family',
    'format_module_score_labels',
    'gradient_end_mid_breaks',
    'gradient_end_limits',
    'get_day_colors',
    'get_msigdb_gene_set',
    'public_gene_set_sources',
    'make_module_specs',
    'make_permutation_summary_table',
    'make_timepoint_correlation_matrix',
    'metadata',
    'output_paths',
    'p_custom_heatmap',
    'p_custom_permutation',
    'p_custom_summary',
    'p_public_summary',
    'p_public_score_correlation',
    'pca_gene_count',
    'pca_gene_fraction',
    'pca_gene_vars',
    'pca_genes',
    'pca_input_mat',
    'p_reference_overview',
    'p_reference_pca',
    'p_reference_timepoint_correlation',
    'pc1_flip',
    'pca',
    'pca_df',
    'required_packages',
    'missing_packages',
    'plot_public_score_correlation',
    'plot_score_heatmap',
    'plot_module_boxplot',
    'plot_module_pairs',
    'plot_pca_score',
    'plot_reference_pca',
    'plot_score_arrow_relationship',
    'plot_timepoint_correlation_heatmap',
    'score_ctrl',
    'score_gradients',
    'score_nbin',
    'score_seed',
    'score_to_long',
    'save_figure',
    'save_report_figure',
    'save_svg',
    'reference_day_palette',
    'style_permutation_plot',
    'timepoint_correlation_mat',
    'var_explained',
    'vst_mat',
    'vst_tbl',
    'x_label',
    'y_label'),
  envir = base::.GlobalEnv)
