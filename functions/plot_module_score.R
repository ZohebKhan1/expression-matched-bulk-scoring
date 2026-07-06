#' Plot bulk module scores as boxplots
#'
#' Creates a publication-style boxplot for one bulk module score column across
#' a sample annotation column such as time point, genotype, treatment, or batch.
#'
#' @param score_df Data frame containing module scores and sample annotations.
#' @param score_col Character scalar naming the module score column.
#' @param x_col Character scalar naming the x-axis annotation column.
#' @param group_col Optional character scalar naming a grouping column used for
#'   grouped boxplots.
#' @param x_order Optional vector specifying the order of x-axis values.
#' @param group_order Optional vector specifying the order of group values.
#' @param fill Single fill color used when `group_col = NULL`.
#' @param palette Optional named vector of colors used when `group_col` is set.
#' @param title Optional plot title.
#' @param tag Optional panel tag.
#' @param x_label Optional x-axis title. Defaults to `x_col`.
#' @param y_label Y-axis title.
#' @param zero_line If `TRUE`, draw a dashed horizontal line at zero.
#' @param show_points If `TRUE`, overlay individual sample points.
#' @param connect_medians If `TRUE`, draw a smooth trend through group medians
#'   across x-axis values.
#' @param include_zero_in_limits If `TRUE`, include zero when calculating
#'   automatic y-axis limits.
#' @param y_limits Optional numeric length-two y-axis limits.
#' @param y_padding_mult Fraction of the observed y-range used as automatic
#'   y-axis padding.
#' @param y_padding_min Minimum automatic y-axis padding.
#' @param base_size Base ggplot font size.
#' @param base_family Base ggplot font family.
#'
#' @return A `ggplot` object.
plot_bulk_module_score_boxplot <- function(score_df,
                                           score_col,
                                           x_col,
                                           group_col = NULL,
                                           x_order = NULL,
                                           group_order = NULL,
                                           fill = '#D95F02',
                                           palette = NULL,
                                           title = NULL,
                                           tag = NULL,
                                           x_label = NULL,
                                           y_label = 'Bulk module score',
                                           zero_line = TRUE,
                                           show_points = TRUE,
                                           connect_medians = TRUE,
                                           include_zero_in_limits = FALSE,
                                           y_limits = NULL,
                                           y_padding_mult = 0.16,
                                           y_padding_min = 0.08,
                                           base_size = 12,
                                           base_family = 'sans') {
  .check_boxplot_args(
    score_df = score_df,
    score_col = score_col,
    x_col = x_col,
    group_col = group_col,
    y_limits = y_limits,
    y_padding_mult = y_padding_mult,
    y_padding_min = y_padding_min)

  plot_data <- .prepare_boxplot_data(
    score_df = score_df,
    score_col = score_col,
    x_col = x_col,
    group_col = group_col,
    x_order = x_order,
    group_order = group_order)
  has_group <- !is.null(group_col)
  y_limits <- .get_boxplot_y_limits(
    score_values = plot_data$.score,
    y_limits = y_limits,
    include_zero = include_zero_in_limits,
    y_padding_mult = y_padding_mult,
    y_padding_min = y_padding_min)
  x_label <- if (is.null(x_label)) x_col else x_label
  x_breaks <- seq_along(levels(plot_data$.x_plot))
  x_labels <- levels(plot_data$.x_plot)
  size_scale <- base_size / 12

  if (has_group) {
    p <- .plot_grouped_module_boxplot(
      plot_data = plot_data,
      group_col = group_col,
      palette = palette,
      zero_line = zero_line,
      show_points = show_points,
      connect_medians = connect_medians,
      size_scale = size_scale)
  } else {
    p <- .plot_single_module_boxplot(
      plot_data = plot_data,
      fill = fill,
      zero_line = zero_line,
      show_points = show_points,
      connect_medians = connect_medians,
      size_scale = size_scale)
  }

  p +
    ggplot2::labs(
      title = title,
      tag = tag,
      x = x_label,
      y = y_label) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'plain', size = 13, lineheight = 1.05),
      plot.tag = ggplot2::element_text(face = 'bold', size = 14),
      plot.tag.position = base::c(0.01, 0.99),
      axis.title = ggplot2::element_text(size = 13),
      axis.title.y = ggplot2::element_text(size = 13, margin = ggplot2::margin(r = 0)),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(12, 8, 32, 8))
}

.check_boxplot_args <- function(score_df,
                                score_col,
                                x_col,
                                group_col,
                                y_limits,
                                y_padding_mult,
                                y_padding_min) {
  if (!base::is.data.frame(score_df)) {
    base::stop('score_df must be a data frame.', call. = FALSE)
  }
  for (arg_name in base::c('score_col', 'x_col')) {
    arg_value <- base::get(arg_name)
    if (!base::is.character(arg_value) || base::length(arg_value) != 1L || base::is.na(arg_value)) {
      base::stop(arg_name, ' must be a non-missing character scalar.', call. = FALSE)
    }
  }
  if (!is.null(group_col) && (!base::is.character(group_col) || base::length(group_col) != 1L || base::is.na(group_col))) {
    base::stop('group_col must be NULL or a non-missing character scalar.', call. = FALSE)
  }
  required_cols <- base::c(score_col, x_col, group_col)
  missing_cols <- base::setdiff(required_cols, base::colnames(score_df))
  if (base::length(missing_cols) > 0L) {
    base::stop(
      'score_df is missing required column(s): ',
      base::paste(missing_cols, collapse = ', '),
      call. = FALSE)
  }
  if (!base::is.numeric(score_df[[score_col]])) {
    base::stop('score_col must name a numeric column.', call. = FALSE)
  }
  if (!is.null(y_limits) && (!base::is.numeric(y_limits) || base::length(y_limits) != 2L || base::any(!base::is.finite(y_limits)))) {
    base::stop('y_limits must be NULL or a finite numeric vector of length two.', call. = FALSE)
  }
  if (!base::is.numeric(y_padding_mult) || base::length(y_padding_mult) != 1L || !base::is.finite(y_padding_mult) || y_padding_mult < 0) {
    base::stop('y_padding_mult must be a non-negative numeric scalar.', call. = FALSE)
  }
  if (!base::is.numeric(y_padding_min) || base::length(y_padding_min) != 1L || !base::is.finite(y_padding_min) || y_padding_min < 0) {
    base::stop('y_padding_min must be a non-negative numeric scalar.', call. = FALSE)
  }
}

.prepare_boxplot_data <- function(score_df,
                                  score_col,
                                  x_col,
                                  group_col,
                                  x_order,
                                  group_order) {
  plot_data <- score_df
  plot_data$.score <- plot_data[[score_col]]
  plot_data$.x_value <- plot_data[[x_col]]
  keep <- base::is.finite(plot_data$.score) & !base::is.na(plot_data$.x_value)
  if (!is.null(group_col)) {
    plot_data$.group_value <- plot_data[[group_col]]
    keep <- keep & !base::is.na(plot_data$.group_value)
  }
  plot_data <- plot_data[keep, , drop = FALSE]
  if (base::nrow(plot_data) == 0L) {
    base::stop('No finite score values remain after filtering missing x/group values.', call. = FALSE)
  }

  if (is.null(x_order)) {
    x_order <- base::unique(plot_data$.x_value)
  }
  plot_data$.x_plot <- base::factor(plot_data$.x_value, levels = x_order)
  if (base::any(base::is.na(plot_data$.x_plot))) {
    base::stop('x_order does not include all observed x_col values.', call. = FALSE)
  }
  plot_data$.x_index <- base::as.integer(plot_data$.x_plot)

  if (!is.null(group_col)) {
    if (is.null(group_order)) {
      group_order <- base::unique(plot_data$.group_value)
    }
    plot_data$.group_plot <- base::factor(plot_data$.group_value, levels = group_order)
    if (base::any(base::is.na(plot_data$.group_plot))) {
      base::stop('group_order does not include all observed group_col values.', call. = FALSE)
    }
  }

  plot_data
}

.get_boxplot_y_limits <- function(score_values,
                                  y_limits,
                                  include_zero,
                                  y_padding_mult,
                                  y_padding_min) {
  if (!is.null(y_limits)) {
    return(y_limits)
  }

  range_values <- score_values
  if (include_zero) {
    range_values <- base::c(range_values, 0)
  }
  y_range <- base::range(range_values, na.rm = TRUE)
  y_span <- base::diff(y_range)
  y_padding <- base::max(y_span * y_padding_mult, y_padding_min)
  y_range + base::c(-y_padding, y_padding)
}

.plot_single_module_boxplot <- function(plot_data,
                                        fill,
                                        zero_line,
                                        show_points,
                                        connect_medians,
                                        size_scale) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$.x_index, y = .data$.score))
  if (zero_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = 'black', linewidth = 0.28 * size_scale, linetype = 'solid')
  }
  p <- p +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = .data$.x_plot),
      fill = fill,
      width = 0.62,
      linewidth = 0.35 * size_scale,
      outlier.shape = NA,
      alpha = 0.94,
      color = 'grey15')
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        position = ggplot2::position_identity(),
        shape = 21,
        size = 1.25 * size_scale,
        stroke = 0.2 * size_scale,
        alpha = 0.92,
        fill = fill,
        color = 'black')
  }
  if (connect_medians) {
    medians <- dplyr::summarise(
      dplyr::group_by(plot_data, .data$.x_plot),
      .x_index = dplyr::first(.data$.x_index),
      .median_score = stats::median(.data$.score, na.rm = TRUE),
      .groups = 'drop')
    p <- p +
      ggplot2::geom_line(
        data = medians,
        ggplot2::aes(x = .data$.x_index, y = .data$.median_score, group = 1),
        inherit.aes = FALSE,
        color = 'black',
        linewidth = 0.55 * size_scale)
  }
  p
}

.plot_grouped_module_boxplot <- function(plot_data,
                                         group_col,
                                         palette,
                                         zero_line,
                                         show_points,
                                         connect_medians,
                                         size_scale) {
  dodge <- ggplot2::position_dodge(width = 0.72)
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$.x_index, y = .data$.score, fill = .data$.group_plot, group = .data$.group_plot))
  if (zero_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = 'black', linewidth = 0.28 * size_scale, linetype = 'solid')
  }
  p <- p +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = base::interaction(.data$.x_plot, .data$.group_plot)),
      width = 0.62,
      linewidth = 0.35 * size_scale,
      outlier.shape = NA,
      alpha = 0.94,
      color = 'grey15',
      position = dodge)
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        position = dodge,
        shape = 21,
        size = 1.25 * size_scale,
        stroke = 0.2 * size_scale,
        alpha = 0.92,
        color = 'black')
  }
  if (connect_medians) {
    medians <- dplyr::summarise(
      dplyr::group_by(plot_data, .data$.x_plot, .data$.group_plot),
      .x_index = dplyr::first(.data$.x_index),
      .median_score = stats::median(.data$.score, na.rm = TRUE),
      .groups = 'drop')
    median_counts <- base::table(medians$.group_plot)
    p <- p +
      ggplot2::geom_line(
        data = medians,
        ggplot2::aes(
          x = .data$.x_index,
          y = .data$.median_score,
          group = .data$.group_plot,
          color = .data$.group_plot),
        inherit.aes = FALSE,
        linewidth = 0.55 * size_scale)
  }
  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_fill_manual(name = group_col, values = palette) +
      ggplot2::scale_color_manual(name = group_col, values = palette)
  } else {
    p <- p +
      ggplot2::labs(fill = group_col, color = group_col)
  }
  p
}
