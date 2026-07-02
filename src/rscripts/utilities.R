library(ggplot)

#' Aesthetic function to keep consistent plots appearance for MDS ggplots
#'
#' @param .x
#' @param group
#' @param shape_var
#' @param shapes
#' @param size_var
#'
#' @returns
#' @export
#' @import ggplot
#' @examples
mds_aes <- function(.x,
                    group = "group",
                    shape_var = "rep",
                    shapes = c(16, 17, 4, 15),
                    size_var = NULL) {
  if (!is.null(size_var)) {
    .x +
      geom_point(aes(
        color = .data[[group]],
        shape = .data[[shape_var]],
        size = .data[[size_var]]
      ), stroke = 1) +
      scale_shape_manual(values = shapes) +
      theme(legend.position = "bottom", legend.box = "horizontal") + scale_color_manual(values = order_colours) +
      guides(
        shape = guide_legend(
          title = shape_var,
          title.position = "top",
          nrow = 2
        ),
        color = guide_legend(
          title = "Group",
          title.position = "top",
          nrow = 2
        )
      )
  } else {
    .x +
      geom_point(aes(color = .data[[group]], shape = .data[[shape_var]]),
                 size = 4,
                 stroke = 1) +
      scale_shape_manual(values = shapes) +
      theme(legend.position = "bottom", legend.box = "horizontal") + scale_color_manual(values = order_colours) +
      guides(
        shape = guide_legend(
          title = "Replicate",
          title.position = "top",
          nrow = 2
        ),
        color = guide_legend(
          title = "Group",
          title.position = "top",
          nrow = 2
        ),
        size = guide_legend(
          title = "library size",
          title.position = "top",
          nrow = 2
        )
      )
  }
}







mds_aes <- function(.x,
                    shapes = c(16, 17, 4, 15),
                    colour_var = "cellType",
                    shape_var = "protocol",
                    size_var = NULL,
                    named_colours = named_cell_line_colours) {
  cell_line_colours <- RColorBrewer::brewer.pal(9, "Oranges")
  # hex <- hue_pal()(6)
  named_cell_line_colours <- list(
    "A549" = cell_line_colours[3],
    "H9" = cell_line_colours[4],
    "HEYA8" = cell_line_colours[5],
    "HepG2" = cell_line_colours[6],
    "Hct116" = cell_line_colours[7],
    "K562" = cell_line_colours[8],
    "MCF7" = cell_line_colours[9]
  )
  
  
  
  if (!is.null(size_var)) {
    .x +
      geom_point(aes(
        color = .data[[colour_var]],
        shape = .data[[shape_var]],
        size = .data[[size_var]]
      ), stroke = 1) +
      scale_shape_manual(values = shapes) +
      theme(legend.position = "bottom", legend.box = "horizontal") + scale_color_manual(values = named_colours) +
      guides(
        shape = guide_legend(title.position = "top", nrow = 1),
        color = guide_legend(
          title = "cell line",
          title.position = "top",
          nrow = 1
        )
      )
  } else {
    .x +
      geom_point(aes(color = .data[[colour_var]], shape = .data[[shape_var]]),
                 size = 4,
                 stroke = 1) +
      scale_shape_manual(values = shapes) +
      theme(legend.position = "bottom", legend.box = "horizontal") + scale_color_manual(values = named_colours) +
      guides(
        shape = guide_legend(title.position = "top", nrow = 1),
        color = guide_legend(
          title = "cell line",
          title.position = "top",
          nrow = 1
        )
      )
  }
}




# Custom function for scatter plots with correlation values
scatter_with_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # Calculate correlation
  corr <- cor.test(x, y, method = "pearson")
  est <- corr$estimate
  
  # Format correlation value
  cor_label <- sprintf("r = %.2f", est)
  
  # Create the scatter plot with correlation text
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5,
               size = 0.1,
               shape = 1) +
    geom_density_2d_filled(alpha = 0.8) +
    scale_fill_brewer(type = "seq",
                      palette = "Blues",
                      name = "density") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0))) +
    theme(legend.position = "bottom") +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "darkblue",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = Inf,
      # Position in the top right corner
      y = Inf,
      label = cor_label,
      hjust = 1.05,
      vjust = 1.05,
      size = 3,
      fontface = "bold",
      color = "darkblue" # ifelse(abs(est) > 0.7, "darkred",
      # ifelse(abs(est) > 0.4, "darkred", "darkred"))
    ) +
    coord_cartesian(expand = FALSE) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    ) +
    labs(x = c("1", "2", "3"), y = c("1", "2", "3"))
  
  return(p)
}


label_diagonal <- function(data, mapping, var_labels, ...) {
  # Create an empty plot
  p <- ggplot(data = data, mapping = mapping) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "gray95", color = "gray80"),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Extract the variable name from the mapping
  var_name <- as.character(mapping$x)[2]
  
  # Get the formatted label or use the original name if not found
  formatted_label <- var_labels[var_name]
  if (is.na(formatted_label))
    formatted_label <- var_name
  
  # Add the variable name text
  p <- p +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = formatted_label,
      size = 3.5,
      fontface = "bold"
    )
  
  return(p)
}

my_diag <- function(data, mapping, ...) {
  label_diagonal(data, mapping, var_labels, ...)
}
