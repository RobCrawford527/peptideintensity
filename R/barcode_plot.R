#' Plot Protein Barcodes For Proteins Of Interest
#'
#' @param data Data frame containing peptide differential abundance results
#' @param protein_id Column containing protein IDs.
#' @param proteins Protein(s) to plot. No limit, but more than 5 or so will look messy.
#' @param comparison Column containing the comparison.
#' @param comparisons Comparison(s) to plot. All by default.
#' @param diff Column containing differential abundance data.
#' @param start Column containing peptide start positions.
#' @param end Column containing peptide end positions.
#' @param length Column containing protein lengths.
#' @param colour Column containing colour indications. Should have three levels: "down", "unchanged" and "up".
#'
#' @return Barcode plots for the selected proteins and comparisons, indicating peptides that have changed in abundance between samples.
#' @export
#'
#' @examples
#'
barcode_plot <- function(data,
                         protein_id,
                         proteins,
                         comparison,
                         comparisons = NULL,
                         diff,
                         start,
                         end,
                         length,
                         colour){

  # define colours to use in plot
  colours <- viridis::viridis(n = 2,
                              begin = 0.15,
                              end = 0.85,
                              direction = -1,
                              option = "inferno")
  colours <- c(colours[1], "grey75", colours[2])

  # filter data to keep only proteins and comparisons of interest
  # remove peptides with NA values for diff
  if (is.null(comparisons)){
    comparisons <- unique(data[[comparison]])
  }
  data <- dplyr::filter(data,
                        {{ protein_id }} %in% proteins & !is.na( {{ diff }} ) & {{ comparison }} %in% comparisons)

  # create plot
  plot <- ggplot2::ggplot(data,
                          mapping = ggplot2::aes(fill = {{ colour }} )) +
    ggplot2::geom_rect(mapping = ggplot2::aes(xmin = ( {{ start }} - 1) / {{ length }} * 100,
                                              xmax = {{ end }} / {{ length }} * 100,
                                              ymin = -15,
                                              ymax = 15),
                       alpha = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 100),
                                expand = c(0, 0),
                                breaks = NULL) +
    ggplot2::scale_y_continuous(limits = NULL,
                                expand = c(0, 0),
                                breaks = NULL) +
    ggplot2::coord_equal(ratio = 1) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::facet_grid(rows = dplyr::vars( {{ comparison }} ),
                        cols = dplyr::vars( {{ protein_id }} )) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA),
                   legend.position = "bottom")

  # return plot
  plot
}
