#' Plot The Midpoint Of A Cumulative Peptide Intensity Distribution
#'
#' @param input A data frame from plot_intensity_distribution, containing
#'     cumulative peptide intensity distributions for the protein, condition(s)
#'     and replicate(s) of interest.
#' @param plot A ggplot2 object created by plot_intensity_distribution.
#' @param shapes A vector of point shapes, defined in
#'     plot_intensity_distribution.
#'
#' @return A ggplot2 object containing points plotted where 50% of the
#'     cumulative intensity is reached for each distribution.
#'
#' @examples
#'
plot_midpoint <- function(input,
                          plot,
                          shapes){

  # define midpoint for each condition
  # midpoint = position in protein at which 50% of cumulative intensity is reached
  midpoint <- data.frame()
  for (e in unique(input[,"experiment"])){
    # select experiment
    input_e <- dplyr::filter(.data = input,
                             experiment == e)
    # determine absolute distance from 0.5 cumulative intensity
    # find minimum distance from 0.5
    input_e[,"midpoint"] <- abs(input_e[,"cumulative"] - 0.5)
    minimum <- min(input_e[,"midpoint"])

    # select row(s) corresponding to minimum distance from midpoint
    # if more than one, select only the first
    input_m <- dplyr::filter(.data = input_e,
                             midpoint == minimum)
    input_m <- dplyr::filter(.data = input_m,
                             distance == min(input_m[,"distance"]))

    # combine with other midpoints
    midpoint <- rbind.data.frame(midpoint,
                                 input_m)
  }

  # add elements to existing plot
  plot <- plot +
    # add scale for shape
    ggplot2::scale_shape_manual(values = shapes) +

    # plot midpoints
    ggplot2::geom_point(data = midpoint,
                        mapping = ggplot2::aes(fill = condition,
                                               shape = replicate,
                                               colour = NULL),
                        size = 1.6)

  # print midpoint data frame
  print(midpoint)

  # return updated plot
  plot
}
