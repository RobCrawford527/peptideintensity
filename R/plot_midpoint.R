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

    # select row corresponding to minimum distance
    # combine with other midpoints
    midpoint <- rbind.data.frame(midpoint,
                                 dplyr::filter(.data = input_e,
                                               midpoint == minimum))
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
