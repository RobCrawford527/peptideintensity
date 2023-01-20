#' Plot Cumulative Peptide Intensity Distributions
#'
#' @param data A multi-level list containing cumulative peptide intensity
#'     distributions, either for individual replicates or means.
#' @param protein A single protein name to plot.
#' @param condition A vector of conditions to plot.
#' @param midpoint Logical indicating whether or not to plot the midpoints of
#'     the cumulative peptide intensity distributions.
#' @param ribbon Logical indicating whether to plot a ribbon showing mean
#'     plus/minus standard deviation. Only applicable if plotting mean
#'     distributions.
#'
#' @return A plot of cumulative peptide intensity against position in protein
#'     (from N terminus to C terminus) for the protein of interest. If
#'     midpoint = TRUE, midpoints are printed to the console.
#' @export
#'
#' @examples
#'
plot_intensity_distribution <- function(data,
                                        protein,
                                        condition,
                                        midpoint = FALSE,
                                        ribbon = FALSE){

  # select appropriate protein and conditions
  input <- data[[protein]][condition]

  # combine distributions into single data frame
  input_df <- data.frame()
  for (i in names(input)){
    # select condition
    input_i <- input[[i]]

    for (j in names(input_i)){
      # select replicate
      # select distribution data frame
      input_ij <- input_i[[j]][["distribution"]]

      # add experiment, condition and replicate information
      input_ij[,"experiment"] <- paste(i, j, sep = "_")
      input_ij[,"condition"] <- i
      input_ij[,"replicate"] <- j

      # combine with input_df
      input_df <- rbind.data.frame(input_df,
                                   input_ij)
    }
  }

  # define colour options
  # define breaks if breaks = NULL
  # define shape options
  colours <- viridis::viridis(n = length(unique(input_df[,"condition"])),
                              begin = 0.25,
                              end = 0.85,
                              option = "inferno")
  breaks <- condition
  shapes <- c(21:25)[1:length(unique(input_df[,"replicate"]))]

  # set up plot structure
  # plot distributions
  plot <- ggplot2::ggplot(data = input_df,
                          mapping = ggplot2::aes(x = distance,
                                                 y = cumulative,
                                                 group = experiment,
                                                 colour = condition)) +
    # set x and y axes
    ggplot2::coord_equal(ratio = 0.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1.02),
                                breaks = seq(0, 1, 0.25),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02),
                                breaks = seq(0, 1, 0.25),
                                expand = c(0, 0)) +

    # set colours
    ggplot2::scale_color_manual(values = colours,
                                breaks = breaks) +

    # set plot and axis titles
    ggplot2::ggtitle(label = protein) +
    ggplot2::xlab(label = "Position (N to C)") +
    ggplot2::ylab(label = "Cumulative intensity") +

    # set theme
    ggplot2::theme_classic() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +

    # plot horizontal marker lines
    ggplot2::geom_hline(yintercept = seq(0, 1, 0.25),
                        lty = 2,
                        colour = "grey90")

  # plot ribbons if appropriate
  if (ribbon == TRUE){
    plot <- plot +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = min,
                                                  ymax = max,
                                                  fill = condition,
                                                  colour = NULL),
                           alpha = 0.5)
  }

  # plot distributions
  plot <- plot +
    ggplot2::geom_line(linewidth = 0.8)

  # plot midpoints if appropriate
  if (midpoint == TRUE){
    plot <- peptideintensity:::plot_midpoint(input = input_df,
                                             plot = plot,
                                             shapes = shapes)
  }

  # add scale for fill if appropriate
  if (midpoint == TRUE | ribbon == TRUE){
    plot <- plot +
      # add scale for fill
      ggplot2::scale_fill_manual(values = colours,
                                 breaks = breaks)
  }

  # return plot
  plot
}
