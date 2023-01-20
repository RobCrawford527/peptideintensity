#' Statistical Testing Of Cumulative Peptide Intensity Distributions
#'
#' @param data A multi-level list (protein, condition, replicate) containing
#'     cumulative peptide intensity distributions and their associated metadata.
#'     These may be individual replicates or means.
#' @param ref_condition The condition to use as reference.
#' @param ref_replicate The replicate to use as reference (default = "mean").
#' @param nboots The number of bootstrap iterations to perform.
#' @param correction Logical indicating whether or not to correct the
#'     p-values (using the BH procedure).
#' @param alpha The significance level to use.
#' @param sep Separating character. Must not be present in experiment names.
#'
#' @return A list of data frames containing the Wasserstein test results for
#'     each comparison.
#' @export
#'
#' @examples
#'
wasserstein_test <- function(data,
                             ref_condition,
                             ref_replicate = "mean",
                             nboots = 5000,
                             correction = TRUE,
                             alpha = 0.01,
                             sep = "/"){

  # create empty output data frame
  # select appropriate data and perform wasserstein test
  # write into data frame
  output <- data.frame()
  for (i in names(data)){
    # select protein
    # define reference condition
    input <- data[[i]]
    reference <- input[[ref_condition]][[ref_replicate]]

    for (j in names(input)[names(input) != ref_condition]){
      # select sample condition
      input_j <- input[[j]]

      for (k in names(input_j)){
        # select sample replicate
        sample <- input_j[[k]]

        # perform wasserstein test
        # write into data frame
        output <- rbind.data.frame(output,
                                   peptideintensity:::wass_test_dist(reference = reference,
                                                                     sample = sample,
                                                                     nboots = nboots,
                                                                     sep = sep))
      }
    }
  }

  # perform p-value correction using BH procedure, if appropriate
  # create empty output list
  # separate comparisons to correct each individually
  corrected <- list()
  for (c in unique(output[,"comparison"])){
    # select comparison of interest
    output_c <- dplyr::filter(.data = output,
                              comparison == c)

    # perform BH procedure (or not, if correction = FALSE)
    # write output into list
    corrected[[c]] <- peptideintensity:::bh_correction(input = output_c,
                                                       correction = correction,
                                                       alpha = alpha)
  }

  # return output list
  corrected
}
