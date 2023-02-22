#' Statistical Testing Of Peptide Intensity Distributions
#'
#' @param data A multi-level list (protein, condition, replicate) containing
#'     cumulative peptide intensity distributions and their associated metadata.
#'     These may be either individual replicates (in which case the
#'     variation between replicates will be taken into account) but may also be
#'     means.
#' @param ref_condition The condition to use as reference.
#' @param nboots The number of bootstrap iterations to perform.
#' @param method P-value correction method to use in p.adjust. Should be one of
#'     p.adjust.methods (default "BH").
#'
#' @return A list of data frames containing the Wasserstein test results for
#'     each comparison.
#' @export
#'
#' @examples
#'
wasserstein_test_2 <- function(data,
                               ref_condition,
                               nboots = 5000,
                               method = "BH"){

  # create empty output list
  # select appropriate data and perform wasserstein test
  # write into data frame
  output <- list()
  for (p in names(data)){
    # select protein
    # define reference condition
    input <- data[[p]]
    reference <- input[[ref_condition]]

    # check that the reference condition is not empty
    if (!is.null(reference)){
      for (s in names(input)[names(input) != ref_condition]){
        # select sample condition
        sample <- input[[s]]

        # check that the sample condition is not empty
        if (!is.null(sample)){

          # perform wasserstein test and write into data frame
          # output is list of data frames, one for each sample condition
          output[[s]] <- rbind.data.frame(output[[s]],
                                          wass_test_2(reference = reference,
                                                      sample = sample,
                                                      nboots = nboots))
        }
      }
    }
  }

  # perform p-value correction using BH procedure
  # create empty output list
  # separate comparisons to correct each individually
  corrected <- lapply(output,
                      function(x) dplyr::mutate(x,
                                                adj_pval = stats::p.adjust(x$pval,
                                                                           method = method)))

  # return corrected output list
  corrected
}
