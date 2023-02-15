#' Statistical Testing Of Cumulative Peptide Intensity Distributions
#'
#' @param data A multi-level list (protein, condition, replicate) containing
#'     cumulative peptide intensity distributions and their associated metadata.
#'     These may be individual replicates or means.
#' @param ref_condition The condition to use as reference.
#' @param ref_replicate The replicate to use as reference (default = "mean").
#' @param nboots The number of bootstrap iterations to perform.
#' @param method P-value correction method to use in p.adjust. Choose from
#'     "holm", "hochberg", "hommel", "bonferroni", "BH" (default), "BY", "fdr",
#'     or "none".
#' @param alpha The significance level to use.
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
                             method = "BH",
                             alpha = 0.01){

  # create empty output list
  # select appropriate data and perform wasserstein test
  # write into data frame
  output <- list()
  for (i in names(data)){
    # select protein
    # define reference condition
    input <- data[[i]]
    reference <- input[[ref_condition]][[ref_replicate]]

    # check that there is a reference distribution
    # if not, cannot perform test
    if (!is.null(reference)){
      for (j in names(input)[names(input) != ref_condition]){
        # select sample condition
        input_j <- input[[j]]

        for (k in names(input_j)){
          # select sample replicate
          sample <- input_j[[k]]

          # check that there is a sample distribution
          # if not, cannot perform test
          if (!is.null(sample)){
            # perform wasserstein test
            # write into data frame for comparison of interest
            output[[j]] <- rbind.data.frame(output[[j]],
                                            wass_test_dist(reference = reference,
                                                           sample = sample,
                                                           nboots = nboots,
                                                           replicate = k))
          }
        }
      }
    }
  }

  # perform p-value correction using BH procedure
  # create empty output list
  # separate comparisons to correct each individually
  corrected <- lapply(output,
                      function(X) dplyr::mutate(X,
                                                adj_pval = stats::p.adjust(X$pval,
                                                                           method = method)))

  # return output list
  #corrected
  output
}
