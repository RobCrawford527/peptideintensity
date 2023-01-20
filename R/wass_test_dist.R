#' Compare Peptide Intensity Distributions Using The Wasserstein Test
#'
#' @param reference A data frame containing the reference cumulative peptide
#'     intensity distribution.
#' @param sample A data frame containing the sample cumulative peptide
#'     intensity distribution to test.
#' @param nboots The number of bootstrap iterations to perform.
#' @param sep Separating character. Must not be present in experiment names.
#'
#' @return A data frame containing the output of the Wasserstein test for the
#'     comparison of interest.
#'
#' @examples
#'
wass_test_dist <- function(reference,
                           sample,
                           nboots,
                           sep){

  # create combined data frame
  combined <- data.frame(residue = reference[["distribution"]]$residue,
                         distance = reference[["distribution"]]$distance,
                         reference = reference[["distribution"]]$cumulative,
                         sample = sample[["distribution"]]$cumulative)

  # perform wasserstein test
  result <- twosamples::wass_test(a = combined[,"sample"],
                                  b = combined[,"reference"],
                                  nboots = nboots)
  names(result) <- NULL

  # determine dominant direction of shift
  # scale from -1 to +1
  # ... negative represents N-shift
  # ... 0 represents identical distributions
  # ... positive represents C-shift
  direction <- ifelse(sum(abs(combined[,"reference"] - combined[,"sample"])) == 0,
                      0,
                      sum(combined[,"reference"] - combined[,"sample"]) / sum(abs(combined[,"reference"] - combined[,"sample"])))

  # create output data frame
  output <- data.frame(protein = reference$name,
                       comparison = paste(sample$experiment,
                                          reference$experiment,
                                          sep = sep),
                       difference = result[1],
                       p_val = result[2],
                       direction = direction,
                       adj_difference = result[1] * direction)

  # return output data frame
  output
}
