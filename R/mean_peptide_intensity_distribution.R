#' Calculate Mean Peptide Intensity Distributions
#'
#' @param data The output from peptide_intensity_distribution. A multi-level
#'     list (protein, condition, replicate) containing cumulative peptide
#'     intensity distributions and their associated metadata.
#'
#' @return A multi-level list in the same format as the input, containing mean
#'     cumulative peptide intensity distributions for each protein and
#'     condition.
#' @export
#'
#' @examples
#'
mean_peptide_intensity_distribution <- function(data){

  # calculate mean intensity distributions
  output <- lapply(data,
                   function(X) lapply(X,
                                      function(Y) peptideintensity:::mean_int_dist_1(Y)))

  # return output data frame
  output
}
