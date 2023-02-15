#' Calculate Mean Peptide Intensity Distributions
#'
#' @param data A multi-level list (protein, condition, replicate) containing
#'     cumulative peptide intensity distributions and their associated metadata.
#'     For example, the output from peptide_intensity_distribution.
#'
#' @return A multi-level list in the same format as the input, containing mean
#'     cumulative peptide intensity distributions for each protein and
#'     condition along with metadata.
#' @export
#'
#' @examples
#'
mean_peptide_intensity_distribution <- function(data){

  # calculate mean intensity distributions
  output <- lapply(data,
                   function(X) lapply(X,
                                      function(Y) mean_int_dist_1(Y)))

  # discard null elements of list
  # remove conditions without mean
  # then remove proteins without any condition means
  output <- lapply(output,
                   function(X) purrr::discard(X,
                                              is.null))
  output <- purrr::discard(output,
                           function(X) length(X) == 0)

  # return output list
  output
}
