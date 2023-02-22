#' Perform Wasserstein Test On Permutations Of A Peptide Intensity Distribution
#'
#' @param vec1 A peptide intensity distribution for the first permutation ("reference").
#' @param vec2 A peptide intensity distribution for the second permutation ("sample").
#'
#' @return The absolute difference between the two permutations.
#' @export
#'
#' @examples
#'
wass_stat_permutation <- function(vec1,
                                  vec2){

  # check that vectors are of equal length
  if (length(vec1) == length(vec2)){

    # define length
    len <- length(vec1)

    # order vectors
    vec1 <- vec1[order(vec1)]
    vec2 <- vec2[order(vec2)]

    # calculate difference metrics - only abs_diff for permutations
    # ... abs_diff: absolute difference between samples
    abs_diff <- sum(abs(vec1 - vec2)) / len

    # return abs_diff
    abs_diff
  }
}
