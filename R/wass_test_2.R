#' Compare Peptide Intensity Distributions Using The Wasserstein Test
#'
#' @param reference A list containing data for the reference condition.
#' @param sample A list containing data for the sample condition.
#' @param nboots The number of bootstrap iterations to perform.
#' @param protein The protein being tested.
#' @param ref_condition The reference condition.
#' @param sam_condition The sample condition being tested.
#'
#' @return A data frame containing the output of the Wasserstein test for the
#'     comparison of interest.
#' @export
#'
#' @examples
#'
wass_test_2 <- function(reference,
                        sample,
                        nboots = 5000,
                        protein,
                        ref_condition,
                        sam_condition){

  # create experiment list for reference and sample
  # create combined experiment list
  ref_exp <- unlist(lapply(reference, function(x) x[["experiment"]]))
  sam_exp <- unlist(lapply(sample, function(x) x[["experiment"]]))
  names(ref_exp) <- NULL
  names(sam_exp) <- NULL
  combined_exp <- c(sam_exp, ref_exp)


  ### calculate true test statistic ###
  # create reference and sample data frames with only cumulative peptide intensity distributions
  reference_df <- sapply(reference, function(x) x[["distribution"]]$cumulative)
  sample_df <- sapply(sample, function(x) x[["distribution"]]$cumulative)

  # determine mean in each position for reference and sample
  mean_ref <- rowMeans(reference_df)
  mean_sam <- rowMeans(sample_df)

  # calculate true test statistic using means
  true_test_stat <- wass_stat_extended(mean_ref,
                                       mean_sam)


  ### permutation test ###
  # create combined data frame
  # create combined vector and order it
  combined_df <- cbind(reference_df,
                       sample_df)
  len <- nrow(combined_df)
  combined <- vector()
  for (i in colnames(combined_df)){
    combined <- c(combined,
                  combined_df[,i])
  }
  combined <- combined[order(combined)]

  # set up for permutation test
  # ... n: number of iterations performed
  # ... larger: number of permutations with greater difference than true test statistic
  n <- 0
  larger <- 0

  # iterate through desired number of permutations
  # check if number of permutations performed is less than number desired
  while (n < nboots){

    # determine permutation of combined sample
    permutation <- combined[sample.int(n = length(combined),
                                       size = length(combined),
                                       replace = FALSE)]

    # create list of individual "samples" from permutation
    permutation_samples <- list()
    for (i in 1:length(combined_exp)){

      # select samples one-by-one
      j <- combined_exp[i]

      # select appropriate indices from vector
      # order new sample
      permutation_samples[[j]] <- permutation[(1 + len*(i-1)):(len + len*(i-1))]
      permutation_samples[[j]] <- permutation_samples[[j]][order(permutation_samples[[j]])]
    }

    # create reference and sample data frames for permutation
    reference_df_p <- as.data.frame(sapply(permutation_samples[ref_exp], function(x) x))
    sample_df_p <- as.data.frame(sapply(permutation_samples[sam_exp], function(x) x))

    # determine mean in each position for reference and sample for permutation
    mean_ref_p <- rowMeans(reference_df_p)
    mean_sam_p <- rowMeans(sample_df_p)

    # calculate test statistic (abs_diff) for permutation
    test_stat_p <- wass_stat_permutation(mean_ref_p,
                                         mean_sam_p)

    # compare test statistic for permutation with true test statistic
    # if equal or greater, increment larger by one
    if (test_stat_p >= true_test_stat["abs_diff"]){
      larger <- larger + 1
    }

    # increment n by one
    n <- n + 1
  }

  # calculate p-value from permutations
  # ... if no permutations equal or larger, assign minimum p-value
  # ... otherwise p-value = larger / nboots
  if (larger == 0){
    pval <- 1 / (nboots * 2)
  } else if (larger >= 1){
    pval <- larger / nboots
  }


  ### create output data frame ###
  # write data into one-line data frame
  output <- data.frame(protein = protein,
                       reference = ref_condition,
                       sample = sam_condition,
                       ref_reps = length(ref_exp),
                       sam_reps = length(sam_exp),
                       comparison = paste(sam_condition,
                                          ref_condition,
                                          sep = "_vs_"),
                       abs_diff = true_test_stat["abs_diff"],
                       net_diff = true_test_stat["net_diff"],
                       direction = true_test_stat["direction"],
                       signed_diff = true_test_stat["signed_diff"],
                       pval = pval)

  # return output data frame
  output
}
