#' Statistical testing of moment for single protein
#'
#' @param moment_data Data frame containing moments (for individual replicates)
#' @param peptide_data Data frame containing all peptides
#' @param protein Single protein to test
#' @param test_sample Sample to test
#' @param reference_sample Sample to use as reference
#' @param replicate_threshold Number of replicates required
#' @param peptide_threshold Number of peptides required
#' @param iterations Number of iterations of simulation to perform
#'
#' @return Data frame with one row, containing test result
#'
moment_test_single <- function(moment_data,
                               peptide_data,
                               protein,
                               test_sample,
                               reference_sample,
                               replicate_threshold = 2,
                               peptide_threshold = 1,
                               iterations = 100){

  # filter moment data
  moment_data_test <- dplyr::filter(moment_data, GENENAME == protein & sample == test_sample)
  moment_data_reference <- dplyr::filter(moment_data, GENENAME == protein & sample == reference_sample)

  # check that at least one replicate is present for each sample
  if (nrow(moment_data_test) >= replicate_threshold & nrow(moment_data_reference) >= replicate_threshold){

    # filter peptide data
    # group by sample and experiment
    peptide_data_reference <- dplyr::filter(peptide_data, GENENAME == protein & sample == reference_sample)
    peptide_data_reference <- dplyr::group_by(peptide_data_reference,
                                              sample, experiment, replicate)

    # determine mean number of test peptides and min number of reference peptides
    # determine number of peptide to use in simulations
    # use reference_peptides - 1 to gain some variation info
    test_peptides <- round(mean(moment_data_test[["n_peptides"]]), 0)
    reference_peptides <- min(moment_data_reference[["n_peptides"]])
    number_of_peptides <- min(c(test_peptides, reference_peptides - 1))

    # check that number of peptides meets threshold
    if (number_of_peptides >= peptide_threshold){

      # perform simulations
      simulated_moments <- moment_simulation(data = peptide_data_reference,
                                             number_of_peptides = number_of_peptides,
                                             iterations = iterations)

      # create results data frame
      # determine mean moments and p-values
      result <- data.frame(GENENAME = protein,
                           test_sample = test_sample,
                           test_peptides = test_peptides,
                           test_mean = mean(moment_data_test[["moment"]]),
                           test_sd = sd(moment_data_test[["moment"]]),
                           reference_sample = reference_sample,
                           reference_peptides = round(mean(moment_data_reference[["n_peptides"]]), 0),
                           reference_mean = mean(moment_data_reference[["moment"]]),
                           reference_sd = sd(moment_data_reference[["moment"]]))
      result <- dplyr::mutate(result,
                              actual_diff = test_mean - reference_mean,
                              simulation_peptides = number_of_peptides,
                              simulation_mean = mean(simulated_moments),
                              simulation_sd = sd(simulated_moments),
                              simulation_diff = test_mean - simulation_mean,
                              pval = t.test(x = moment_data_test[["moment"]],
                                            y = simulated_moments,
                                            alternative = "two.sided")[["p.value"]])
    } else {
      # create empty result data frame
      result <- data.frame()
    }
  } else {
    # create empty result data frame
    result <- data.frame()
  }

  # return result data frame
  result
}
