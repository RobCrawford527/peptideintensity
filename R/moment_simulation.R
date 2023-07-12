moment_simulation <- function(data,
                              number_of_peptides,
                              iterations = 100){

  # create empty data frame for results
  sim_result <- vector()

  # set up number of iterations
  # randomly sample and calculate simulated moments
  for (i in 1:iterations){

    # randomly sample from data
    data_sliced <- dplyr::slice_sample(data, n = number_of_peptides, replace = FALSE)

    # calculate moment based on random sample
    sim <- peptideintensity::moment(data = data_sliced,
                                    protein_id = GENENAME,
                                    experiment = experiment,
                                    value = intensity,
                                    position = Mid.position,
                                    length = protein_length,
                                    sep = "~",
                                    anchor = "centre",
                                    normalise = TRUE)

    # combine iteration into result data frame
    sim_result[i] <- mean(sim$moment)
  }

  # return overall simulated mean moment
  sim_result
}
