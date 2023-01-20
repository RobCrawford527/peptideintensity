mean_peptide_intensity_distribution <- function(data){

  # calculate mean intensity distributions
  output <- lapply(data,
                   function(X) lapply(X,
                                      function(Y) peptideintensity:::mean_int_dist_1(Y)))

  # return output data frame
  output
}

example_dist_mean <- mean_peptide_intensity_distribution(example_dist)
