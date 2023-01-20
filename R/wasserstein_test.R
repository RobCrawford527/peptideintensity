wasserstein_test <- function(data,
                             ref_condition,
                             ref_replicate = "mean",
                             nboots = 5000,
                             correction = TRUE,
                             alpha = 0.01){
  
  # create empty output data frame
  # select appropriate data and perform wasserstein test
  # write into data frame
  output <- data.frame()
  for (i in names(data)){
    # select protein
    # define reference condition
    input <- data[[i]]
    reference <- input[[ref_condition]][[ref_replicate]]
    
    for (j in names(input)[names(input) != ref_condition]){
      # select sample condition
      input_j <- input[[j]]
      
      for (k in names(input_j)){
        # select sample replicate
        sample <- input_j[[k]]
        
        # perform wasserstein test
        # write into data frame
        output <- rbind.data.frame(output,
                                   peptideintensity:::wass_test_dist(reference = reference,
                                                                     sample = sample,
                                                                     nboots = nboots))
      }
    }
  }
  
  # perform p-value correction using BH procedure, if appropriate
  # create empty output list
  # separate comparisons to correct each individually
  corrected <- list()
  for (c in unique(output[,"comparison"])){
    # select comparison of interest
    output_c <- dplyr::filter(.data = output,
                              comparison == c)
    
    # perform BH procedure (or not, if correction = FALSE)
    # write output into list
    corrected[[c]] <- peptideintensity:::bh_correction(input = output_c,
                                                       correction = correction,
                                                       alpha = alpha)
  }
  
  # return output list
  corrected
}
