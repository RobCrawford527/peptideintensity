wasserstein_test_2 <- function(data,
                               ref_condition,
                               nboots = 5000,
                               method = "BH"){
  
  # create empty output list
  # select appropriate data and perform wasserstein test
  # write into data frame
  output <- list()
  for (p in names(data)){
    # select protein
    # define reference condition
    input <- data[[p]]
    reference <- input[[ref_condition]]
    
    # check that the reference condition is not empty
    if (!is.null(reference)){
      for (s in names(input)[names(input) != ref_condition]){
        # select sample condition
        sample <- input[[s]]
        
        # check that the sample condition is not empty
        if (!is.null(sample)){
          
          # perform wasserstein test and write into data frame
          # output is list of data frames, one for each sample condition
          output[[s]] <- rbind.data.frame(output[[s]],
                                          wass_test_2(reference = reference,
                                                      sample = sample,
                                                      nboots = nboots)
        }
      }
    }
  }
  
  # perform p-value correction using BH procedure
  # create empty output list
  # separate comparisons to correct each individually
  corrected <- lapply(output,
                      function(x) dplyr::mutate(x,
                                                adj_pval = stats::p.adjust(x$pval,
                                                                           method = method)))
  
  # return corrected output list
  corrected
}
