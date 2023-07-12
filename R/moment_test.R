moment_test <- function(moment_data,
                        peptide_data,
                        proteins = NULL,
                        test_sample,
                        reference_sample,
                        replicate_threshold = 2,
                        peptide_threshold = 1,
                        iterations = 100){

  # create result data frame
  result <- data.frame()

  # determine list of proteins to test if not specified
  if (is.null(proteins)){
    proteins <- unique(moment_data[["GENENAME"]])
  }

  # test each protein
  for (p in proteins){
    result_single <- moment_test_single(moment_data = moment_data,
                                        peptide_data = peptide_data,
                                        protein = p,
                                        test_sample = test_sample,
                                        reference_sample = reference_sample,
                                        replicate_threshold = replicate_threshold,
                                        peptide_threshold = peptide_threshold,
                                        iterations = iterations)

    # combine protein result with overall result
    result <- rbind.data.frame(result, result_single)
  }

  # adjust p-values
  if (nrow(result) > 0){
    result <- dplyr::mutate(result,
                            adj_pval = p.adjust(pval, method = "BH"))
  }

  # return overall result
  result
}
