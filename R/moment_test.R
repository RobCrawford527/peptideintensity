#' Statistical testing of protein moments
#'
#' @param moment_data Data frame containing moments (for individual replicates)
#' @param peptide_data Data frame containing all peptides
#' @param proteins Vector of protein(s) to analyse (default is all proteins in data)
#' @param test_sample Sample to test
#' @param reference_sample Sample to use as reference
#' @param replicate_threshold Number of replicates required
#' @param peptide_threshold Number of peptides required
#' @param iterations Number of iterations of simulation to perform
#'
#' @return Data frame containing test results
#' @export
#'
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
