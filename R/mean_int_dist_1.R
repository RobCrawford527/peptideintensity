mean_int_dist_1 <- function(input){

  # check there is data for at least two reps
  if (length(input) > 1){

    # calculate mean intensity distribution
    distribution <- peptideintensity:::mean_int_dist_2(input)

    # write output list
    # find metadata and distribution
    output <- list()
    output[["mean"]] <- list(name = input[[1]][["name"]],
                             experiment = paste(input[[1]][["condition"]],
                                                "mean",
                                                sep = "_"),
                             condition = input[[1]][["condition"]],
                             replicate = "mean",
                             length = input[[1]][["length"]],
                             nonzero_peptides = sapply(input,
                                                       function(X) X[["nonzero_peptides"]]),
                             coverage = sapply(input,
                                               function(X) X[["coverage"]]),
                             distribution = distribution)

    # return output list
    output
  }
}
