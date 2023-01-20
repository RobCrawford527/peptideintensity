int_dist_1 <- function(input){

  # check that at least one non-zero peptide is present
  if (nrow(input) > 0){

    # calculate cumulative peptide intensity distribution
    distribution <- peptideintensity:::int_dist_2(input)

    # write output list
    # find metadata and distribution
    output <- list(name = input[1, "name"],
                   experiment = input[1, "experiment"],
                   condition = input[1, "condition"],
                   replicate = input[1, "replicate"],
                   length = nrow(distribution),
                   nonzero_peptides = nrow(input),
                   coverage = sum(distribution[,"intensity"] != 0) / nrow(distribution),
                   distribution = distribution)

    # return output list
    output
  }
}
