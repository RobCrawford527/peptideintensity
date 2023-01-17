peptide_intensity_distribution <- function(data,
                                           names_col,
                                           sam_col,
                                           seq_col,
                                           val_col,
                                           start_col,
                                           end_col,
                                           len_col,
                                           sep = "/"){
  
  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == sam_col] <- "sample"
  colnames(data)[colnames(data) == seq_col] <- "sequence"
  colnames(data)[colnames(data) == val_col] <- "value"
  colnames(data)[colnames(data) == start_col] <- "start"
  colnames(data)[colnames(data) == end_col] <- "end"
  colnames(data)[colnames(data) == len_col] <- "length"
  
  # unite name and sample columns
  data <- tidyr::unite(data = data,
                       col = "name_sample",
                       name,
                       sample,
                       sep = sep)
  
  # keep only rows and columns of interest
  data_1 <- dplyr::filter(.data = data[,c("name_sample",
                                          "sequence",
                                          "value",
                                          "start",
                                          "end",
                                          "length")],
                          value != 0)
  
  # separate into list, one item for each name/sample combination
  data_2 <- list()
  for (i in unique(data_1[,"name_sample"])){
    data_2[[i]] <- dplyr::filter(.data = data_1,
                                 name_sample == i)
  }
  
  ### intensity distributions
  # calculate cumulative peptide intensity distributions
  distribution <- lapply(data_2,
                         function(X) int_dist(X))
  
  ### metadata
  # find protein length
  # calculate numbers of non-zero peptides
  # calculate coverage (proportion of residues detected)
  length <- lapply(data_2,
                   function(X) X[1, "length"])
  nonzero <- lapply(data_2,
                    function(X) nrow(X))
  coverage <- lapply(distribution,
                     function(X) sum(X[,"intensity"] != 0))
  
  # create empty list for output
  # write in metadata and distribution for each name/sample combination
  output <- list()
  for (i in unique(names(data_2))){
    # split name and sample
    name_sample <- strsplit(i, split = "/")[[1]]
    
    # write output list
    # contains metadata vector and distribution data frame
    output[[i]] <- list(metadata = c(name = name_sample[1],
                                     sample = name_sample[2],
                                     length = length[[i]],
                                     nonzero_peptides = nonzero[[i]],
                                     coverage = coverage[[i]] / length[[i]]),
                        distribution = distribution[[i]])
  }
  
  # return output list
  output
}
