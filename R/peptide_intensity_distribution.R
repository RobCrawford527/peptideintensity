peptide_intensity_distribution <- function(data,
                                           names_col,
                                           exp_col,
                                           cond_col,
                                           rep_col,
                                           seq_col,
                                           val_col,
                                           start_col,
                                           end_col,
                                           len_col,
                                           sep = "/"){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == exp_col] <- "experiment"
  colnames(data)[colnames(data) == cond_col] <- "condition"
  colnames(data)[colnames(data) == rep_col] <- "replicate"
  colnames(data)[colnames(data) == seq_col] <- "sequence"
  colnames(data)[colnames(data) == val_col] <- "value"
  colnames(data)[colnames(data) == start_col] <- "start"
  colnames(data)[colnames(data) == end_col] <- "end"
  colnames(data)[colnames(data) == len_col] <- "length"

  # keep only rows and columns of interest
  data_1 <- dplyr::filter(.data = data[,c("name",
                                          "experiment",
                                          "condition",
                                          "replicate",
                                          "sequence",
                                          "value",
                                          "start",
                                          "end",
                                          "length")],
                          value != 0)

  # separate into list by protein name
  data_2 <- list()
  for (i in unique(data_1[,"name"])){
    data_1_filtered <- dplyr::filter(.data = data_1,
                                     name == i)
    for (j in unique(data_1_filtered[,"experiment"])){
      data_2[[i]][[j]] <- dplyr::filter(.data = data_1_filtered,
                                        experiment == j)
    }
  }

  # calculate cumulative peptide intensity distributions
  output <- lapply(data_2,
                   function(X) lapply(X,
                                      function(Y) int_dist_1(Y)))

  # return output list
  output
}

example_dist <- peptide_intensity_distribution(data = example,
                                               names_col = "GENENAME",
                                               exp_col = "experiment",
                                               cond_col = "condition",
                                               rep_col = "replicate",
                                               seq_col = "Sequence",
                                               val_col = "intensity",
                                               start_col = "Start.position",
                                               end_col = "End.position",
                                               len_col = "protein_length",
                                               sep = "/")
