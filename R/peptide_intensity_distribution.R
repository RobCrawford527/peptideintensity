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
    data_1a <- dplyr::filter(.data = data_1,
                             name == i)
    for (j in unique(data_1a[,"condition"])){
      data_1b <- dplyr::filter(.data = data_1a,
                               condition == j)
      for (k in unique(data_1b[,"replicate"])){
        data_2[[i]][[j]][[k]] <- dplyr::filter(.data = data_1b,
                                               replicate == k)
      }
    }
  }

  # calculate cumulative peptide intensity distributions
  output <- lapply(data_2,
                   function(X) lapply(X,
                                      function(Y) lapply(Y,
                                                         function(Z) peptideintensity:::int_dist_1(Z))))

  # return output list
  output
}
