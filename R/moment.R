moment <- function(data,
                   names_col,
                   sam_col,
                   val_col,
                   pos_col,
                   len_col,
                   sep = "/"){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == sam_col] <- "sample"
  colnames(data)[colnames(data) == val_col] <- "value"

  # calculate position of peptide relative to centre of protein
  # calculate weighted intensity (intensity multipled by relative position)
  data[,"relative_position"] <- (data[,pos_col] - (data[,len_col]/2)) / data[,len_col]
  data[,"weight"] <- data[,"value"] * data[,"relative_position"]

  # unite name and sample columns
  data <- tidyr::unite(data = data,
                       col = "name_sample",
                       name,
                       sample,
                       sep = sep)

  # keep only rows and columns of interest
  data_1 <- dplyr::filter(.data = data[,c("name_sample",
                                          "value",
                                          "weight")],
                          value != 0)

  # separate into list, one item for each name/sample combination
  data_2 <- list()
  for (i in unique(data_1[,"name_sample"])){
    data_2[[i]] <- dplyr::filter(.data = data_1,
                                 name_sample == i)
  }

  # calculate numbers of non-zero peptides
  # calculate raw and weighted intensity sums
  # calculate relative moments (weighted intensity sum / raw intensity sum)
  nonzero <- sapply(data_2,
                    function(X) nrow(X))
  sums <- lapply(data_2,
                 function(X) colSums(X[,c("intensity", "weight")],
                                     na.rm = TRUE))
  moment <- sapply(sums,
                   function(X) X["weight"] / X["intensity"])

  # write moments and number of peptides into data frame
  # separate names and samples
  output <- data.frame(name_sample = gsub(".weight",
                                          "",
                                          names(moment)),
                       moment = moment,
                       nonzero = nonzero)
  row.names(output) <- NULL
  output <- tidyr::separate(data = moment,
                            col = name_sample,
                            into = c("name", "sample"),
                            sep = sep)

  # revert column names
  colnames(output)[colnames(output) == "name"] <- names_col
  colnames(output)[colnames(output) == "sample"] <- sam_col
  colnames(output)[colnames(output) == "value"] <- val_col

  # return output data frame
  output
}
