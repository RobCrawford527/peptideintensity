moment <- function(data,
                   names_col,
                   sam_col,
                   val_col,
                   pos_col,
                   len_col){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_col] <- "sam"

  # define proteins and samples
  proteins <- unique(data[,"names"])
  samples <- unique(data[,"sam"])

  # calculate relative position for each peptide
  # calculate weighted intensity for each peptide
  data[,"rel_pos"] <- (data[,pos_col] - (data[,len_col]/2)) / data[,len_col]
  data[,"weighted_intensity"] <- data[,"rel_pos"] * data[,val_col]

  # create output data frame
  output <- data.frame(names = NA,
                       sam = NA,
                       peptides = NA,
                       nonzero_peptides = NA,
                       n_peptides = NA,
                       c_peptides = NA,
                       moment = NA,
                       relative_moment = NA)[0,]

  # calculate moments and write into data frame
  for (i in proteins){
    # filter dataset
    data_i <- dplyr::filter(.data = data,
                            names == i)
    for (j in samples){
      # filter dataset
      data_ij <- dplyr::filter(.data = data_i,
                               sam == j)

      # calculate moment for protein/sample combination
      # remove combinations with no non-zero peptides
      moment <- rbind.data.frame(output,
                                 data.frame(names = i,
                                            sam = j,
                                            peptides = nrow(data_ij),
                                            nonzero_peptides = sum(data_ij[,val_col] > 0),
                                            n_peptides = sum(data_ij[,"rel_pos"] <= 0 & data_ij[,val_col] > 0),
                                            c_peptides = sum(data_ij[,"rel_pos"] > 0 & data_ij[,val_col] > 0),
                                            moment = sum(data_ij[,"weighted_intensity"]),
                                            relative_moment = sum(data_ij[,"weighted_intensity"]) / sum(data_ij[,val_col])))
      moment <- na.omit(moment)

      # write in to output data frame
      output <- rbind.data.frame(output,
                                 moment)
    }
  }

  # return output data frame
  output
}
