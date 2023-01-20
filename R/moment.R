#' Calculate Peptide Intensity Moments
#'
#' @param data A data frame containing peptide-level proteomics data.
#' @param names_col Column containing protein names.
#' @param exp_col Column containing experiment.
#' @param val_col Column containing peptide intensity data.
#' @param pos_col Column containing peptide midpoint.
#' @param len_col Column containing protein length.
#' @param sep Separating character. Must not be present in either protein names
#'     or experiment.
#'
#' @return A data frame containing peptide intensity moments for each protein
#'     and sample.
#' @export
#'
#' @examples
#'
moment <- function(data,
                   names_col,
                   exp_col,
                   val_col,
                   pos_col,
                   len_col,
                   sep = "/"){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == exp_col] <- "experiment"
  colnames(data)[colnames(data) == val_col] <- "value"

  # calculate position of peptide relative to centre of protein
  # calculate weighted intensity (intensity multipled by relative position)
  data[,"relative_position"] <- (data[,pos_col] - (data[,len_col]/2)) / data[,len_col]
  data[,"weight"] <- data[,"value"] * data[,"relative_position"]

  # unite name and experiment columns
  data <- tidyr::unite(data = data,
                       col = "name_experiment",
                       name,
                       experiment,
                       sep = sep)

  # keep only rows and columns of interest
  data_1 <- dplyr::filter(.data = data[,c("name_experiment",
                                          "value",
                                          "weight")],
                          value != 0)

  # separate into list, one item for each name/experiment combination
  data_2 <- list()
  for (i in unique(data_1[,"name_experiment"])){
    data_2[[i]] <- dplyr::filter(.data = data_1,
                                 name_experiment == i)
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
  # separate names and experiments
  output <- data.frame(name_experiment = gsub(".weight",
                                              "",
                                              names(moment)),
                       moment = moment,
                       nonzero = nonzero)
  row.names(output) <- NULL
  output <- tidyr::separate(data = moment,
                            col = name_experiment,
                            into = c("name", "experiment"),
                            sep = sep)

  # revert column names
  colnames(output)[colnames(output) == "name"] <- names_col
  colnames(output)[colnames(output) == "experiment"] <- exp_col
  colnames(output)[colnames(output) == "value"] <- val_col

  # return output data frame
  output
}
