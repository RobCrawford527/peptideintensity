#' Calculate Peptide Intensity Moments
#'
#' @param data A data frame containing peptide-level proteomics data.
#' @param names_col Column containing protein names.
#' @param exp_col Column containing experiment.
#' @param val_col Column containing peptide intensity data.
#' @param pos_col Column containing peptide midpoint.
#' @param len_col Column containing protein length.
#' @param sep Separating character. Must not be present in either protein names
#'     or experiment. Default = "/".
#' @param anchor Must be set to "centre" (default), "N" or "C". Defines which
#'     position in the protein the moment is relative to.
#' @param normalise Logical indicating whether the position should be
#'     normalised to the protein length (default = TRUE).
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
                   sep = "/",
                   anchor = "centre",
                   normalise = TRUE){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == exp_col] <- "experiment"
  colnames(data)[colnames(data) == val_col] <- "value"


  # calculate position of peptide relative to the appropriate anchor point (N-terminus, C-terminus or centre)
  if (anchor == "centre"){
    data[,"position"] <- data[,pos_col] - (data[,len_col]/2)
  } else if (anchor == "N"){
    data[,"position"] <- data[,pos_col]
  } else if (anchor == "C"){
    data[,"position"] <- data[,len_col] - data[,pos_col]
  } else {
    stop("no anchor point defined")
  }

  # normalise position if appropriate
  if (normalise == TRUE){
    data[,"position"] <- data[,"position"] / data[,len_col]
  }

  # calculate weighted intensity (intensity multipled by relative position)
  data[,"weight"] <- data[,"value"] * data[,"position"]

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
                 function(X) colSums(X[,c("value", "weight")],
                                     na.rm = TRUE))
  moment <- sapply(sums,
                   function(X) X["weight"] / X["value"])
  names(moment) <- gsub(".weight",
                        "",
                        names(moment))

  # write moments and number of peptides into data frame
  # separate names and experiments
  output <- data.frame(name_experiment = names(moment),
                       moment = moment,
                       nonzero = nonzero)
  row.names(output) <- NULL
  output <- tidyr::separate(data = output,
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
