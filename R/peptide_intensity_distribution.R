#' Generate Cumulative Peptide Intensity Distributions
#'
#' @param data A data frame containing peptide-level proteomics data.
#' @param names_col Column containing protein names.
#' @param exp_col Column containing experiment.
#' @param cond_col Column containing condition.
#' @param rep_col Column containing replicate.
#' @param seq_col Column containing peptide sequence.
#' @param val_col Column containing peptide intensity data.
#' @param start_col Column containing start residue of peptide.
#' @param end_col Column containing end residue of peptide.
#' @param len_col Column containing protein length.
#' @param sep Separating character. Must not be present in either protein names
#'     or experiment.
#' @param threshold Threshold for number of peptides. Peptide intensity
#'     distribution is not calculated if below threshold. Default is 1.
#'
#' @return A multi-level list containing cumulative peptide intensity
#'     distributions (cumulative proportion of intensity observed against
#'     cumulative proportion of N-to-C distance). The list is separated
#'     first by protein, then by condition, then by replicate.
#' @export
#'
#' @examples
#'
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
                                           sep = "/",
                                           threshold = 1){

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

  # separate into list by protein name, condition and replicate
  data_2 <- list()
  for (i in unique(data_1[,"name"])){
    data_1a <- dplyr::filter(.data = data_1,
                             name == i)
    for (j in unique(data_1a[,"condition"])){
      data_1b <- dplyr::filter(.data = data_1a,
                               condition == j)
      for (k in unique(data_1b[,"replicate"])){
        data_1c <- dplyr::filter(.data = data_1b,
                                 replicate == k)

        # check if number of peptides exceeds threshold
        # write into list if it does
        if (nrow(data_1c) >= threshold){
          data_2[[i]][[j]][[k]] <- data_1c
        }
      }
    }
  }

  # calculate cumulative peptide intensity distributions
  output <- lapply(data_2,
                   function(X) lapply(X,
                                      function(Y) lapply(Y,
                                                         function(Z) int_dist_1(Z))))

  # return output list
  output
}
