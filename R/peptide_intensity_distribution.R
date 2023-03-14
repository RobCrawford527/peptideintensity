#' Generate Cumulative Peptide Intensity Distributions
#'
#' @param data A data frame containing peptide-level proteomics data.
#' @param protein_id Column containing protein names.
#' @param experiment Column containing experiment (i.e. sample + replicate).
#' @param sample Column containing sample.
#' @param replicate Column containing replicate.
#' @param sequence Column containing peptide sequence.
#' @param value Column containing peptide intensity data.
#' @param start Column containing start residue of peptide.
#' @param end Column containing end residue of peptide.
#' @param length Column containing protein length.
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
                                           protein_id,
                                           experiment,
                                           sample,
                                           replicate,
                                           sequence,
                                           value,
                                           start,
                                           end,
                                           length,
                                           sep = "/",
                                           threshold = 1){

  # keep only rows and columns of interest
  data_1 <- dplyr::filter(dplyr::select(data,
                                        {{ protein_id }},
                                        {{ experiment }},
                                        {{ sample }},
                                        {{ replicate }},
                                        {{ sequence }},
                                        {{ value }},
                                        {{ start }},
                                        {{ end }},
                                        {{ length }} ),
                          {{ value }} != 0)

  # separate into list by protein name, condition and replicate
  data_2 <- list()
  for (i in unique(dplyr::pull(data_1, {{ protein_id }} ))){
    data_1a <- dplyr::filter(data_1,
                             {{ protein_id }} == i)
    for (j in unique(dplyr::pull(data_1a, {{ sample }} ))){
      data_1b <- dplyr::filter(data_1a,
                               {{ sample }} == j)
      for (k in unique(dplyr::pull(data_1b, {{ replicate }} ))){
        data_1c <- dplyr::filter(data_1b,
                                 {{ replicate }} == k)

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
