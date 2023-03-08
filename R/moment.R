#' Calculate Peptide Intensity Moments
#'
#' @param data A data frame containing peptide-level proteomics data.
#' @param protein_id Column containing protein names.
#' @param experiment Column containing experiment.
#' @param value Column containing peptide intensity data.
#' @param position Column containing peptide midpoint.
#' @param length Column containing protein length.
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
                   protein_id,
                   experiment,
                   value,
                   position,
                   length,
                   sep = "/",
                   anchor = "centre",
                   normalise = TRUE){

  # calculate position of peptide relative to the appropriate anchor point (N-terminus, C-terminus or centre)
  if (anchor == "centre"){
    data <- dplyr::mutate(data, relative_position = {{ position }} - ( {{ length }} / 2 ))
  } else if (anchor == "N"){
    data <- dplyr::mutate(data, relative_position = {{ position }})
  } else if (anchor == "C"){
    data <- dplyr::mutate(data, relative_position = {{ length }} - {{ position }})
  } else {
    stop("no anchor point defined: must be 'centre', 'N' or 'C'")
  }

  # normalise position if appropriate
  if (normalise == TRUE){
    data <- dplyr::mutate(data, relative_position = relative_position / {{ length }})
  }

  # calculate weighted intensity (intensity multipled by relative position)
  data <- dplyr::mutate(data, weight = {{ value }} * relative_position)

  # keep only rows of interest
  # group by name and experiment
  output <- dplyr::filter(data,
                          {{ value }} != 0)
  output <- dplyr::group_by(output,
                            {{ protein_id }}, {{ experiment }})

  # count number of non-zero peptides
  # calculate protein moment
  output <- dplyr::summarise(output,
                             n_peptides = dplyr::n(),
                             moment = sum(weight, na.rm = TRUE) / sum( {{ value }}, na.rm = TRUE))
  output <- dplyr::select(output,
                          {{ protein_id }}, {{ experiment }}, n_peptides, moment)
  output <- dplyr::distinct(output)

  # return output data frame
  output
}
