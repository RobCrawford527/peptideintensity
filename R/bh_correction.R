bh_correction <- function(input,
                          correction,
                          alpha){
  
  # order data frame
  # ... first by difference (largest to smallest)
  # ... then by p-value (smallest to largest)
  output <- input[order(input[,"difference"],
                        decreasing = TRUE),]
  output <- output[order(output[,"p_val"],
                         decreasing = FALSE),]  
  
  # rank proteins
  # assign critical value depending on correction term
  # assign significance
  output[,"rank"] <- 1:nrow(output)
  if (correction == TRUE){
    output[,"critical_val"] <- (output[,"rank"] / nrow(output)) * alpha
  } else {
    output[,"critical_val"] <- alpha
  }
  output[,"significant"] <- ifelse(output[,"p_val"] < output[,"critical_val"],
                                   TRUE,
                                   FALSE)
  
  # find lowest ranked significant protein 
  # adjust so that all higher ranked proteins are significant, regardless of critical value
  if (sum(output[,"significant"]) > 0){
    max <- max(dplyr::filter(.data = output,
                             significant == TRUE)[,"rank"])
    output[,"significant"] <- ifelse(output[,"rank"] <= max,
                                     TRUE,
                                     FALSE)
  }
  
  # return output data frame
  output
}
