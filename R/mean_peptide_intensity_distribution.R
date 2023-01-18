mean_peptide_intensity_distribution <- function(data,
                                                samples){


  for (i in names(data)){

    data[[i]]


  }


  temp <- example_dist[sapply(example_dist,
                              function(X) grepl("BY4741_Total_Untreated", X$condition) & grepl("RPS3", X$name))]
  temp_combined <- merge(temp[[1]]$distribution,
                         merge(temp[[2]]$distribution[,c(2,4)],
                               temp[[3]]$distribution[,c(2,4)],
                               by = "distance",
                               all = TRUE),
                         by = "distance",
                         all = TRUE)
  temp_combined$mean <- rowMeans(temp_combined[,4:6])
  temp_combined$min <- temp_combined$mean - rowSds(as.matrix(temp_combined[,4:6]))
  temp_combined$max <- temp_combined$mean + rowSds(as.matrix(temp_combined[,4:6]))





  combined <- as.matrix(data.frame(rep1 = temp[[1]]$distribution[,"cumulative"],
                                   rep2 = temp[[2]]$distribution[,"cumulative"],
                                   rep3 = temp[[3]]$distribution[,"cumulative"]))
  temp_combined <- data.frame(residue = temp[[1]]$distribution[,"residue"],
                              distance = temp[[1]]$distribution[,"distance"],
                              cumulative = rowMeans(combined, na.rm = TRUE),
                              min = rowMeans(combined, na.rm = TRUE) - rowSds(combined, na.rm = TRUE),
                              max = rowMeans(combined, na.rm = TRUE) + rowSds(combined, na.rm = TRUE))






  temp_combined <- temp_combined[,c("residue",
                                    "distance",
                                    "mean",
                                    "min",
                                    "max")]
  colnames(temp_combined)[colnames(temp_combined) == "mean"] <- "cumulative"



}
