mean_int_dist_2 <- function(input){
  
  # check there is data for at least two reps
  if (length(input) > 1){
    
    # create output data frame
    mean <- data.frame(residue = input[[1]][["distribution"]][,"residue"],
                       distance = input[[1]][["distribution"]][,"distance"])
    
    # write individual replicate distributions into combined data frame
    for (i in names(input)){
      mean[,i] <- input[[i]][["distribution"]][,"cumulative"]
    }
    
    # calculate row means and SDs
    rowmeans <- rowMeans(mean[,names(input)], na.rm = TRUE)
    rowsds <- matrixStats::rowSds(as.matrix(mean[,names(input)]), na.rm = TRUE)
    
    # write into data frame
    # calculate minimum and maximum for each point
    mean[,"cumulative"] <- rowmeans
    mean[,"min"] <- rowmeans - rowsds
    mean[,"max"] <- rowmeans + rowsds
    
    # keep only columns of interest
    mean <- mean[,c("residue",
                    "distance",
                    "cumulative",
                    "min",
                    "max")]
    
    # return mean distribution
    mean
  }
}
