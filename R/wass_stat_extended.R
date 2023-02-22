wass_stat_extended <- function(vec1,
                               vec2){
  
  # check that vectors are of equal length
  if (length(vec1) == length(vec2)){
    
    # define length
    len <- length(vec1)
    
    # order vectors
    vec1 <- vec1[order(vec1)]
    vec2 <- vec2[order(vec2)]
    
    # calculate difference metrics
    # ... abs_diff: absolute difference between samples
    # ... net_diff: net difference between samples
    # ... direction: net_diff / abs_diff
    # ... signed_diff: abs_diff * sign of direction
    abs_diff <- sum(abs(vec1 - vec2)) / len
    net_diff <- sum(vec1 - vec2) / len
    direction <- ifelse(abs_diff == 0, 0, net_diff / abs_diff)
    signed_diff <- abs_diff * ifelse(direction >= 0, 1, -1)
    
    # create output vector containing distance metrics
    output <- c(abs_diff = abs_diff,
                net_diff = net_diff,
                direction = direction,
                signed_diff = signed_diff)
    
    # return output vector
    output
  }
}
