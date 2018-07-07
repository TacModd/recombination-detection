### OUTDATED: no longer intended for use; here if I need the old output format ###
# takes returned value of partitions function plus desired threshold
# iteratively calculates probability of k + 1 or more sampled in n
# for all possible k with first and last k indices forming bounds of n
# returns results as a sparse list of lists
spt1 = function(ptns, sig){
  results = list()
  for (i in 1:length(ptns$pattern.IDs)){
    indices = which(ptns$pattern.indices == i)
    p = length(indices)/length(ptns$pattern.indices)
    if (ptns$pattern.counts[i] > 1){
      tempcount = 0
      tempvector = list()
      for (j in 1:(length(indices) - 1)){
        for (k in 1:(length(indices) - j)){
          n = indices[j+k] - indices[j] + 1
          r = 1 - pbinom(k, n, p)
          if (r <= sig){
            tempcount = tempcount + 1
            tempvector[[tempcount]] = c(indices[j], indices[j + k], k + 1, n, r)
          }
          
        }
      }
      results[[i]] = tempvector
    }
  }
  results
}
