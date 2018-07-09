# takes a partition object and a desired significance threshold
# for each partition generates all possible permutations of
# windows using the ith partition's sites as bounds, then
# calculates probability of k + 1 sites or more sampled within
# those bounds, assuming a binomial distribution
# returns a result object (a list of matrices)
rbn.permute = function(partitions, sig){
  results = list()
  outertempcount = 0
  for (i in 1:length(partitions$pattern.IDs)){
    if (partitions$pattern.counts[i] > 1){
      indices = which(partitions$pattern.indices == i)
      p = length(indices)/length(partitions$pattern.indices)
      innertempcount = 0
      tempvector = list()
      for (j in 1:(length(indices) - 1)){
        for (k in 1:(length(indices) - j)){
          n = indices[j+k] - indices[j] + 1
          r = 1 - pbinom(k, n, p)
          if (r <= sig){
            innertempcount = innertempcount + 1
            tempvector[[innertempcount]] = c(i, indices[j], indices[j + k], k + 1, n, r)
          }
        }
      }
      tempmatrix = matrix(0, nrow=innertempcount, ncol=6)
      innertempcount = 1
      if (length(tempvector) > 0){
        for (j in 1:length(tempvector)){
          if (!is.null(tempvector[[j]])){
            tempmatrix[innertempcount, ] = tempvector[[j]]
            innertempcount = innertempcount + 1
          }
        }
        outertempcount = outertempcount + 1
        results[[outertempcount]] = tempmatrix
      }
    }
  }
  results
}
