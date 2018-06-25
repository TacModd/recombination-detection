# takes returned value of partitions function plus desired threshold
# iteratively calculates probability of k + 1 or more sampled in n
# for all possible k with first and last k indices forming bounds of n
# returns results as a list of matrices
res.testall = function(ptns, sig){
  results = list()
  outertempcount = 0
  for (i in 1:length(ptns$pattern.IDs)){
    if (ptns$pattern.counts[i] > 1){
      indices = which(ptns$pattern.indices == i)
      p = length(indices)/length(ptns$pattern.indices)
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
