# implements (basic) user set window method

rbn.userbasic = function(partitions, sig, n){
  results = list()
  outertempcount = 0
  for (i in 1:length(partitions$pattern.IDs)){
    if (partitions$pattern.counts[i] > 2){
      indices = which(partitions$pattern.indices == i)
      p = length(indices)/length(partitions$pattern.indices)
      innertempcount = 0
      tempvector = list()
      j = indices[1]
      while (j < indices[length(indices)]){
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        tempindices = tempindices + j - 1
        q = length(tempindices) - 1
        r1 = 1 - pbinom(q, n, p)
        if (r1 <= sig & q > 0){
          n = tempindices[length(tempindices)] - j
          sigval = 1 - pbinom(q, n, p)
          tempvector[[j]] = c(i, j, tempindices[length(tempindices)], q+1, n, sigval)
          innertempcount = innertempcount + 1
          j = tempindices[q + 1] + 1
        } else {
          if (length(tempindices) > 0){
            if (j == tempindices[1]){
              j = j + 1
            } else {
              j = tempindices[1]
            }
          } else {
            j = j + n + 1
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
