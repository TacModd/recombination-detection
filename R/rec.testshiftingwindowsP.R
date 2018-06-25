# just implements the same as rec.testshifingwindows but with a poisson instead of a binomial model

rec.testshiftingwindowsP = function(ptns, sig){
  results = list()
  outertempcount = 0
  for (i in 1:length(ptns$pattern.IDs)){
    if (ptns$pattern.counts[i] > 2){
      indices = which(ptns$pattern.indices == i)
      p = length(indices)/length(ptns$pattern.indices)
      innertempcount = 0
      tempvector = list()
      j = 1
      while (j < (length(indices)-1)){
        n1 = indices[j+2]-indices[j]+1
        r1 = 1 - pbinom(2, n1, p)
        if (r1 <= sig){
          k = j+2
          while (k < length(indices)){
            n1 = indices[k]-indices[j]
            r1 = 1 - ppois((k-j-1), p*n1)
            n2 = indices[k+1]-indices[j]
            r2 = 1 - ppois((k-j), p*n2)
            
            if (r2 <= r1){
              k = k + 1
            } else {
              break
            }
          }
          n = indices[k]-indices[j]
          sigval = 1 - ppois((k-j-1), p*n)
          tempvector[[j]] = c(indices[j], indices[k], k-j, n, sigval)
          innertempcount = innertempcount + 1
          j = j + k
        } else {
          j = j + 1
        }
      }
      tempmatrix = matrix(0, nrow=innertempcount, ncol=5)
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
