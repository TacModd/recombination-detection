# alternative implementation of test 1

spt1 = function(ptns, sig){
  results = list()
  for (i in 1:length(ptns$pattern.IDs)){
    if (ptns$pattern.counts[i] > 2){
      indices = which(ptns$pattern.indices == i)
      p = length(indices)/length(ptns$pattern.indices)
      tempvector = list()
      j = 1
      while (j < (length(ptns)-1)){
        n1 = ptns[j+2]-ptns[j]+1
        r1 = 1 - pbinom(2, n1, p)
        if (r1 <= sig){
          k = j+2
          while (k < length(ptns)){
            n1 = ptns[k]-ptns[j]+1
            r1 = 1 - pbinom((k-j), n1, p)
            n2 = ptns[k+1]-ptns[j]+1
            r2 = 1 - pbinom((k-j+1), n2, p)
            
            if (r2 <= r1){
              k = k + 1
            } else {
              break
            }
          }
          n = ptns[k]-ptns[j]+1
          lp = p^(k-j+1) * (ptns[k]-ptns[j]+1)
          tempvector[[j]] = c(j, k, k-j+1, n, lp)
          j = j + k
        } else {
          j = j + 1
        }
      }
      results[[i]] = tempvector
    }
  }
  results
}
