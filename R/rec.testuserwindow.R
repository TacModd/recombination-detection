# implements (advanced) user set window method
# note does not work well in practise - perhaps consider
# changing implementation to testing genome length/n times
# and shift bounds from there?

rec.testuserwindow = function(ptns, sig, n){
  results = list()
  outertempcount = 0
  for (i in 1:length(ptns$pattern.IDs)){
    if (ptns$pattern.counts[i] > 2){
      indices = which(ptns$pattern.indices == i)
      p = length(indices)/length(ptns$pattern.indices)
      innertempcount = 0
      tempvector = list()
      j = indices[1]
      while (j < indices[length(indices)]){
        tempindices = which(ptns$pattern.indices[j:(j+n-1)] == i)
        tempindices = tempindices + j - 1
        q = length(tempindices) - 1
        r1 = 1 - pbinom(q, n, p)
        if (r1 <= sig & q > 0){
          k = 0
          q1 = length(tempindices) - k - 1
          while (q1 > 2){
            q1 = length(tempindices) - k - 1
            n1 = tempindices[length(tempindices) - k] - j
            r1 = 1 - pbinom(q1, n1, p)
            q2 = length(tempindices) - k - 2
            n2 = tempindices[length(tempindices) - k - 1] - j
            r2 = 1 - pbinom(q2, n2, p)

            if (r2 <= r1){
              k = k + 1
            } else {
              break
            }
          }
          q = length(tempindices) - k - 1
          n = tempindices[length(tempindices) - k] - j
          sigval = 1 - pbinom(q, n, p)
          tempvector[[j]] = c(i, j, tempindices[length(tempindices) - k], q+1, n, sigval)
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
