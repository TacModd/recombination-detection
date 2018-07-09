# implements (basic) user set window method

rbn.userbasic = function(partitions, sig, n){
  # initialise results object
  results = list()
  # initialise a count value to keep track of recombined partitions
  outertempcount = 0
  # for each unique partition ID:
  for (i in 1:length(partitions$pattern.IDs)){
    # if there are at least 3 partitions belonging to said ID:
    if (partitions$pattern.counts[i] > 2){
      # get the partition indices
      indices = which(partitions$pattern.indices == i)
      # get the partition probability
      p = length(indices)/length(partitions$pattern.indices)
      # initialise a count value to keep track of detected rbn events
      innertempcount = 0
      # initalise a vector to store event details
      tempvector = list()
      # initialise a marker equal to the 1st partition index
      j = indices[1]
      # while we haven't iterated over all members of the current ptn:
      while (j < indices[length(indices)]){
        #
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        #
        tempindices = tempindices + j - 1
        #
        q = length(tempindices) - 1
        #
        r1 = 1 - pbinom(q, n, p) # fix
        #
        if (r1 <= sig & q > 0){ # fix?
          n = tempindices[length(tempindices)] - j
          sigval = 1 - pbinom(q, n, p) # doesn't need fixing
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
