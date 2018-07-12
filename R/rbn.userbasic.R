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
      # initalise a vector to temporarily store event details
      tempvector = list()
      # initialise a left bound marker equal to the 1st partition index
      # j = indices[1] # tintroduces a bias on the first test
      j = 1
      # while we haven't iterated past the last ptn index:
      while (j < indices[length(indices)]){
        # get indices between j and j + user window size
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        # correct index values (above values from 1)
        tempindices = tempindices + j - 1
        # get the number of 'events' (partitions) within window
        k = length(tempindices)
        # calculate the probability of k or more events
        r = pbinom(k-1, n, p, log=TRUE) # 1 - pbinom(k-1, n, p)
        # if probability below significance threshold:
        if (r > log(1 - sig) & k > 1){ # (r1 <= sig & k > 0); k > 1 necessary?
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(i, j, tempindices[length(tempindices)], k, n, 1 - exp(r)) # log(1-r)
          # update left bound marker
          #j = tempindices[q + 1] + 1
          j = j + n
          
        # if not below significance threshold 
        } else {
          # add exception if !(k > 1)?
          # just update left bound marker
          j = j + n
        }
      }
      # initialise matrix to store results for ith partition
      tempmatrix = matrix(0, nrow=innertempcount, ncol=6)
      # if at least 1 rbn event was found:
      if (length(tempvector) > 0){
        # reset rbn event count
        innertempcount = 1
        # for each event:
        for (j in 1:length(tempvector)){
          # if the event is not null (is this check necessary??):
          if (!is.null(tempvector[[j]])){
            # add event details to matrix
            tempmatrix[innertempcount, ] = tempvector[[j]]
            # update rbn event count
            innertempcount = innertempcount + 1
          }
        }
        # update recombined ptns count
        outertempcount = outertempcount + 1
        # add matrix for ith ptn to result list
        results[[outertempcount]] = tempmatrix
      }
    }
  }
  # return result object
  results
}
