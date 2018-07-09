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
      j = indices[1] # this introduces a bias on the first test - change?
      # while we haven't iterated past the last ptn index:
      while (j < indices[length(indices)]){
        # get indices between j and j + user window size
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        # correct index values (above values from 1)
        tempindices = tempindices + j - 1
        # get the number of 'events' (partitions) within window
        k = length(tempindices)
        # calculate the probability of k or more events
        r1 = 1 - pbinom(k-1, n, p) # fix
        # if probability below significance threshold:
        if (r1 <= sig & k > 0){ # fix?
          # no point in below line, should probably be removed
          #n = tempindices[length(tempindices)] - j
          # recalculate probability
          sigval = 1 - pbinom(k-1, n, p) # needs fixing to avoid underflow (convert to log and back)
          # store event details (needs fixing)
          tempvector[[j]] = c(i, j, tempindices[length(tempindices)], k, n, sigval)
          # update rbn event count
          innertempcount = innertempcount + 1
          # update left bound marker
          j = tempindices[q + 1] + 1 # not wrong but should be changed (simplified)
          
        # if not below significance threshold (should also be simplified probably)
        } else {
          # if there were actually partitions in window to test
          if (length(tempindices) > 0){
            # if marker equals the first ptn index
            if (j == tempindices[1]){
              # update it to the next
              j = j + 1
            # otherwise update it to the first
            } else {
              j = tempindices[1]
            }
          # otherwise add window size to marker
          } else {
            j = j + n
          }
        }
      }
      # initialise matrix to store results for ith partition
      tempmatrix = matrix(0, nrow=innertempcount, ncol=6)
      # reset rbn event count (line should go below check)
      innertempcount = 1
      # if at least 1 rbn event was found:
      if (length(tempvector) > 0){
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
