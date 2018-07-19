# implements (advanced) user set window method
# note does not appear to work well in practise 
# - requires further testing

rbn.usercomplex = function(partitions, sig, n){
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
      j = indices[1]
      # while we haven't iterated past the last ptn index:
      while (j < indices[length(indices)]){
        # get indices between j and j + user window size
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        # correct index values (above values start from 1 but should start from j)
        tempindices = tempindices + j - 1
        # get the number of 'events' (partitions) within window
        k = length(tempindices)
        # calculate the probability of k or more events
        r1 = pbinom(k-1, n, p, log=TRUE) # r1 = 1 - pbinom(1, n, p) #
        # if probability below significance threshold:
        if (r1 > log(1 - sig) & k > 1){ # (r1 <= sig & k > 1) #
          # set a counter
          q = 0
          # while there are at least 3 events
          while (k > 2){
            # update the number of events
            k = length(tempindices) - q
            # calculate probability of k-q events
            n1 = tempindices[length(tempindices) - q] - j
            r1 = pbinom(k-1, n1, p, log=TRUE) # r1 = 1 - pbinom(k-1, n1, p) #
            # calculate probability of k-q-1 events
            n2 = tempindices[length(tempindices) - q - 1] - j
            r2 = pbinom(k-2, n2, p, log=TRUE) # r2 = 1 - pbinom(k-2, n2, p) #
            # if shrinking the window reduces probability
            if (r2 > r1){ # (r2 <= r1) #
              # keep shrinking window (via right bound)
              q = q + 1
            # otherwise stop with current values
            } else {
              break
            }
          }
          
          # get number of events
          k = length(tempindices) - q
          # get size of window
          n = tempindices[length(tempindices) - q] - j
          # get log probability
          logsigval = pbinom(k-1, n, p, log=TRUE) # to reduce underflow
          # update result count
          innertempcount = innertempcount + 1
          # temporarily store results
          tempvector[[innertempcount]] = c(i, j, tempindices[length(tempindices) - q], k, n, 1 - exp(logsigval)) # log(1-logsigval)
          # update left bound marker to just after former right bound
          j = tempindices[k] + 1
        
        # if not below significance threshold
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
            j = j + n + 1
          }
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
          #if (!is.null(tempvector[[j]])){
            # add event details to matrix
          tempmatrix[innertempcount, ] = tempvector[[j]]
            # update rbn event count
          innertempcount = innertempcount + 1
          #}
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
