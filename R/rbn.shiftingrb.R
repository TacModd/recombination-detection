# takes a partition object and a desired significance threshold
# for each partition generates an initial window, then if window
# significant attempts to expand window to reduce likelihood
# further. otherwise starts new window with prior right bound
# as the new left bound (achieves linear complexity).

rbn.shiftingrb = function(partitions, sig, correction='neither'){
  # initialise result object
  results = list()
  # initialise a count value to keep track of recombined partitions
  outertempcount = 0
  # initialise a count value to keep track of all significant events
  m = 0
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
      # initialise a left bound marker
      j = 1
      # while we haven't iterated past the 2nd last ptn:
      while (j < (length(indices)-1)){
        # calculate window size to include 3 partitions
        n1 = indices[j+2]-indices[j]+1
        # calculate probability of 3 or more events
        r1 = pbinom(2, n1, p, log=TRUE) # r1 = 1 - pbinom(2, n1, p) # 
        # probability below significant threshold:
        if (r1 > log(1 - sig)){ # (r1 <= sig) # 
          # set a right bound marker
          k = j+2
          # while we haven't iterated past the last ptn:
          while (k < length(indices)){
            # calculate the probability of k-j+1 or more events
            n1 = indices[k]-indices[j]+1
            r1 = pbinom((k-j), n1, p, log=TRUE) # r1 = 1 - pbinom((k-j), n1, p) # 
            # calculate the probability of k-j+2 or more events
            n2 = indices[k+1]-indices[j]+1
            r2 = pbinom((k-j+1), n2, p, log=TRUE) # r2 = 1 - pbinom((k-j+1), n2, p) # 
            # if expanding the window reduces probability
            if (r2 >= r1){ # (r2 <= r1)
              # keep expanding window (via right bound)
              k = k + 1
            # otherwise stop and start a new window
            } else {
              break
            }
          }
          # get size of window
          n = indices[k]-indices[j]+1
          # get log probability
          logsigval = pbinom((k-j), n, p, log=TRUE)
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(i, indices[j], indices[k], k-j+1, n, 1 - exp(logsigval))
          # update left bound marker to partition following right bound marker
          j = j + k
        # if starting window not significant:
        } else {
          # update left bound to next partition
          j = j + 1
        }
      }
      # update all significant events count
      m = m + innertempcount
      # if at least 1 rbn event was found:
      if (length(tempvector) > 0){
        # run (optional) local correction
        if (correction == 'local') {
          tempvector = local.bonferroni(tempvector, sig, innertempcount)
          # skip rest of loop if 0 results after correction (otherwise continue)
          if (length(tempvector) == 0){
            next
          }
        }
        # initialise matrix to store results for ith partition
        tempmatrix = matrix(0, nrow=innertempcount, ncol=6)
        # reset rbn event count
        innertempcount = 1
        # for each event:
        for (j in 1:length(tempvector)){
          # add event details to matrix
          tempmatrix[innertempcount, ] = tempvector[[j]]
          # update rbn event count
          innertempcount = innertempcount + 1
        }
        # update recombined ptns count
        outertempcount = outertempcount + 1
        # add matrix for ith ptn to result list
        results[[outertempcount]] = tempmatrix
      }
    }
  }
  # run (optional) global correction
  if (correction == 'global') {
    # correct significance threshold
    results = global.bonferroni(results, sig, m)
  }
  # return result object
  results
}
