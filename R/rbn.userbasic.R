# implements (basic) user set window method

rbn.userbasic = function(partitions, sig, n, correction='neither'){
  
  ### initialise global variables
  # initialise results object
  results = list()
  # initialise a count value to keep track of recombined partitions
  outertempcount = 0
  # initialise a count value to keep track of all significant events
  m = 0
  
  ### for each unique partition ID:
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
      j = 1
      
      ### initialise variables for partition
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
        if (r > log(1 - sig)){ # (r1 <= sig & k > 0); k > 0/1 not necessary
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(i, j, tempindices[length(tempindices)], k, n, 1 - exp(r)) # log(1-r)
          # update left bound marker
          j = j + n
        # if not below significance threshold 
        } else {
          # just update left bound marker
          j = j + n
        }
      }
      # update all significant events count
      m = m + innertempcount
      
      # if at least 1 rbn event was found:
      if (length(tempvector) > 0){
        
        ### run (optional) local correction
        if (correction == 'local') {
          tempvector = local.bonferroni(tempvector, sig, innertempcount)
          # skip rest of loop if 0 results after correction (otherwise continue)
          if (length(tempvector) == 0){
            next
          }
        }
        
        ### store results
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
  
  ### run (optional) global correction
  if (correction == 'global') {
    # correct significance threshold
    results = global.bonferroni(results, sig, m)
  }
  
  ### return result object
  results
}
