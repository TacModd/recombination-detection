# takes a partition object and a desired significance threshold
# for each partition generates an initial window, then if window
# significant attempts to expand window to reduce likelihood
# further. otherwise starts new window with prior right bound
# as the new left bound (achieves linear complexity).
# assumes a poisson (not binomial) distribution.

rbn.shiftingrb.P = function(partitions, sig, correction='neither'){
  
  ### initialise global variables
  # initialise result object
  results = list()
  # initialise a count value to keep track of recombined partitions
  outertempcount = 0
  # initialise a count value to keep track of all significant events
  m = 0
  
  ### for each unique partition ID:
  for (i in 1:length(partitions$pattern.IDs)){
    # if there are at least 3 partitions belonging to said ID:
    if (partitions$pattern.counts[i] > 2){
      
      ### initialise variables for partition
      # get the partition indices
      indices = which(partitions$pattern.indices == i)
      # get the partition probability
      p = length(indices)/length(partitions$pattern.indices)
      # initialise a count value to keep track of detected rbn events
      innertempcount = 0
      # initialise a vector to temporarily store event details
      tempvector = list()
      # initialise a left bound marker
      j = 1
      
      ### identify recombination events
      # while we haven't iterated past the 2nd last ptn:
      while (j < (length(indices)-1)){
        # calculate window size to include 2 partitions
        n1 = indices[j+1]-indices[j]+0 # +1
        # calculate probability of 1 or more events # 3
        r1 = ppois(0, n1*p, log=TRUE) # r1 = 1 - ppois(0, n1*p) # 2
        if (r1 > log(1 - sig)){ # (r1 <= sig) #
          # set a right bound marker
          k = j+1 # +2
          # while we haven't iterated past the last ptn:
          while (k < length(indices)){
            # calculate the probability of k-j+0 events # +1
            n1 = indices[k]-indices[j]+0 # +1
            r1 = ppois((k-j-1), n1*p, log=TRUE) # r1 = 1 - ppois((k-j-1), n1*p) # (k-j)
            # calculate the probability of k-j+1 events # +2
            n2 = indices[k+1]-indices[j]+0 # +1
            r2 = ppois((k-j), n2*p, log=TRUE) # r2 = 1 - ppois((k-j), n2*p) # (k-j+1)
            # if expanding the window size reduces probability:
            if (r2 >= r1){ # (r2 <= r1) #
              # keep expanding window (via right bound)
              k = k + 1
            # otherwise stop and start a new window
            } else {
              break
            }
          }
          # get the size of the window
          n = indices[k]-indices[j]+0 # +1
          # get the (log) probability 
          logsigval = ppois((k-j-1), n*p, log=TRUE) # sigval = 1 - ppois((k-j), n*p)  # (k-j)
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(i, indices[j]+1, indices[k], k-j-1, n, 1 - exp(logsigval)) # indices[j]; k-j
          # again, specification can never detect 1st partition
          # update left bound marker to right bound marker # partition following
          j = j + k + 0 # + 1
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
