# takes a partition object and a desired significance threshold
# for each partition generates an initial window, then if window
# significant attempts to expand window to reduce likelihood
# further. otherwise starts new window with prior right bound
# as the new left bound (achieves linear complexity).
# assumes a poisson (not binomial) distribution.

rbn.shiftingrbP = function(partitions, sig){
  # initialise result object
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
      # initialise a vector to temporarily store event details
      tempvector = list()
      # initialise a left bound marker
      j = 1
      # while we haven't iterated past the 2nd last ptn:
      while (j < (length(indices)-1)){
        # calculate window size to include 3 partitions
        n1 = indices[j+2]-indices[j]+1
        # calculate probability of 3 or more events
        r1 = 1 - ppois(2, n1*p) # r1 = ppois(2, n1*p, log=TRUE) #
        if (r1 <= sig){ # (r1 > log(1 - sig)) #
          # set a right bound marker
          k = j+2
          # while we haven't iterated past the last ptn:
          while (k < length(indices)){
            # calculate the probability of k-j+1 events
            n1 = indices[k]-indices[j]
            r1 = 1 - ppois((k-j), n1*p) # r1 = ppois((k-j), n1*p, log=TRUE) #
            # calculate the probability of k-j+2 events
            n2 = indices[k+1]-indices[j]
            r2 = 1 - ppois((k-j+1), n2*p) # r1 = ppois((k-j+1), n1*p, log=TRUE) #
            # if expanding the window size reduces probability:
            if (r2 <= r1){ # (r2 >= r1) #
              # keep expanding window (via right bound)
              k = k + 1
            # otherwise stop and start a new window
            } else {
              break
            }
          }
          # get the size of the window
          n = indices[k]-indices[j]
          # get the (log) probability 
          sigval = 1 - ppois((k-j), n*p) # logsigval = ppois((k-j), n*p, log=TRUE)
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(indices[j], indices[k], k-j, n, sigval)
          # update left bound marker to partition following right bound marker
          j = j + k
        # if starting window not significant:
        } else {
          # update left bound to next partition
          j = j + 1
        }
      }
      # initialise matrix to store results for ith partition
      tempmatrix = matrix(0, nrow=innertempcount, ncol=5)
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
