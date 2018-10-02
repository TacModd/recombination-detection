### explores modifications to improve user basic function
# includes options to modify window bounds after discovery of event
# bounds left unmodified by default
# option to reduce left and right bound to outside partitions
# option to keep reducing bounds if it further reduces probability


rbn.usercomplex = function(partitions, sig, n, correction='neither', minsize=3, boundedit=0){
  
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
    if (partitions$pattern.counts[i] >= minsize){
      # get the partition indices
      indices = which(partitions$pattern.indices == i)
      # get the partition probability
      p = length(indices)/length(partitions$pattern.indices)
      # initialise a count value to keep track of detected rbn events
      innertempcount = 0
      # initalise a vector to temporarily store event details
      tempvector = list()
      # initialise a left bound marker equal to the 1st partition index
      # j = indices[1] # introduces a bias on the first test
      j = 1
      
      ### initialise variables for partition
      # while we haven't iterated past the SECOND last ptn index:
      while (j <= indices[length(indices) - 1]){
        # get indices between j and j + user window size
        tempindices = which(partitions$pattern.indices[j:(j+n-1)] == i)
        # correct index values (above values from 1)
        tempindices = tempindices + j - 1
        # get the number of 'events' (partitions) within window
        k = length(tempindices)
        # calculate the probability of k or more events
        r = pbinom(k-1, n, p, log=TRUE) # 1 - pbinom(k-1, n, p)
        
        # if probability not below significance threshold (or not enough events); k < 2 not necessary but useful
        if (r <= log(1 - sig) | k < 2){
          # just update left bound marker
          j = j + n
        }
        
        ### otherwise record event according to boundedit parameter # was if (r > log(1 - sig) & k > 1 & boundedit = 0){
        # if 0, hey! teacher! leave those bounds alone!
        else if (boundedit == 0){ # (r1 <= sig & k > 0)
          # update rbn event count
          innertempcount = innertempcount + 1
          # store event details
          tempvector[[innertempcount]] = c(i, j, j + n - 1, k, n, 1 - exp(r)) # log(1-r)
          # update left bound marker
          j = j + n
        }
        
        # if 1, edit events bounds to nearest partitions
        else if (boundedit == 1){
          # update rbn event count
          innertempcount = innertempcount + 1
          # update event size
          n = tempindices[length(tempindices)] - tempindices[1] + 1
          # update event probability
          r = pbinom(k-1, n, p, log=TRUE)
          # store event details, reducing bounds to nearest partitions
          tempvector[[innertempcount]] = c(i, tempindices[1], tempindices[length(tempindices)], k, n, 1 - exp(r)) # log(1-r)
          # update left bound marker appropriately 
          j = j + n
        }
        
        # if 2, iteratively edit bounds to reduce probability of event further
        else if (boundedit == 2){
          # set a counter
          y = 0
          # while there are at least 3 events, attempt to shrink right bound
          while (k > 2){
            # calculate probability of k-y events
            n1 = tempindices[length(tempindices) - y] - j + 1
            r1 = pbinom(k-1, n1, p, log=TRUE) # r1 = 1 - pbinom(k-1, n1, p) #
            # calculate probability of k-y-1 events
            n2 = tempindices[length(tempindices) - y - 1] - j + 1
            r2 = pbinom(k-2, n2, p, log=TRUE) # r2 = 1 - pbinom(k-2, n2, p) #
            # if shrinking the window reduces probability
            if (r2 > r1){ # (r2 <= r1) #
              # keep shrinking window
              y = y + 1
              # update the number of events
              k = length(tempindices) - y
            # otherwise stop with current values
            } else {
              break
            }
          }
          
          # set a counter (this needs to be 1 rather than 0)
          x = 1
          # while there are at least 3 events, attempt to shrink left bound
          while (k > 2){
            # calculate probability of k-y-x events
            n1 = tempindices[length(tempindices) - y] - tempindices[x] + 1
            r1 = pbinom(k-1, n1, p, log=TRUE) # r1 = 1 - pbinom(k-1, n1, p) #
            # calculate probability of k-y-x-1 events
            n2 = tempindices[length(tempindices) - y] - tempindices[x + 1] + 1
            r2 = pbinom(k-2, n2, p, log=TRUE) # r2 = 1 - pbinom(k-2, n2, p) #
            # if shrinking the window reduces probability
            if (r2 > r1){ # (r2 <= r1) #
              # keep shrinking window
              x = x + 1
              # update the number of events
              k = length(tempindices) - y - x
            # otherwise stop with current values
            } else {
              break
            }
          }
          
          # update rbn event count
          innertempcount = innertempcount + 1
          # update number of ptns
          k = length(tempindices) - y - x
          # update event size
          n = tempindices[length(tempindices) - y] - tempindices[x] + 1
          # update event probability
          r = pbinom(k-1, n, p, log=TRUE)
          # temporarily store results
          tempvector[[innertempcount]] = c(i, tempindices[x], tempindices[length(tempindices) - y], k, n, 1 - exp(r)) # log(1-r)
          # update left bound marker to just after former right bound
          j = tempindices[length(tempindices) - y] + 1
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
