# takes a partition object and a desired significance threshold
# for each partition generates all possible permutations of
# windows using the ith partition's sites as bounds, then
# calculates probability of k + 1 sites or more sampled within
# those bounds, assuming a binomial distribution
# returns a result object (a list of matrices)

# adding familywise error correction (bonferroni)
##########   implementation strategy:   ##########
# 'global' correction
# after running all tests:
#   m = count significant tests
#   divide initial significance by m to get new threshold
#   how to eliminate results that don't pass new threshold?
# (due to number of results this strategy may be too conservative)

#########   alternative specification:   #########
# 'local' correction
# after running tests for each partition:
#   correct significance by innertempcount
#   eliminate results before recording them
# (not implemented)

rbn.permute.v2 = function(partitions, sig, correction='neither'){
  # initialise result object
  results = list()
  # initialise a count value to keep track of recombined partitions
  outertempcount = 0
  # initialise a count value to keep track of all significant events
  m = 0
  # for each unique partition ID:
  for (i in 1:length(partitions$pattern.IDs)){
    # if there are at least 2 partitions belonging to said ID:
    if (partitions$pattern.counts[i] > 1){
      # get the partition indices
      indices = which(partitions$pattern.indices == i)
      # get the partition probability
      p = length(indices)/length(partitions$pattern.indices)
      # initialise a count value to keep track of detected rbn events
      innertempcount = 0
      # initalise a vector to temporarily store event details
      tempvector = list()
      # for each possible ptn index as left bound:
      for (j in 1:(length(indices) - 1)){
        # for each possible ptn index as right bound:
        for (k in 1:(length(indices) - j)){
          # calculate n
          n = indices[j+k] - indices[j] + 1
          # calculate the probability of k+1 or more events
          r = pbinom(k, n, p, log=TRUE) # r = 1 - pbinom(k, n, p) #
          # if probability below significance threshold:
          if (r > log(1 - sig)){ # (r1 <= sig) #
            # update rbn event count
            innertempcount = innertempcount + 1
            # store event details
            tempvector[[innertempcount]] = c(i, indices[j], indices[j + k], k + 1, n, 1 - exp(r)) # log(1-r)
          }
        }
      }
      # initialise matrix to store results for ith partition
      tempmatrix = matrix(0, nrow=innertempcount, ncol=6)
      # update all significant events count
      m = m + innertempcount
      # if at least 1 rbn event was found:
      if (length(tempvector) > 0){
        # hypothetical (optional) local correction goes here
        if (correction == 'local') {
          tempvector = local.bonferroni(tempvector, sig, innertempcount)
        }
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
  # run (optional) correction
  if (correction == 'global') {
    # correct significance threshold
    results = global.bonferroni(results, sig, m)
  }
  # return result object
  results
}
