# file adds helper functions for test correction

# adding familywise error correction (analogous to bonferroni)
##########   implementation strategy:   ##########
# 'global' correction
# after running all tests:
#   m = count significant tests
#   divide initial significance by m to get new threshold
#   eliminate results that don't pass new threshold
# (due to number of results this strategy may be too conservative)

#########   alternative specification:   #########
# 'local' correction
# after running tests for each partition:
#   correct significance by innertempcount
#   eliminate results before recording them
### may return empty list 
### ---> need to build test cases
### confirmed returns empty list if all results eliminated
### ---> add test cases to demonstrate
### confirmed returns subscript out of bounds if fed an empty list
### ---> best way to reduce unnecessary checks is to alter calling function

# call: tempvector = local.bonferroni(tempvector, sig, innertempcount)
# hypothetical (optional) local correction goes here
local.bonferroni = function(tempvector, sig, innertempcount){
  # for each result
  for (x in 1:length(tempvector)){
    # if result > corrected threshold
    if (tempvector[[x]][6] > sig/innertempcount){
      # remove result
      tempvector[[x]] = NA
    }
  }
  # clean empty results
  tempvector = tempvector[lapply(tempvector, length) > 1]
  tempvector
}

# call: results = global.bonferroni(results, sig, m)
# run (optional) correction
global.bonferroni = function(results, sig, m){
  # correct significance threshold
  sig = sig / m
  # for each partition
  for (x in 1:length(results)){
    # for each result
    for (y in 1:(length(results[[x]])/6)){
      # if result > new threshold
      if (results[[x]][y, 6] > sig){
        # remove result
        results [[x]][y, ] = NA
      }
    }
    # clean empty results
    results[[x]] = na.omit(results[[x]])
    # remove empty matrices
    if (all(is.na(results[[x]][, 1]))){
      results[[x]] = NA
    }
  }
  # clean empty partitions (matrices)
  results = results[lapply(results, length) > 1]
  results
}
