# takes a multiple alignment object, a partition object and a result object.
# masks sites in the multiple alignment object that are recombined 
# according to result object
# need to test that it's working on toy data
mask = function(sequences, partitions, results){
  # this line might be introducing a bug -- need to confirm
  # sequences = as.matrix(sequences)
  # for each recombined partition:
  for (i in 1:length(results)){
    # find the sites that accord with that partition
    indices = which(partitions$pattern.indices == results[[i]][1, 1])
    # then for each recombination event
    for (j in 1:nrow(results[[i]])){
      # get the bounds of the recombination event (as positions within indices)
      tempm = c(which(indices == results[[i]][j, 2]), which(indices == results[[i]][j, 3]))
      # and mask each partition site that falls within those bounds (inclusive)
      sequences[, indices[tempm[1]]:indices[tempm[length(tempm)]]] = as.DNAbin('n')
    }
  }
  # return the masked alignment object
  sequences
}
