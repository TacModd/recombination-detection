# takes a multiple alignment (DNAbin) object, a partition object and a result object.
# masks sites in the multiple alignment object that are recombined 
# according to result object

mask = function(sequences, partitions, results){
  # convert sequences to matrix in case sequences in list format
  sequences = as.matrix(sequences)
  # for each recombined partition:
  for (i in 1:length(results)){
    # find the sites that accord with that partition
    indices = which(partitions$pattern.indices == results[[i]][1, 1])
    # then for each recombination event
    for (j in 1:nrow(results[[i]])){
      # get the bounds of the recombination event (as positions within indices)
      tempm = c(which(indices == results[[i]][j, 2]), which(indices == results[[i]][j, 3]))
      # and mask each partition site that falls within those bounds (inclusive)
      sequences[, indices[tempm[1]:tempm[length(tempm)]]] = as.DNAbin('n')
    }
  }
  # could convert back to list but not necessary (better to leave up to user)
  # sequences = as.list(sequences)
  # return the masked alignment object
  sequences
}
