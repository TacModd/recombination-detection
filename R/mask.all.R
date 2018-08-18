# takes a multiple alignment (DNAbin) object, a partition object and a result object.
# masks ALL sites in the multiple alignment object that are within recombined bounds
# according to result object

mask.all = function(sequences, partitions, results, aslist=T){
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
      # mask everything within bounds
      sequences[, indices[tempm[1]]:indices[tempm[length(tempm)]]] = as.DNAbin('n')
    }
  }
  # convert back to list useful for writing as .phy file but the user may not want this
  if (aslist){
    sequences = as.list(sequences)
  }
  # return the masked alignment object
  sequences
}