# takes a multiple alignment (DNAbin) object, a partition object and a result object.
# masks ALL sites in the multiple alignment object that are within recombined bounds
# according to result object

# the prior mask.all function is severely inefficient. this function rewrites mask.all
# to remove needless overwriting of sites.

mask.all = function(sequences, partitions, results, aslist=T){
  # convert sequences to matrix in case sequences in list format
  sequences = as.matrix(sequences)
  
  # initialise vector to keep track of sites to mask
  mask.sites = 1:length(partitions$pattern.indices)
  
  # for each recombined partition:
  for (i in 1:length(results)){
    # find the sites that accord with that partition
    indices = which(partitions$pattern.indices == results[[i]][1, 1])
    # then for each recombination event
    for (j in 1:nrow(results[[i]])){
      ## get the bounds of the recombination event (as positions within indices)
      #tempm = c(which(indices == results[[i]][j, 2]), which(indices == results[[i]][j, 3]))
      
      # get the bounds of the recombination event
      site.bounds = c(results[[i]][j, 2], results[[i]][j, 3])
      
      ## mask everything within bounds
      #sequences[, indices[tempm[1]]:indices[tempm[length(tempm)]]] = as.DNAbin('n')
      
      # record everything within bounds for masking
      mask.sites[site.bounds[1]:site.bounds[2]] = 0
    }
  }
  # mask actual sites
  sequences[, which(mask.sites == 0)] = as.DNAbin('n')
  
  # convert back to list useful for writing as .phy file but the user may not want this
  if (aslist){
    sequences = as.list(sequences)
  }
  # return the masked alignment object
  sequences
}
