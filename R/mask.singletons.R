# removes partition members and any singleton partitions between event bounds
# but not ALL partitions
# logic being that we don't have the power to detect singletons but otherwise
# we want to reduce the likelihood of false positives

mask.singletons = function(sequences, partitions, results, aslist=T){
  # convert sequences to matrix in case sequences in list format
  sequences = as.matrix(sequences)
  
  # get singleton partition indices
  singletons = match(which(partitions$pattern.counts == 1), partitions$pattern.indices)
  
  # initialise vector to keep track of sites to mask
  mask.sites = 1:length(partitions$pattern.indices)
  
  # for each recombined partition:
  for (i in 1:length(results)){
    # find the sites that accord with that partition
    indices = which(partitions$pattern.indices == results[[i]][1, 1])
    
    # add the singleton sites so that they can be detected if they occur within events
    indices = c(indices, singletons)
    
    # then for each recombination event
    for (j in 1:nrow(results[[i]])){
      ## get the bounds of the recombination event (as positions within indices)
      #tempm = c(which(indices == results[[i]][j, 2]), which(indices == results[[i]][j, 3]))
      
      # get the partition sites within the bounds of the recombination event
      sites = indices[indices >= results[[i]][j, 2] & indices <= results[[i]][j, 3]]
      
      ## mask each partition site
      #sequences[, indices[tempm[1]:tempm[length(tempm)]]] = as.DNAbin('n')
      
      # record each site for masking (inclusive)
      mask.sites[sites] = 0
    }
  }
  # mask the recorded sites
  sequences[, which(mask.sites == 0)] = as.DNAbin('n')
  
  # convert back to list useful for writing as .phy file but the user may not want this
  if (aslist){
    sequences = as.list(sequences)
  }
  # return the masked alignment object
  sequences
}
