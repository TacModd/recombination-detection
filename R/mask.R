# masks partitions in DNAbin found to be recombined by detection method
# need to test that it's working on toy data
mask = function(results, sequences, partitions){
  sequences = as.matrix(sequences)
  for (i in 1:length(results)){
    indices = which(partitions$pattern.indices == results[[i]][1, 1])
    for (j in 1:nrow(results[[i]])){
      tempm = c(which(indices == results[[i]][j, 2]), which(indices == results[[i]][j, 3]))
      sequences[, indices[tempm[1]]:indices[tempm[length(tempm)]]] = as.DNAbin('n')
    }
  }
  sequences
}
