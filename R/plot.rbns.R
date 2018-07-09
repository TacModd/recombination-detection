# takes a partition object and a result object.
# plots all identified recombined regions and bounds
# according to result object

plot.rbns = function(partitions, results){
  # plot an empty frame the length of the genome
  plot(1, type='n', ylim = c(-1, 1), xlim = c(1, length(partitions$pattern.indices)), 
      yaxt='n', xlab='Recombination positions in genome', ylab='')
  # then for each recombined partition according to result object
  for (k in 1:length(results)){
    #i = res[[k]][1, 1]
    #indices = which(ptns$pattern.indices == i)
    # overlap each recombined site according to kth result
    for (j in 1:nrow(results[[k]])){
      points(results[[k]][j, 2:3], c(0, 0), col='red', pch = 19)
      lines(results[[k]][j, 2:3], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
    }
  }
}
