# takes a partition object, a result object and an integer 
# for the kth result. determines the ith partition from k and
# plots the distribution of sites from the ith partition AND 
# the recombined sites from the kth result across the length 
# of the genome. recombined regions are highlighted in pink

plot.rbn = function(partitions, results, k){
  # get ith partition from k
  i = res[[k]][1, 1]
  # get indices of sites for ith partition
  indices = which(partitions$pattern.indices == i)
  # plot the partition sites
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(partitions$pattern.indices)), 
      yaxt='n', xlab='Partition positions in genome', ylab='')
  # overlap recombined sites according to kth result
  for (j in 1:nrow(results[[k]])){
    points(results[[k]][j, 2:3], c(0, 0), col='red', pch = 19)
    lines(results[[k]][j, 2:3], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
  }
}
