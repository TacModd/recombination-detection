# takes a partition object, a result object and an integer 
# for the jth result. determines the ith partition from j and
# plots the distribution of sites from the ith partition AND 
# the recombined sites from the jth result across the length 
# of the genome. recombined regions are highlighted in pink

plot.rbn = function(partitions, results, j){
  i = res[[j]][1, 1]
  indices = which(partitions$pattern.indices == i)
  #fit = cmdscale(dist(indices), eig=FALSE, k=1)
  #plot(fit[, 1], rep(0, length(indices)), pch=19)
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(partitions$pattern.indices)))
  for (k in 1:nrow(results[[j]])){
    points(results[[j]][k, 2:3], c(0, 0), col='red', pch = 19)
    lines(results[[j]][k, 2:3], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
  }
}
