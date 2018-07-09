# takes a partition object. plots the distribution of sites 
# from the ith partition across the length of the genome

plot.ptn = function(partitions, i){
  # get indices of sites for ith partition
  indices = which(partitions$pattern.indices == i)
  #fit = cmdscale(dist(indices), eig=FALSE, k=1)
  #plot(fit[, 1], rep(0, length(indices)), pch=19)
  # plot sites
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(0, length(partitions$pattern.indices)))
}
