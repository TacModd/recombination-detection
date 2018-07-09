# 'test' based on clustering and height cutoff
# takes only one partition at a time for now

rbn.cluster = function(partitions, i){
  # get the indices for the ith partition
  indices = which(partitions$pattern.indices == i)
  # 
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(partitions$pattern.indices)))
  # 
  fit = hclust(dist(indices), method = 'ward.D2')
  # plot the fitted tree
  plot(fit)
  # return fit so user can cutree
  fit
}
