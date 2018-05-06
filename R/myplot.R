# plots the distribution of sites i from ptns across the length of the genome
# need to add functionality to show hypothetical recombined regions
myplot = function(ptns, i){
  indices = which(ptns$pattern.indices == i)
  fit = cmdscale(dist(indices), eig=FALSE, k=1)
  plot(fit[, 1] + (max(indices)-min(indices))/2, rep(0, length(indices)), pch=19, xlim = c(0, length(ptns$pattern.indices)))
}
