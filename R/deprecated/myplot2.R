### OUTDATED: works only with old output format ###
# plots the distribution of sites i from ptns AND recombined sites j from res 
# across the length of the genome
# recombined regions are highlighted in pink
myplot2 = function(ptns, i, res, j){
  indices = which(ptns$pattern.indices == i)
  #fit = cmdscale(dist(indices), eig=FALSE, k=1)
  #plot(fit[, 1], rep(0, length(indices)), pch=19)
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(ptns$pattern.indices)))
  for (k in 1:length(res[[j]])){
    points(res[[j]][[k]][1:2], c(0, 0), col='red', pch = 19)
    lines(res[[j]][[k]][1:2], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
  }
}
