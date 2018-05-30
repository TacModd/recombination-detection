# an updated version of myplot2 that works with the new output format

myplot2v2 = function(ptns, res, j){
  i = res[[j]][1, 1]
  indices = which(ptns$pattern.indices == i)
  #fit = cmdscale(dist(indices), eig=FALSE, k=1)
  #plot(fit[, 1], rep(0, length(indices)), pch=19)
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(ptns$pattern.indices)))
  for (k in 1:nrow(res[[j]])){
    points(res[[j]][k, 2:3], c(0, 0), col='red', pch = 19)
    lines(res[[j]][k, 2:3], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
  }
}
