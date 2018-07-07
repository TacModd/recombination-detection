# plots all identified recombined regions and bounds

myplot3 = function(ptns, res){
  plot(1, type='n', ylim = c(-1, 1), xlim = c(1, length(ptns$pattern.indices)))
  for (i in 1:length(res)){
    #j = res[[i]][1, 1]
    #indices = which(ptns$pattern.indices == j)
    #points(indices, rep(0, length(indices)), pch=19)
    for (k in 1:nrow(res[[i]])){
      points(res[[i]][k, 2:3], c(0, 0), col='red', pch = 19)
      lines(res[[i]][k, 2:3], c(0, 0), lwd = c(100, 1), lend = 1, col='pink')
    }
  }
}
