# 'test' based on clustering and height cutoff
# takes only one partition at a time for now

rec.testcluster = function(ptns, i){
  indices = which(ptns$pattern.indices == i)
  #fit = cmdscale(dist(indices), eig=FALSE, k=1)
  #plot(fit[, 1], rep(0, length(indices)), pch = 19)
  plot(indices, rep(0, length(indices)), pch=19, xlim = c(1, length(ptns$pattern.indices)))
  fit = hclust(dist(indices), method = 'ward.D2')
  plot(fit)
  # return fit so user can cutree
  fit
}
