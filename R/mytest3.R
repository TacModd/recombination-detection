# 'test' based on clustering and height cutoff

spt3 = function(ptns, i){
  indices = which(ptns$pattern.indices == i)
  fit = cmdscale(dist(indices), eig=FALSE, k=1)
  plot(fit[, 1], rep(0, length(indices)), pch = 19)
  fit2 = hclust(dist(indices), method = 'ward.D2')
  plot(fit2)
}
