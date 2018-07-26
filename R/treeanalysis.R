library(NELSI)

library(phangorn)



### Simulate a time tree and attach dates to the names of the tips

time_tree<- rtree(10)

# allnode is from NELSI
time_tree$tip.label <- paste0(time_tree$tip.label, '_', round(allnode.times(time_tree, tipsonly = T, rev =T), 2)+2010)

plot(time_tree)

add.scale.bar()



#### Simulate rates to get a phylogram

phylogram <- time_tree

phylogram$edge.length <- phylogram$edge.length * abs(rnorm(length(phylogram$edge.length), 0.001, 0.1))

plot(phylogram)

add.scale.bar()

write.tree(unroot(phylogram), file = 'phylogram.tree')



####################

#### Comparing trees

tr1 <- unroot(rtree(10))

tr2 <- rNNI(tr1)

par(mfrow = c(1, 2))

plot(tr1)

plot(tr2, direction = 'left')

# unroot tree before dist.topo (ph95) distance measure
# ph85 is the default
dist.topo(tr1, tr2, method = 'PH85')

dist.topo(tr1, tr2, method = 'score')

dist.topo(unroot(rtree(10)), unroot(rtree(10)))

tr_alt <- tr1

tr_alt$edge.length <- rexp(length(tr_alt$edge.length), rate = 10)

plot(tr1)

plot(tr_alt)

dist.topo(tr1, tr_alt, method = 'score')



# compare total tree length

diff(c(sum(tr1$edge.length), sum(tr_alt$edge.length)))

par(mfrow = c(1, 1))

diff_br <- ladderize(tr1)$edge.length - ladderize(tr_alt)$edge.length

plot(tr1)

edgelabels(round(diff_br, 2))
