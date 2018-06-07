# script to generate toy sequence data with recombination for testing
# (plots images of data as well)

library(phangorn)

set.seed(12345678)
tr_core <- rtree(20)
set.seed(12345678)
tr_core$edge.length <- rlnorm(n = length(tr_core$edge.length), 
                              meanlog = -6.3, sdlog = 1)
plot(tr_core)
edgelabels(round(tr_core$edge.length, 2))

set.seed(12345678)
tr_recombine = rNNI(tr_core, moves = 3)
par(mfrow = c(1, 2))
plot(tr_core)
edgelabels()
plot(tr_recombine)
edgelabels()

tr_recombine$edge.length[8] = tr_recombine$edge.length[8] * 50
plot(tr_recombine)
edgelabels()


core_genome = as.DNAbin(simSeq(tr_core, l = 1000))
recombination = as.DNAbin(simSeq(tr_recombine, l = 100, rate = 10))

concatenated = cbind(core_genome[, 1:500], recombination, core_genome[, 501:1000])
par(mfrow = c(1, 2))
image(concatenated)

pttns = patterns(concatenated)
pttns = sapply(1:ncol(pttns), function(x) paste(pttns[, x], collapse = ''))
ptns = partitions(pttns)

m = matrix(ptns$pattern.indices, nrow = length(ptns$pattern.indices), ncol = nchar(ptns$patterns[1]), byrow = FALSE)
image(x=1:length(ptns$pattern.indices), y=1:nchar(ptns$patterns[1]), z = m, col=rainbow(length(ptns$pattern.IDs)))
