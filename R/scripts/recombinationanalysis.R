### this script extracts recombined sites as integer arrays for comparison
### and analysis of how they distribute across partitions


### first, the recombined sites of the GC2 dataset

# extracting sites from the relevant file
dat = read.table('GC2subset_dates.fasta.vcf', head=T, comment='&', skip=3)
dat$POS # <-- non recombined sites

# getting the sites of all partitions
allSNPs = which(GC2ptns$pattern.indices > 0)
length(allSNPs) # <-- total number of SNPs = 28959
length(allSNPs) - length(dat$POS) # <-- no. recombined SNPs = 9926

# gets the recombined sites
GC2_recombined_sites = allSNPs[!allSNPs %in% dat$POS]
length(GC2_recombined_sites) # <-- sanity check = 9926


### second, the recombined sites from my test results

# gets my test total number of recombined sites
count = 0
for (i in 1:length(my_test_result)){
  count = count + sum(my_test_result[[i]][, 4])
}
count # no. recombined SNPs = eg. 7124

# gets my test actual recombined sites (as integer array)
my_recombined_sites = c()
# for each partition
for (i in 1:length(my_test_result)){
  # get the relevant partition sites (indices)
  indices = which(GC2ptns$pattern.indices == my_test_result[[i]][1, 1])
  # then for each separate result
  for (j in 1:nrow(my_test_result[[i]])){
    # store sites within the listed bounds in a temporary array
    tempa = indices[which(indices == my_test_result[[i]][j, 2]):which(indices == my_test_result[[i]][j, 3])]
    # and add them to the overall array
    my_recombined_sites = c(my_recombined_sites, tempa)
  }
}

length(my_recombined_sites) # <-- sanity check should = count from before eg. 7124

# sorts sites into ascending order
my_recombined_sites = sort(my_recombined_sites)


### example comparisons

# if recombination detection has been successful then a molecular clock structure 
# should be sucessfully recovered by tempest. this is the case for the previously
# obtained results and checking how our resulting tree performs should be the first
# port of call.

# the fact that two other programs (which use different methods for detection) agree 
# on recombined sites and recovered molecular clock data affords us a fair degree of 
# confidence that - even if not exactly correct - those results are generally in the 
# right direction (i.e. the identified sites are mostly correct).

# thus, ideally, our attempts at recombination detection should yield a similar suite
# of sites, and we can attempt to compare our obtained sites against the previously
# obtained sites as of way of determining - if we are going wrong - how badly and why.


# with that in mind, the following compares the overlap between two given sets:
length(my_recombined_sites[my_recombined_sites %in% GC2_recombined_sites])


# eg. comparing the 7124 set of sites against the immutable 9926 set (where the 
# difference in length already implies some false negatives) we should obtain close 
# to 7124 if we are on the right track. on the other hand a result of 2258 suggests a
# high degree of inaccuracy (a high number of false positives AND negatives).

# so then we might want to look at how both sets distribute across the genomes and 
# across partitions.

# the following function is a basic means to plot the distribution of recombination across
# the genome, but it isn't very informative

plot.arr = function(partitions, int.arr){
  # plot an empty frame the length of the genome
  plot(1, type='n', ylim = c(-1, 1), xlim = c(1, length(partitions$pattern.indices)), 
       yaxt='n', xlab='Recombination positions in genome', ylab='')
  # plot sites
  points(int.arr, c(rep(0, length(int.arr))), col='red', pch = 19)
}

plot.arr(GC2ptns, GC2_recombined_sites)
plot.arr(GC2ptns, my_recombined_sites)

# to look at the distribution across partitions, we can try using a histogram

# get the partitions for each recombined site
GC2_rbn_ptns = GC2ptns$pattern.indices[GC2_recombined_sites]

# get counts for the number of each recombined partition
GC2_rbn_ptns_counts = sort(table(GC2_rbn_ptns), decreasing = T)
GC2_rbn_ptns_counts[1:20] # <-- most common partitions

# my recombined sites
my_rbn_ptns = GC2ptns$pattern.indices[hmm]

# get counts again
my_rbn_ptns_counts = sort(table(my_rbn_ptns), decreasing = T)
my_rbn_ptns_counts[1:20] # <-- most common partitions

length(my_rbn_ptns_counts) # <-- 646
length(GC2_rbn_ptns_counts) # <-- 771

length(GC2_rbn_ptns_counts[GC2_rbn_ptns_counts == 1]) # <-- 451
length(my_rbn_ptns_counts[my_rbn_ptns_counts == 1]) # <-- 0

# comparing the above tables reveals the degree of overlap between
# the most common partitions for each set. it seems that the baseline
# results are much more densely concentrated within just a few common
# partitions, whereas my results are much more evenly spread amongst
# partitions

# it is also interesting to note that quite a lot of partitions for the
# baseline results have only one recombining member - a physical impossibility
# using at least some of my methods. assuming these partitions are important 
# (which can be tested) this may represent an upper bound on what accuracy is 
# achievable using binomial tests of spatial distribution, and perhaps for 
# poisson tests as well. 451/9926 is a fair chunk (about 4.5%).

# so to test the relevance of these sites we need to remove them from the list
# and then mask the remaining sites

# partitions with only 1 recombined member
ptns_to_remove = GC2_rbn_ptns_counts[GC2_rbn_ptns_counts == 1]
# removing them
GC2_rbn_ptns_revised = 
  GC2_rbn_ptns[!GC2_rbn_ptns %in% 
                 as.integer(names(ptns_to_remove))]
# sites
sites_to_remove = which(!GC2_rbn_ptns %in% as.integer(names(ptns_to_remove)))

GC2_rbn_sites_revised = GC2_recombined_sites[sites_to_remove]

GC2m = as.matrix(GC2d)

GC2m[, GC2_rbn_sites_revised] = as.DNAbin('n')

# write and we are ready to phyml

