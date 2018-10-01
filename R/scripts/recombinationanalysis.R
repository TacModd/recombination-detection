### this script extracts recombined sites as integer arrays for comparison
### and analysis of how they distribute across partitions

# my_test_result = 

##########################################################
##### first, the recombined sites of the GC2 dataset #####

# extract sites from the relevant file
G.data = read.table('GC2subset_dates.summary_of_snp_distribution.vcf', head=T, comment='&', skip=3)
G.data$POS # <-- non recombined sites

# get the sites of all partitions
allSNPs = which(!is.na(GC2ptns$pattern.indices))
length(allSNPs) # <-- total number of SNPs = 28959
sum(GC2ptns$pattern.counts) # <-- sanity check = 28959

length(allSNPs) - length(dat$POS) # <-- no. recombined SNPs = 28517

# get the recombined sites
G_recombined_sites = allSNPs[!allSNPs %in% G.data$POS]
length(GC2_recombined_sites) # <-- sanity check = 28517


#############################################################
##### second, the recombined sites from my test results #####

# get my test number of recombined sites
count = 0
for (i in 1:length(my_test_result)){
  count = count + sum(my_test_result[[i]][, 4])
}
count # no. recombined SNPs = eg. 24421

# function obtains test recombined sites (as integer array)
get.rbn.sites = function(testResult, partitions=GC2ptns, SNPs=allSNPs){
  # initialise vector to keep track of sites
  mask.sites = 1:length(partitions$pattern.indices)
  
  # gets my test actual recombined sites (as integer array)
  # for each partition
  for (i in 1:length(testResult)){
    # get the relevant partition sites (indices)
    indices = which(partitions$pattern.indices == testResult[[i]][1, 1])
    # then for each separate result
    for (j in 1:nrow(testResult[[i]])){
      # get the partition sites within the bounds of the recombination event
      sites = SNPs[SNPs >= testResult[[i]][j, 2] & SNPs <= testResult[[i]][j, 3]]
      # and add them to the overall array
      mask.sites[sites] = 0
    }
  }
  my_recombined_sites = which(mask.sites == 0)
  
  # sorts sites into ascending order
  my_recombined_sites = sort(my_recombined_sites)
  
  # return sites
  my_recombined_sites
}

# get my test recombined sites
my_recombined_sites = get.rbn.sites(my_test_result)

# a much simpler way to get the rbn site count using the above function
count = length(my_recombined_sites) 
count # <-- sanity check should = count from before eg. 24421


####################################################
##### comparing recombined sites between tests #####

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
overlap = length(my_recombined_sites[my_recombined_sites %in% GC2_recombined_sites])

# alternatively, the following gets count of partitions exclusive to my test result
exclusive = length(my_recombined_sites[!my_recombined_sites %in% GC2_recombined_sites])


# eg. comparing the 24421 set of sites against the "immutable" 28517 set (where the 
# difference in length already implies quite a few false negatives) we should obtain 
# close to a 24421 overlap  and 0 exclusive if we are on the right track. on the other 
# hand a lower and higher result respectively suggests a high degree of inaccuracy
# (a high number of false positives in addition to false negatives).

# so then we may want to look at how both sets distribute across the genomes and 
# across partitions.

# the following function is a basic means to plot the distribution of recombination across
# the genome, but it isn't very informative in R due to the number of sites

plot.arr = function(partitions, int.arr){
  # plot an empty frame the length of the genome
  plot(1, type='n', ylim = c(-1, 1), xlim = c(1, length(partitions$pattern.indices)), 
       yaxt='n', xlab='Recombination positions in genome', ylab='')
  # plot sites
  points(int.arr, c(rep(0, length(int.arr))), col='red', pch = 19)
}

plot.arr(GC2ptns, GC2_recombined_sites)
plot.arr(GC2ptns, my_recombined_sites)

# to look at the distribution across partitions, we could try using a histogram?


#############################################################
##### examining the recombined partitions between tests #####

# the following attempts to determine how different partitions
# may be discovered according to different tests - even ignoring
# differences in counts, differences in identified partitions
# could result in significantly different trees

# get the partitions for recombined sites according to Gubbins
G_rbn_ptns = GC2ptns$pattern.indices[G_recombined_sites]

# get counts for each recombined partition
G_rbn_ptns_counts = sort(table(G_rbn_ptns), decreasing = T)
sum(G_rbn_ptns_counts) # <-- sanity check = 28517
G_rbn_ptns_counts[1:20] # <-- most common partitions

length(G_rbn_ptns_counts[G_rbn_ptns_counts == 6300]) # <-- most common ptn id has 6300 rbn sites
length(G_rbn_ptns_counts[G_rbn_ptns_counts == 244]) # <-- 2 ptn ids have 244 rbn sites

length(G_rbn_ptns_counts[G_rbn_ptns_counts == 2]) # <-- 355 ptns have 2 rbn sites
length(G_rbn_ptns_counts[G_rbn_ptns_counts == 1]) # <-- 1514 ptns have 1

# get the partitions for recombined sites according to my results
my_rbn_ptns = GC2ptns$pattern.indices[my_recombined_sites]

# get the counts for each partition
my_rbn_ptns_counts = sort(table(my_rbn_ptns), decreasing = T)
my_rbn_ptns_counts[1:20] # <-- most common partitions

length(my_rbn_ptns_counts[my_rbn_ptns_counts == 1]) # <-- 0 ptns

# so looking at recombined partition counts, it is interesting to note that 
# quite a lot of partitions for the baseline results (1514 )have only one 
# recombining member - a physical impossibility using at least some of my methods. 
# moreover because Gubbins had only 442 non recombined sites, most of them must
# actually come from partitions with only one member (singleton partitions).
# that's important, because there is no way I can design a spatial test with the
# power to detect singleton partitions specifically: if these partitions are 
# important, then this represents an upper limit to how well spatial tests perform
# unless I amend the masking function appropriately. shall have to test how they
# affect the resulting tree


##################################################################
##### scripts to test the importance of singleton partitions #####

# the following script removes all the rbn sites according to Gubbins but 
# the sites for which only one ptn member was identified (=1514)

# table of partitions with only 1 recombining member according to baseline results
ptns_to_keep = G_rbn_ptns_counts[G_rbn_ptns_counts == 1]
# table of all the recombining partitions EXCLUDING those partitions
G_ptns_to_remove = G_rbn_ptns[!G_rbn_ptns %in% as.integer(names(ptns_to_keep))]
# vector of partitions with only 1 recombining member according to baseline results
ptns_to_keep = which(G_rbn_ptns %in% as.integer(names(ptns_to_keep)))
# vector of partitions with more than 1 recombining member according to baseline results
G_ptns_to_remove = which(!G_rbn_ptns %in% as.integer(names(ptns_to_keep)))
# all the recombining sites EXCLUDING the sites with only 1 recombining member
indices_to_remove = G_recombined_sites[G_ptns_to_remove]

GC2m_Gs_k = as.matrix(GC2d)

GC2m_Gs_k[, indices_to_remove] = as.DNAbin('n')

GC2m_Gs_k = as.list(GC2m_Gs_k)

# the following script removes all the rbn sites according to Gubbins but 
# the singleton partitions - the hypothetical best result using spatial tests
# to test for each partition specifically (=1504)

GC2_ptns_counts = sort(table(GC2ptns$pattern.indices), decreasing=T)

GC2_singleton_ptns = GC2_ptns_counts[GC2_ptns_counts == 1]

rbn_non_singleton_ptns = which(!G_rbn_ptns %in% as.integer(names(GC2_singleton_ptns)))

singleton_indices_to_keep = G_recombined_sites[rbn_non_singleton_ptns]

GC2m_s_k = as.matrix(GC2d)

GC2m_s_k[, singleton_indices_to_keep] = as.DNAbin('n')

GC2m_s_k = as.list(GC2m_s_k)

write.dna(GC2m_s_k, 'Gub_sk.phy')

# the following script removes only those singleton partitions that
# Gubbins identified

GC2_ptns_counts = sort(table(GC2ptns$pattern.indices), decreasing=T)

GC2_singleton_ptns = GC2_ptns_counts[GC2_ptns_counts == 1]

rbn_singleton_ptns = which(G_rbn_ptns %in% as.integer(names(GC2_singleton_ptns)))

singleton_indices_to_remove = G_recombined_sites[rbn_singleton_ptns]

GC2m_s_r = as.matrix(GC2d)

GC2m_s_r[, singleton_indices_to_remove] = as.DNAbin('n')

GC2m_s_r = as.list(GC2m_s_r)

write.dna(GC2m_s_r, 'Gub_sr.phy')


###########################################################
##### checking partition patterns (1/2s vs 1/2/3/4s) ######

# one way to deal with singleton partition problem is to consider
# merging singleton partitions somehow

# first I would like to get a broad overview of partition pattern similarity 
# between singleton partitions and other partitions. to do that, we want to 
# using hamming distance.

# first, the singleton partition patterns:
grp1 = GC2ptns$patterns[as.integer(names(ptns_to_remove))]
length(gr1) # = 1514
# second, all other partitions
other_partitions = G_rbn_ptns_counts[G_rbn_ptns_counts != 1]
grp2 = unique(GC2ptns$patterns[as.integer(names(other_partitions))])
length(grp2) # = 2506 - 1514 = 992

# now finding the best hamming distances for each singleton partition
answers = c()

for (i in 1:length(grp1)){
  answer = length(strsplit(grp1[1], '')[[1]])
  for (j in 1:length(grp2)){
    answer = min(sum(strsplit(grp1[i], '')[[1]] != strsplit(grp2[j], '')[[1]]), answer)
  }
  answers = c(answers, answer)
}

# and seeing what the results look like
sort(table(answers), decreasing = T)

# a look at the resulting table suggests many singleton partitions are not especially
# similar to non-singleton partitions

# now I'm curious to see: if changing all triallelic and quadallelic partition
# patterns to biallelic patterns makes a big difference. I'm not sure how I'd use this
# information just now, but it might still be useful to know

bi.answers = c()

for (i in 1:length(grp1)){
  str1 = strsplit(grp1[i], '')[[1]]
  str1[str1 == '3'] = '2'
  str1[str1 == '4'] = '2'
  answer = length(strsplit(grp1[1], '')[[1]])
  for (j in 1:length(grp2)){
    str2 = strsplit(grp2[j], '')[[1]]
    str2[str2 == '3'] = '2'
    str2[str2 == '4'] = '2'
    answer = min(sum(str1 != str2), answer)
  }
  bi.answers = c(bi.answers, answer)
}

# and seeing what the results look like
sort(table(bi.answers), decreasing = T)

length(answers) # 1514
length(bi.answers) # 1514

bi.answers
  1   0   2   3   4   5   6   7   9   8  11  10  12  13  14  15  19 
456 368 166 104  81  53  47  41  39  33  25  23  18  17  11   6   6 
 16  22  17  18  20  25  21  23 
  4   4   3   3   2   2   1   1 

answers
  1   2   3   4   5   6   9   7   8  10  11  12  13  14  15  17  16 
565 250 146  97  68  60  46  41  31  31  25  25  24  21  19  13  11 
 18  19  21  22  23  25  26  24  20 
 10   5   5   4   4   4   4   3   2

# a fair amount of partitions are very similar, as expected, but there are
# also plenty that are not. changing 3 and 4 alleles to 2 alleles appears
# to do only a moderate amount to resolve singleton partitions
