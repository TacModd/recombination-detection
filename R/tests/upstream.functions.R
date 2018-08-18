# tests get.patterns and get.partitions

# initialise a matrix of of sequences with constant sites
sequences = matrix('a', nrow=5, ncol=5)
# introduce 3 SNPs (2 partitions since columns 1 & 4 match)
sequences[2, 1] = 'g'
sequences[2, 4] = 'g'
sequences[3, 3] = 'g'

# test code body of get.patterns
for (i in 1:ncol(sequences)){
  # get the unique values
  k = unique(sequences[, i])
  # generate the patterns (by matching against the vector of
  # unique values)
  for (j in 1:nrow(sequences)){
    patterns[j, i] = match(sequences[j, i], k)
  }
}
# should produce a matrix of integers corresponding to pattern of variation for each column
View(patterns)
# should produce a character array of integer columns/partitions
output = sapply(1:ncol(patterns), function(x) paste(patterns[, x], collapse = ''))
output
# [1] "12111" "11111" "11211" "12111" "11111"

# should produce a partition object
partitions = get.partitions(output)
# with the following values:
# partitions.patterns ---------> '12111' '11211'
# partitions.pattern.IDs ------> 1 2
# partitions.pattern.indices --> 1 NA 2 1 NA
# partitions.pattern.counts ---> 2 1

View(partitions)

# now testing compatability with DNAbin class

library(phangorn)

genome = as.DNAbin(sequences) 
genome # should be a DNAbin matrix
as.list(genome) # should be a DNAbin list

TESTING = get.patterns(genome) # should succeed
TESTING = get.patterns(as.list(genome)) # should fail
TESTING = get.patterns(as.matrix(genome)) # should succeed
# successful run should produce same output as before
TESTING
# [1] "12111" "11111" "11211" "12111" "11111"
