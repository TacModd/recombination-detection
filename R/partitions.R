# this function returns a partition object, which includes the number of unique partition
# patterns, their ids, their indices, and their counts in a list of lists
partitions = function(patterns){
  # obtain partitions from pattern object (ignore constant sites)
  partitions = unique(patterns[patterns != as.character(paste(rep(1, nchar(patterns[1])), collapse=''))])
  # return partition object
  list(patterns = partitions, pattern.IDs = match(partitions, partitions), pattern.indices = match(patterns, partitions),
       pattern.counts = c(table(match(patterns, partitions))))
  
}
