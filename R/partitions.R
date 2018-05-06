# this function returns a partition object, which includes the number of unique partition
# patterns, their ids, their indices, and their counts
partitions = function(patterns){
  #constant.sites = unique(grep(as.character(paste('1{', nchar(patterns[1]), '}', sep='')), patterns, value = TRUE))
  #sites = [patterns[patterns != constant.sites]
  # this way is faster
  sites = unique(patterns[patterns != as.character(paste(rep(1, nchar(patterns[1])), collapse=''))])
  list(patterns = sites, pattern.IDs = match(sites, sites), pattern.indices = match(patterns, sites),
       pattern.counts = c(table(match(patterns, sites))))
  
}
