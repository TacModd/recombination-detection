# this function returns a partition object, which includes the number of unique partition
# patterns, their ids, their indices, and their counts
partitions = function(patterns){
  #constant.sites = unique(grep('1{66}', GC2pttns, value = TRUE))
  #sites = GC2pttns[GC2pttns != constant.sites]
  # this way is faster
  sites = unique(patterns[patterns != as.character(paste(rep(1, nchar(patterns[1])), collapse=''))])
  list(patterns = sites, pattern.IDs = match(sites, sites), pattern.indices = match(patterns, sites),
       pattern.counts = c(table(match(patterns, sites))))
  
}
