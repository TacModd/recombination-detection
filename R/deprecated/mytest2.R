### OUTDATED: test no longer intended to be completed ###
# incomplete
# window based on set n
# for i in range pattern.IDs
spt2 = function(ptns, i){
  indices = which(ptns$pattern.indices == i)
  p = length(indices)/length(ptns$pattern.indices)
  n = 2
  
  if (ptns$pattern.counts[i] > 1){
    while (1 - (1-p)**n - n*(p)*(1-p)**(n-1) < 0.05){
      n = n + 1
    }
    
    for (j in range(1:length(indices) - 1)){
      q = length(which(ptns$pattern.indices[j:j+n] == i)) - 1
      if (length(q) > 1){
        r = 1 - pbinom(q, n, p)
      }
    }
    r
  }
}
