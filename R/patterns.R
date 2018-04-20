# patterns generates integers for columns of sites: the order of these integers 
# reflects the pattern of site base variation across sequences
patterns = function(sequences){
  pttns = matrix(nrow = nrow(sequences), ncol = ncol(sequences))
  for (i in 1:ncol(sequences)){
    k = unique(sequences[, i])
    for (j in 1:nrow(sequences)){
      pttns[j, i] = match(sequences[j, i], k)
    }
  }
  pttns
}
