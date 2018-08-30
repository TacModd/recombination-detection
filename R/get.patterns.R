# takes a multiple aligned sequence object (as a matrix) and generates 
# integer sequences (patterns) representing columns of sites: these 
# patterns reflect the pattern of site base variation across DNA 
# sequences. returns a pattern object.
get.patterns = function(sequences){
  # initialise matrix to store patterns
  patterns = matrix(nrow = nrow(sequences), ncol = ncol(sequences))
  # for each column in the aligned sequences:
  for (i in 1:ncol(sequences)){
    
    
    # get the unique values in order of frequency (so major allele is assigned 1)
    k = as.raw(names(sort(table(as.numeric(sequences[, i])), decreasing = T)))
    
    ## get the unique values
    #k = unique(sequences[, i])
    
    
    # generate the patterns (by matching against the vector of
    # unique values)
    for (j in 1:nrow(sequences)){
      patterns[j, i] = match(sequences[j, i], k)
    }
  }
  # create the final patterns object (single-column matrix of strings)
  patterns = sapply(1:ncol(patterns), function(x) paste(patterns[, x], collapse = ''))
  # return the patterns object
  patterns
}
