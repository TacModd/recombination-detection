# takes a multiple aligned sequence object and generates integer 
# sequences (patterns) representing columns of sites: these patterns 
# reflect the pattern of site base variation across DNA sequences. 
# returns a pattern object.
get.patterns = function(sequences){
  # initialise matrix
  patterns = matrix(nrow = nrow(sequences), ncol = ncol(sequences))
  # for each column in the aligned sequences:
  for (i in 1:ncol(sequences)){
    # get the unique values
    k = unique(sequences[, i])
    # generate the patterns (by matching against the vector of
    # unique values)
    for (j in 1:nrow(sequences)){
      patterns[j, i] = match(sequences[j, i], k)
    }
  }
  # return the patterns object
  patterns
}
