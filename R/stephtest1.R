# first of stephens' tests - corresponds to equation 5 in paper
lengthdtest = function(n, d, s){
  numerator = n*choose(d, s - 1) - (s - 1)*choose(d + 1, s)
  denominator = choose(n, s)
  numerator/denominator
}
