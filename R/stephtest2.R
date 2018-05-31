# 2nd of Stephens' tests - corresponds to equations 8, 9 and 10 in paper
lengthrtest = function(s, r, k){
  pm = 0
  for (i in k:r){
    numerator = choose(s + r - i - 3, r - i)
    denominator = choose(s + r - 2, r)
    pm = pm + numerator/denominator
  }
  1 - (1 - pm)**(s-1)
}
