update_stats = function(agg, x) {
  n1 = agg$n
  agg$n = agg$n + 1
  delta = x - agg$M1
  delta_n = delta / agg$n
  delta_n2 = delta_n * delta_n
  term1 = delta * delta_n * n1
  agg$M1 = agg$M1 + delta_n
  agg$M4 = agg$M4 + term1 * delta_n2 * (agg$n * agg$n - 3 * agg$n + 3) + 6 * delta_n2 * agg$M2 - 4 * delta_n * agg$M3
  agg$M3 = agg$M3 + term1 * delta_n * (agg$n - 2) - 3 * delta_n * agg$M2
  agg$M2 = agg$M2 + term1
  return(agg)
}

combine_stats = function(agg1, agg2) {
  combined = data.table()
  combined$n = agg1$n + agg2$n

  delta = agg2$M1 - agg1$M1
  delta2 = delta * delta
  delta3 = delta * delta2
  delta4 = delta2 * delta2

  combined$M1 = (agg1$n * agg1$M1 + agg2$n * agg2$M1) / combined$n

  combined$M2 = agg1$M2 + agg2$M2 +
    delta2 * agg1$n * agg2$n / combined$n

  combined$M3 = agg1$M3 + agg2$M3 +
    delta3 * agg1$n * agg2$n * (agg1$n - agg2$n) / (combined$n * combined$n)
  combined$M3 = combined$M3 + 3.0 * delta * (agg1$n * agg2$M2 - agg2$n * agg1$M2) / combined$n

  combined$M4 = agg1$M4 + agg2$M4 + delta4 * agg1$n * agg2$n * (agg1$n * agg1$n - agg1$n * agg2$n + agg2$n * agg2$n) /
    (combined$n * combined$n * combined$n)
  combined$M4 = combined$M4 + 6.0 * delta2 * (agg1$n * agg1$n * agg2$M2 + agg2$n * agg2$n * agg1$M2) / (combined$n * combined$n) +
    4.0 * delta * (agg1$n * agg2$M3 - agg2$n * agg1$M3) / combined$n

  return (combined)
}

# agg_all = combine_stats(agg, agg2)
#
# e1071::moment(all, 2, center = T)*length(all)
# e1071::moment(all, 3, center = T)*length(all)
# e1071::moment(all, 4, center = T)*length(all)
# mean(x)
# agg_kur(agg_all)
# agg_skew(agg_all)

agg_mean = function(agg) agg$M1
agg_variance = function(agg) agg$M2 / (agg$n - 1.0)
agg_sd = function(agg) sqrt(agg_variance(agg))

agg_skew = function(agg) {
  g = (sqrt(agg$n) * agg$M3) / (agg$M2^1.5)

  sqrt((agg$n * (agg$n - 1))) * g / (agg$n - 2)
}


agg_kur = function(agg) {
  g = (agg$n * agg$M4) / (agg$M2 * agg$M2) - 3.0
  ((agg$n - 1) / ((agg$n - 2) * (agg$n - 3))) * ((agg$n + 1) * g + 6)
}

agg_n = function(agg) agg$n
