permutan <- function(chrlist, attempts = 10000, ncor=3){
  data.fc <- do.call("FeatureCalls.to.GFF.like", chrlist)
  data.rb <- do.call("rbind", chrlist)
}