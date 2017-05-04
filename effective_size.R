#number of effective samples per iteration

library(coda)
s <- readRDS("paschold/output_stick_paschold.rds")
sapply(1:4, function(chain) effectiveSize(s[[chain]][[1]]$alpha)/5000)

