rm(list=ls()) 
library(ggplot2)

setwd("../data/")
adr_dp = read.csv("adr_dp.csv", row.names = 1)
colnames(adr_dp) = gsub("\\.", " ", colnames(adr_dp))
I = intersect(rownames(adr_dp), colnames(adr_dp))
adr_dp = adr_dp[I, I]

adr_dp_mat = as.matrix(adr_dp)
known_pairs = diag(adr_dp_mat)

unknown_pairs = c(adr_dp_mat[upper.tri(adr_dp_mat)], adr_dp_mat[lower.tri(adr_dp_mat)])
wilcox.test(known_pairs, unknown_pairs, alternative = "less")$p.value