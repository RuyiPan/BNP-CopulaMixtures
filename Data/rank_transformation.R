temp <- read.csv("../Data/Wine/winequality-red2.csv")
colnames(temp)
wine_red_rank <- apply(temp[, c(1,3,8)], 2, function(x) {
  rank(x) / length(x)
})
saveRDS(wine_red_rank, file = "../Data/Wine/wine_red_rank.rds")
cor(wine_red_rank, method = "kendall")
