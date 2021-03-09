# Small example for Manly distance
# make up some test data; two grouping features (2-3 traits)
# calculate via hand:



# calculate via ade4 package:
vec <- sub("\\_.*", "\\1", names(data))
blocks <- rle(vec)$lengths
data <- prep.fuzzy(data, blocks)
data <- ktab.list.df(list(data))
dist.ktab(data, type = "F")
