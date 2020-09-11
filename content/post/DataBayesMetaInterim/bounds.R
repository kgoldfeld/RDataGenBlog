library(ldbounds)
t <- c(40/80, 60/80, 80/80)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
threshholds <- obf.bd$diff.pr

summary(obf.bd)
plot(obf.bd)
