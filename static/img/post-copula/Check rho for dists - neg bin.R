set.seed(58)

dt <- genCorData(250, mu = c(0, 0), sigma = 1, 
                 rho = 0.8, corstr = "cs", 
                 cnames = c("X1","X2"))

dt[, U1 := pnorm(X1)]
dt[, U2 := pnorm(X2)]

dt[, Y1 := qpois(U1, 3)]
dt[, Y2 := qpois(U2, 3)]

dt[, B1 := qbinom(U1, 1, .5)]
dt[, B2 := qbinom(U2, 1, .5)]

dt[, N1 := qnbinom(U1, 5, .5)]
dt[, N2 := qnbinom(U2, 5, .5)]


dtCor <- dt[, .(X = cor(X1, X2), U = cor(U1, U2),
                Y = cor(Y1, Y2), B = cor(B1, B2),
                N = cor(N1, N2)
)]

dtS <- dtCor[,.(X = paste("Est. corr:", sprintf("%1.2f",X)),
         U = paste("Est. corr:", sprintf("%1.2f",U)),
         Y = paste("Est. corr:", sprintf("%1.2f",Y)),
         B = paste("Est. corr:", sprintf("%1.2f",B)),
         N = paste("Est. corr:", sprintf("%1.2f",N))
         )]

p1 <- ggplot(data = dt, aes(x=X1, y=X2)) +
  geom_point(color="#6285BA", size = 1) +
  theme_ksg("grey95") +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  ggtitle("Normal") +
  annotate(geom = "text", label =  dtS$X, x = -2, y = 2, fontface = 2)

p2 <- ggplot(data = dt, aes(x=U1, y=U2)) +
  geom_point(color="#6285BA", size = 1) +
  theme_ksg("grey95") +
  scale_x_continuous(limits = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  theme(legend.position = "none") +
  ggtitle("Uniform")  +
  annotate(geom = "text", label =  dtS$U, x = .10, y = 0.92, fontface = 2)


p3 <- ggplot(data = dt, aes(x=Y1, y=Y2)) +
  geom_jitter(color="#6285BA", height = .15, width = .15, size = 1) +
  theme_ksg("grey95") +
  scale_x_continuous(limits = c(-.2, 7.2), breaks = c(0:7)) +
  scale_y_continuous(limits = c(-.2, 7.2), breaks = c(0:7)) +
  theme(legend.position = "none") +
  ggtitle("Poisson") +
  annotate(geom = "text", label =  dtS$Y, x = 0.6, y = 6.6, fontface = 2)

p4 <- ggplot(data = dt, aes(x=B1, y=B2), size = 1) +
  geom_jitter(color="#6285BA", height = .1, width = .1, size = 1) +
  theme_ksg("grey95") +
  scale_x_continuous(limits = c(-.2, 1.2), breaks = c(0:1)) +
  scale_y_continuous(limits = c(-.2, 1.2), breaks = c(0:1)) +
  theme(legend.position = "none") +
  ggtitle("Binary")+
  annotate(geom = "text", label =  dtS$B, x = .5, y = 0.5, fontface = 2)

p5 <- ggplot(data = dt, aes(x=N1, y=N2), size = 1) +
  geom_jitter(color="#6285BA", height = .1, width = .1, size = 1) +
  theme_ksg("grey95") +
 # scale_x_continuous(limits = c(-.2, 1.2), breaks = c(0:1)) +
#  scale_y_continuous(limits = c(-.2, 1.2), breaks = c(0:1)) +
  theme(legend.position = "none") +
  ggtitle("Negative Binomial")+
  annotate(geom = "text", label =  dtS$N, x = 2.5, y = 20, fontface = 2)

gridExtra::grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
