library(simstudy)
library(ggplot2)
library(data.table)

x <- seq(-2, 8, by = .05)

lo <- -2 + 1.5*x
p1 <- 1/(1 + exp(-lo))

lo <- -2 + 1.3*x
p2 <- 1/(1 + exp(-lo))

lo <- -2 + 1.1*x
p3 <- 1/(1 + exp(-lo))

lo <- -2 + 0.9*x
p4 <- 1/(1 + exp(-lo))

lo <- -2 + 0.8*x
p5 <- 1/(1 + exp(-lo))

dt <- data.table(x = rep(x,5), p = c(p1,p2,p3,p4,p5), group = rep(1:5, each=length(x)))

p1 <- ggplot(data=dt, aes(x=x, y=p, group=group)) +
  geom_line(aes(color = factor(group)), size = 4) +
  theme(  panel.background = element_rect(fill="grey40", color = "grey25", size = 1),
          panel.grid.minor= element_blank(),
          panel.grid.major.x = element_line(size=.5, color = "white"),
          panel.grid.major.y = element_line(size=.5, color = "white"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.background = element_rect(fill="grey95")
  )



jpeg("rdatagen.jpg", height = 100, width = 100)
p1

dev.off()
