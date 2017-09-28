library(simstudy)
library(ggplot2)
library(data.table)

x <- seq(-2, 5, by = .05)

lo <- -2 + 1.2*x
p <- 1/(1 + exp(-lo))

dt <- data.table(x, p)

p1 <- ggplot(data=dt, aes(x=x, y=p)) +
  geom_line(size=6) +
  geom_area(fill = "red") +
  theme(  panel.background = element_rect(fill="grey75", color = "grey5", size = 1),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.background = element_rect(fill="grey95")
  )



jpeg("rdatagen.jpg", height = 100, width = 100)
p1

dev.off()
