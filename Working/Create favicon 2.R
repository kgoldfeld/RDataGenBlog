library(simstudy)
library(ggplot2)

defI <- defData(varname = "x", formula = "0;1", dist = "uniform")
defI <- defData(defI, varname = "y", formula = "2 + .5*x + 2.5*x^2", dist = "nonrandom")

set.seed(450)

dt <- genData(100, defI)

p1 <- ggplot(data=dt, aes(x=x, y=y)) +
  geom_smooth(color = "#fd9700", se = FALSE, size = 10) +
  theme(  panel.background = element_rect(fill="grey25", color = "grey5", size = 1),
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
