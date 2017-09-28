library(simstudy)
library(ggplot2)


defC <- defData(varname = "rInt", formula = 0, variance = .75, dist = "normal", id = "cluster")
defC <- defData(defC, varname = "rSlp", formula = 0, variance = .5, dist = "normal")
defC <- defData(defC, varname = "num", formula = 100)

defI <- defDataAdd(varname = "x", formula = "0;1", dist = "uniform")
defI <- defDataAdd(defI, varname = "y", formula = "rInt + (2+rSlp) * x")


set.seed(450)

dt <- genData(4, defC)
dt <- genCluster(dtClust = dt, "cluster", "num", "id")
dt <- addColumns(dtDefs = defI, dtOld = dt)


p1 <- ggplot(data=dt, aes(x=x, y=y, group = cluster))+
  geom_smooth(aes(color = factor(cluster)), se = FALSE, size = 4) +
  theme(  panel.background = element_rect(fill="grey30", color = "grey5", size = 1),
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
