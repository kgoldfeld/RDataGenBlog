probreject[, ICC := cvar/(1+cvar)]

# prob6 <- probreject
# prob12 <- probreject
# prob24 <- probreject
# prob6re <- probreject
# prob24re <- probreject
# prob6cr <- probreject

probreject <- prob24re
nsites <- 24

p24f <- ggplot(data = probreject, aes(x=ICC, y=preject)) +
  geom_jitter(width=.0025, height = .01, size = 1, color = "#6e9009") +
  scale_y_continuous(name = "P(Type I error)", 
                     breaks = c(.05, .25, .50, .75, 1), 
                     limits = c(-0.025,1.025)
  ) +
  scale_x_continuous(name = "ICC",
                     breaks = ICCs) +
  ggtitle(paste("Fixed effect:", nsites, "sites")) +
  # scale_color_manual(values = c("#80a58c", "#a58099"), labels = c("`Fixed`","`Random`")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(.8,.25)) 

p24r <- ggplot(data = probreject, aes(x=ICC, y=prejectre)) +
  geom_jitter(width=.0025, height = .01, size = 1, color = "#2b0990") +
  scale_y_continuous(name = "P(Type I error)", 
                     breaks = c(.05, .25, .50, .75, 1), 
                     limits = c(-0.025,1.025)
  ) +
  scale_x_continuous(name = "ICC",
                     breaks = ICCs) +
  ggtitle(paste("Random effect:", nsites, "sites")) +
  # scale_color_manual(values = c("#80a58c", "#a58099"), labels = c("`Fixed`","`Random`")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(.8,.25)) 

p6cr <- ggplot(data = probreject, aes(x=ICC, y=preject)) +
  geom_jitter(width=.0025, height = .01, size = 1, color = "#902b09") +
  scale_y_continuous(name = "P(Type I error)", 
                     breaks = c(.05, .25, .50, .75, 1), 
                     limits = c(-0.025,1.025)
  ) +
  scale_x_continuous(name = "ICC",
                     breaks = ICCs) +
  ggtitle(paste("Randomized within cluster:", nsites, "sites")) +
  # scale_color_manual(values = c("#80a58c", "#a58099"), labels = c("`Fixed`","`Random`")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(.8,.25)) 


png("/Users/goldfk01/RDataGenBlog/static/img/post-smallcluster/RwithinC6.png", 
    width = 750, height = 400)
p6cr
dev.off()
