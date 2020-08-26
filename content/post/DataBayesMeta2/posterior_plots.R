library(data.table)
library(ggplot2)

load("~/Local R Projects/RDataGenBlog/Working/divergence_c.rda")

p <- merge(parsamp, diverge[n_divergent >= 50], by = "iternum")

p50 <- ggplot(data = p, aes( x = sample, group = iternum)) +
  geom_density(fill = "blue", alpha = .2) +
  theme(panel.grid = element_blank()) +
  ylim(0, 1.6) +
  xlim(-1.5, 2.5) +
  ggtitle("50+ divergent transitions")

p <- merge(parsamp, diverge[n_divergent < 50], by = "iternum")

p49 <- ggplot(data = p, aes( x = sample, group = iternum)) +
  geom_density(fill = "blue", alpha = .2) +
  theme(panel.grid = element_blank()) +
  ylim(0, 1.6) +
  xlim(-1.5, 2.5) +
  ggtitle("< 50 divergent transitions")


post_plot_c <- ggplot(data = parsamp, aes( x = sample, group = iternum)) +
  geom_density(fill = "blue", alpha = .2) +
  theme(panel.grid = element_blank()) +
  ylim(0, 1.6) +
  xlim(-1.5, 2.5) +
  ggtitle("Original parameterization")


load("~/Local R Projects/RDataGenBlog/Working/divergence_nc.rda")

post_plot_nc <- ggplot(data = parsamp, aes( x = sample, group = iternum)) +
  geom_density(fill = "blue", alpha = .2) +
  theme(panel.grid = element_blank()) +
  ylim(0, 1.6) +
  xlim(-1.5, 2.5) +
  ggtitle("Non-centered parameterization")

post_plot <- gridExtra::grid.arrange(post_plot_c, post_plot_nc, nrow = 1)
post_div <- gridExtra::grid.arrange(p49, p50, nrow = 1)

ggsave(post_plot, filename = "static/img/post-bayesdiag/post_plot.png", 
       width = 6, height = 3, scale = 1.3)

ggsave(post_div, filename = "static/img/post-bayesdiag/post_div.png", 
       width = 6, height = 3, scale = 1.3)
