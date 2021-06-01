library(ggplot2)
library(paletteer)
library(ggpubr)

load("data/bp.rda")

pdata <- final_tab[t_sigma == 5 & mu.lor == -1]

p1 <- ggplot(data = pdata, aes(x = nobs, y = p_95)) +
  geom_hline(yintercept = c(0.80), color = "white", lty = 1) +
  geom_line(size = 2, color = "#2561dd") +
  scale_y_continuous(limits = c(0, 1), breaks = c(.2, .4, .6, .8, 1), name = "proportion") +
  xlab("total sample size") +
  # facet_grid(sigma.lor ~ t_sigma) +
  # scale_color_paletteer_d(name = "log(OR)", palette="wesanderson::Moonrise2") +
  ggtitle("proportion of posterior distributions with P(log(OR) < 0) > 0.95") +
  # facet_grid(~t_sigma) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12))

ggsave(file = "img/power_curve.png",
       p1, width = 4, height = 2, scale=1.5)



