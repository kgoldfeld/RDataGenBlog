x = c(- 2.6, -1.7, -.5, .9)
dx <- data.table(i = c(1:4),
        y = c("inconclusive", "equivalent", "non-inferior", "superior"),
        x = x, xend = x + 3, est = x + 1.5)

ggplot(data=dx) + 
  geom_vline(xintercept = c(-2, 0, 2), color = c("grey70", "white", "grey70")) +
  geom_segment(aes(x = x, xend = xend, y = i, yend = i), size = .8) +
  geom_point(aes(x=est, y = i), size = 2) +
  geom_point(aes(x=x, y = i), shape = 22, fill = "red") +
  geom_point(aes(x=xend, y = i), shape = 22, fill = "red") +
  geom_text(aes(x = est, y = i + .15, label = y)) +
  scale_x_continuous(limits = c(-4, 4), breaks = c(-2,0,2), 
                     labels = c(expression(-Delta), 0, expression(Delta))) +
  scale_y_continuous(expand=c(0.2,0)) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) 

x = c(- 2.6, -1.7, -.5, .9)
dx <- data.table(i = c(1:4),
                 y = c("inconclusive", "inconclusive", "inconclusive", "superior"),
                 x = x, xend = x + 3, est = x + 1.5)

ggplot(data=dx) + 
  geom_vline(xintercept = c(-.35, 0, .35), color = c("grey70", "white", "grey70")) +
  geom_segment(aes(x = x, xend = xend, y = i, yend = i), size = .8) +
  geom_point(aes(x=est, y = i), size = 2) +
  geom_point(aes(x=x, y = i), shape = 22, fill = "red") +
  geom_point(aes(x=xend, y = i), shape = 22, fill = "red") +
  geom_text(aes(x = est, y = i + .15, label = y)) +
  scale_x_continuous(limits = c(-4, 4), breaks = c(-.35,0,.35), 
                     labels = c(expression(-Delta), 0, expression(Delta))) +
  scale_y_continuous(expand=c(0.2,0)) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) 
  
