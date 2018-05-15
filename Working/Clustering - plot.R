###

dv <- data.table()

for (i in seq(6, 60, by = 6)) {
  for (j in seq(0, 1, by = .1)) {
    v <- var(seq(-j, j, length = i )) 
    
    dv <- rbind(dv, data.table(i=i, j=j, v=v))
  }
  
}

dv[, s := sqrt(v)]

dcast(dv, i ~ j, value.var = "v")
dcast(dv, i ~ j, value.var = "s")

sds <- seq(0,.42, by = .06)
vars <- sds^2
vars
ICC <- vars/(1+vars)

###


probs <- cbind(probf, probc)
setnames(probs, c("j", "pf","pfre", "var", "pc", "pcre"))
probs[, ICC := round(ICC,3)]

pmelt <- melt(probs, id.vars = "ICC", measure.vars = c("pf", "pc"))

ggplot(data = pmelt, aes(x=ICC, y=value, group = variable)) +
  geom_point(aes(color=variable)) +
  geom_line(aes(color=variable)) +
  scale_y_continuous(name = "P(Type I error)", 
                     breaks = c(.05, .25, .50, .75, 1), 
                     limits = c(0,1)
  ) +
  xlab("Estimated ICC") +
  ggtitle(paste(nsites, "sites")) +
  scale_color_manual(values = c("#80a58c", "#a58099"), labels = c("`Fixed`","`Random`")) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(.8,.25)) 


