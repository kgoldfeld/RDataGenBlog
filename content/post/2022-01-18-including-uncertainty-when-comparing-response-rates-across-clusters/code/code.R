library(ggplot2)
library(simstudy)
library(data.table)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggdist)

def <- defData(varname = "n", formula = 80, variance = .3, dist="negBinomial")
def <- defData(def, varname = "y", formula = 1.35, 
               variance = "n", link = "logit", dist = "binomial")
def <- defData(def, varname = "p", formula = "y/n", dist = "nonrandom")

set.seed(4601)
dd <- genData(30, def, id = "site")

ggplot(data = dd, aes(x = site, y = p)) +
  geom_col(aes(fill = factor((p< 0.75)))) +
  geom_hline(yintercept = 0.75, color = "white", size = .5) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "response rate") +
  scale_fill_manual(values = c("grey70", "#e78787"))

#--- Simple model with no predictors - just intercept

fit1 <- glm(cbind(y, n - y) ~ 1, 
            data = dd, family = "binomial")

summary(fit1)

# compare estimated probability with observed probability

lOR <- coef(summary(fit1))[1]
1/(1+exp(-lOR)) # convert log odds to probability
with(dd, sum(y)/sum(n)) # calculate prob from data

#--- Model with site-specific intercepts

fit2 <- glm(cbind(y, n - y) ~ factor(site) - 1, 
            data = dd, family = "binomial")

summary(fit2)

# create data.table for plotting

sites <- rownames(coef(summary(fit2)))
p_est <- 1/(1 + exp(-coef(summary(fit2))[,"Estimate"]))

ci <- data.table(confint(fit2)) # ci on log-odds scale
setnames(ci, c("l95_lo", "u95_lo"))
ci[, `:=`(l95_p = 1/(1+exp(-l95_lo)), u95_p = 1/(1+exp(-u95_lo)))] # convert to prob

dp <- data.table(sites, p_est, ci[, .(l95_p, u95_p)])

setkey(dp, p_est) # print order is by probability
dp[, index := .I]

ggplot(data = dp, aes(x = p_est, y = index)) +
  geom_vline(xintercept = 0.75, color = "white", size = .8) +
  geom_point(size = .5) +
  geom_segment(aes(x = l95_p, xend = u95_p, yend = index),
               color = "grey40", size = .2) +
  theme(panel.grid = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 11, face = "bold")) +
  scale_x_continuous(limits = c(0.2,1), breaks = seq(0.2,1, by = 0.2),
                     name = "response rate") +
  ylab("site")

####

mod <- cmdstan_model("code/binom.stan")

data_list <- list(S = nrow(dd), y = dd$y, n = dd$n)

fit <- mod$sample(
  data = data_list,
  refresh = 500,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 500,
  iter_sampling = 2500,
  show_messages = FALSE
)

fit$summary(c("alpha", "beta", "mu", "s2", "disp"))

post_array <- fit$draws()
# mcmc_trace(post_array, pars = c("alpha", "beta", "mu","s2", "disp"), 
#            facet_args = list(nrow = 3))

df <- data.frame(as_draws_rvars(fit$draws(variables = "theta")))
df$index <- rank(median(df$theta))

ggplot(data = df, aes(dist = theta, y = index)) +
  geom_vline(xintercept = 0.75, color = "white", size = .8) +
  stat_dist_pointinterval(.width = c(.95), 
                          interval_color = "grey80",
                          interval_size = 1,
                          point_size = .2) +
  theme(panel.grid = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 11, face = "bold")) +
  scale_x_continuous(limits = c(0.2,1), breaks = seq(0.2,1, by = 0.2),
                     name = "response rate") +
  ylab("site")
