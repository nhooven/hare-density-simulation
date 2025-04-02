# EXTRA

# How similar are Gaussian speeds vs. conditional speeds?
ggplot(data = all.speeds %>% mutate(scenario = factor(scenario,
                                                      labels = c("1 hr, 100%, 4 wk",
                                                                 "2 hr, 100%, 4 wk",
                                                                 "4 hr, 100%, 4 wk",
                                                                 "1 hr, 80%, 4 wk",
                                                                 "2 hr, 80%, 4 wk",
                                                                 "4 hr, 80%, 4 wk",
                                                                 "1 hr, 100%, 2 wk",
                                                                 "2 hr, 100%, 2 wk",
                                                                 "4 hr, 100%, 2 wk",
                                                                 "1 hr, 80%, 2 wk",
                                                                 "2 hr, 80%, 2 wk",
                                                                 "4 hr, 80%, 2 wk"))),
       aes(x = speed.mean.est,
           y = speed.cond.est)) +
  
  facet_wrap(~ scenario,
             nrow = 4) +
  
  theme_bw() +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_errorbar(aes(y = speed.cond.est,
                    ymin = speed.cond.lo,
                    ymax = speed.cond.hi),
                color = "gray") +
  
  geom_errorbarh(aes(y = speed.cond.est,
                     xmin = speed.mean.lo,
                     xmax = speed.mean.hi),
                 color = "gray") +
  
  geom_point(color = "black",
             fill = "#CC33FF",
             size = 1.5,
             shape = 21) +
  
  xlab("Mean speed (m/s) from Gaussian movement process") +
  ylab("Mean simulated speed (m/s)") +
  
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5))
