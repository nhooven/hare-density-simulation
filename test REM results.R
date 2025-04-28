# sanity check - did it work?

ggplot(data = boot.rem.test) +
  
  theme_bw() +
  
  facet_wrap(~ scenario) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_errorbarh(aes(xmin = l95.REM.D,
                     xmax = u95.REM.D,
                     y = true.D),
                 height = 0) +
  
  geom_point(aes(x = mean.REM.D,
                 y = true.D),
             alpha = 0.15) +
  
  coord_cartesian(xlim = c(0, 7),
                  ylim = c(0, 7))

# calculate bias
boot.rem.test.1 <- boot.rem.test %>%
  
  mutate(perc.bias = ((mean.REM.D - true.D) / true.D * 100),
         abs.perc.bias = abs((mean.REM.D - true.D) / true.D * 100)) %>%
  
  # rename factor levels
  mutate(scenario.fac = factor(scenario,
                               labels = c("0.5 hr, 100% , 4 wk",
                                          "0.5 hr, 100% , 2 wk",
                                          "0.5 hr, 60% , 4 wk",
                                          "0.5 hr, 60% , 2 wk",
                                          "1 hr, 100% , 4 wk",
                                          "1 hr, 100% , 2 wk",
                                          "1 hr, 60% , 4 wk",
                                          "1 hr, 60% , 2 wk",
                                          "4 hr, 100% , 4 wk",
                                          "4 hr, 100% , 2 wk",
                                          "4 hr, 60% , 4 wk",
                                          "4 hr, 60% , 2 wk")))


# plot bias as a histogram
ggplot(data = boot.rem.test.1,
       aes(x = perc.bias)) +
  
  theme_bw() +
  
  facet_wrap(~ scenario.fac) +
  
  geom_vline(xintercept = 0) +
  
  geom_histogram(binwidth = 5,
                 fill = "white",
                 color = "black") +
  
  theme(panel.grid = element_blank()) +
  
  scale_x_continuous(breaks = c(-50, 0, 50, 100, 150)) +
  
  xlab("Percent bias")

# plot CVs
ggplot(data = boot.rem.test.1,
       aes(x = cv.REM.D)) +
  
  theme_bw() +
  
  facet_wrap(~ scenario.fac) +
  
  geom_histogram(binwidth = 0.05,
                 fill = "white",
                 color = "black") +
  
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5)) +
  
  xlab("Coefficient of variation")

# summary of percent bias and CV
boot.rem.test.summary <- boot.rem.test.1 %>% 
  
  group_by(scenario, abund) %>%
  
  summarize(mean.abs.bias = mean(abs.perc.bias, na.rm = T),
            sd.abs.bias = sd(abs.perc.bias, na.rm = T),
            mean.cv = mean(cv.REM.D, na.rm = T),
            sd.cv = sd(cv.REM.D, na.rm = T)) %>%
  
  ungroup() %>%
  
  mutate(scenario.fac = factor(scenario,
                               labels = c("0.5 hr, 100% , 4 wk",
                                          "0.5 hr, 100% , 2 wk",
                                          "0.5 hr, 60% , 4 wk",
                                          "0.5 hr, 60% , 2 wk",
                                          "1 hr, 100% , 4 wk",
                                          "1 hr, 100% , 2 wk",
                                          "1 hr, 60% , 4 wk",
                                          "1 hr, 60% , 2 wk",
                                          "4 hr, 100% , 4 wk",
                                          "4 hr, 100% , 2 wk",
                                          "4 hr, 60% , 4 wk",
                                          "4 hr, 60% , 2 wk")),
         density = factor(abund,
                          labels = c("0.4", "0.8", "1.6", "3.2")))

# column plots
ggplot(data = boot.rem.test.summary,
       aes(y = scenario.fac,
           x = mean.abs.bias,
           fill = density)) +
  
  theme_bw() +
  
  geom_col(position = "dodge",
           color = "black") +
  
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title.y = element_blank()) +
  
  xlab("Mean absolute percent bias") +
  
  scale_fill_viridis_d(option = "magma",
                       begin = 0.25)

# CV
ggplot(data = boot.rem.test.summary,
       aes(y = scenario.fac,
           x = mean.cv,
           fill = density)) +
  
  theme_bw() +
  
  geom_col(position = "dodge",
           color = "black") +
  
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title.y = element_blank()) +
  
  xlab("Mean coefficient of variation") +
  
  scale_fill_viridis_d(option = "magma",
                       begin = 0.25)
