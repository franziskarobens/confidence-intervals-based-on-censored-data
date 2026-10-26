# replicate Figure 1 from Jeng 

library(ggplot2)

ggplot(data.frame(t = ball_bearing_data), aes(t)) +
  stat_ecdf(geom = "point") +
  geom_vline(xintercept = c(40, 60), linetype = "dashed", col = "red") +
  scale_x_log10() +
  coord_cartesian(xlim = c(5, 200), ylim = c(0, 1)) +
  labs(x = "t (log scale)", y = "Proportion failing")

