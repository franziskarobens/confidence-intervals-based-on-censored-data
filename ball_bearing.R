

ball_bearing_data <- c(
  17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 
  48.40, 51.84, 51.96, 54.12, 55.56, 67.80, 
  68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 
  105.12, 105.84, 127.92, 128.04, 173.40
)

delta_40 <- ball_bearing_data <= 40
sum(delta_40)

delta_60 <- ball_bearing_data <= 60
sum(delta_60)



mle_sev(t = ball_bearing_data, t_c = 60)
mle_sev(t = ball_bearing_data, t_c = 40)


prob_plot <- function() {
  x <- seq(from = 10, to = 200, by = 1)
  y <- sapply(x, function(val) {
    mean(ball_bearing_data <= val)
  })
  
  plot(x, y, 
       type = "p",
       log = "x",
       xlim = c(5, 200),
       xlab = "Ball Bearing Life",
       ylab = "Proportion failing",
       main = "Probability Plot (Log X-Axis)")
  
  abline(v = 40, lty = 2)
  abline(v = 60, lty = 2)
}

prob_plot()
