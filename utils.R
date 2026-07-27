
get_delta <- function(data, t_c) {
  # delta = 0 if data is censored
  return(ifelse(data <= t_c, 1, 0))
}

