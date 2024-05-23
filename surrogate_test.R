library(tidyverse)
library(fpp3)

# Function to perform surrogate test for independence
surrogate.test = function(data, lag, N = 1000, test.stat = "ljung-box") {
  
  # data: a tsibble or numeric vector object
  # lag: number of lags in portmanteau test statistic
  # N: number of permutations to perform
  # test.stat: either "ljung-box" or "box-pierce"
  
  if (is_tsibble(data)) {
    if (length(measures(data)) != 1) {  
      stop("data must be a tsibble with one measurement variable")
    }
    # Extract time series 
    data = data %>% 
      pull(as.character(measures(data)[[1]])) 
  }

  n = length(data)

  Q.null = rep(NA, N)  # Open test statistic vectors
  
  if (test.stat == "ljung-box") {
    
    # Observed test statistic
    r.obs = acf(data, plot = FALSE)$acf[2:(lag + 1)]
    Q.obs = n * (n + 2) * sum(r.obs ^ 2 / (n - 1:lag))
    
    # Null distribution
    for (i in 1:N) {
      surrogate = sample(data, n)  # Permute data (kill autocorrelation, maintain amplitude)
      r = acf(surrogate, plot = FALSE)$acf[2:(lag + 1)]   # Estimate autocorrelation
      Q.null[i] = n * (n + 2) * sum(r ^ 2 / (n - 1:lag))  # Ljung-Box test statistic
    }
    
  }
  
  if (test.stat == "box-pierce") {
    
    # Observed test statistic
    r.obs = acf(data, plot = FALSE)$acf[2:(lag + 1)]
    Q.obs = n * sum(r.obs ^ 2)
    
    # Null distribution
    for (i in 1:N) {
      surrogate = sample(data, n)  # Permute data (kill autocorrelation, maintain amplitude)
      r = acf(surrogate, plot = FALSE)$acf[2:(lag + 1)]  # Estimate autocorrelation
      Q.null[i] = n * sum(r ^ 2)                         # Box-Pierce test statistic
    }
    
  }
  
  # Compute p-value
  p.value = mean(Q.null >= Q.obs)  # p-value

  # Output
  output = list(Q.null = Q.null,
                Q.obs = Q.obs,
                test.stat = test.stat,
                p.value = p.value)
  
  class(output) = "surrogate"
  
  return(output)
  
}


# Function to plot surrogate null distribution and observed test statistic
# Requires ggplot2
plot.surrogate = function(obj) {
  
  # obj: Object of class "surrogate"

  ggplot(data = data.frame(Q = obj$Q.null),
         mapping = aes(x = Q)) +
    geom_histogram(fill = "navy", colour = "black") +
    geom_vline(xintercept = obj$Q.obs,
               linetype = "dashed") +
    labs(x = "Test statistic",
         y = "Count") 
  
}


# Example 1: Data set is a numeric vector
#            Simulated AR(1): Should not be independent
set.seed(1)
data = arima.sim(n = 256, 
                 model = list(ar = 0.35))  # Simulate data
s = surrogate.test(data, lag = 10)         # Run surrogate test
s$p.value                                  # Extract p-value
s %>% 
  plot.surrogate() + 
  theme_bw() +
  labs(title = "Surrogate data test with Ljung-Box statistic")
  
# Example 2: Data set is a tsibble object
#            Simulated white noise: Should be independent
data = tsibble(Value = rnorm(256),
               Time = 1:256,
               index = Time)
s = surrogate.test(data, lag = 8, N = 5000, 
                   test.stat = "box-pierce")  # Run surrogate test
s$p.value                                     # Extract p-value
s %>% 
  plot.surrogate() + 
  theme_bw() +
  labs(title = "Surrogate data test with Box-Pierce statistic")
