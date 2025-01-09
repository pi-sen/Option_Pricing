install.packages("ggplot2")
install.packages("knitr")

# Load required libraries
library(ggplot2)
library(knitr)

# Parameters and Simulation Setup
S0 <- 100  # Initial asset price
r <- 0.04  # Risk-free interest rate (per year)
q <- 0.02  # Dividend yield
sigma <- 0.2  # Volatility (per year)
K <- 100  # Strike price
T <- 0.5  # Maturity (in years)

# Black-Scholes Exact Price
black_scholes_put <- function(S0, K, T, r, q, sigma) {
  d1 <- (log(S0 / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  put_price <- K * exp(-r * T) * pnorm(-d2) - S0 * exp(-q * T) * pnorm(-d1)
  return(put_price)
}

exact_price <- black_scholes_put(S0, K, T, r, q, sigma)
print(paste("Exact Black-Scholes Price:", round(exact_price, 4)))

# Monte Carlo Simulation
monte_carlo_put <- function(S0, K, T, r, q, sigma, N) {
  set.seed(123)
  Z <- rnorm(N)
  ST <- S0 * exp((r - q - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
  payoff <- pmax(K - ST, 0)
  option_price <- exp(-r * T) * mean(payoff)
  standard_error <- sd(payoff) / sqrt(N)
  CI <- c(option_price - 1.96 * standard_error, option_price + 1.96 * standard_error)
  return(list(price = option_price, se = standard_error, CI = CI))
}

sample_sizes <- c(100, 1000, 5000, 10000, 50000, 100000)
prices <- numeric(length(sample_sizes))
standard_errors <- numeric(length(sample_sizes))

for (i in 1:length(sample_sizes)) {
  N <- sample_sizes[i]
  result <- monte_carlo_put(S0, K, T, r, q, sigma, N)
  prices[i] <- result$price
  standard_errors[i] <- result$se
}

plot(sample_sizes, prices, type = "b", col = "blue", xlab = "Sample Size", ylab = "Option Price Estimate",
     main = "Convergence of Monte Carlo Estimate")
grid()

# Monte Carlo with Antithetic Variates
monte_carlo_antithetic <- function(S0, K, T, r, q, sigma, N) {
  set.seed(123)
  Z <- rnorm(N / 2)
  ST1 <- S0 * exp((r - q - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
  ST2 <- S0 * exp((r - q - 0.5 * sigma^2) * T - sigma * sqrt(T) * Z)
  payoff1 <- pmax(K - ST1, 0)
  payoff2 <- pmax(K - ST2, 0)
  payoff <- (payoff1 + payoff2) / 2
  option_price <- exp(-r * T) * mean(payoff)
  standard_error <- sd(payoff) / sqrt(N)
  CI <- c(option_price - 1.96 * standard_error, option_price + 1.96 * standard_error)
  return(list(price = option_price, se = standard_error, CI = CI))
}

result_antithetic <- monte_carlo_antithetic(S0, K, T, r, q, sigma, 100000)
print(result_antithetic)

# Antithetic Variates Effectiveness
methods <- c("Standard Monte Carlo", "Antithetic Variates")
standard_errors <- c(0.0232, 0.0116)

comparison_data <- data.frame(Method = methods, StandardError = standard_errors)

ggplot(comparison_data, aes(x = Method, y = StandardError, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Reduction in Standard Error with Antithetic Variates",
    x = "Method",
    y = "Standard Error"
  ) +
  theme_minimal() +
  geom_text(aes(label = round(StandardError, 4)), vjust = -0.5, size = 5)

# Final Comparison Table
results <- data.frame(SampleSize = sample_sizes,
                      OptionPrice_Std = NA, SE_Std = NA, Time_Std = NA,
                      OptionPrice_Antithetic = NA, SE_Antithetic = NA, Time_Antithetic = NA)

for (i in 1:length(sample_sizes)) {
  N <- sample_sizes[i]
  
  # Standard Monte Carlo
  time_std <- system.time({
    res_std <- monte_carlo_put(S0, K, T, r, q, sigma, N)
  })
  results$OptionPrice_Std[i] <- res_std$price
  results$SE_Std[i] <- res_std$se
  results$Time_Std[i] <- time_std["elapsed"]
  
  # Antithetic Variates
  time_antithetic <- system.time({
    res_antithetic <- monte_carlo_antithetic(S0, K, T, r, q, sigma, N)
  })
  results$OptionPrice_Antithetic[i] <- res_antithetic$price
  results$SE_Antithetic[i] <- res_antithetic$se
  results$Time_Antithetic[i] <- time_antithetic["elapsed"]
}

print(results)
