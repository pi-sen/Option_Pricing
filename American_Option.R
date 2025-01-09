# Required Libraries
library(ggplot2)

# Parameters
S0 <- 100     # Initial stock price
K <- 100      # Strike price
r <- 0.05     # Risk-free rate
sigma <- 0.2  # Volatility
T <- 1        # Time to maturity
N <- 10000    # Number of paths for LSM
M <- 50       # Number of time steps for LSM
N1 <- 500     # Steps for Binomial Tree
N2 <- 1000    # Steps for Richardson Extrapolation in BBSR
grid_points <- 100  # FDM grid resolution

# LSM Pre-Optimized
lsm_pre_optimized <- function(S0, K, r, sigma, T, N, M) {
  dt <- T / M
  discount <- exp(-r * dt)
  set.seed(123)
  paths <- matrix(0, nrow = N, ncol = M + 1)
  paths[, 1] <- S0
  for (i in 2:(M + 1)) {
    z <- rnorm(N)
    paths[, i] <- paths[, i - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * z)
  }
  cashflows <- matrix(0, nrow = N, ncol = M + 1)
  cashflows[, M + 1] <- pmax(K - paths[, M + 1], 0)
  for (t in M:1) {
    itm <- which(paths[, t] < K)
    if (length(itm) > 0) {
      X <- cbind(1, paths[itm, t], paths[itm, t]^2)
      Y <- discount * cashflows[itm, t + 1]
      beta <- solve(t(X) %*% X, t(X) %*% Y)
      continuation_value <- X %*% beta
      cashflows[itm, t] <- pmax(K - paths[itm, t], continuation_value)
    }
  }
  price <- mean(cashflows[, 2] * discount)
  return(price)
}

# LSM Post-Optimized
lsm_post_optimized <- function(S0, K, r, sigma, T, N, M) {
  dt <- T / M
  discount <- exp(-r * dt)
  set.seed(123)
  paths <- matrix(S0 * exp((r - 0.5 * sigma^2) * (0:M) * dt + sigma * sqrt(dt) * cbind(0, matrix(rnorm(N * M), N))), nrow = N)
  cashflows <- matrix(pmax(K - paths[, M + 1], 0), nrow = N)
  for (t in M:1) {
    itm <- paths[, t] < K
    X <- cbind(1, paths[itm, t], paths[itm, t]^2)
    Y <- discount * cashflows[itm]
    if (sum(itm) > 0) {
      beta <- solve(t(X) %*% X, t(X) %*% Y)
      continuation_value <- X %*% beta
      cashflows[itm] <- pmax(K - paths[itm, t], continuation_value)
    }
  }
  price <- mean(cashflows * discount)
  return(price)
}

# Binomial Tree
binomial_tree <- function(S0, K, r, sigma, T, N) {
  dt <- T / N
  u <- exp(sigma * sqrt(dt))
  d <- 1 / u
  p <- (exp(r * dt) - d) / (u - d)
  stock_prices <- matrix(0, nrow = N + 1, ncol = N + 1)
  stock_prices[1, 1] <- S0
  for (i in 2:(N + 1)) {
    stock_prices[1:i, i] <- S0 * u^(0:(i - 1)) * d^((i - 1):0)
  }
  option_values <- matrix(0, nrow = N + 1, ncol = N + 1)
  option_values[, N + 1] <- pmax(K - stock_prices[, N + 1], 0)
  for (i in N:1) {
    option_values[1:i, i] <- exp(-r * dt) * (p * option_values[1:i, i + 1] + (1 - p) * option_values[2:(i + 1), i + 1])
  }
  return(option_values[1, 1])
}

# BBSR Method
bbsr <- function(S0, K, r, sigma, T, N1, N2) {
  price1 <- binomial_tree(S0, K, r, sigma, T, N1)
  price2 <- binomial_tree(S0, K, r, sigma, T, N2)
  extrapolated <- 2 * price2 - price1
  return(list(price1 = price1, price2 = price2, extrapolated = extrapolated))
}

fdm_implementation <- function(S0, K, r, sigma, T, M, grid_points) {
  dt <- T / M
  ds <- 2 * S0 / grid_points
  S <- seq(0, 2 * S0, by = ds)
  V <- matrix(0, nrow = grid_points + 1, ncol = M + 1)
  
  # Boundary Conditions
  V[, M + 1] <- pmax(K - S, 0)
  V[1, ] <- K * exp(-r * seq(0, T, by = dt))
  V[grid_points + 1, ] <- 0
  
  # Coefficients for Implicit Method
  i <- 1:grid_points
  alpha <- -0.5 * dt * (sigma^2 * i^2 - r * i)
  beta <- 1 + dt * (sigma^2 * i^2 + r)
  gamma <- -0.5 * dt * (sigma^2 * i^2 + r * i)
  
  # Construct Tridiagonal Matrix
  A <- diag(beta[2:grid_points])
  for(i in 1:(grid_points-2)) {
    A[i, i+1] <- gamma[i+1]
    A[i+1, i] <- alpha[i+2]
  }
  
  # Solve Backward in Time
  for (j in M:1) {
    rhs <- V[2:grid_points, j + 1]
    rhs[1] <- rhs[1] - alpha[2] * V[1, j]
    rhs[grid_points - 1] <- rhs[grid_points - 1] - gamma[grid_points] * V[grid_points + 1, j]
    V[2:grid_points, j] <- solve(A, rhs)
    V[, j] <- pmax(V[, j], K - S)  # Apply early exercise condition
  }
  
  # Interpolate to Get Price at S0
  price <- approx(S, V[, 1], xout = S0)$y
  return(price)
}


# Price Calculations and Timings
start_time <- Sys.time()
price_lsm_pre <- lsm_pre_optimized(S0, K, r, sigma, T, N, M)
time_lsm_pre <- Sys.time() - start_time

start_time <- Sys.time()
price_lsm_post <- lsm_post_optimized(S0, K, r, sigma, T, N, M)
time_lsm_post <- Sys.time() - start_time

start_time <- Sys.time()
price_binomial <- binomial_tree(S0, K, r, sigma, T, N1)
time_binomial <- Sys.time() - start_time

start_time <- Sys.time()
price_bbsr <- bbsr(S0, K, r, sigma, T, N1, N2)$extrapolated
time_bbsr <- Sys.time() - start_time

start_time <- Sys.time()
price_fdm <- fdm_implementation(S0, K, r, sigma, T, M, grid_points)
time_fdm <- Sys.time() - start_time

# Results DataFrame
results <- data.frame(
  Method = c("LSM Pre-Optimized", "LSM Post-Optimized", "Binomial Tree", "BBSR", "FDM"),
  Price = c(price_lsm_pre, price_lsm_post, price_binomial, price_bbsr, price_fdm),
  Time = c(time_lsm_pre, time_lsm_post, time_binomial, time_bbsr, time_fdm)
)

print(results)

# Visualization
ggplot(results, aes(x = Method, y = Price, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Option Pricing Comparison", x = "Method", y = "Option Price") +
  theme_minimal() +
  geom_text(aes(label = round(Price, 2)), vjust = -0.5, size = 3)

ggplot(results, aes(x = Method, y = Time, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Computation Time Comparison", x = "Method", y = "Time (seconds)") +
  theme_minimal() +
  geom_text(aes(label = round(Time, 3)), vjust = -0.5, size = 3)

#Convergence using FDM

# Convergence Analysis for FDM
fdm_convergence <- function(S0, K, r, sigma, T, M_values, grid_points) {
  prices <- numeric(length(M_values))
  for (i in seq_along(M_values)) {
    prices[i] <- fdm_implementation(S0, K, r, sigma, T, M_values[i], grid_points)
  }
  return(prices)
}

# Define M values for testing
M_values <- seq(10, 200, by = 10)

# Calculate prices
prices_fdm <- fdm_convergence(S0, K, r, sigma, T, M_values, grid_points)

# Convergence Plot
library(ggplot2)
data <- data.frame(M_values = M_values, prices = prices_fdm)
ggplot(data, aes(x = M_values, y = prices)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(title = "Convergence of FDM Option Price",
       x = "Number of Time Steps (M)",
       y = "Option Price") +
  theme_minimal()

#other convergence plots

# Parameters for convergence analysis
M_values <- seq(10, 200, by = 10)  # Number of time steps for methods
N_values <- seq(1000, 20000, by = 2000)  # Number of paths for LSM

# Convergence for LSM (Pre- and Post-Optimization)
lsm_pre_convergence <- sapply(N_values, function(N) {
  lsm_pre_optimized(S0, K, r, sigma, T, N, M)
})

lsm_post_convergence <- sapply(N_values, function(N) {
  lsm_post_optimized(S0, K, r, sigma, T, N, M)
})

# Convergence for Binomial Tree and BBSR
binomial_convergence <- sapply(M_values, function(M) {
  binomial_tree(S0, K, r, sigma, T, M)  # Binomial tree should use M values
})

# BBSR convergence uses N1 and N2, let's keep it simple for now with M for N1
bbsr_convergence <- sapply(M_values, function(M) {
  bbsr(S0, K, r, sigma, T, M, 10000)$extrapolated  # Only use the extrapolated price
})

# Data Preparation for Plotting
convergence_data <- data.frame(
  N = rep(N_values, 2),
  M = rep(M_values, 2),
  Method = c(rep("LSM (Pre-Optimized)", length(N_values)),
             rep("LSM (Post-Optimized)", length(N_values))),
  Price = c(lsm_pre_convergence, lsm_post_convergence)
)

# Data Preparation for Plotting
binomial_data <- data.frame(
  M = rep(M_values, 2),
  Method = c(rep("Binomial Tree", length(M_values)),
             rep("BBSR", length(M_values))),
  Price = c(binomial_convergence, bbsr_convergence)
)

# Plot Convergence for LSM
ggplot(convergence_data, aes(x = N, y = Price, color = Method)) +
  geom_line(size = 1) +
  labs(title = "Convergence of LSM Method",
       x = "Number of Paths (N)",
       y = "Option Price") +
  theme_minimal()

# Plot Convergence for Binomial Tree and BBSR
ggplot(binomial_data, aes(x = M, y = Price, color = Method)) +
  geom_line(size = 1) +
  labs(title = "Convergence of Binomial Tree and BBSR",
       x = "Number of Time Steps (M)",
       y = "Option Price") +
  theme_minimal()