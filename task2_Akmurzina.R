#install.packages("quantmod")

library(quantmod)

# Загружаем акции BYD за 4 года
byd <- getSymbols(
  "002594.SZ",
  src = "yahoo",
  from = "2022-10-19",
  to = "2026-10-20",
  auto.assign = FALSE
)

# Получаем скорректированные цены
prices <- Ad(byd)

# Создаём таблицу
data <- data.frame(
  Date = as.Date(index(prices)),
  AdjClose = as.numeric(prices)
)

# Сохраняем в CSV
write.csv(data, "BYD_stock_prices.csv", row.names = FALSE)

View(data)

getwd()

#Data preparation
nrow(data)
head(data)
tail(data)
colSums(is.na(data))
sum(duplicated(data$Date))
range(data$Date)

#Calculate daily logarithmic returns
data$Return <- c(NA, diff(log(data$AdjClose)))
summary(data$Return)


# Task 2: Exploratory data analysis

returns <- na.omit(data$Return)

plot(
  data$Date, data$AdjClose,
  type = "l",
  main = "BYD Adjusted Stock Prices",
  xlab = "Date",
  ylab = "Price (CNY)"
)

plot(
  data$Date[-1], returns,
  type = "l",
  main = "BYD Daily Logarithmic Returns",
  xlab = "Date",
  ylab = "Log Return"
)

hist(
  returns,
  breaks = 40,
  probability = TRUE,
  main = "Distribution of BYD Returns",
  xlab = "Log Return"
)

curve(
  dnorm(x, mean(returns), sd(returns)),
  add = TRUE,
  col = "red",
  lwd = 2
)


# Task 2: Descriptive statistics

#install.packages("moments")
library(moments)

stats <- data.frame(
  Mean = mean(returns),
  Median = median(returns),
  SD = sd(returns),
  Variance = var(returns),
  Min = min(returns),
  Max = max(returns),
  Skewness = skewness(returns),
  Kurtosis = kurtosis(returns)
)

round(stats, 5)

sum(returns == 0)

data$Date[which.max(data$AdjClose)]; max(data$AdjClose)
data$Date[which.min(data$AdjClose)]; min(data$AdjClose)
data[order(-abs(data$Return)), c("Date","Return")][1:10, ]
data[which.max(data$AdjClose), c("Date","AdjClose")]
data[which.min(data$AdjClose), c("Date","AdjClose")]
# Task 2: Five largest daily losses

data$ReturnPct <- 100 * (exp(data$Return) - 1)

data$LossCNY <- c(
  NA,
  -diff(data$AdjClose)
)

largest_losses <- head(
  data[order(data$ReturnPct), 
       c("Date", "ReturnPct", "LossCNY")],
  5
)

largest_losses$ReturnPct <- round(largest_losses$ReturnPct, 2)
largest_losses$LossCNY <- round(largest_losses$LossCNY, 2)

largest_losses

# 3. Даты максимума и минимума цены
data[which.max(data$AdjClose), c("Date","AdjClose")]
data[which.min(data$AdjClose), c("Date","AdjClose")]

# 4. Топ-5 убытков с тремя знаками (порядок по логдоходности, не по округлению)
top5 <- head(data[order(data$Return), c("Date","Return","ReturnPct","LossCNY")], 5)
top5$ReturnPct <- round(top5$ReturnPct, 3)
top5

# 1. Убыток в тенге (подставьте курсы Нацбанка РК)
fx <- data.frame(
  Date   = as.Date(c("2025-04-07","2025-03-21","2024-10-09","2025-05-26","2022-10-28")),
  CNYKZT = c(69.27, NA, NA, NA, NA)
)
top5 <- merge(top5, fx, by = "Date")
top5$LossKZT <- round(top5$LossCNY * top5$CNYKZT, 2)
top5[order(top5$Return), ]



# Task 3: Probability Distributions and Risk Models

packages <- c("fitdistrplus", "sn")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(fitdistrplus)
library(sn)

if (!"Return" %in% names(data)) {
  data$Return <- c(NA, diff(log(data$AdjClose)))
}

data$Date <- as.Date(data$Date)
returns <- as.numeric(na.omit(data$Return))


# Task 3.1: Daily Returns

normal_fit <- fitdist(returns, "norm")

dstudent <- function(x, mu, sigma, nu, log = FALSE) {
  z <- (x - mu) / sigma
  d <- dt(z, df = nu, log = TRUE) - log(sigma)
  if (log) d else exp(d)
}

student_nll <- function(par) {
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- 2 + exp(par[3])
  
  -sum(dstudent(
    returns, mu, sigma, nu, log = TRUE
  ))
}

student_starts <- c(3, 5, 8, 12, 20)

student_fits <- lapply(student_starts, function(nu0) {
  optim(
    par = c(
      median(returns),
      log(sd(returns)),
      log(nu0 - 2)
    ),
    fn = student_nll,
    method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-12)
  )
})

student_fits <- Filter(
  function(x) x$convergence == 0 && is.finite(x$value),
  student_fits
)

if (length(student_fits) == 0) {
  stop("Student t fitting failed")
}

student_result <- student_fits[[
  which.min(sapply(student_fits, function(x) x$value))
]]

student_params <- c(
  mu = student_result$par[1],
  sigma = exp(student_result$par[2]),
  nu = 2 + exp(student_result$par[3])
)

student_aic <- 2 * 3 + 2 * student_result$value

skewt_nll <- function(par) {
  xi <- par[1]
  omega <- exp(par[2])
  alpha <- par[3]
  nu <- 2 + exp(par[4])
  
  log_density <- tryCatch(
    dst(
      returns,
      xi = xi,
      omega = omega,
      alpha = alpha,
      nu = nu,
      log = TRUE
    ),
    error = function(e) rep(-Inf, length(returns))
  )
  
  if (any(!is.finite(log_density))) {
    return(1e100)
  }
  
  -sum(log_density)
}

skewt_starts <- expand.grid(
  alpha = c(-2, 0, 2),
  nu = c(4, 8, 15)
)

skewt_fits <- lapply(
  seq_len(nrow(skewt_starts)),
  function(i) {
    optim(
      par = c(
        median(returns),
        log(sd(returns)),
        skewt_starts$alpha[i],
        log(skewt_starts$nu[i] - 2)
      ),
      fn = skewt_nll,
      method = "L-BFGS-B",
      lower = c(-0.2, log(1e-5), -20, log(0.01)),
      upper = c(0.2, log(1), 20, log(198)),
      control = list(maxit = 10000)
    )
  }
)

skewt_fits <- Filter(
  function(x) x$convergence == 0 && is.finite(x$value),
  skewt_fits
)

if (length(skewt_fits) == 0) {
  stop("Skewed Student t fitting failed")
}

skewt_result <- skewt_fits[[
  which.min(sapply(skewt_fits, function(x) x$value))
]]

skewt_params <- c(
  xi = skewt_result$par[1],
  omega = exp(skewt_result$par[2]),
  alpha = skewt_result$par[3],
  nu = 2 + exp(skewt_result$par[4])
)

skewt_aic <- 2 * 4 + 2 * skewt_result$value

returns_aic <- data.frame(
  Model = c("Normal", "Student t", "Skewed Student t"),
  AIC = c(normal_fit$aic, student_aic, skewt_aic)
)

returns_parameters <- data.frame(
  Model = c("Normal", "Student t", "Skewed Student t"),
  Location = c(
    normal_fit$estimate["mean"],
    student_params["mu"],
    skewt_params["xi"]
  ),
  Scale = c(
    normal_fit$estimate["sd"],
    student_params["sigma"],
    skewt_params["omega"]
  ),
  Alpha = c(NA, NA, skewt_params["alpha"]),
  Nu = c(NA, student_params["nu"], skewt_params["nu"])
)

rownames(returns_parameters) <- NULL

print(returns_parameters, digits = 6)
print(returns_aic, row.names = FALSE)

hist(
  returns,
  breaks = 40,
  probability = TRUE,
  main = "BYD Daily Returns: Fitted Distributions",
  xlab = "Log Return",
  col = "lightgray",
  border = "white"
)

x <- seq(min(returns), max(returns), length.out = 500)

lines(
  x,
  dnorm(
    x,
    mean = normal_fit$estimate["mean"],
    sd = normal_fit$estimate["sd"]
  ),
  col = "red",
  lwd = 2
)

lines(
  x,
  dstudent(
    x,
    mu = student_params["mu"],
    sigma = student_params["sigma"],
    nu = student_params["nu"]
  ),
  col = "blue",
  lwd = 2
)

lines(
  x,
  dst(
    x,
    xi = skewt_params["xi"],
    omega = skewt_params["omega"],
    alpha = skewt_params["alpha"],
    nu = skewt_params["nu"]
  ),
  col = "darkgreen",
  lwd = 2
)

legend(
  "topright",
  legend = c("Normal", "Student t", "Skewed Student t"),
  col = c("red", "blue", "darkgreen"),
  lwd = 2,
  bty = "n"
)


# Task 3.2: Loss Size Distribution

losses <- -100 * returns[returns < 0]

loss_count <- length(losses)
loss_proportion <- loss_count / length(returns)

cat("\nLoss days:", loss_count, "\n")
cat("Loss proportion:", round(loss_proportion * 100, 2), "%\n")

gamma_fit <- fitdist(losses, "gamma")
lognormal_fit <- fitdist(losses, "lnorm")

gamma_shape <- unname(gamma_fit$estimate["shape"])
gamma_rate <- unname(gamma_fit$estimate["rate"])

ln_meanlog <- unname(lognormal_fit$estimate["meanlog"])
ln_sdlog <- unname(lognormal_fit$estimate["sdlog"])

loss_parameters <- data.frame(
  Model = c("Gamma", "Lognormal"),
  Parameter1 = c(gamma_shape, ln_meanlog),
  Parameter2 = c(gamma_rate, ln_sdlog)
)

loss_aic <- data.frame(
  Model = c("Gamma", "Lognormal"),
  AIC = c(gamma_fit$aic, lognormal_fit$aic)
)

loss_comparison <- data.frame(
  Model = c("Empirical", "Gamma", "Lognormal"),
  Mean = c(
    mean(losses),
    gamma_shape / gamma_rate,
    exp(ln_meanlog + ln_sdlog^2 / 2)
  ),
  Variance = c(
    var(losses),
    gamma_shape / gamma_rate^2,
    (exp(ln_sdlog^2) - 1) *
      exp(2 * ln_meanlog + ln_sdlog^2)
  )
)

print(loss_parameters, row.names = FALSE, digits = 6)
print(loss_aic, row.names = FALSE, digits = 6)
print(loss_comparison, row.names = FALSE, digits = 6)

hist(
  losses,
  breaks = 35,
  probability = TRUE,
  main = "BYD Loss Size Distribution",
  xlab = "Loss Size (% log return)",
  col = "lightgray",
  border = "white"
)

x_loss <- seq(
  0.0001,
  max(losses) * 1.2,
  length.out = 500
)

lines(
  x_loss,
  dgamma(
    x_loss,
    shape = gamma_shape,
    rate = gamma_rate
  ),
  col = "blue",
  lwd = 2
)

lines(
  x_loss,
  dlnorm(
    x_loss,
    meanlog = ln_meanlog,
    sdlog = ln_sdlog
  ),
  col = "red",
  lwd = 2
)

legend(
  "topright",
  legend = c("Gamma", "Lognormal"),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)


# Task 3.3: Monthly Loss Frequency

trading_data <- data[!is.na(data$Return), ]

trading_data$Month <- format(
  trading_data$Date,
  "%Y-%m"
)

trading_data$Decline2 <- as.integer(
  100 * expm1(trading_data$Return) < -2
)

monthly_events <- aggregate(
  Decline2 ~ Month,
  data = trading_data,
  FUN = sum
)

monthly_days <- aggregate(
  Decline2 ~ Month,
  data = trading_data,
  FUN = length
)

monthly <- data.frame(
  Month = monthly_events$Month,
  Events = monthly_events$Decline2,
  TradingDays = monthly_days$Decline2
)

print(monthly, row.names = FALSE)

poisson_rate <- sum(monthly$Events) /
  sum(monthly$TradingDays)

poisson_mu <- poisson_rate * monthly$TradingDays

poisson_nll <- -sum(dpois(
  monthly$Events,
  lambda = poisson_mu,
  log = TRUE
))

poisson_aic <- 2 + 2 * poisson_nll

nb_nll <- function(par) {
  rate <- exp(par[1])
  theta <- exp(par[2])
  mu <- rate * monthly$TradingDays
  
  -sum(dnbinom(
    monthly$Events,
    size = theta,
    mu = mu,
    log = TRUE
  ))
}

nb_result <- optim(
  par = c(log(poisson_rate), log(10)),
  fn = nb_nll,
  method = "L-BFGS-B",
  lower = c(log(0.0001), log(0.01)),
  upper = c(log(10), log(10000)),
  control = list(maxit = 10000)
)

if (nb_result$convergence != 0) {
  stop("Negative Binomial optimization failed")
}

nb_rate <- exp(nb_result$par[1])
nb_theta <- exp(nb_result$par[2])

nb_aic <- 4 + 2 * nb_result$value

frequency_aic <- data.frame(
  Model = c("Poisson", "Negative Binomial"),
  AIC = c(poisson_aic, nb_aic)
)

frequency_parameters <- data.frame(
  Model = c("Poisson", "Negative Binomial"),
  DailyRate = c(poisson_rate, nb_rate),
  Theta = c(NA, nb_theta)
)

print(frequency_parameters, row.names = FALSE, digits = 6)
print(frequency_aic, row.names = FALSE, digits = 6)

if (nb_theta > 9900) {
  cat(
    "Negative Binomial theta reached the upper bound.\n",
    "The model is close to the Poisson limit.\n"
  )
}

monthly$PoissonExpected <- poisson_mu
monthly$NBExpected <- nb_rate * monthly$TradingDays

plot(
  seq_len(nrow(monthly)),
  monthly$Events,
  type = "b",
  pch = 16,
  main = "Monthly Frequency of Losses Above 2%",
  xlab = "Month Index",
  ylab = "Number of Loss Days",
  col = "black"
)

lines(
  seq_len(nrow(monthly)),
  monthly$PoissonExpected,
  col = "red",
  lwd = 2
)

lines(
  seq_len(nrow(monthly)),
  monthly$NBExpected,
  col = "blue",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c("Observed", "Poisson", "Negative Binomial"),
  col = c("black", "red", "blue"),
  lty = c(1, 1, 2),
  lwd = 2,
  bty = "n"
)


# Task 3.4: Extreme Loss Probability

largest_loss <- max(losses)
threshold <- -largest_loss / 100

p_normal <- pnorm(
  threshold,
  mean = normal_fit$estimate["mean"],
  sd = normal_fit$estimate["sd"]
)

p_student <- pt(
  (threshold - student_params["mu"]) /
    student_params["sigma"],
  df = student_params["nu"]
)

p_skewt <- pst(
  threshold,
  xi = skewt_params["xi"],
  omega = skewt_params["omega"],
  alpha = skewt_params["alpha"],
  nu = skewt_params["nu"]
)

p_gamma <- loss_proportion * pgamma(
  largest_loss,
  shape = gamma_shape,
  rate = gamma_rate,
  lower.tail = FALSE
)

p_lognormal <- loss_proportion * plnorm(
  largest_loss,
  meanlog = ln_meanlog,
  sdlog = ln_sdlog,
  lower.tail = FALSE
)

extreme_results <- data.frame(
  Model = c(
    "Normal",
    "Student t",
    "Skewed Student t",
    "Gamma",
    "Lognormal"
  ),
  Probability = c(
    p_normal,
    p_student,
    p_skewt,
    p_gamma,
    p_lognormal
  )
)

extreme_results$OnceEveryNDays <- ifelse(
  extreme_results$Probability > 0,
  1 / extreme_results$Probability,
  Inf
)

observed_exceedances <- sum(returns < threshold)
observed_at_least <- sum(returns <= threshold)

cat("\nLargest logarithmic loss (%):",
    round(largest_loss, 4), "\n")

cat("Equivalent simple loss (%):",
    round(100 * (1 - exp(-largest_loss / 100)), 2), "\n")

cat("Strict exceedances:", observed_exceedances, "\n")
cat("Observations at least this large:", observed_at_least, "\n")

print(extreme_results, row.names = FALSE, digits = 8)


# Task 3: Final Model Comparison

all_aic <- rbind(
  data.frame(Variable = "Returns", returns_aic),
  data.frame(Variable = "Loss Size", loss_aic),
  data.frame(Variable = "Loss Frequency", frequency_aic)
)

best_models <- do.call(
  rbind,
  lapply(
    split(all_aic, all_aic$Variable),
    function(df) df[which.min(df$AIC), ]
  )
)

rownames(best_models) <- NULL

cat("\nAll model AIC values:\n")
print(all_aic, row.names = FALSE, digits = 6)

cat("\nBest models by AIC:\n")
print(best_models, row.names = FALSE, digits = 6)



# Task 4: Model Selection and Goodness-of-Fit

library(sn)
library(fitdistrplus)

set.seed(123)
B <- 199
alpha_level <- 0.05

dstudent <- function(x, mu, sigma, nu, log = FALSE) {
  result <- dt((x - mu) / sigma, df = nu, log = TRUE) -
    log(sigma)
  if (log) result else exp(result)
}

pstudent <- function(x, mu, sigma, nu) {
  pt((x - mu) / sigma, df = nu)
}

fit_student_sample <- function(x) {
  objective <- function(par) {
    -sum(dstudent(
      x,
      mu = par[1],
      sigma = exp(par[2]),
      nu = 2 + exp(par[3]),
      log = TRUE
    ))
  }
  
  fit <- optim(
    c(
      median(x),
      log(sd(x)),
      log(max(student_params["nu"] - 2, 0.01))
    ),
    objective,
    method = "BFGS",
    control = list(maxit = 3000, reltol = 1e-9)
  )
  
  if (fit$convergence != 0 || !is.finite(fit$value)) {
    stop("Student t fitting failed")
  }
  
  c(
    mu = fit$par[1],
    sigma = exp(fit$par[2]),
    nu = 2 + exp(fit$par[3])
  )
}

fit_skewt_sample <- function(x) {
  objective <- function(par) {
    result <- tryCatch(
      dst(
        x,
        xi = par[1],
        omega = exp(par[2]),
        alpha = par[3],
        nu = 2 + exp(par[4]),
        log = TRUE
      ),
      error = function(e) rep(-Inf, length(x))
    )
    
    if (any(!is.finite(result))) return(1e100)
    -sum(result)
  }
  
  fit <- optim(
    c(
      unname(skewt_params["xi"]),
      log(unname(skewt_params["omega"])),
      unname(skewt_params["alpha"]),
      log(unname(skewt_params["nu"]) - 2)
    ),
    objective,
    method = "L-BFGS-B",
    lower = c(-0.2, log(1e-5), -20, log(0.01)),
    upper = c(0.2, log(1), 20, log(198)),
    control = list(maxit = 3000)
  )
  
  if (fit$convergence != 0 || !is.finite(fit$value)) {
    stop("Skewed Student t fitting failed")
  }
  
  c(
    xi = fit$par[1],
    omega = exp(fit$par[2]),
    alpha = fit$par[3],
    nu = 2 + exp(fit$par[4])
  )
}

ks_statistic <- function(x, cdf, pars) {
  z <- sort(x)
  n <- length(z)
  
  fitted_cdf <- do.call(
    cdf,
    c(list(z), as.list(pars))
  )
  
  max(
    seq_len(n) / n - fitted_cdf,
    fitted_cdf - (seq_len(n) - 1) / n
  )
}

bootstrap_ks <- function(x, pars, cdf, generator,
                         refit, B = 199) {
  observed <- ks_statistic(x, cdf, pars)
  simulated <- numeric(0)
  attempts <- 0
  
  while (length(simulated) < B && attempts < 3 * B) {
    attempts <- attempts + 1
    
    result <- tryCatch({
      sample_b <- generator(length(x), pars)
      pars_b <- refit(sample_b)
      ks_statistic(sample_b, cdf, pars_b)
    }, error = function(e) NA_real_)
    
    if (is.finite(result)) {
      simulated <- c(simulated, result)
    }
  }
  
  if (length(simulated) < B) {
    stop("Too many bootstrap fitting failures")
  }
  
  c(
    Statistic = observed,
    PValue = (1 + sum(simulated >= observed)) / (B + 1)
  )
}


# Task 4.1: Returns Q-Q Plots

normal_pars <- c(
  mean = unname(normal_fit$estimate["mean"]),
  sd = unname(normal_fit$estimate["sd"])
)

p <- ppoints(length(returns))
observed_returns <- sort(returns)

q_normal <- qnorm(
  p,
  mean = normal_pars["mean"],
  sd = normal_pars["sd"]
)

q_student <- student_params["mu"] +
  student_params["sigma"] *
  qt(p, df = student_params["nu"])

q_skewt <- qst(p, dp = skewt_params)

par(mfrow = c(1, 3))

plot(
  q_normal, observed_returns,
  pch = 20,
  col = "gray40",
  main = "Q-Q: Normal",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles"
)
abline(0, 1, col = "red", lwd = 2)

plot(
  q_student, observed_returns,
  pch = 20,
  col = "gray40",
  main = "Q-Q: Student t",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles"
)
abline(0, 1, col = "blue", lwd = 2)

plot(
  q_skewt, observed_returns,
  pch = 20,
  col = "gray40",
  main = "Q-Q: Skewed Student t",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles"
)
abline(0, 1, col = "darkgreen", lwd = 2)

par(mfrow = c(1, 1))


# Task 4.2: Returns Goodness-of-Fit

normal_test <- bootstrap_ks(
  returns,
  normal_pars,
  pnorm,
  function(n, pars) {
    rnorm(n, pars["mean"], pars["sd"])
  },
  function(x) {
    c(
      mean = mean(x),
      sd = sqrt(mean((x - mean(x))^2))
    )
  },
  B
)

student_test <- bootstrap_ks(
  returns,
  student_params,
  pstudent,
  function(n, pars) {
    pars["mu"] +
      pars["sigma"] * rt(n, df = pars["nu"])
  },
  fit_student_sample,
  B
)

skewt_test <- bootstrap_ks(
  returns,
  skewt_params,
  function(x, xi, omega, alpha, nu) {
    pst(
      x,
      xi = xi,
      omega = omega,
      alpha = alpha,
      nu = nu
    )
  },
  function(n, pars) {
    rst(n, dp = pars)
  },
  fit_skewt_sample,
  B
)

returns_tests <- data.frame(
  Model = c(
    "Normal",
    "Student t",
    "Skewed Student t"
  ),
  KS = unname(c(
    normal_test["Statistic"],
    student_test["Statistic"],
    skewt_test["Statistic"]
  )),
  PValue = unname(c(
    normal_test["PValue"],
    student_test["PValue"],
    skewt_test["PValue"]
  )),
  AIC = returns_aic$AIC
)

print(returns_tests, row.names = FALSE, digits = 6)


# Task 4.3: Loss Size Q-Q Plots

p_loss <- ppoints(length(losses))
observed_losses <- sort(losses)

q_gamma <- qgamma(
  p_loss,
  shape = gamma_shape,
  rate = gamma_rate
)

q_lognormal <- qlnorm(
  p_loss,
  meanlog = ln_meanlog,
  sdlog = ln_sdlog
)

par(mfrow = c(1, 2))

plot(
  q_gamma, observed_losses,
  pch = 20,
  col = "gray40",
  main = "Q-Q: Gamma",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles"
)
abline(0, 1, col = "blue", lwd = 2)

plot(
  q_lognormal, observed_losses,
  pch = 20,
  col = "gray40",
  main = "Q-Q: Lognormal",
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles"
)
abline(0, 1, col = "red", lwd = 2)

par(mfrow = c(1, 1))


# Task 4.4: Loss Size Goodness-of-Fit

gamma_pars <- c(
  shape = gamma_shape,
  rate = gamma_rate
)

lognormal_pars <- c(
  meanlog = ln_meanlog,
  sdlog = ln_sdlog
)

gamma_test <- bootstrap_ks(
  losses,
  gamma_pars,
  pgamma,
  function(n, pars) {
    rgamma(
      n,
      shape = pars["shape"],
      rate = pars["rate"]
    )
  },
  function(x) {
    fitdist(x, "gamma")$estimate
  },
  B
)

lognormal_test <- bootstrap_ks(
  losses,
  lognormal_pars,
  plnorm,
  function(n, pars) {
    rlnorm(
      n,
      meanlog = pars["meanlog"],
      sdlog = pars["sdlog"]
    )
  },
  function(x) {
    lx <- log(x)
    c(
      meanlog = mean(lx),
      sdlog = sqrt(mean((lx - mean(lx))^2))
    )
  },
  B
)

loss_tests <- data.frame(
  Model = c("Gamma", "Lognormal"),
  KS = unname(c(
    gamma_test["Statistic"],
    lognormal_test["Statistic"]
  )),
  PValue = unname(c(
    gamma_test["PValue"],
    lognormal_test["PValue"]
  )),
  AIC = loss_aic$AIC
)

print(loss_tests, row.names = FALSE, digits = 6)


# Task 4.5: Loss Frequency Goodness-of-Fit

frequency_statistic <- function(y, mu, theta = Inf) {
  variance <- if (is.infinite(theta)) {
    mu
  } else {
    mu + mu^2 / theta
  }
  
  sum((y - mu)^2 / variance)
}

fit_nb_sample <- function(y) {
  objective <- function(par) {
    rate <- exp(par[1])
    theta <- exp(par[2])
    mu <- rate * monthly$TradingDays
    
    -sum(dnbinom(
      y,
      size = theta,
      mu = mu,
      log = TRUE
    ))
  }
  
  fit <- optim(
    c(
      log(max(sum(y) / sum(monthly$TradingDays), 1e-5)),
      log(max(min(nb_theta, 10000), 0.01))
    ),
    objective,
    method = "L-BFGS-B",
    lower = c(log(0.0001), log(0.01)),
    upper = c(log(10), log(10000)),
    control = list(maxit = 3000)
  )
  
  if (fit$convergence != 0 || !is.finite(fit$value)) {
    stop("Negative Binomial refitting failed")
  }
  
  c(
    rate = exp(fit$par[1]),
    theta = exp(fit$par[2])
  )
}

poisson_mu <- poisson_rate * monthly$TradingDays
nb_mu <- nb_rate * monthly$TradingDays

poisson_stat <- frequency_statistic(
  monthly$Events,
  poisson_mu
)

nb_stat <- frequency_statistic(
  monthly$Events,
  nb_mu,
  nb_theta
)

poisson_sim <- numeric(B)
nb_sim <- numeric(B)

for (b in seq_len(B)) {
  y_p <- rpois(
    nrow(monthly),
    lambda = poisson_mu
  )
  
  fitted_rate_p <- sum(y_p) / sum(monthly$TradingDays)
  fitted_mu_p <- fitted_rate_p * monthly$TradingDays
  
  poisson_sim[b] <- frequency_statistic(
    y_p,
    fitted_mu_p
  )
  
  repeat {
    y_nb <- rnbinom(
      nrow(monthly),
      mu = nb_mu,
      size = nb_theta
    )
    
    fitted_nb <- tryCatch(
      fit_nb_sample(y_nb),
      error = function(e) NULL
    )
    
    if (!is.null(fitted_nb)) break
  }
  
  fitted_mu_nb <- fitted_nb["rate"] * monthly$TradingDays
  
  nb_sim[b] <- frequency_statistic(
    y_nb,
    fitted_mu_nb,
    fitted_nb["theta"]
  )
}

frequency_tests <- data.frame(
  Model = c("Poisson", "Negative Binomial"),
  Pearson = c(poisson_stat, nb_stat),
  PValue = c(
    (1 + sum(poisson_sim >= poisson_stat)) / (B + 1),
    (1 + sum(nb_sim >= nb_stat)) / (B + 1)
  ),
  AIC = frequency_aic$AIC
)

print(frequency_tests, row.names = FALSE, digits = 6)


# Task 4.6: AIC, BIC and Final Comparison

results <- rbind(
  data.frame(
    Variable = "Returns",
    Model = returns_tests$Model,
    Test = "KS",
    Statistic = returns_tests$KS,
    PValue = returns_tests$PValue,
    AIC = returns_tests$AIC,
    k = c(2, 3, 4),
    n = length(returns)
  ),
  data.frame(
    Variable = "Loss Size",
    Model = loss_tests$Model,
    Test = "KS",
    Statistic = loss_tests$KS,
    PValue = loss_tests$PValue,
    AIC = loss_tests$AIC,
    k = c(2, 2),
    n = length(losses)
  ),
  data.frame(
    Variable = "Loss Frequency",
    Model = frequency_tests$Model,
    Test = "Pearson",
    Statistic = frequency_tests$Pearson,
    PValue = frequency_tests$PValue,
    AIC = frequency_tests$AIC,
    k = c(1, 2),
    n = nrow(monthly)
  )
)

results$BIC <- results$AIC +
  results$k * (log(results$n) - 2)

results$Decision <- ifelse(
  results$PValue < alpha_level,
  "Reject H0",
  "Do not reject H0"
)

final_table <- results[
  ,
  c(
    "Variable",
    "Model",
    "Test",
    "Statistic",
    "PValue",
    "AIC",
    "BIC",
    "Decision"
  )
]

print(
  final_table,
  row.names = FALSE,
  digits = 5
)

selected_models <- do.call(
  rbind,
  lapply(
    split(final_table, final_table$Variable),
    function(df) {
      valid <- df[df$PValue >= alpha_level, ]
      
      if (nrow(valid) > 0) {
        valid[which.min(valid$AIC), ]
      } else {
        chosen <- df[which.min(df$AIC), ]
        chosen$Decision <- "Best AIC; GOF rejected"
        chosen
      }
    }
  )
)

rownames(selected_models) <- NULL

cat("\nFinal Selected Models:\n")

print(
  selected_models,
  row.names = FALSE,
  digits = 5
)

