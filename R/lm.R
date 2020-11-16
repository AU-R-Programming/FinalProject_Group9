

my_lm = function(response, covariates, alpha = 0.05) {

  # Make sure data formats are appropriate
  response <- as.vector(response)
  covariates <- as.matrix(covariates)

  # Define parameters
  n <- length(response)
  p <- dim(covariates)[2]
  df <- n - p

  # Estimate beta through Eq. (6.1)
  beta.hat <- solve(t(covariates)%*%covariates)%*%t(covariates)%*%response

  # Estimate of the residual variance (sigma2) from Eq. (6.3)
  # Compute residuals
  resid <- response - covariates%*%as.matrix(beta.hat)
  sigma2.hat <- (1/df)*t(resid)%*%resid

  # Estimate of the variance of the estimated beta from Eq. (6.2)
  var.beta <- sigma2.hat*solve(t(covariates)%*%covariates)

  # Estimate of the confidence interval based on alpha
  quant <- 1 - alpha/2
  ci.beta <- c(beta.hat - qnorm(p = quant)*sqrt(var.beta), beta.hat +
                 qnorm(p = quant)*sqrt(var.beta))

  # Return all estimated values
  # return(list(beta = beta.hat, sigma2 = sigma2.hat,
  #             variance_beta = var.beta, ci = ci.beta))

  # Create myLm class object to hold output
  values = list(beta = beta.hat, sigma2 = sigma2.hat, variance_beta = var.beta,
       ci = ci.beta)
  class(values) <- "myLm"

  return(values)
}

print.myLm = function(x) {
  cat("Beta: ", x$beta, "\n")
  cat("Sigma^2: ", x$sigma2, "\n")
  cat("Variance Beta: ", x$variance_beta, "\n")
  cat("Confidence Interval: ", x$ci[1], "-", x$ci[2], "\n")
}

# Residuals vs Fitted Plot
plot.myLm = function(x) {

}

# qqPlot
qqPlot.myLm = function(x) {

}

# Plot histogram with residuals
hist.myLm = function(x) {

}



library(gamair)
data(hubble)
fit = my_lm(hubble$y,hubble$x)


fit







