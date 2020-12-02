

myLm = function(response, covariates) {

  stopifnot(length(response) == nrow(covariates))


  ### Much of this code was modified from the textbook for the course
  # Make sure data formats are appropriate
  response <- as.vector(response)
  covariates <- as.matrix(covariates)



  # Define sample size
  n <- length(response)

  #add column of 1s for intercept
  covariates=cbind(rep(1,n),covariates)

  #Define parameters
  p <- dim(covariates)[2]
  df <- n - p

  # Estimate beta through Eq. (6.1)
  beta.hat <- solve(t(covariates)%*%covariates)%*%t(covariates)%*%response
  rownames(beta.hat)[1] = "intercept"

  # Estimate of the residual variance (sigma2) from Eq. (6.3)
  # Compute residuals
  yHat=covariates%*%as.matrix(beta.hat)
  resid <- response - covariates%*%as.matrix(beta.hat)
  sigma2.hat <- (1/df)*t(resid)%*%resid
  sigma.hat = sigma2.hat^0.5

  # Estimate of the variance of the estimated beta from Eq. (6.2)
  #A=(solve(t(covariates)%*%covariates))%*%t(covariates)
  #var.beta=A%*%t(A)%*%sigma2.hat
  var.beta <- as.numeric(sigma2.hat)*(solve(t(covariates)%*%covariates))

  #Compute the Mean Square Prediction Error
  MSPE=1/n*sum(resid^2)

  #Compute the F-statistic
  yBar=mean(response)
  SSM=sum((covariates%*%as.matrix(beta.hat)-yBar)^2)
  SSE=sum(resid^2)
  DFM=p-1
  DFE=n-p
  MSM=SSM/DFM
  MSE=SSE/DFE
  FStat=MSM/MSE
  Pvalue=pf(FStat,df1=DFM,df2=DFE,lower.tail=FALSE)

  # Create myLm class object
  values = list(
    betas=beta.hat,
    sigma2 = sigma2.hat,
    variance_beta = var.beta,
    residuals = resid,
    predictors = covariates,
    y=response,
    df=df,
    yHat=yHat,
    MSPE=MSPE,
    SSM=SSM,
    SSE=SSE,
    FStatistic=FStat,
    Pvalue=Pvalue)
  class(values) <- "myLm"

  return(values)
}

print.myLm = function(x) {
  cat("Coefficients:\n")
  cat("-------------\n")
  ids = rownames(x$betas)
  for (i in 1:length(ids)) {
    cat(ids[i], ": ", x$betas[i], "\n")
  }
}

#calculate confidence intervals using asymptotic or bootstrap methods
confint.myLm=function(x,alpha=0.05,approach="asymp"){
  stopifnot(alpha > 0 && alpha < 1)
  stopifnot(is.element(approach, c("asymp", "boot")))

  quant <- 1 - alpha/2
  if(approach=="asymp"){
    #ci.beta <- c(beta.hat - qnorm(p = quant)*sqrt(var.beta), beta.hat +
    #               qnorm(p = quant)*sqrt(var.beta))
    #modified to used t-statistics rather than z-statistic
    ci.betas=data.frame(lwr=x$betas-qt(p=quant,x$df)*sqrt(diag(x$variance_beta)),upr=
              x$betas+qt(p=quant,x$df)*sqrt(diag(x$variance_beta)))
  } else if(approach=="boot"){
    betaMatrix=matrix(NA,nrow=1000,ncol=length(x$betas))
    for(i in 1:1000){
      bootData=cbind(x$y,x$predictors[,-1])
      tempVector=1:length(bootData[,1])
      tempSample=sample(tempVector,length(tempVector),replace=TRUE)
      tempData=bootData[tempSample,]
      tempResults=myLm(tempData[,1],tempData[,-1])$betas
      betaMatrix[i,]=tempResults
    }
    ci.betas=matrix(NA,length(x$predictors[1,]),2)
    for(i in 1:length(x$predictors[1,])){
      ci.betas[i,]=quantile(betaMatrix[,i],probs=c(1-quant,quant))
    }
  }
  rownames(ci.betas) <- rownames(x$betas)
  return(ci.betas)
}


# Residuals vs Fitted Plot
plot.myLm = function(x) {
  fittedvalues=x$yHat
  plot(fittedvalues, x$residuals,ylab="residuals",xlab="Fitted Y Values")
}

# qqPlot
qqPlot = function(x) {
  qqnorm(x$residuals)
  qqline(x$residuals)
}

# Plot histogram with residuals
hist.myLm = function(x) {
  hist(x$residuals,xlab="residuals")
}

# library(MASS)
# library(Group9LinearModel)
# data(Boston)
# fit = myLm(Boston$crim,Boston[c("age", "medv")])
# fit
#
# confint(fit)
# confint(fit, alpha=.1)
# confint(fit, alpha=.1, approach="boot")
#
# plot(fit)
# qqPlot(fit)
# hist(fit)

