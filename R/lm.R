

myLm = function(response, covariates) {

  
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
    #ci = ci.beta,
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
  cat("Beta: ", x$beta, "\n")
  cat("Sigma^2: ", x$sigma2, "\n")
  cat("Variance Beta: ", x$variance_beta, "\n")
  cat("Confidence Interval: ", x$ci[1], "-", x$ci[2], "\n")
}

#calculate confindence intervals using asymptotic or bootstrap methods
confint.myLm=function(x,alpha=0.05,approach="asymp"){
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
    ci.betas[1,]=quantile(betaMatrix[,1],probs=c(1-quant,quant))
    ci.betas[2,]=quantile(betaMatrix[,2],probs=c(1-quant,quant))
    ci.betas[3,]=quantile(betaMatrix[,3],probs=c(1-quant,quant))
  }
  return(ci.betas)
}
  

# Residuals vs Fitted Plot
plot.myLm = function(x) {
  #fittedValues = x$predictors * x$beta
  fittedvalues=x$yHat
  p = plot(fittedvalues, x$residuals)
  return(p)
}

# qqPlot
qqPlot = function(x) {
  stopifnot(class(x) == "myLm")
  p = qqnorm(x$residuals)
  qqline(x$residuals)
  return(p)
}

# Plot histogram with residuals
hist.myLm = function(x) {
  hist(x$residuals)
}



library(gamair)
data(hubble)
fit = myLm(hubble$y,hubble$x)

fit
plot(fit)
qqPlot(fit)
hist(fit)


