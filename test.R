library(devtools)
install_github("AU-R-Programming/FinalProject_Group9", build_vignettes=TRUE)

# install.packages("/Users/kerry/Dropbox/repos/r-class/FinalProject_Group9",
                 # repos = NULL,
                 # type = "source")

library(MASS)
library(Group9LinearModel)
data(Boston)
fit = myLm(Boston$crim,Boston[c("age", "medv")])
fit

confint(fit)
confint(fit, alpha=.1)
confint(fit, alpha=.1, approach="boot")

plot(fit)
qqPlot(fit)
hist(fit)

?myLm
browseVignettes("Group9LinearModel")
