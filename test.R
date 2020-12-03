library(devtools)
install_github("AU-R-Programming/FinalProject_Group9")

install.packages("/Users/kerry/Dropbox/repos/r-class/FinalProject_Group9",
                 repos = NULL,
                 type = "source")


library(MASS)
library(Group9LinearModel)
data(Boston)
fit = Group9LinearModel::myLm(Boston$crim,Boston[c("age", "medv")])
print.myLm(fit)

confint.myLm(fit)
confint.myLm(fit, alpha=.1)
confint.myLm(fit, alpha=.1, approach="boot")

plot.myLm(fit)
qqPlot(fit)
Group9LinearModel::hist.myLm(fit)


?Group9LinearModel
