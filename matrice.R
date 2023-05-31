rm(list=ls())

library(MASS)

# number of observations per group
n <- 50

# create vector mu and Sigma variance-covariance matrix
mu <- c(-1.5, 1.5)
Sigma <- matrix(0.8, nrow = 2, ncol = 2); diag(Sigma) <- 1; print(Sigma)

# generate the data from a bivariate normal distribution
biv_normal <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

# Pearson correlation coefficient
cor <- cor(x = biv_normal[,1], y = biv_normal[,2]); print(cor)

# linear model 
Model1 <- lm(biv_normal[,1] ~ biv_normal[,2]); summary(Model1)

# plot the model
plot(biv_normal[,1] ~ biv_normal[,2]); abline(Model1)

## play with matrices computation 
# create the X matrix
# start
X <- as.matrix(biv_normal[,2], ncol = 1)
X <- cbind(rep(1, 50), biv_normal[,2])
X <- as.matrix(X, ncol = 2); print(X)
# end

# X' matrix
Xt <- t(X)
# X'X matrix
XtX <- Xt%*%X
# (X'X)^-1
XtXmin1 <- solve(XtX)
# X'Y
XtY <- Xt%*%biv_normal[,1]
# display the (k+1)*1 vector concatenating the beta coefficients 
beta <- XtXmin1%*%XtY; print(beta)

# perform the same "more manually" 
# X'X
XprimeX <- matrix(c(n, sum(biv_normal[,2]), sum(biv_normal[,2]), sum(biv_normal[,2]^2)),
                  nrow = 2, byrow = T); print(XprimeX)
# X'Y
XprimeY <- matrix(c(sum(biv_normal[,1]), sum(biv_normal[,1]*biv_normal[,2])),
                  nrow = 2); print(XprimeY)
# (X'X)^-1
XprimeXmin1 <- solve(XprimeX)

# display the (k+1)*1 vector concatenating the beta coefficients 
betaVec <- matrix(c(sum(XprimeXmin1[1,1]*XprimeY[1,1], XprimeXmin1[1,2]*XprimeY[2,1]),
                    sum(XprimeXmin1[2,1]*XprimeY[1,1], XprimeXmin1[2,2]*XprimeY[2,1])),
                  nrow = 2); print(betaVec)

# find the hat matrix 
H <- X%*%ginv(t(X)%*%X)%*%t(X); print(H)
# H <- X%*%solve(t(X)%*%X)%*%t(X)
# put the hat on Y
Yhat <- H%*%biv_normal[,1]; print(Yhat)
# residuals
e <- biv_normal[,1] - Yhat
# total variance
eq <- (biv_normal[,1] - mean(biv_normal[,1]))^2
# R^2
Rsquare <- sum(eq)/sum(sum(e^2), sum(eq)); print(Rsquare)
# Adjusted-R^2
AdjRsquare <- 1-((1-Rsquare)*(nrow(biv_normal)-1)/(nrow(biv_normal)- (ncol(biv_normal)-1)-1)); print(AdjRsquare)
# Residual standard error
ResSE <- sqrt(sum(e^2)/48)

## find beta 0 (i.e., the intercept) and beta 1 (i.e., the slope) by Covariance/Variance
beta1 <- cov(x = biv_normal[,1], y = biv_normal[,2])/var(biv_normal[,2]); print(beta1)
beta0 <- mean(biv_normal[,1] - beta1*mean(biv_normal[,2])); print(beta0)

