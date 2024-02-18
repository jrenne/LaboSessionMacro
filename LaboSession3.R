# =============================================================
# Macroeconometrics
# Laboratory session 2 - Estimation of ARMA Processes
# Jean-Paul Renne, 2016
# =============================================================
# This session  requires a preliminary run of various_proc_TS.R
# =============================================================


# Exercise 1
# --------------------------------

# Load data:
data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession3a.csv")
y <- data$V1
plot(y,type="l")


# Function that computes the density of vector eps, where the elements
# of vector eps are Gaussian i.i.d. variables N(mu,sigma^2):
log.l.gaussian <- function(eps,mu,sigma2){
  vec.log.f <- - 1/2*log(2*pi*sigma2) - (eps - mu)^2/(2*sigma2)
  return(vec.log.f)
}

eps <- c(1,0,-1.2,3)
log.l.gaussian(eps,.5,10)


# Question 1.a:
# Replace the "1111111" in the line beginning by "sigma.eps"
exact.log.L <- function(param,Y){
  c <- param[1]
  phi <- param[2]
  sigma2 <- abs(param[3])
  T <- length(Y)
  sigma.eps <- Y[2:T] - c - phi*Y[1:(T-1)]
  log.Lik <- log.l.gaussian(sigma.eps,0,sigma2)
  log.L <- sum(log.Lik) - 1/2*log(2*pi*sigma2) + 1/2*log(1-phi^2) -
    (y[1] - c/(1-phi))^2 / (2* (sigma2/(1 - phi^2)))
  return(-log.L) # we return -log-lik because the optim procedure of R minimises functions
}

exact.log.L(c(15,0,2.5^2),y)



# Question 1.b:
# Optimization of the exact log-likelihood.

param0 <- c(15,0,2.5^2)

res.opt <- optim(param0,exact.log.L,Y=y,
                 # method="Nelder-Mead",
                 method="BFGS",
                 control=list(trace=FALSE,maxit=5000),
                 hessian=TRUE)


# Question 1.c:
# Compute 95% confidence intervals for the parameters:
M <- res.opt$hessian
V <- solve(M)

theta.MLE <- res.opt$par

lower.bounds <- theta.MLE - 1.96*sqrt(diag(V))
upper.bounds <- theta.MLE + 1.96*sqrt(diag(V))
print(cbind(lower.bounds,upper.bounds))


# Question 1.d:
T <- length(y)
eq.ols <- lm(y[2:T]~y[1:(T-1)])
summary(eq.ols)
summary.eq <- summary(eq.ols)
summary.eq$sigma^2

# Question 1.e:
print(
  cbind(eq.ols$coefficients,res.opt$par[1:2])
)


# Question 1.f:

fit.arima <- arima(y,order=c(2,0,0),method = "ML")
theta <- fit.arima$coef


# Wald test:
R <- matrix(0,2,4)
R[1,2] <- 1
R[2,3] <- 1
Var.coef <- fit.arima$var.coef
xi.W <- 1111111111111111111

# LR test:
111111111111111
xi.LR <- 2*(log.L.unrestr - log.L.restr)
pchisq(xi.LR,df=2)


# Question 1.g:
T <- length(y)
vec.AIC <- NULL
vec.BIC <- NULL
vec.HQ  <- NULL
for(p in 1:10){
  fit.arima <- arima(y,order=c(p,0,0),method="ML")
  vec.AIC <- c(vec.AIC,111111111111111111)
  vec.BIC <- c(vec.BIC,111111111111111111)
  vec.HQ <- c(vec.HQ,111111111111111111)
}
plot(vec.AIC,lwd=2,type="l")
plot(vec.BIC,lwd=2,type="l")
plot(vec.HQ,lwd=2,type="l")


# Question 1.h:
max.p <- 5
max.q <- 5
matrix.AIC <- matrix(NaN,max.p+1,max.q+1)
matrix.BIC <- matrix(NaN,max.p+1,max.q+1)
matrix.HQ  <- matrix(NaN,max.p+1,max.q+1)
for(p in 0:max.p){
  for(q in 0:max.q){
    fit.arima <- arima(y,order=c(p,0,q),method="ML")
    matrix.AIC[p+1,q+1] <- 111111111111111
    matrix.BIC[p+1,q+1] <- 111111111111111
    matrix.HQ[p+1,q+1]  <- 111111111111111
  }
}




# Exercise 2
# --------------------------------

# Load data:
data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession3b.csv")
y <- data$V1
plot(y,type="l")
T <- length(y)



# Exercise 3
# --------------------------------


#swiss_cpi <- read.csv("~/Google Drive/cours/UNIL/Rcode/data/swiss_cpi.csv")
swiss_cpi <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/swiss_cpi.csv")
dates <- as.Date(
  paste("01/",swiss_cpi$Date,sep=""),
  "%d/%m/%Y")
cpi <- log(swiss_cpi$CPI)
T <- length(cpi)
y <- 100 * (cpi[13:T] - cpi[1:(T-12)])
dates <- dates[13:T]

y <- y[300:length(y)]
dates <- dates[300:length(dates)]
plot(dates,y,type="l")


