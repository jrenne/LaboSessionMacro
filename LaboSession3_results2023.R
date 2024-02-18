# =============================================================
# Macroeconometrics
# Laboratory session 3 - Estimation of ARMA Processes
# Jean-Paul Renne, 2023
# =============================================================
# This session  requires a preliminary run of various_proc_TS.R
# =============================================================


# Exercise 1
# --------------------------------

# Load data:
data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession3a.csv")
#data <- read.csv("~/Google Drive/cours/UNIL/Rcode/data/dataLaboSession3a.csv", sep="")
y <- data$V1
plot(y,type="l")


# Function that computes the density of vector eps, where the elements
# of vector eps are Gaussian i.i.d. variables N(mu,sigma^2):
log.l.gaussian <- function(eps,mu,sigma2){
  vec.log.f <- - 1/2*log(2*pi*sigma2) - (eps - mu)^2/(2*sigma2)
  return(vec.log.f)
}

x <- rnorm(10)
log.l.gaussian(x,0,2)

# Question 1.a:

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

exact.log.L(c(1,.6,1),y)


# Question 1.b:
# Optimization of the exact log-likelihood.

param0 <- c(1,.5,1)

res.opt <- optim(param0,exact.log.L,Y=y,
                 # method="Nelder-Mead",
                 method="BFGS",
                 # method="CG",
                 control=list(trace=FALSE,maxit=5000),
                 hessian=TRUE)
print(res.opt$par)

# Question 1.c:
# Compute 95% confidence intervals for the parameters:
I_1 <- solve(res.opt$hessian)
eigen(I_1)$values
std.dev <- sqrt(diag(I_1))
CI <- cbind(res.opt$par - 2.58*std.dev,res.opt$par + 2.58*std.dev)


# Question 1.d:
# The conditional log-likelihood is maximised by the OLS estimator:
T <- length(y)
eq.ols <- lm(y[2:T]~y[1:(T-1)])
summary(eq.ols)
summary.eq <- summary(eq.ols)
summary.eq$sigma

# Question 1.e:
print(
  cbind(eq.ols$coefficients,res.opt$par[1:2])
)



# Question 1.f:

fit.arima <- arima(y,order=c(3,0,0),method = "ML")
theta <- fit.arima$coef

# Wald test:
R <- matrix(0,2,4)
R[1,2] <- 1
R[2,3] <- 1
Var.coef <- fit.arima$var.coef
xi.W <- t(R%*%theta) %*% solve(R%*%Var.coef%*%t(R)) %*% (R%*%theta)
# Critical value for 95%:
qchisq(.95,df=2)
# p-value of the test:
pchisq(xi.W,df=2)

# LR test:
fit.arima <- arima(y,order=c(1,0,0),method = "ML")
log.L.restr <- fit.arima$loglik
fit.arima <- arima(y,order=c(3,0,0),method = "ML")
log.L.unrestr <- fit.arima$loglik
xi.LR <- 2*(log.L.unrestr - log.L.restr)
pchisq(xi.LR,df=2)


# Question 1.g:
T <- length(y)
vec.AIC <- NULL
vec.BIC <- NULL
vec.HQ  <- NULL
for(p in 1:10){
  fit.arima <- arima(y,order=c(p,0,0),method="ML")
  vec.AIC <- c(vec.AIC,
               -2*fit.arima$loglik/T + (p+2)*2/T)
  vec.BIC <- c(vec.BIC,
               -2*fit.arima$loglik/T + (p+2)*log(T)/T)
  vec.HQ <- c(vec.HQ,
              -2*fit.arima$loglik/T + (p+2)*2*log(log(T))/T)
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
    matrix.AIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*2/T
    matrix.BIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*log(T)/T
    matrix.HQ[p+1,q+1]  <- -2*fit.arima$loglik/T + (p+q+2)*2*log(log(T))/T
  }
}






# Exercise 2
# --------------------------------

# Load data:
data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession3b.csv")
#data <- read.csv("~/Google Drive/cours/UNIL/Rcode/data/dataLaboSession3b.csv", sep="")

y <- data$V1
plot(y,type="l")

T <- length(y)
max.p <- 5
max.q <- 5
matrix.AIC <- matrix(NaN,max.p+1,max.q+1)
matrix.BIC <- matrix(NaN,max.p+1,max.q+1)
matrix.HQ  <- matrix(NaN,max.p+1,max.q+1)
for(p in 0:max.p){
  for(q in 0:max.q){
    fit.arima <- arima(y,order=c(p,0,q),method="ML")
    matrix.AIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*2/T
    matrix.BIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*log(T)/T
    matrix.HQ[p+1,q+1]  <- -2*fit.arima$loglik/T + (p+q+2)*2*log(log(T))/T
  }
}
print(min(matrix.AIC))
print(min(matrix.BIC))
print(min(matrix.HQ))

fit.arima <- arima(y,order=c(3,0,1),method = "ML")



# Exercise 3
# --------------------------------


#swiss_cpi <- read.csv("~/Google Drive/cours/UNIL/Rcode/data/swiss_cpi.csv")
swiss_cpi <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/swiss_cpi.csv")
dates <- as.Date(
  paste("01/",swiss_cpi$Date,sep=""),
  "%d/%m/%Y")
cpi <- log(swiss_cpi$CPI)
T <- length(cpi) # update sample length
y <- 100 * (cpi[13:T] - cpi[1:(T-12)])
dates <- dates[13:T]

first.date <- 300
y <- y[first.date:length(y)]
dates <- dates[first.date:length(dates)]
plot(dates,y,type="l")

T <- length(y) # update sample length


# Look for appropriate specification:

max.p <- 5
max.q <- 5
matrix.AIC <- matrix(NaN,max.p+1,max.q+1)
matrix.BIC <- matrix(NaN,max.p+1,max.q+1)
matrix.HQ  <- matrix(NaN,max.p+1,max.q+1)
for(p in 0:max.p){
  for(q in 0:max.q){
    fit.arima <- arima(y,order=c(p,0,q),method="CSS")
    matrix.AIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*2/T
    matrix.BIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*log(T)/T
    matrix.HQ[p+1,q+1]  <- -2*fit.arima$loglik/T + (p+q+2)*2*log(log(T))/T
  }
}

rownames(matrix.AIC) <- paste("p =",0:max.p)
rownames(matrix.BIC) <- paste("p =",0:max.p)
rownames(matrix.HQ)  <- paste("p =",0:max.p)

colnames(matrix.AIC) <- paste("q =",0:max.q)
colnames(matrix.BIC) <- paste("q =",0:max.q)
colnames(matrix.HQ)  <- paste("q =",0:max.q)

print(matrix.AIC)
print(matrix.BIC)
print(matrix.HQ)

print(min(matrix.AIC))
print(min(matrix.BIC))
print(min(matrix.HQ))

# Based on AIC, best specification:
p <- 5
q <- 3

fit.arima <- arima(y,order=c(p,0,q)) # based on AIC
abs(eigen(make.F(fit.arima$coef[1:5]))$value)

phi <- fit.arima$coef[1:p]
theta <- c(1,fit.arima$coef[(p+1):(p+q)])
mu <- fit.arima$coef[p+q+1]
c <- mu * (1 - sum(phi))
sim.y <- sim.arma(c = c,phi = phi,theta = theta,
                  sigma = sqrt(fit.arima$sigma2),
                  T = 300, # length of simulated sample
                  y.0 = rep(mu,p), # starting value
                  nb.sim = 1, # simulate only one sample
                  make.IRF = 0
                  )

plot(sim.y,type="l")


