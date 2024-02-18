# =============================================================
# Macroeconometrics
# Laboratory session 2 - VAR models
# Jean-Paul Renne, 2016
# =============================================================
# This session  requires a preliminary run of various_proc_TS.R
# =============================================================

library(vars)



# Exercise 1 - VAR simulation
# --------------------------------

nb.sim <- 100
c <- c(1,0)
B <- matrix(c(1,.5,0,1),2,2)
Phi <- list(
  matrix(c(.3,0,.3,.5),2,2),
  matrix(c(.3,0,0,.3),2,2)
)

# Question 1.a:
PHI <- make.PHI(Phi)
abs(eigen(PHI)$values)


# Question 1.b:
c.star <- c(c,rep(0,2))
mu <- solve(diag(4) - PHI) %*% c.star
Omega.star <- matrix(0,4,4)
Omega.star[1:2,1:2] <- B %*% t(B)
unc.var <- Omega.star
PHI.k <- PHI
for(i in 1:200){
  unc.var <- unc.var +  PHI.k %*% Omega.star %*% t(PHI.k)
  PHI.k <- PHI.k %*% PHI
}

y0.star <- mu
nb.sim <- 1000
Y <- simul.VAR(c,Phi,B,nb.sim,y0.star)
plot(Y[,1],type="l")

print(c(mean(Y[,1]),mean(Y[,2])))
print(var(Y))

# Question 1.c:
print(unc.var[3:4,1:2])

# Question 1.d:
Theta.0 <- diag(2)
Theta.1 <- PHI[1:2,1:2]
Theta.2 <- (PHI %*% PHI)[1:2,1:2]
Theta.3 <- (PHI %*% PHI %*% PHI)[1:2,1:2]
Theta.4 <- (PHI %*% PHI %*% PHI %*% PHI)[1:2,1:2]
Omega <- B %*% t(B)
cond.var <- Omega +
  Theta.1 %*% Omega %*% t(Theta.1) +
  Theta.2 %*% Omega %*% t(Theta.2) +
  Theta.3 %*% Omega %*% t(Theta.3) +
  Theta.4 %*% Omega %*% t(Theta.4)
  
print(sqrt(cond.var[1,1]))


# Question 1.e:
nb.sim <- 40
par(mfrow=c(2,2))
par(plt=c(.1,.9,.1,.7))
Y <- simul.VAR(c,Phi,B,nb.sim,y0.star,indic.IRF = 1,u.shock = c(1,0))
plot(Y[,1],type="l",main="Response of y1 to u1")
plot(Y[,2],type="l",main="Response of y2 to u1")
Y <- simul.VAR(c,Phi,B,nb.sim,y0.star,indic.IRF = 1,u.shock = c(0,1))
plot(Y[,1],type="l",main="Response of y1 to u2")
plot(Y[,2],type="l",main="Response of y2 to u2")


# Exercise 2 - Blanchard and Quah
# --------------------------------

# Load data:
dataBQ <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession4_BQ.csv")
dataBQ$dates <- as.Date(dataBQ$dates,"%m/%d/%Y")

First.date <- "1948-04-01"
Last.date  <- "1988-01-01"
Mid.date   <- "1974-01-01"

# Find position of the previous dates in the dataframe:
fst.date.indic <- which(dataBQ$dates==First.date)
lst.date.indic <- which(dataBQ$dates==Last.date)
mid.date.indic <- which(dataBQ$dates==Mid.date)

first.sample.part <- 1:(mid.date.indic-fst.date.indic)
secnd.sample.part <- (mid.date.indic-fst.date.indic + 1):(lst.date.indic-fst.date.indic+1)

# Reduce data.frame according to selection of dates:
dataBQ <- dataBQ[fst.date.indic:lst.date.indic,]

# Take log of GDP:
y <- log(dataBQ$gdp)
# Lagged values of log of GDP:
y_1 <- c(NaN,log(dataBQ$gdp[1:(length(dataBQ$gdp)-1)]))
# Compute growth rate:
delta.Y <- (y - y_1)*100


# ---------------------------------------------
# Prepare graphics showing data:
par(mfrow=c(2,1))

# ------ GDP growth:
plot(dataBQ$dates,delta.Y,type="l",main="GDP quarterly growth rate",xlab="",ylab="",lwd=2)
# Remove means of GDP growth on two subsamples:
delta.Y[first.sample.part] <- delta.Y[first.sample.part] - mean(delta.Y[first.sample.part],na.rm=TRUE)
delta.Y[secnd.sample.part] <- delta.Y[secnd.sample.part] - mean(delta.Y[secnd.sample.part],na.rm=TRUE)
lines(dataBQ$dates,delta.Y,col="red",lwd=2)


plot(dataBQ$dates,dataBQ$u,type="l",ylim=c(-3,12),main="Unemployment rate",xlab="",ylab="",lwd=2)
# Remove trend in unemployment series:
u <- dataBQ$u
trend <- 1:length(u)
eq <- lm(u ~ trend)
u <- eq$residuals
lines(dataBQ$dates,u,col="red",lwd=2)
abline(h=0,col="grey")


# Question 2.a:

y <- cbind(delta.Y,u)
y <- y[2:dim(y)[1],]

est.VAR <- VAR(y,p=8)



# get estimated Phi matrices:
Phi <- Acoef(est.VAR)

# Check stationarity of the estimated model:
PHI <- make.PHI(Phi)
print(
  abs(eigen(PHI)$values)
)


# Question 2.b:

sum.PHI.k <- diag(dim(PHI)[1])
PHI.k <- PHI
for(i in 1:500){
  sum.PHI.k <- sum.PHI.k +  PHI.k
  PHI.k <- PHI.k %*% PHI
}

# Compute the covariance matrix of VAR residuals:
resids <- residuals(est.VAR)
Omega <- var(resids)


# Question 2.c:

f <- function(param){
  B <- matrix(param,2,2)
  X <- Omega - B %*% t(B)
  Theta <- sum.PHI.k[1:2,1:2] %*% B
  loss <- 10000 * ( X[1,1]^2 + X[2,1]^2 + X[2,2]^2 + Theta[1,1]^2 )
  return(loss)
}

param0 <- c(1,0,0,1)
res.opt <- optim(param0,f,
                 #method="Nelder-Mead",
                 method="BFGS",
                 hessian=FALSE)
print(res.opt$par)
B.hat <- matrix(res.opt$par,2,2)
print(cbind(Omega,B.hat %*% t(B.hat)))

# Question 2.d:

nb.sim <- 40
par(mfrow=c(2,2))
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(1,0))
plot(cumsum(Y[,1]),type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
plot(Y[,2],type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(0,1))
plot(cumsum(Y[,1]),type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
plot(Y[,2],type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")



