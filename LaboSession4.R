# ==============================================================================
# Macroeconometrics
# Jean-Paul Renne, 2024
# Laboratory session 4 - VAR models
# ==============================================================================
# This session  requires a preliminary run of various_proc_TS.R
# ==============================================================================

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


# Question 1.b:

mu <- XXXX

Omega.star <- matrix(0,4,4)
Omega.star[1:2,1:2] <- B %*% t(B)

unc.var <- XXX

y0.star <- mu
nb.sim <- 1000
Y <- simul.VAR(c,Phi,B,nb.sim,y0.star)
plot(Y[,1],type="l")

print(c(mean(Y[,1]),mean(Y[,2])))
print(var(Y))


# Question 1.c:


# Question 1.d:
Theta.0 <- XXXX
Theta.1 <- XXXX
Theta.2 <- XXXX
Theta.3 <- XXXX
Theta.4 <- XXXX
Omega <- B %*% t(B)
cond.var <- XXXX
  

# Question 1.e:
nb.sim <- 40
par(mfrow=c(2,2))

Y <- XXXX
plot(Y[,1],type="l",main="Response of y1 to u1")
plot(Y[,2],type="l",main="Response of y2 to u1")

Y <- XXXX
plot(Y[,1],type="l",main="Response of y1 to u2")
plot(Y[,2],type="l",main="Response of y2 to u2")



# Exercise 2 - Blanchard and Quah
# --------------------------------

# Load data:
library (readr)
urlfile="https://github.com/jrenne/LaboSessionMacro/blob/main/Data/dataLaboSession4_BQ.csv"
mydata<-read_csv(url(urlfile))

dataBQ <- read.csv("https://github.com/jrenne/LaboSessionMacro/blob/main/Data/dataLaboSession4_BQ.csv")
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

y <- XXXX

est.VAR <- VAR(y,p=8)

# get estimated Phi matrices:
Phi <- Acoef(est.VAR)

# Check stationarity of the estimated model:
XXXX

# Question 2.b:

sum.PHI.k <- diag(dim(PHI)[1])
PHI.k <- PHI
for(i in 1:500){
  sum.PHI.k <- XXXX
  PHI.k <- XXXX
}

# Compute the covariance matrix of VAR residuals:
resids <- residuals(est.VAR)
Omega <- var(resids)


# Question 2.c:

f <- function(param){
  B <- matrix(param,2,2)
  X <- Omega - B %*% t(B)
  Theta <- sum.PHI.k[1:2,1:2] %*% B
  loss <- XXXX
  return(loss)
}

XXXXX

# Question 2.d:

nb.sim <- 40
par(mfrow=c(2,2))
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(1,0))
plot(cumsum(Y[,1]),type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
plot(Y[,2],type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(0,1))
plot(cumsum(Y[,1]),type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")
plot(Y[,2],type="l",lwd=2,xlab="",ylab="",main="Effect of XXX on XXX")



