# =============================================================
# Macroeconometrics
# Laboratory session 1 - Long run variance
# Jean-Paul Renne, 2016
# =============================================================

# Exercise 1
# --------------------------------















# Question 1.a:

autocov <- function(X,n){
  T <- length(X)
  X.1 <- X[1:(T-n)] - mean(X)
  X.2 <- X[(n+1):T] - mean(X)
  return(1/T * sum(X.1 * X.2))
}

# test:
X <- rnorm(100)
print(autocov(X,2))


# Question 1.b:

NW.Weights <- function(q){
  weights <- 1 - (1:q)/(q+1)
  return(weights)
}

# test:
print(NW.Weights(5))


# Question 1.c:

NW.LongRunVariance <- function(X,q){
  gamma0 <- autocov(X,0)
  LRV <- gamma0
  weights <- NW.Weights(q)
  for(i in 1:q){
    LRV <- LRV + 2 * weights[i] * autocov(X,i)
  }
  return(LRV)
}

# test:
T <- 10000
X <- rnorm(T)
q <- 3
print(NW.LongRunVariance(X,q))

# Question 1.d:

data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession1a.csv")
x <- data$x
plot(data$x,type="l")
abline(h=mean(x),col="red")
T <- length(data$x)

q <- 10
NW.LongRunVariance(data$x,q)
q <- 20
NW.LongRunVariance(data$x,q)
q <- 50
NW.LongRunVariance(data$x,q)
q <- 100
NW.LongRunVariance(data$x,q)


# With too small a q:
q <- 1
lower.bound <- mean(x) - 1.96*sqrt(NW.LongRunVariance(data$x,q))/sqrt(T)
upper.bound <- mean(x) + 1.96*sqrt(NW.LongRunVariance(data$x,q))/sqrt(T)
abline(h=lower.bound,col="blue",lty=2)
abline(h=upper.bound,col="blue",lty=2)

q <- 100
lower.bound <- mean(x) - 1.96*sqrt(NW.LongRunVariance(data$x,q))/sqrt(T)
upper.bound <- mean(x) + 1.96*sqrt(NW.LongRunVariance(data$x,q))/sqrt(T)
abline(h=lower.bound,col="red",lty=2)
abline(h=upper.bound,col="red",lty=2)



# Exercise 2
# --------------------------------

# These functions will be used in that exercise:
VAR.sum.AR <- function(rho,sigma,T){
  sum.var <- 0
  sum.rho <- 0
  for(i in 0:(T-1)){
    sum.rho <- sum.rho + rho^i # for iteration i, this is equal to 1 + rho + rho^2 + ... + rho^{i}
    sum.var <- sum.var + sum.rho^2*sigma^2
  }
  return(T * 1/T^2 *(sum.var + sum.rho^2*sigma^2*rho^2/(1-rho^2)))
}

simul.ar <- function(c,rho,sigma,T,x.0){
  x <- x.0
  for(t in 2:T){
    x <- c(x,c+rho*x[length(x)]+sigma*rnorm(1))
  }
  return(x)
}

# Question 2.d:
c <- 1
rho <- .5
sigma <- 1
T <- 400
VAR.sum.AR(rho,sigma,T)

y.0 <- c/(1-rho) + sqrt(sigma^2/(1-rho^2))*rnorm(1)
y <- simul.ar(c,rho,sigma,T,y.0)
plot(y,type="l")
NW.LongRunVariance(y,q=30)




# Exercise 3 - Debt over GDP
# --------------------------------

data_fred <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/data_fred_quarterly.csv")

# Compute appropriate date format:
data_fred$DATE <- as.Date(data_fred$DATE,"%Y-%m-%d")

# Remove missing variables:
D <- data_fred$GFDEGDQ188S[77:275]
dates <- data_fred$DATE[77:275]

# Length of the sample:
T <- length(D)

# Question 3.a:
var.question.a <- var(D)/T
stdv <- sd(D)/sqrt(T)

# Question 3.b:
plot(dates,D,type="l")
mu <- mean(D,na.rm=TRUE)
abline(h=mu,col="red")
abline(h=mu+2*stdv,col="green")
abline(h=mu-2*stdv,col="green")

# Question 3.c:
NW.LR <- NW.LongRunVariance(D,10)
stdv <- sqrt(NW.LR/T)
abline(h=mu+2*stdv,col="blue")
abline(h=mu-2*stdv,col="blue")
NW.LR <- NW.LongRunVariance(D,100)
stdv <- sqrt(NW.LR/T)
abline(h=mu+2*stdv,col="blue",lty=2)
abline(h=mu-2*stdv,col="blue",lty=2)


# Question 3.e:
gdp <- log(data_fred$GDPC1)
T <- length(gdp)
trend <- 1:T
eq <- lm(gdp~trend)
summary(eq)

plot(data_fred$DATE,eq$residuals,type="l")
abline(h=0,col="red")


