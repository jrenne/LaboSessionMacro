# ==============================================================================
# Macroeconometrics
# Jean-Paul Renne, 2024
# Laboratory session 1 - Long run variance
# ==============================================================================

# Exercise 1
# --------------------------------

# Question 1.a:
# (replace the "111111111")

autocov <- function(X,n){
  T <- length(X)
  X.1 <- X[1:(T-n)] - mean(X)
  X.2 <- 111111111
  return(111111111)
}

# test:
X <- rnorm(100)
print(autocov(X,2))


# Question 1.b:
# (replace the "111111111")

NW.Weights <- function(q){
  weights <- 111111111
  return(weights)
}

# test:
print(NW.Weights(5))


# Question 1.c:
# (replace the "111111111")

NW.LongRunVariance <- function(X,q){
  111111111
  return(LRV)
}

# test:
T <- 10000
X <- rnorm(T)
q <- 3
print(NW.LongRunVariance(X,q))



# Question 1.d:

data <- read.csv("https://raw.githubusercontent.com/jrenne/LaboSessionMacro/main/Data/dataLaboSession1a.csv")
plot(data$x,type="l")




# Exercise 2
# --------------------------------

# These functions will be used in that exercise:
VAR.sum.AR <- function(rho,sigma,T){
  sum.var <- sigma^2
  rho.k <- 1
  sum.rho <- 0
  for(i in 1:T){
    sum.var <- sum.var + sum.rho^2*sigma^2
    sum.rho <- sum.rho + rho.k
    rho.k <- rho.k * rho
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




# Exercise 3 - Debt over GDP
# --------------------------------

data_fred <- read.csv("https://raw.githubusercontent.com/jrenne/LaboSessionMacro/main/Data/dataFredQuarterlyMScELabSess1.csv")

# Compute appropriate date format:
data_fred$DATE <- as.Date(data_fred$DATE,"%Y-%m-%d")

# Remove missing variables:
indic_available_data <- complete.cases(data_fred)
D     <- data_fred$GFDEGDQ188S[indic_available_data]
dates <- data_fred$DATE[indic_available_data]

# Length of the sample:
T <- length(D)

# Question 3.a:


# Question 3.b:


# Question 3.c:


# Question 3.e:



