# =============================================================
# Macroeconometrics
# Laboratory session 2 - ARMA Processes
# Jean-Paul Renne, 2016
# =============================================================
# This session  requires a preliminary run of various_proc_TS.R
# =============================================================


# Exercise 1
# --------------------------------

# Question 1.b:

phi1 <- 0.9
phi2 <- -0.2
sigma2 <- 1

M <- matrix(
  c(1,phi1,phi2,-phi1,phi2-1,phi1,-phi2,0,-1),3,3)

vector.gamma012 <- solve(M) %*% matrix(c(sigma2,0,0),3,1)

vector.gamma <- vector.gamma012

gamma.j_1 <- vector.gamma012[3] # when j=3, gamma.j_1 has
# to be equal to gamma.2
gamma.j_2 <- vector.gamma012[2]
for(j in 3:10){
  gamma.j <- phi1*gamma.j_1 + phi2*gamma.j_2
  gamma.j_2 <- gamma.j_1
  gamma.j_1 <- gamma.j
  vector.gamma <- c(vector.gamma,gamma.j)
}

plot(0:10,vector.gamma)

# Question 1.c:
# replace the "111111"
c <- 0
mu <- c/(1-phi1-phi2)
N <- 1000
y.sim <- sim.arma(c=c,phi=c(phi1,phi2),theta=1,
                  sigma=sqrt(sigma2),T=400,
                  y.0=c(mu,mu),nb.sim=N)

plot(y.sim[,2],type="l")

n <- 2 # order of the autocovariance
all.gamma.n <- NULL
for(i in 1:N){
  all.gamma.n <- c(all.gamma.n,autocov(y.sim[,i],n))
}

# The next command plots an histogram with 20 vertical bars:
hist(all.gamma.n,breaks = 20)
abline(v=vector.gamma[n+1],col="red",lwd=3)

# Question 1.d:




# Exercise 2
# --------------------------------

# Question 2.a:

# Question 2.b:

# Question 2.c:



# Exercise 3
# --------------------------------

data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession2a.csv")
y <- data$V1
plot(y,type="l")



# Exercise 4
# --------------------------------

data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession2b.csv")
w <- data$V1
x <- data$V2
y <- data$V3
z <- data$V4

