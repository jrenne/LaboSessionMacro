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

M <- rbind(
  c(1,-phi1,-phi2),
  c(-phi1,1-phi2,0),
  c(-phi2,-phi1,1)
)

gamma <- solve(M) %*% c(sigma2,0,0)




for(i in 3:10){
  gamma <- c(gamma,phi1*gamma[i]+phi2*gamma[i-1])
}
plot(0:10,gamma)



# Question 1.c:
c <- 1
mu <- c/(1-phi1-phi2)
N <- 10000
y.sim <- sim.arma(c=0,phi=c(phi1,phi2),theta=1,
                  sigma=sqrt(sigma2),T=100,
                  y.0=c(mu,mu),nb.sim=N)

plot(y.sim[,1],type="l")


n <- 0 # 1 for gamma_1
all.gamma.n <- NULL
for(i in 1:N){
  all.gamma.n <- c(all.gamma.n,autocov(y.sim[,i],n))
}

# The next command plots an histogram with 20 vertical bars:
hist(all.gamma.n,breaks = 20)
abline(v=gamma[n+1],col="red",lwd=4)
abline(v=mean(all.gamma.n),col="blue",lwd=4)


# Question 1.d:

c <- 1
mu <- c/(1-phi1-phi2)
N <- 1000
y.sim <- sim.arma(c=0,phi=c(phi1,phi2),theta=1,
                  sigma=sqrt(sigma2),T=500,
                  y.0=c(mu,mu),nb.sim=N)
n <- 1
all.gamma.n <- NULL
for(i in 1:N){
  all.gamma.n <- c(all.gamma.n,autocov(y.sim[,i],n))
}

# The next command plots an histogram with 20 vertical bars:
hist(all.gamma.n,breaks = 20)
abline(v=gamma[n+1],col="red")


# Exercise 2
# --------------------------------

# Question 2.a:
phi <- c(0.5,0.6,-0.2)
F <- make.F(phi)
print(eigen(F)$values)
print(abs(eigen(F)$values))

# Check spectral decomp:
V <- eigen(F)$vectors
V %*% diag(eigen(F)$values) %*% solve(V)
# Compute the dyn mult of order 10:
V[1,] %*% diag(eigen(F)$values)^10 %*% solve(V)[,1]


phi <- c(0.8,-0.4,0.4)
F <- make.F(phi)
print(eigen(F)$values)
print(abs(eigen(F)$values))

phi <- c(1.3,-1.3,0.2)
F <- make.F(phi)
print(eigen(F)$values)
print(abs(eigen(F)$values))

# Question 2.b:

phi <- c(0.5,0.6,-0.2)
y.sim <- sim.arma(c=0,phi,theta=1,
                  sigma=sqrt(sigma2),T=100,
                  y.0=c(0,0,0),nb.sim=1)
plot(y.sim,type="l")

phi <- c(0.8,-0.4,0.4)
y.sim <- sim.arma(c=0,phi,theta=1,
                  sigma=sqrt(sigma2),T=100,
                  y.0=c(0,0,0),nb.sim=1)
plot(y.sim,type="l")

phi <- c(1.3,-1.3,0.2)
y.sim <- sim.arma(c=0,phi,theta=1,
                  sigma=sqrt(sigma2),T=100,
                  y.0=c(0,0,0),nb.sim=1)
plot(y.sim,type="l")

# Question 2.c:

phi <- c(0.5,0.6,-0.2)
dyn.mult <- make.dyn.mult(phi,max.dyn.mult=40)
plot(dyn.mult,type="l")

phi <- c(0.8,-0.4,-0.4)
dyn.mult <- make.dyn.mult(phi,max.dyn.mult=40)
plot(dyn.mult,type="l")

phi <- c(0.5,0.6,-0.2)
dyn.mult <- make.dyn.mult(phi,max.dyn.mult=40)
plot(dyn.mult,type="l")



# Exercise 3
# --------------------------------

data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession2a.csv")
#data <- read.csv("~/Google Drive/cours/UNIL/Rcode/data/dataLaboSession2a.csv", sep="")

y <- data$V1
plot(y,type="l")

mean(y)
gamma0 <- autocov(y,0)
gamma1 <- autocov(y,1)
gamma2 <- autocov(y,2)
rho <- gamma1/gamma0
discri <- (1/rho)^2 - 4
solutions <- (1/rho + c(1,-1)*sqrt(discri))/2
theta <- solutions[2]
sigma2 <- gamma0/(1+theta^2)


# Exercise 4
# --------------------------------

data <- read.csv("http://jeanpaul.renne.pagesperso-orange.fr/UNIL/dataLaboSession2b.csv")

w <- data$V1
x <- data$V2
y <- data$V3
z <- data$V4
par(mfrow=c(1,2))
acf(w)
pacf(w)
acf(x)
pacf(x)
acf(y)
pacf(y)
acf(z)
pacf(z)


