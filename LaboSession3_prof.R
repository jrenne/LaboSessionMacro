
setwd("~/Google Drive/cours/UNIL/Rcode")

# Session on the estimation of ARMA processes.


# For Exercise 1:

y <- sim.arma(c=3,phi=0.8,theta=c(1),
              sigma=2,T=300,
              y.0=3/(1-.8),nb.sim=1)
plot(y,type="l")

mean(y)
autocov(y,0)
autocov(y,1)
autocov(y,2)

write.csv(y,file="data/dataLaboSession3a.csv",row.names = FALSE)

y <- sim.arma(c=1,phi=c(0.4,0,0.4),theta=c(1,1),
              sigma=2,T=500,
              y.0=rep(1/(1-.8),3),nb.sim=1)
plot(y,type="l")

mean(y)
autocov(y,0)
autocov(y,1)
autocov(y,2)

write.csv(y,file="data/dataLaboSession3b.csv",row.names = FALSE)

