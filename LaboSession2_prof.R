
setwd("~/Google Drive/cours/UNIL/Rcode")

# For Exercise 3:

y <- sim.arma(c=3,phi=0,theta=c(1,.5),
                  sigma=2,T=200,
                  y.0=0,nb.sim=1)
plot(y,type="l")

mean(y)
autocov(y,0)
autocov(y,1)
autocov(y,2)

write.csv(y,file="data/dataLaboSession2a.csv",row.names = FALSE)


w <- sim.arma(c=3,phi=0,theta=c(1,1),
                   sigma=2,T=200,
                   y.0=0,nb.sim=1)

x <- sim.arma(c=0,phi=c(.3,.4,.2),theta=1,
                   sigma=2,T=200,
                   y.0=c(0,0,0),nb.sim=1)

y <- sim.arma(c=0,phi=c(.3,.4,.2),theta=c(1,1,1),
                   sigma=2,T=200,
                   y.0=c(0,0,0),nb.sim=1)

z <- sim.arma(c=3,phi=0,theta=c(1,1,1,1,1,1),
                   sigma=2,T=200,
                   y.0=0,nb.sim=1)

data <- cbind(w,x,y,z)
write.csv(data,file="data/dataLaboSession2b.csv",row.names = FALSE)

