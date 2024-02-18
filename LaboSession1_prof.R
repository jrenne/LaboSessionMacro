
setwd("~/Google Drive/cours/UNIL/Rcode")

rho <- .85
c <- 1
sigma <- 1
x <- c/(1-rho) + sqrt(sigma^2/(1-rho^2))*rnorm(1)
T <- 1000
for(i in 2:T){
  x <- c(x,c + rho*x[length(x)] + sigma * rnorm(1))
}
plot(x,type="l")
write.csv(x,file="data/dataLaboSession1a.csv",row.names = FALSE)

print(VAR.sum.AR(rho,sigma,T))
print(NW.LongRunVariance(x,100))
print(NW.LongRunVariance(x,1))

