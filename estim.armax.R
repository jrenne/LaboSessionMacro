# =============================================================
# Macroeconometrics
# Jean-Paul Renne, 2016
# =============================================================
# This session  requires a preliminary run of various_proc_TS.R
# =============================================================

setwd("~/Google Drive/Teaching/UNIL")


c <- 1
phi <- c(0.8)
theta <- c(1,0.5)
sigma <- 2
beta <- .5

T <- 300

X.sim <- (runif(T)<.2)*rgamma(T,shape=.1,scale=20)
plot(X.sim,type="l")

Y <- sim.arma(c,phi,theta,sigma,T,y.0=rep(0,length(phi)),nb.sim=1,make.IRF=0,
              X=X.sim,beta=beta)

plot(Y,type="l")
lines(X.sim,col="red")

res <- armax.log.L(Y,c,phi,theta,sigma,X=X.sim,beta=beta)
print(sum(res$vec.log.l))

x <- estim.armax(Y,p=length(phi),q=length(theta)-1,X=X.sim)



Ramey <- read.csv("~/Google Drive/Teaching/UNIL/Rcode/data/MonetaryShocks/Ramey2.csv")
Ramey$DATES <- as.Date(Ramey$DATES,"%m/%d/%Y")

T <- dim(Ramey)[1]

# Construct inflation:
Ramey$infl <- Ramey$LCPI - c(rep(NaN,3),Ramey$LCPI[1:(T-3)])

# Construct growth series:
Ramey$growth <- Ramey$LIP - c(rep(NaN,12),Ramey$LIP[1:(length(Ramey$LIP)-12)])


# Prepare matrix of exogenous variables:
vec.lags <- c(9,12,18)
Matrix.of.Exog <- NULL

shocks <- Ramey$ED2_TC

for(i in 1:length(vec.lags)){
  Matrix.of.Exog <- cbind(Matrix.of.Exog,
                          c(rep(NaN,vec.lags[i]),shocks[1:(T-vec.lags[i])]))
}


# Look for dates where data are available:
aux <- apply(Matrix.of.Exog,1,sum)
indic.good.dates <- which(!is.na(aux))


# Estimate ARMAX:
p <- 1
q <- 1
x <- estim.armax(Ramey$growth[indic.good.dates],p,q,X=Matrix.of.Exog[indic.good.dates,])




FILE = "/figures/Figure_Ramey.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=4)

par(mfrow=c(1,1),plt=c(.15,.95,.15,.95))

plot(Ramey$DATES[indic.good.dates],100*Ramey$growth[indic.good.dates],type="l",
     ylim=c(-15,10),xlab="",ylab="in percent",lwd=2)
abline(h=0,col="grey")
lines(Ramey$DATES[indic.good.dates],40*Ramey$ED3_TC[indic.good.dates],col="blue",lwd=2)
#lines(Ramey$DATES[indic.good.dates],1*(Ramey$FFR[indic.good.dates]-Ramey$FFR[indic.good.dates-1]),col="red")

legend("bottomleft", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c("12-month growth rate of IP","Gertler-Karadi shocks (x 40)"),
       lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       col=c("black","blue"), # gives the legend lines the correct color and width
       pt.bg=c(1),
       #text.width = 2,
       #cex=1.0,# size of the text
       #pch = c(1,2,3,4),#symbols,
       pt.cex = c(1),
       bg="white",
       seg.len = 4
)


dev.off()



irf <- sim.arma(0,x$phi,x$beta,x$sigma,T=60,y.0=rep(0,length(x$phi)),nb.sim=1,make.IRF=1,
                X=NaN,beta=NaN)



irf.function <- function(THETA){
  c <- THETA[1]
  phi <- THETA[2:(p+1)]
  if(q>0){
    theta <- c(1,THETA[(1+p+1):(1+p+q)])
  }else{
    theta <- 1
  }
  sigma <- THETA[1+p+q+1]
  r <- dim(Matrix.of.Exog)[2] - 1
  beta <- THETA[(1+p+q+1+1):(1+p+q+1+(r+1))]
  
  irf <- sim.arma(0,phi,beta,sigma=sd(Ramey$ED3_TC,na.rm=TRUE),T=60,y.0=rep(0,length(x$phi)),nb.sim=1,make.IRF=1,
                  X=NaN,beta=NaN)
  return(irf)
}

IRF.0 <- 100*irf.function(x$THETA)
eps <- .00000001
d.IRF <- NULL
for(i in 1:length(x$THETA)){
  THETA.i <- x$THETA
  THETA.i[i] <- THETA.i[i] + eps
  IRF.i <- 100*irf.function(THETA.i)
  d.IRF <- cbind(d.IRF,
                 (IRF.i - IRF.0)/eps
                 )
}

mat.var.cov.IRF <- d.IRF %*% x$I %*% t(d.IRF)


FILE = "/figures/Figure_Ramey2.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=4)

par(mfrow=c(1,1),plt=c(.15,.95,.2,.95))
plot(IRF.0,type="l",
     ylim=c(min(IRF.0-2*sqrt(diag(mat.var.cov.IRF))),max(IRF.0+2*sqrt(diag(mat.var.cov.IRF)))),
     xlab="Number of months after shock",lwd=2,ylab="in percent")
abline(h=0,col="black",lty=3,lwd=2)
lines(IRF.0+2*sqrt(diag(mat.var.cov.IRF)),col="blue",lty=2,lwd=2)
lines(IRF.0-2*sqrt(diag(mat.var.cov.IRF)),col="blue",lty=2,lwd=2)

dev.off()
