
setwd("~/Google Drive/cours/UNIL/Rcode")

# Session on the estimation of VAR models.

First.date <- "1950-01-01"
Last.date  <- "2015-01-01"

# Load data from txt files (extracted fromn the FFED database)

QVAR_Monthly <- read.delim("~/Google Drive/MyPapers/ICA/Rcodes/data/QVAR/QVAR_Monthly.txt")
QVAR_Monthly[3:dim(QVAR_Monthly)[1],2:dim(QVAR_Monthly)[2]] <- 
  1/3 *(
    QVAR_Monthly[1:(dim(QVAR_Monthly)[1]-2),2:dim(QVAR_Monthly)[2]] +
      QVAR_Monthly[2:(dim(QVAR_Monthly)[1]-1),2:dim(QVAR_Monthly)[2]] +
      QVAR_Monthly[3:(dim(QVAR_Monthly)[1]-0),2:dim(QVAR_Monthly)[2]]
  )

QVAR_Quarterly <- read.delim("~/Google Drive/MyPapers/ICA/Rcodes/data/QVAR/QVAR_Quarterly.txt")

QVAR_Monthly$DATE <- as.Date(QVAR_Monthly$DATE,"%Y-%m-%d")
QVAR_Quarterly$DATE <- as.Date(QVAR_Quarterly$DATE,"%Y-%m-%d")

QVAR_Monthly$month <- as.integer(format(QVAR_Monthly$DATE,"%m"))
QVAR_Monthly <- subset(QVAR_Monthly,month %in% c(1,4,7,10))

DATA <- merge(QVAR_Monthly,QVAR_Quarterly,by="DATE",all=TRUE)


# Select first and last dates of the estimation sample:

first.date <- which(DATA$DATE==as.Date(First.date,"%Y-%m-%d"))
last.date <- which(DATA$DATE==as.Date(Last.date,"%Y-%m-%d"))

# first.date <- which(DATA$DATE==as.Date("1979-10-01","%Y-%m-%d"))
# last.date <- which(DATA$DATE==as.Date("2002-4-01","%Y-%m-%d"))

# first.date <- which(DATA$DATE==as.Date("1980-01-01","%Y-%m-%d"))
# last.date <- which(DATA$DATE==as.Date("2007-4-01","%Y-%m-%d"))

# first.date <- 10
# last.date <- 250


# Prepare data:

# Output:
gdp <- log(DATA$GDPC1[first.date:last.date])
trend <- 1:length(gdp)
trend2 <- trend^2
trend3 <- trend^3
eq <- lm(gdp ~ trend)
y.gdp <- 100*eq$residuals

gdp.pot <- log(DATA$GDPPOT[first.date:last.date])
y.gdp.gap <- 100*(gdp - gdp.pot)

plot(y.gdp.gap,type="l")

# Unemployment
u.gap <- DATA$UNRATE[first.date:last.date] - DATA$NROU[first.date:last.date]
y.u.gap <- u.gap

# Choice of real activity variable:
y <- y.gdp.gap

# Inflation:
p <- log(DATA$GDPDEF)
infl <- 400 * (p - c(NaN,p[1:(length(p)-1)]))
#infl <- 100 * (p - c(rep(NaN,3),p[1:(length(p)-3)]))
infl <- infl[first.date:last.date]

# Commodities:
aux <- log(DATA$OILPRICE)
commo <- 400 * (aux - c(rep(NaN,1),aux[1:(length(aux)-1)]))
commo <- 100 * (aux - c(rep(NaN,4),aux[1:(length(aux)-4)]))
#commo <- aux
commo_1 <- commo[(first.date-1):(last.date-1)]
commo_2 <- commo[(first.date-2):(last.date-2)]
commo <- commo[first.date:last.date]

# Fed fund rate:
r <- DATA$FEDFUNDS[first.date:last.date]

vec.dates <- DATA$DATE[first.date:last.date]

par(mfrow=c(2,2))
par(plt=c(.1,.9,.15,.7))
plot(vec.dates,y,type="l",main="Real activity")
lines(vec.dates,y.gdp.gap,col="red")
abline(h=0,col="grey")
plot(vec.dates,infl,type="l",main="Inflation")
abline(h=0,col="grey")
plot(vec.dates,r,type="l",main="Federal fund rate")
plot(vec.dates,commo,type="l",main="Commodity prices")
abline(h=0,col="grey")


# Prepare database for Blanchard-Quah paper:
dataBQ <- data.frame(dates=DATA$DATE,u=DATA$UNRATE, gdp=DATA$GDPC1)
write.csv(dataBQ,file="data/dataLaboSession4_BQ.csv",row.names = TRUE)

library(RCurl)
#ftpUpload("data/dataLaboSession4_BQ.csv", "ftp://renne.jean-paul@wanadoo.fr:rgfd7tu@perso-ftp.orange.fr/UNIL/dataLaboSession4_BQ.csv")
#jeanpaul.renne.pagesperso-orange.fr/UNIL/

First.date <- "1948-01-01"
Last.date  <- "1988-01-01"
Mid.date <- "1974-01-01"

fst.date.indic <- which(dataBQ$dates==First.date)
lst.date.indic <- which(dataBQ$dates==Last.date)
mid.date.indic <- which(dataBQ$dates==Mid.date)

first.sample.part <- 1:(mid.date.indic-fst.date.indic)
secnd.sample.part <- (mid.date.indic-fst.date.indic + 1):(lst.date.indic-fst.date.indic+1)

dataBQ <- dataBQ[fst.date.indic:lst.date.indic,]

#attach(dataBQ)
y <- log(dataBQ$gdp)
y_1 <- c(NaN,log(dataBQ$gdp[1:(length(dataBQ$gdp)-1)]))
delta.Y <- (y - y_1)*100
delta.Y[first.sample.part] <- delta.Y[first.sample.part] - mean(delta.Y[first.sample.part],na.rm=TRUE)
delta.Y[secnd.sample.part] <- delta.Y[secnd.sample.part] - mean(delta.Y[secnd.sample.part],na.rm=TRUE)
par(mfrow=c(2,1))
plot(dataBQ$dates,delta.Y,type="l")
plot(dataBQ$dates,dataBQ$u,type="l")

u <- dataBQ$u
trend <- 1:length(u)
eq <- lm(u ~ trend)
u <- eq$residuals

y <- cbind(delta.Y,u)
y <- y[2:dim(y)[1],]

est.VAR <- VAR(y,p=8)

Phi <- Acoef(est.VAR)
PHI <- make.PHI(Phi)

sum.PHI.k <- diag(dim(PHI)[1])
PHI.k <- PHI
for(i in 1:200){
  sum.PHI.k <- sum.PHI.k +  PHI.k
  PHI.k <- PHI.k %*% PHI
}

resids <- residuals(est.VAR)
Omega <- var(resids)

f <- function(param){
  B <- matrix(param,2,2)
  X <- Omega - B %*% t(B)
  Theta <- sum.PHI.k[1:2,1:2] %*% B
  loss <- 10000 *( X[1,1]^2 + X[2,1]^2 + X[2,2]^2 + Theta[1,1]^2 )
  return(loss)
}

param0 <- c(1,0,0,1)
res.opt <- optim(param0,f,
                 #method="Nelder-Mead",
                 method="BFGS",
                 # method="CG",
                 control=list(trace=TRUE,maxit=10000),
                 hessian=FALSE)
print(res.opt$par)
B.hat <- matrix(res.opt$par,2,2)
print(cbind(Omega,B.hat %*% t(B.hat)))

nb.sim <- 40
par(mfrow=c(2,2))
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(1,0))
plot(cumsum(Y[,1]),type="l")
plot(Y[,2],type="l")
Y <- simul.VAR(c,Phi,B.hat,nb.sim,y0.star,indic.IRF = 1,u.shock = c(0,1))
plot(cumsum(Y[,1]),type="l")
plot(Y[,2],type="l")








