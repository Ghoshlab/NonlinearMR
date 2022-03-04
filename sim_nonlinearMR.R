rm(list=ls())
library(MASS)
library(sisVIVE)
library(randomForestSRC)
library(hdm)
library(NonpModelCheck)
library(boot)


est.fun <- function(dat,ind){
dat <- dat[ind,]
y <- dat[,1]
D <- dat[,2]
X.o <- dat[,3:22]


dat.x <- data.frame(D,X.o)

# 2SLS

fit.x.lm <- lm(D ~., data = dat.x)
pred.x.lm <- predict(fit.x.lm)
fit2.lm <- lm(y ~ pred.x.lm)
ceff.lm <- as.numeric(summary(fit2.lm)$coefficients[2,1])


tune.x <- tune(D ~., data=dat.x, samptype = "swr")
fit.var.sel.rf <- var.select(D ~., data=dat.x,nodesize=tune.x$optimal[1],mtry = tune.x$optimal[2],conservative = "medium", ntree=500, samptype = "swr")
# use randomforest

fit.x.rf <- rfsrc(D ~., data = dat.x, mtry = tune.x$optimal[2], nodesize=tune.x$optimal[1], ntree=500, samptype = "swr")
num.list <- match(fit.var.sel.rf$topvars,names(dat.x)[2:21])
num.list.s<-sort(num.list)

fit.x.rf.sel.lm <- lm(D ~ as.matrix(X.o[,num.list.s]))
pred.x.rf.sel.lm <- predict(fit.x.rf.sel.lm)
fit.y.rf.sel.lm <- lm(y ~ pred.x.rf.sel.lm)
ceff.rf.lm <- as.numeric(summary(fit.y.rf.sel.lm)$coefficients[2,1])


pred.x.rf.list <-predict(fit.x.rf)
pred.x.rf <- pred.x.rf.list$predicted.oob
fit2.rf <- lm(y ~ pred.x.rf)
ceff.rf <- as.numeric(summary(fit2.rf)$coefficients[2,1])




fit.rlassoIV <- rlassoIV(x = NULL, d = D, y = y, z = X.o, select.X = FALSE, select.Z = TRUE)
ceff.rlasso <- as.numeric(fit.rlassoIV$coefficients)
ceff.rlasso.se <- as.numeric(fit.rlassoIV$se)

npvarIV.sel <- npvarselec(as.matrix(X.o), as.vector(D), method = "forward", p = 5, degree.pol = 0, kernel.type = "trun.normal", bandwidth = "CV")	
fit.npvarIV.lm <- lm(D ~ as.matrix(X.o[,npvarIV.sel$selected]))
pred.npvarIV.lm <- predict(fit.npvarIV.lm)
fit.y.npvarIV.lm <- lm(y ~ pred.npvarIV.lm)
ceff.npvarIV <- as.numeric(summary(fit.y.npvarIV.lm)$coefficients[2,1])

num.list.s.rf <- rep(0,p)
num.list.s.rf[num.list.s] <- 1
num.list.rlasso <- rep(0,p)
if(sum(fit.rlassoIV$selected)!=0){
num.list.rlasso <- as.numeric(fit.rlassoIV$selected)
}
num.list.npar <- rep(0,p)
num.list.npar[npvarIV.sel$selected] <- 1 

# sisVIVE (Kang et al. 2016)

fit.sisVIVE <- cv.sisVIVE(y,D,as.matrix(X.o),K=10)
ceff.sisVIVE <- fit.sisVIVE$beta
num.list.sisVIVE <- rep(0,p)
num.list.sisVIVE[fit.sisVIVE$alpha!=0] <- 1
return(c(ceff.lm,ceff.rf.lm,ceff.rf,ceff.rlasso,ceff.npvarIV,num.list.s.rf,num.list.rlasso,num.list.npar,ceff.rlasso.se,ceff.sisVIVE,num.list.sisVIVE))
}

est.fun.boot <- function(dat,ind){
dat <- dat[ind,]
y <- dat[,1]
D <- dat[,2]
X.o <- dat[,3:22]


dat.x <- data.frame(D,X.o)

# 2SLS

fit.x.lm <- lm(D ~., data = dat.x)
pred.x.lm <- predict(fit.x.lm)
fit2.lm <- lm(y ~ pred.x.lm)
ceff.lm <- as.numeric(summary(fit2.lm)$coefficients[2,1])


tune.x <- tune(D ~., data=dat.x, samptype = "swr")
fit.var.sel.rf <- var.select(D ~., data=dat.x,nodesize=tune.x$optimal[1],mtry = tune.x$optimal[2],conservative = "medium", ntree=500, samptype = "swr")
# use randomforest

fit.x.rf <- rfsrc(D ~., data = dat.x, mtry = tune.x$optimal[2], nodesize=tune.x$optimal[1], ntree=500, samptype = "swr")
num.list <- match(fit.var.sel.rf$topvars,names(dat.x)[2:21])
num.list.s<-sort(num.list)
fit.x.rf.sel.lm <- lm(D ~ as.matrix(X.o[,num.list.s]))
pred.x.rf.sel.lm <- predict(fit.x.rf.sel.lm)
fit.y.rf.sel.lm <- lm(y ~ pred.x.rf.sel.lm)
ceff.rf.lm <- as.numeric(summary(fit.y.rf.sel.lm)$coefficients[2,1])


pred.x.rf.list <-predict(fit.x.rf)
pred.x.rf <- pred.x.rf.list$predicted.oob
fit2.rf <- lm(y ~ pred.x.rf)
ceff.rf <- as.numeric(summary(fit2.rf)$coefficients[2,1])


fit.rlassoIV <- rlassoIV(x = NULL, d = D, y = y, z = X.o, select.X = FALSE, select.Z = TRUE)
ceff.rlasso <- as.numeric(fit.rlassoIV$coefficients)
ceff.rlasso.se <- as.numeric(fit.rlassoIV$se)

npvarIV.sel <- npvarselec(as.matrix(X.o), as.vector(D), method = "forward", p = 5, degree.pol = 0, kernel.type = "trun.normal", bandwidth = "CV")	
fit.npvarIV.lm <- lm(D ~ as.matrix(X.o[,npvarIV.sel$selected]))
pred.npvarIV.lm <- predict(fit.npvarIV.lm)
fit.y.npvarIV.lm <- lm(y ~ pred.npvarIV.lm)
ceff.npvarIV <- as.numeric(summary(fit.y.npvarIV.lm)$coefficients[2,1])

num.list.s.rf <- rep(0,p)
num.list.s.rf[num.list.s] <- 1
num.list.rlasso <- rep(0,p)
if(sum(fit.rlassoIV$selected)!=0){
num.list.rlasso <- as.numeric(fit.rlassoIV$selected)
}
num.list.npar <- rep(0,p)
num.list.npar[npvarIV.sel$selected] <- 1 
fit.sisVIVE <- cv.sisVIVE(y,D,as.matrix(X.o),K=10)
ceff.sisVIVE <- fit.sisVIVE$beta
num.list.sisVIVE <- rep(0,p)
num.list.sisVIVE[fit.sisVIVE$alpha!=0] <- 1

return(c(ceff.lm,ceff.rf.lm,ceff.rf,ceff.rlasso,ceff.npvarIV,ceff.sisVIVE))
}



seednum <- 1
seednum

set.seed(seednum)

p <- 20
rhox <- 0.75
n <- 200
Sigx <- matrix(0,p,p)
for(j in 1:p){
   for(k in 1:p){
	if(j != k){
      	Sigx[j,k] <- rhox^(abs(j-k))
	}else{
		Sigx[j,k] <- 1
	}
   }
}
rho = 0.5
X.o <- mvrnorm(n,mu = rep(0,p), Sigma=Sigx)
Sig <- matrix(c(1,rho,rho,1),2,2)
eps <- mvrnorm(n, mu = c(0,0), Sigma = Sig)
gamma <- c(0.2,-0.4,0.3)
eta <- 0
zeta <- 0.5
beta <- 0.5
xi1 <- 0.2
xi2 <- -0.3
U <- rnorm(n)
X <- cbind(X.o[,3]^3, X.o[,6]^2, sin(1/2*pi*X.o[,9]))

D = zeta + X %*% gamma  + xi1*abs(U) + eps[,2]
y = eta + beta * D + + xi2*U^2 + eps[,1]


dat <- data.frame(y=y, D=D, X1=X.o[,1],X2=X.o[,2],X3=X.o[,3],X4=X.o[,4],X5=X.o[,5],X6=X.o[,6],X7=X.o[,7],X8=X.o[,8],X9=X.o[,9],X10=X.o[,10],X11=X.o[,11],X12=X.o[,12],X13=X.o[,13],X14=X.o[,14],X15=X.o[,15],X16=X.o[,16],X17=X.o[,17],X18=X.o[,18],X19=X.o[,19],X20=X.o[,20])
ind <- 1:n
est.all <- est.fun(dat,ind)

ceff.lm <- est.all[1]
ceff.rf.lm <- est.all[2]
ceff.rf <- est.all[3]
ceff.rlasso <- est.all[4]
ceff.npvarIV <- est.all[5]
ceff.sisVIVE <- est.all[67]

rf.sel <- est.all[6:25]
lasso.sel <- est.all[26:45]
npar.sel <- est.all[46:65]
rlasso.se <- est.all[66]
sisVIVE.sel <- est.all[68:87]
est.all.boot <- boot(dat,est.fun.boot,R=100)
se.all.boot <- apply(est.all.boot[[2]],2,sd)

res.all <- cbind(t(est.all[1:5]),t(se.all.boot),t(est.all[6:25]),t(est.all[26:45]),t(est.all[46:65]),est.all[66],est.all[67],t(est.all[68:87]))
res.seed <- paste('res.',seednum,'.csv',sep = '')
write.csv(res.all,file = res.seed)