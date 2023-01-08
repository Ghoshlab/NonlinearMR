rm(list=ls())
library(MASS)
library(sisVIVE)
library(randomForestSRC)
library(hdm)
library(gbm)
library(NonpModelCheck)
library(boot)
library(naivereg)


est.fun <- function(dat,ind){
  dat <- dat[ind,]
  y <- dat[,1]
  D <- dat[,2]
  X.o <- dat[,3:(p+2)]
  
  
  dat.x <- data.frame(D,X.o)
  
  # naive OLS
  
  naive.lm <- lm(y ~D, data = dat)
  ceff.na.lm <- as.numeric(summary(naive.lm)$coefficients[2,1])
  # 2SLS
  
  fit.x.lm <- lm(D ~., data = dat.x)
  pred.x.lm <- predict(fit.x.lm)
  fit2.lm <- lm(y ~ pred.x.lm)
  ceff.lm <- as.numeric(summary(fit2.lm)$coefficients[2,1])
  
  # use randomforest
  
  fit.x.rf <- rfsrc(D ~., data = dat.x, ntree=1000, samptype = "swr")
  pred.x.rf.list <-predict(fit.x.rf)
  pred.x.rf <- pred.x.rf.list$predicted.oob
  fit.rf <- lm(y ~ pred.x.rf)
  ceff.rf <- as.numeric(summary(fit.rf)$coefficients[2,1])
  
  # variable selection with random forest
  
  fit.x.rf.vs <- var.select(D ~., data = dat.x, refit = TRUE, nodesize=20)
  pred.x.rf.list.vs <-predict(fit.x.rf.vs$rfsrc.refit.obj)
  pred.x.rf.vs <- pred.x.rf.list.vs$predicted.oob
  fit.rf.vs <- lm(y ~ pred.x.rf.vs)
  ceff.rf.vs <- as.numeric(summary(fit.rf.vs)$coefficients[2,1])
  # use gbm
  
  fit.x.xgb <- gbm(D~., data = dat.x, distribution = "gaussian", interaction.depth=10,bag.fraction=0.5,train.fraction=0.5,n.trees = 1000)
  pred.x.xgb <- predict(fit.x.xgb, n.trees = 1000)
  fit.xgb <- lm(y ~ pred.x.xgb)
  ceff.xgb <- as.numeric(summary(fit.xgb)$coefficients[2,1])
  
  
  fit.rlassoIV <- rlassoIV(x = NULL, d = D, y = y, z = X.o, select.X = FALSE, select.Z = TRUE)
  ceff.rlasso <- as.numeric(fit.rlassoIV$coefficients)
  ceff.rlasso.se <- as.numeric(fit.rlassoIV$se)
  
  fit.zh = naivereg(y=y,x=D,z=as.matrix(X.o),intercept=FALSE,criterion='BIC') 
  ceff.zh <- fit.zh$beta.endogenous
  ceff.zh.se <- fit.zh$std.endogenous
  num.list.zh <- rep(0,p)
  num.list.zh[fit.zh$ind] <- 1 


  npvarIV.sel <- npvarselec(as.matrix(X.o), as.vector(D), method = "forward", p = 5, degree.pol = 0, kernel.type = "trun.normal", bandwidth = "CV")	
  fit.npvarIV.lm <- lm(D ~ as.matrix(X.o[,npvarIV.sel$selected]))
  pred.npvarIV.lm <- predict(fit.npvarIV.lm)
  fit.y.npvarIV.lm <- lm(y ~ pred.npvarIV.lm)
  ceff.npvarIV <- as.numeric(summary(fit.y.npvarIV.lm)$coefficients[2,1])
  
  num.list <- match(fit.x.rf.vs$topvars,names(dat.x)[2:(p+1)])
  num.list.s<-sort(num.list)
  num.list.rf <- rep(0,p)
  num.list.rf[num.list.s] <- 1
  
  num.list.rlasso <- rep(0,p)
  if(sum(fit.rlassoIV$selected)!=0){
    num.list.rlasso <- as.numeric(fit.rlassoIV$selected)
  }
  num.list.npar <- rep(0,p)
  num.list.npar[npvarIV.sel$selected] <- 1 
  
  # sisVIVE (Kang et al. 2016)
  
  #fit.sisVIVE <- cv.sisVIVE(y,D,as.matrix(X.o),K=10)
  #ceff.sisVIVE <- fit.sisVIVE$beta
  #num.list.sisVIVE <- rep(0,p)
  #num.list.sisVIVE[fit.sisVIVE$alpha!=0] <- 1
  #return(c(ceff.lm,ceff.rf,ceff.xgb,ceff.rlasso,ceff.npvarIV,ceff.sisVIVE,num.list.rlasso,num.list.npar,num.list.sisVIVE,ceff.rlasso.se))
  return(c(ceff.na.lm,ceff.lm,ceff.rf,ceff.rf.vs,ceff.xgb,ceff.rlasso,ceff.zh,ceff.npvarIV,num.list.rf,num.list.rlasso,num.list.zh,num.list.npar,ceff.rlasso.se,ceff.zh.se))
}

est.fun.boot <- function(dat,ind){
  dat <- dat[ind,]
  y <- dat[,1]
  D <- dat[,2]
  X.o <- dat[,3:(p+2)]
  
  
  dat.x <- data.frame(D,X.o)
  
  # naive OLS
  
  naive.lm <- lm(y ~D, data = dat)
  ceff.na.lm <- as.numeric(summary(naive.lm)$coefficients[2,1])
  
  # 2SLS
  
  fit.x.lm <- lm(D ~., data = dat.x)
  pred.x.lm <- predict(fit.x.lm)
  fit2.lm <- lm(y ~ pred.x.lm)
  ceff.lm <- as.numeric(summary(fit2.lm)$coefficients[2,1])
  
  # use randomforest
  
  fit.x.rf <- rfsrc(D ~., data = dat.x, ntree=1000, samptype = "swr")
  pred.x.rf.list <-predict(fit.x.rf)
  pred.x.rf <- pred.x.rf.list$predicted.oob
  fit.rf <- lm(y ~ pred.x.rf)
  ceff.rf <- as.numeric(summary(fit.rf)$coefficients[2,1])
  
  # variable selection with random forest
  
  fit.x.rf.vs <- var.select(D ~., data = dat.x, refit = TRUE, nodesize=20)
  pred.x.rf.list.vs <-predict(fit.x.rf.vs$rfsrc.refit.obj)
  pred.x.rf.vs <- pred.x.rf.list.vs$predicted.oob
  fit.rf.vs <- lm(y ~ pred.x.rf.vs)
  ceff.rf.vs <- as.numeric(summary(fit.rf.vs)$coefficients[2,1])
  
  # use gbm
  
  fit.x.xgb <- gbm(D~., data = dat.x, distribution = "gaussian", interaction.depth=10,bag.fraction=0.5,train.fraction=0.5,n.trees = 1000)
  pred.x.xgb <- predict(fit.x.xgb, n.trees = 1000)
  fit.xgb <- lm(y ~ pred.x.xgb)
  ceff.xgb <- as.numeric(summary(fit.xgb)$coefficients[2,1])
  
  
  #fit.rlassoIV <- rlassoIV(x = NULL, d = D, y = y, z = X.o, select.X = FALSE, select.Z = TRUE)
  #ceff.rlasso <- as.numeric(fit.rlassoIV$coefficients)
  #ceff.rlasso.se <- as.numeric(fit.rlassoIV$se)
  
  # naivereg (Fan and Zhong, 2018) 
  #fit.naivereg <- naivereg(y=y,x = D,z = X.o,intercept=FALSE,criterion='BIC',IV.intercept = FALSE)


  npvarIV.sel <- npvarselec(as.matrix(X.o), as.vector(D), method = "forward", p = 5, degree.pol = 0, kernel.type = "trun.normal", bandwidth = "CV")	
  fit.npvarIV.lm <- lm(D ~ as.matrix(X.o[,npvarIV.sel$selected]))
  pred.npvarIV.lm <- predict(fit.npvarIV.lm)
  fit.y.npvarIV.lm <- lm(y ~ pred.npvarIV.lm)
  ceff.npvarIV <- as.numeric(summary(fit.y.npvarIV.lm)$coefficients[2,1])
  
    
  return(c(ceff.na.lm,ceff.lm,ceff.rf,ceff.rf.vs,ceff.xgb,ceff.npvarIV))
}

set.seed(1)

p <- 100
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
gamma <- c(0.2,-0.4,1)

eta <- 0
zeta <- 0.5
beta <- 0.5
xi1 <- 0.2
xi2 <- -0.3
U <- rnorm(n)

X <- cbind(sqrt(abs(X.o[,1])),exp(X.o[,7]),sin(0.5*pi*X.o[,8]*X.o[,9]))


D = zeta + X %*% gamma  + xi1*abs(U) + eps[,2]
y = eta + beta * D + xi2*U^2 + eps[,1]
dat <- data.frame(y=y, D=D, X.o)
ind <- 1:n


est.all <- est.fun(dat,ind)

ceff.naive <- est.all[1]
ceff.lm <- est.all[2]
ceff.rf <- est.all[3]
ceff.rf.vs <- est.all[4]
ceff.xgb <- est.all[5]
ceff.rlasso <- est.all[6]
ceff.zh <- est.all[7]
ceff.npvarIV <- est.all[8]


rf.sel <- est.all[9:(8+p)]
lasso.sel <- est.all[(p+9):(8+2*p)]
zh.sel <- est.all[(2*p+9):(8+3*p)]
npar.sel <- est.all[(3*p+9):(8+4*p)]
rlasso.se <- est.all[(9+4*p)]
zh.se <- est.all[(10+4*p)]
est.all.boot <- boot(dat,est.fun.boot,R=100)
se.all.boot <- apply(est.all.boot[[2]],2,sd)

res.all <- cbind(t(est.all[1:8]),t(c(se.all.boot[1:5],rlasso.se,zh.se,se.all.boot[6])),t(rf.sel),t(lasso.sel),t(zh.sel),t(npar.sel))
res.seed <- paste('res.',seednum,'.csv',sep = '')
write.csv(res.all,file = res.seed)
