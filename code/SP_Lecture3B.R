library(ggplot2)

#
####
####### example kernels
####
#

# brownian motion
x <- 1:100
K <- matrix(0,100,100)
for(i in 1:100) {
  for(j in 1:100) {
    K[i,j] <- min(i,j)
  }
}
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:20){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}

df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Brownian Motions") +
  geom_line() + theme_bw() + theme(legend.position = "none")
gg

ggsave(gg,filename = "~/Desktop/brownianMotions.pdf",width =5,height=3)


# linear kernel
x <- 1:100
K <- x%*%t(x)
eig.obj <- eigen(K)
vals <- c(sqrt(eig.obj$values[1]),rep(0,99))
sqrtK <- eig.obj$vectors %*% diag(vals) %*% t(eig.obj$vectors)
bm <- sqrtK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:20){
  bm <- sqrtK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}

df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Linear Kernel") +
  geom_line() + theme_bw() + theme(legend.position = "none")
gg

ggsave(gg,filename = "~/Desktop/linearKernel.pdf",width =5,height=3)


# linear kernel 2
x <- 1:100
X <- cbind(10,x)
K <- X%*%t(X)
eig.obj <- eigen(K)
vals <- c(sqrt(eig.obj$values[1:2]),rep(0,98))
sqrtK <- eig.obj$vectors %*% diag(vals) %*% t(eig.obj$vectors)
bm <- sqrtK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:20){
  bm <- sqrtK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}

df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Linear Kernel with Intercept") +
  geom_line() + theme_bw() + theme(legend.position = "none")
gg

ggsave(gg,filename = "~/Desktop/linearKernel2.pdf",width =5,height=3)



# linear kernel 2
x <- 1:100
X <- cbind(10,x)
K <- X%*%t(X) + 100*diag(100)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:20){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Linear Kernel with Intercept and Noise Term") +
  geom_line() + theme_bw() + theme(legend.position = "none")
gg

ggsave(gg,filename = "~/Desktop/linearKernel3.pdf",width =5,height=3)



# expo kernel 2
x <- 1:100
X <- matrix(x,100,100)
K <- exp(-abs(X-t(X)))
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Exponential Kernel") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-6,6)) 
gg

# sigma^2 = 5
x <- 1:100
X <- matrix(x,100,100)
K <- 16*exp(-abs(X-t(X)))
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg2 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Exponential Kernel (bigger sigma)") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-6,6)) 
gg2

# rho = 20
x <- 1:100
X <- matrix(x,100,100)
K <- exp(-abs(X-t(X))/20)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg3 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("Exponential Kernel (bigger rho)") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-6,6)) 
gg3


library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,nrow=1),filename = "~/Desktop/expoKernel.pdf",width =10,height=3)

#############################################################################
# let nu change

# expo kernel 
x <- 1:100
X <- matrix(x,100,100)
K <- exp(-abs(X-t(X))/3)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("nu=1/2") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-3,3)) 
gg

# sigma^2 = 5
x <- 1:100
X <- matrix(x,100,100)
K <- (1+ sqrt(3)*abs(X-t(X))/3) * exp(-sqrt(3)*abs(X-t(X))/3)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg2 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("nu=3/2") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-3,3)) 
gg2

x <- 1:100
X <- matrix(x,100,100)
K <- exp(-(X-t(X))^2/(2*3^2))
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
df <- cbind(bm,1,x)

for(i in 2:5){
  bm <- cholK%*%rnorm(100)
  df <- rbind(df, cbind(bm,i,x))
}


df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

gg3 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + ggtitle("nu=infinity") +
  geom_line() + theme_bw() + theme(legend.position = "none") + ylim(c(-3,3)) 
gg3


library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,nrow=1),filename = "~/Desktop/maternKernels.pdf",width =10,height=3)




#############################################################################

# prediction

x <- 1:100
X <- matrix(x,100,100)
K <- exp(-(X-t(X))^2/(10^2*2)) + 0.00001*diag(100)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
bm <- bm[1:50]
df <- cbind(bm,1,1:50)

for(i in 1:50) {
  pred <- mvtnorm::rmvnorm(n=1,
                           mean= K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm,
                           sigma = K[51:100,51:100] - K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% t(K[51:100,1:50]))
  df <- rbind(df,cbind(t(pred),rep(i,50),51:100))
}



df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

df2 <- data.frame(x=51:100,bm=K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm)

gg <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") +
  geom_line() + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept = 51) +
  geom_line(data = df2, aes(x=x,y=bm),size=1.1,inherit.aes = FALSE) + ylim(c(-4,4))
gg

x <- 1:100
X <- matrix(x,100,100)
K <- exp(-(X-t(X))^2/(10^2*2)) + 0.00001*diag(100)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
bm <- bm[1:50]
df <- cbind(bm,1,1:50)

for(i in 1:50) {
  pred <- mvtnorm::rmvnorm(n=1,
                           mean= K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm,
                           sigma = K[51:100,51:100] - K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% t(K[51:100,1:50]))
  df <- rbind(df,cbind(t(pred),rep(i,50),51:100))
}



df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

df2 <- data.frame(x=51:100,bm=K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm)

gg2 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + 
  geom_line() + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept = 51) +
  geom_line(data = df2, aes(x=x,y=bm),size=1.1,inherit.aes = FALSE) + ylim(c(-4,4))
gg2

x <- 1:100
X <- matrix(x,100,100)
K <- exp(-(X-t(X))^2/(10^2*2)) + 0.00001*diag(100)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(100)
bm <- bm[1:50]
df <- cbind(bm,1,1:50)

for(i in 1:50) {
  pred <- mvtnorm::rmvnorm(n=1,
                           mean= K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm,
                           sigma = K[51:100,51:100] - K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% t(K[51:100,1:50]))
  df <- rbind(df,cbind(t(pred),rep(i,50),51:100))
}



df <- as.data.frame(df)
colnames(df) <- c("bm","id","x")
df$id <- factor(df$id)

df2 <- data.frame(x=51:100,bm=K[51:100,1:50] %*% solve(K[1:50,1:50]) %*% bm)

gg3 <- ggplot(df,aes(x=x,y=bm,color=id)) + ylab("y") + 
  geom_line() + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept = 51) +
  geom_line(data = df2, aes(x=x,y=bm),size=1.1,inherit.aes = FALSE) + ylim(c(-4,4))
gg3

library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,nrow=1),filename = "~/Desktop/predictions.pdf",width =10,height=3)

################################################################################

#
### inference using MH
#

covMat <- function(x,sigma2,rho2,tau2) {
  N <- length(x)
  K <- matrix(0,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      K[i,j] <- sigma2*exp(-0.5*(x[i]-x[j])^2/rho2)
      if(i==j) K[i,j] <- K[i,j] + tau2
    }
  }
  return(K)
}

logLik <- function(f,x,sigma2,rho2,tau2) {
  K <- covMat(x,sigma2,rho2,tau2)
  N <- length(x)
  
  ll <- -1/2*log(det(K)) - t(f)%*%solve(K)%*%f/2
  return(ll)
}

log_priors <- function(sigma2,rho2) {
  # standard exponential priors (or whatever)
  return((-sigma2-rho2)/10)
}

target <- function(f,x,theta) {
  sigma2 <- theta[1]
  rho2   <- theta[2]
  tau2   <- 0.000001
  lp <- logLik(f,x,sigma2,rho2,tau2) + log_priors(sigma2,rho2)
  return(lp)
}

delta <- function(n) {
  return( min(0.01,n^(-0.5)) )
}

metropolis_hastings <- function(maxIts,D,f,x,targetAccept=0.4,propSD=1) {
  chain <- matrix(0,maxIts,D) 
  chain[1,] <- 5

  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  
  for(s in 2:maxIts) {
    thetaStar <- truncnorm::rtruncnorm(n=D, a=0, b=Inf, mean = chain[s-1,], sd = propSD)
    u         <- runif(1)
    
    logA      <- target(f,x,thetaStar) - target(f,x,chain[s-1,]) + # targets
      sum(log(truncnorm::dtruncnorm(x=chain[s-1,],a=0, mean=thetaStar,sd=propSD))) -
      sum(log(truncnorm::dtruncnorm(x=thetaStar,a=0, mean=chain[s-1,],sd=propSD)))
    
    if(log(u) < logA) {
      chain[s,] <- thetaStar
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    } else {
      chain[s,] <- chain[s-1,]
    }
    SampCount <- SampCount + 1
    
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      cat(AcceptRatio,"\n")
      if ( AcceptRatio > targetAccept ) {
        propSD <- propSD * (1 + delta(i-1))
      } else {
        propSD <- propSD * (1 - delta(i-1))
      }
      
      SampCount <- 0
      Acceptances <- 0
    }
    
    if(s %% 1000 == 0) cat(s,"\n")
  }
  
  return(chain)
}

# simulate data (sigma2=0.1,rho2=16,tau2=0.01)
N <- 10
x <- runif(N)*100
X <- matrix(x,N,N)
K <- 2*exp(-(X-t(X))^2/(4^2*2)) + 0.000001*diag(N)
cholK <- t(chol(K))
bm <- cholK%*%rnorm(N)
plot(x,bm)

# get posterior samples
samps <- metropolis_hastings(maxIts = 100000, D=2,f=bm,x=x,propSD = 1)
plot(samps[10001:100000,1],type="l") # 2
abline(h=2,col="red")
quantile(samps[10001:100000,1],p=c(0.025,0.975))
plot(samps[10001:100000,2],type="l") # 16
abline(h=16,col="red")
quantile(samps[10001:100000,2],p=c(0.025,0.975))

#
### plot 10 posterior predictive curves
#
L <- 300
thin <- ceiling(seq(from=10001,to=100000,length.out = 5))
xStar <- seq(from=0,to=100,length.out=L) # unobserved grid points

K <- covMat(c(x,xStar),sigma2=samps[thin[1],1], rho2=samps[thin[1],2], tau2=0.000001)
pred <- mvtnorm::rmvnorm(n=1,
                         mean= K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% bm,
                         sigma = K[(N+1):(L+N),(N+1):(L+N)] - K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% t(K[(N+1):(L+N),1:N]))
df <- cbind(as.numeric(pred),1,xStar)

for(i in 2:5) {
  K <- covMat(c(x,xStar),sigma2=samps[thin[i],1], rho2=samps[thin[i],2], tau2=0.000001)
  pred <- mvtnorm::rmvnorm(n=1,
                           mean= K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% bm,
                           sigma = K[(N+1):(L+N),(N+1):(L+N)] - K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% t(K[(N+1):(L+N),1:N]))
  df <- rbind(df,cbind(as.numeric(pred),i,xStar))
}
df <- as.data.frame(df)
colnames(df) <- c("y","id","x")
df$id <- factor(df$id)

df2 <- data.frame(y=bm,x=x)

gg <- ggplot(df,aes(x=x,y=y,color=id)) +
  geom_line() + theme_bw() + theme(legend.position = "none") +
  geom_point(data=df2,aes(x=x,y=y),inherit.aes = FALSE,size=2)

gg




L <- 300
thin <- ceiling(seq(from=10001,to=100000,length.out = 100))
xStar <- seq(from=0,to=100,length.out=L) # unobserved grid points

K <- covMat(c(x,xStar),sigma2=samps[thin[1],1], rho2=samps[thin[1],2], tau2=0.000001)
pred <- mvtnorm::rmvnorm(n=1,
                         mean= K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% bm,
                         sigma = K[(N+1):(L+N),(N+1):(L+N)] - K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% t(K[(N+1):(L+N),1:N]))
df <- cbind(as.numeric(pred),1,xStar)

for(i in 2:100) {
  K <- covMat(c(x,xStar),sigma2=samps[thin[i],1], rho2=samps[thin[i],2], tau2=0.000001)
  pred <- mvtnorm::rmvnorm(n=1,
                           mean= K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% bm,
                           sigma = K[(N+1):(L+N),(N+1):(L+N)] - K[(N+1):(L+N),1:N] %*% solve(K[1:N,1:N]) %*% t(K[(N+1):(L+N),1:N]))
  df <- rbind(df,cbind(as.numeric(pred),i,xStar))
}
df <- as.data.frame(df)
colnames(df) <- c("y","id","x")
df$id <- factor(df$id)

df2 <- data.frame(y=bm,x=x)

gg2 <- ggplot(df,aes(x=x,y=y,color=id)) +
  geom_line() + theme_bw() + theme(legend.position = "none") +
  geom_point(data=df2,aes(x=x,y=y),inherit.aes = FALSE,size=2)
  
gg2

ggsave(ggarrange(gg,gg2,nrow=1),filename = "~/Desktop/postPreds.pdf",width =8,height=3)


  