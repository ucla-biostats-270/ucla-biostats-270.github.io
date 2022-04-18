#
####
####### A sample Brownian path
####
#
library(ggplot2)

kicks <- rnorm(n=999999)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df <- data.frame(Displacement=bm,Time = (0:999999)/1000000)

gg <- ggplot(data=df,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() + ggtitle("A sample Brownian path")

ggsave("brownianPath.png",height=3,width = 5)

#
####
####### Piecewise linear functions converge to BM
####
#
n = 2
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg


n=5
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df1 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df1,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg

n=10
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df2 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df2,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg

n=20
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df3 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df3,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg

n=100
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df4 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df4,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg

n=1000
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df5 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

gg <- ggplot(data=df5,aes(x=Time,y=Displacement)) +
  geom_line() + theme_bw() 
gg




#
####
####### reflection principle
####
#
set.seed(2)
n=10000
kicks <- rnorm(n=n)
kicks <- c(0,kicks)
bm    <- cumsum(kicks)
df5 <- data.frame(Displacement=bm,Time = (0:n)/(n+1))

# reflected brownian motion
bmr <- c(bm[1:2500] , 2*bm[2501]-bm[2501:10001])
df6 <- data.frame(Displacement=bmr,Time = (0:n)/(n+1))

df <- rbind(df5,df6)
df$Reflected <- c(rep(c("original","reflection"),each=10001))
df$Reflected <- factor(df$Reflected)
df$` ` <- relevel(df$Reflected,ref = "reflection")

gg <- ggplot(data=df,aes(x=Time,y=Displacement,color=` `))  +
  geom_hline(yintercept = bm[2501]) +
  geom_line() + theme_bw() 
gg


ggsave("reflectionPrinciple.png",height=3,width = 5)


#
####
####### Prob max > 10 at time t
####
#
probMax <- vector()
for(i in 1:10000) {
  probMax[i] <- 2-2*pnorm(10/sqrt(i))
}
probMax <- c(0,probMax)

gg <- ggplot(data.frame(Probability=probMax,Time=0:10000),
             aes(x=Time,y=Probability)) +
  geom_line() + theme_bw() + ggtitle("Probability max > 10")
gg
ggsave("~/Desktop/probabilityGreater.png",height=3,width = 5)

#
####
####### Brownian motion with drift
####
#

#
mu <- 0.1
epsilon <- 0.01

n <- 1000
bm <- c(0,cumsum(rnorm(n))) + mu*0:n
df <- data.frame(Displacement=bm,Time=0:n,n=n)
gg <- ggplot(df,aes(x=Time,y=Displacement)) +
  geom_line() + geom_abline(slope = mu+epsilon,color="red") +
  geom_abline(slope = mu-epsilon,color="red") +
  theme_bw() 
gg

n <- 10000
bm <- c(0,cumsum(rnorm(n))) + mu*0:n
df2 <- data.frame(Displacement=bm,Time=0:n,n=n)
gg2 <- ggplot(df2,aes(x=Time,y=Displacement)) +
  geom_line() + geom_abline(slope = mu+epsilon,color="red") +
  geom_abline(slope = mu-epsilon,color="red") +
  theme_bw() 
gg2

n <- 100000
bm <- c(0,cumsum(rnorm(n))) + mu*0:n
df3 <- data.frame(Displacement=bm,Time=0:n,n=n)
gg3 <- ggplot(df3,aes(x=Time,y=Displacement)) +
  geom_line() + geom_abline(slope = mu+epsilon,color="red") +
  geom_abline(slope = mu-epsilon,color="red") +
  theme_bw() 
gg3

library(gridExtra)
ggsave("~/Desktop/bmWithDrift.png",grid.arrange(gg,gg2,gg3,ncol=3),width=12,height=3)

n <- 1000000
bm <- c(0,cumsum(rnorm(n))) + mu*0:n
plot(bm)
lines((mu+epsilon)*0:n,col="red")
lines((mu-epsilon)*0:n,col="red")


#
###
####### Finite differences for ODE u'(x)=5u(x)+2, u(0)=0
###
#
h <- 0.5
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg

h <- 0.1
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg2 <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg2

h <- 0.01
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg3 <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg3

h <- 0.001
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg4 <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg4

h <- 0.0001
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg5 <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg5

h <- 0.00001
u <- rep(0,1/h+1)
for(i in 1:(1/h)) {
  u[i+1] <- u[i] + h*(5*u[i]+2)
}
df <- data.frame(Values=u,Time=seq(from=0,to=1,by=h))
gg6 <- ggplot(df,aes(x=Time,y=Values)) +
  geom_line() + theme_bw() + ggtitle(paste0("h=",h)) +
  geom_hline(yintercept = 58.96519,col="red")
gg6

library(gridExtra)
ggsave("~/Desktop/finiteDifferences.png",grid.arrange(gg,gg2,gg3,gg4,gg5,gg6,ncol=3),width=12,height=6)

#
###
####### Euler-Maruyama for OU process
###
#
sigma <- 1
alpha <- 0
U <- vector()
h <- 1/100000
for(i in 1:20) {
  u <- rep(0,1/h+1)
  for(i in 1:(1/h)) {
    u[i+1] <- (1 - h*alpha)*u[i] + sigma*sqrt(h)*rnorm(1)
  }
  U <- c(U,u)
}
Time <- rep(seq(from=0,to=1,by=h),20)
Values <- U
df <- data.frame(Values,Time,Sample=rep(1:20,each=length(Time)/20))
df$Sample <- factor(df$Sample)
gg <- ggplot(df,aes(x=Time,y=Values,col=Sample)) +
  geom_line()  + ylim(-3,3) +
  theme_bw() + theme(legend.position = "none") + ggtitle(paste0("alpha=",alpha))
gg

sigma <- 1
alpha <- 1
U <- vector()
h <- 1/100000
for(i in 1:20) {
  u <- rep(0,1/h+1)
  for(i in 1:(1/h)) {
    u[i+1] <- (1 - h*alpha)*u[i] + sigma*sqrt(h)*rnorm(1)
  }
  U <- c(U,u)
}
Time <- rep(seq(from=0,to=1,by=h),20)
Values <- U
df <- data.frame(Values,Time,Sample=rep(1:20,each=length(Time)/20))
df$Sample <- factor(df$Sample)
gg2 <- ggplot(df,aes(x=Time,y=Values,col=Sample)) +
  geom_line()  + ylim(-3,3) +
  theme_bw() + theme(legend.position = "none") + ggtitle(paste0("alpha=",alpha))
gg2

sigma <- 1
alpha <- 10
U <- vector()
h <- 1/100000
for(i in 1:20) {
  u <- rep(0,1/h+1)
  for(i in 1:(1/h)) {
    u[i+1] <- (1 - h*alpha)*u[i] + sigma*sqrt(h)*rnorm(1)
  }
  U <- c(U,u)
}
Time <- rep(seq(from=0,to=1,by=h),20)
Values <- U
df <- data.frame(Values,Time,Sample=rep(1:20,each=length(Time)/20))
df$Sample <- factor(df$Sample)
gg3 <- ggplot(df,aes(x=Time,y=Values,col=Sample)) +
  geom_line()  + ylim(-3,3) +
  theme_bw() + theme(legend.position = "none") + ggtitle(paste0("alpha=",alpha))
gg3

ggsave("~/Desktop/ouEulerMaruyama.png",grid.arrange(gg,gg2,gg3,ncol=3),width=12,height=3)

#
###
####### Langevin Monte Carlo
###
#
library(ggplot2)

# gaussian target
gradLogPi <- function(theta) {
  return(-theta)
}

h <- 1/10
length <- 1000
vals <- rep(0,(length/h)+1)
for(i in 1:(length/h)) {
  vals[i+1] <- vals[i] + h*gradLogPi(vals[i]) + sqrt(2*h)*rnorm(1)
}
hist(vals)
qqnorm(vals)
qqline(vals,col="red")


# Langevin Monte Carlo
LMC <- function(h=1/10,timeLength=1000,D=1,initialPos=0) {
  if(D==1) {
    vals <- rep(0,(timeLength/h)+1)
    vals[1] <- initialPos
    for(i in 1:(timeLength/h)) {
      vals[i+1] <- vals[i] + h*gradLogPi(vals[i]) + sqrt(2*h)*rnorm(1)
    }
  } else {
    vals <- matrix(0,(timeLength/h)+1,2)
    vals[1,] <- initialPos
    for(i in 1:(timeLength/h)) {
      vals[i+1,] <- vals[i,] + h*gradLogPi(vals[i,]) + sqrt(2*h)*rnorm(D)
    }
  }
  return(vals)
}

# spherical 2D target
vals <- LMC(h=1/10,timeLength=1000,D=2,initialPos = 100)
colnames(vals) <- c("X","Y")
vals <- as.data.frame(vals)
gg <- ggplot(vals[1:100,],aes(x=X,y=Y)) +
  geom_path() +theme_bw()
gg
gg2 <- ggplot(vals,aes(x=X,y=Y)) +
  geom_path() + geom_density_2d() + theme_bw()
gg2

# spherical 2D target
vals <- LMC(h=1/10,timeLength=10000,D=2,initialPos = 1)
colnames(vals) <- c("X","Y")
vals <- as.data.frame(vals)
gg <- ggplot(vals[1:100,],aes(x=X,y=Y)) +
  geom_path() +theme_bw()
gg
gg2 <- ggplot(vals,aes(x=X,y=Y)) +
  geom_path() + geom_density_2d() +theme_bw()
gg2

# correlated 2D target
gradLogPi <- function(theta) {
  Rinv <- matrix(c(2.78,-2.22,-2.22,2.78),2,2)
  return(-Rinv%*%theta)
}
vals <- LMC(h=1/10,timeLength=10000,D=2,initialPos = 1)
colnames(vals) <- c("X","Y")
vals <- as.data.frame(vals)
gg <- ggplot(vals[1:100,],aes(x=X,y=Y)) +
  geom_path() +theme_bw() + ggtitle("100 steps: h=0.1, time=10")
gg
gg2 <- ggplot(vals,aes(x=X,y=Y)) +
  geom_path() + geom_density_2d() +theme_bw() + ggtitle("100k steps: h=0.1, time=10k")
gg2

ggsave("~/Desktop/LMC.png",gridExtra::grid.arrange(gg,gg2,ncol=2),width=12,height=4)

