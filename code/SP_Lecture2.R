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

