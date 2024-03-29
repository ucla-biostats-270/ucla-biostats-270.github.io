#
###
####### exponential spacings for homogeneous Poisson process
###
#
lambda <- 1
maxTime <- 40
nPoints <- 2*lambda*maxTime
n.t <- cumsum(rexp(nPoints,rate=lambda))
n.t <- n.t[n.t<40]
x   <- 0:length(n.t)
cat(length(x)-1,"\n")

plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE)
points(n.t, x[-length(x)], pch = 1)


#
###
####### another method for sampling homogeneous Poisson process
###
#
lambda <- 1
maxTime <- 40
nPoints <- rpois(n=1,lambda = lambda*maxTime)
nPoints
n.t <- sort(runif(n=nPoints)*maxTime)
x   <- 0:length(n.t)
plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE)
points(n.t, x[-length(x)], pch = 1)

#
###
#######
###
#
lambda <- function(x) {
  y <- 1/x 
  return(y)
}

Lambda <- function(x) {
  y <- log(x)
  return(y)
}

Lambda_inv <- function(y) {
  x <- exp(y) 
  return(x)
}

maxTime <- 40
nPoints <- 100
y.t <- cumsum(rexp(nPoints,rate=1))
n.t <- Lambda_inv(y.t)
n.t <- n.t[n.t<40]
x   <- 0:length(n.t)
cat(length(x)-1,"\n")

plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE,
     xlim = c(0,40))
points(n.t, x[-length(x)], pch = 1)

#
###
####### thinning
###
#

lambda_fun <- function(x) {
  y <- sin(x) + 1 
  return(y)
}

# sample homogeneous Poisson process with lambda>lambda(x)
lambda <- 2
maxTime <- 40
nPoints <- rpois(n=1,lambda = lambda*maxTime)
nPoints
n.t <- sort(runif(n=nPoints)*maxTime)
x   <- 0:length(n.t)
plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE)
points(n.t, x[-length(x)], pch = 1)

# rejection portion
n.t2 <- vector()
for(i in 1:nPoints) {
  u <- runif(1)
  if(u<lambda_fun(n.t[i])/lambda) {
    n.t2 <- c(n.t2,n.t[i])
  }
}
n.t <- n.t2
x   <- 0:length(n.t)
plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE)
points(n.t, x[-length(x)], pch = 1)


# 
#####
########## hawkes processes ####################################################
#####
#

#
### thinning
#

rate <- function(t,mu,a,b,preceding=NULL,arePreceding=TRUE) {
  if(arePreceding) {
    out <- mu + a/sum((0.0000001+(t-preceding)))#a*sum(exp(-b*(t-preceding)))
  } else {
    out <- mu
  }
  return(out)
}

mu <- 0.1 # background rate
a  <- 1 # self-exciting rate coefficient
b  <- 2 # exponential lengthscale
epsilon <- 1E-6
endTime <- 100
t <- 0
events <- vector()
arePreceding <- FALSE

while (t<endTime) {
  upBound <- rate(t=t+epsilon, mu=mu, a=a, b=b,
                  preceding = events, arePreceding = arePreceding)
  t <- t + rexp(n=1,rate = upBound)
  u <- runif(1)
  ratio <- rate(t=t, mu=mu, a=a, b=b,
                preceding = events, arePreceding = arePreceding) / upBound
  if (u<ratio & t< 100) {
    events       <- c(events,t)
    arePreceding <- TRUE
  }
}

n.t <- events
x   <- 0:length(n.t)
plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE,
     xlim=c(0,100))
points(n.t, x[-length(x)], pch = 1)


#
### immigrant/children representation
#
mu <- 0.1 # background rate
a  <- 2 # self-exciting rate coefficient
b  <- 4 # exponential lengthscale
endTime <- 100
t <- 0

backgroundEvents <- vector()
while(t<endTime) {
  t <- t + rexp(n=1,rate=mu)
  if(t < 100)   backgroundEvents <- c(backgroundEvents,t)

}
backgroundEvents
events     <- backgroundEvents
currentGen <- backgroundEvents

while(length(currentGen)>0) {
  
  nextGen    <- vector()
  for(k in 1:length(currentGen)) {
    nChilds <- rpois(n=1,lambda=a/b)
    if (nChilds>0) {
      childs <- currentGen[k] + rexp(n=nChilds,rate=b)
      nextGen <- c(nextGen,childs)
    }
  }
  currentGen <- nextGen
  events <- c(events,currentGen)
  
}
events <- sort(events)
events <- events[events<100]

n.t <- events
x   <- 0:length(n.t)
plot(stepfun(n.t, x), xlab="t", ylab="N",
     do.points = TRUE,pch = 16,col.points = "blue",verticals = FALSE)
points(n.t, x[-length(x)], pch = 1)


#
####
######## spatiotemporal HP with cluster algorithm
####
#

# t less than 100
# space is 100 by 100 square
# background is homogeneous

mu <- 1 # background rate
a  <- 2.5 # self-exciting rate coefficient
b  <- 3 # exponential lengthscale
endTime <- 100
t <- 0

backgroundEvents <- matrix(0,nrow=0,ncol=3)
while(t<endTime) {
  t <- t + rexp(n=1,rate=mu)
  vnt <- c(t,100*runif(2))
  if(t < 100)   backgroundEvents <- rbind(backgroundEvents,vnt)
  
}
backgroundEvents
events     <- backgroundEvents
currentGen <- backgroundEvents

while(length(currentGen)>0) {
  
  nextGen    <- matrix(0,nrow=0,ncol=3)
  for(k in 1:length(currentGen[,1])) {
    nChilds <- rpois(n=1,lambda=a/b)
    if (nChilds>0) {
      childs <- currentGen[k,1] + rexp(n=nChilds,rate=b)
      vnt    <- cbind(childs,matrix(rnorm(2*length(childs))+currentGen[k,2:3],length(childs),2))
      nextGen <- rbind(nextGen,vnt)
    }
  }
  currentGen <- nextGen
  events <- rbind(events,currentGen)
  
}
events <- events[order(events[,1]),]
colnames(events) <- c("t","x","y")
rownames(events) <- NULL

plot(NULL,xlim=c(0,100),ylim=c(0,100))
currentTime <- 0
for(i in 1:(dim(events)[1])){
  Sys.sleep((events[i,1]-currentTime)/10)
  currentTime <- events[i,1]/10
  plot(events[1:i,2],events[1:i,3],xlim=c(0,100),ylim=c(0,100))
}
