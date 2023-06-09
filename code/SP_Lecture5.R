#
#####
######### sampling L-ensemble
#####
#

n <- 100
L <- matrix(rnorm(n^2),n,n)
L <- L %*% t(L) / n

eig.obj <- eigen(L)
vecs    <- eig.obj$vectors # columns are eigenvectors
vals    <- eig.obj$values

states <- vector()
for(i in 1:n) {
  u <- rbinom(n=1,size=1,prob=vals[i]/(1+vals[i]))
  if (u==1) {
    states <- c(states,i)
  }
}
  V <- vecs[,states]
  ncolV <- length(states)
  Y <- vector()
  
  while (ncolV>0) {
    probs <- rep(0,n)
    for(i in 1:n) {
      if(is.matrix(V)) {
        probs[i] <- sum(V[i,]^2)
      } else {
        probs[i] <- sum(V[i]^2)
      }
    }
    item <- sample(x=1:n,size=1,prob=probs)

    Y <- c(Y,item)

    if(is.matrix(V)) {
    Vproj <- V
    Vproj[item,] <- 0
    V <- pracma::gramSchmidt(Vproj)$Q[,1:(ncolV-1)]
    } else {
      Vproj <- V
      Vproj[item] <- 0
      Vproj <- Vproj / sqrt(sum(Vproj^2))
    }
    ncolV <- ncolV - 1
  }
states
Y

#
######
############ 1D grid with kernel function
######
#
n <- 100
x <- 1:n
X <- matrix(x,n,n)
L <- 5*exp(-(X-t(X))^2/(1^2*2))

eig.obj <- eigen(L)
vecs    <- eig.obj$vectors # columns are eigenvectors
vals    <- eig.obj$values

states <- vector()
for(i in 1:n) {
  u <- rbinom(n=1,size=1,prob=vals[i]/(1+vals[i]))
  if (u==1) {
    states <- c(states,i)
  }
}
V <- vecs[,states]
ncolV <- length(states)
Y <- vector()

while (ncolV>0) {
  probs <- rep(0,n)
  for(i in 1:n) {
    if(is.matrix(V)) {
      probs[i] <- sum(V[i,]^2)
    } else {
      probs[i] <- sum(V[i]^2)
    }
  }
  item <- sample(x=1:n,size=1,prob=probs)
  
  Y <- c(Y,item)
  
  if(is.matrix(V)) {
    Vproj <- V
    Vproj[item,] <- 0
    V <- pracma::gramSchmidt(Vproj)$Q[,1:(ncolV-1)]
  } else {
    Vproj <- V
    Vproj[item] <- 0
    Vproj <- Vproj / sqrt(sum(Vproj^2))
  }
  ncolV <- ncolV - 1
}
states
#Y <- sort(Y)
Y
plot(Y,ylim = c(0,101))


#
######
############ 2D grid with kernel function
######
#
N <- 15
x <- expand.grid(1:N,1:N)
#plot(x)
n <- N^2
L <- matrix(0,n,n)
for (i in 1:(n)) {
  for(j in i:(n)) {
    L[i,j] <- 1*exp(-sum((x[i,]-x[j,])^2)/(20^2*2))
  }
}
L <- 0.5*(L + t(L))

eig.obj <- eigen(L)
vecs    <- eig.obj$vectors # columns are eigenvectors
vals    <- eig.obj$values

states <- vector()
for(i in 1:n) {
  u <- rbinom(n=1,size=1,prob=vals[i]/(1+vals[i]))
  if (u==1) {
    states <- c(states,i)
  }
}
V <- vecs[,states]
ncolV <- length(states)
Y <- vector()

while (ncolV>0) {
  probs <- rep(0,n)
  for(i in 1:n) {
    if(is.matrix(V)) {
      probs[i] <- sum(V[i,]^2)
    } else {
      probs[i] <- sum(V[i]^2)
    }
  }
  item <- sample(x=1:n,size=1,prob=probs)
  
  Y <- c(Y,item)
  
  if(is.matrix(V)) {
    Vproj <- V
    Vproj[item,] <- 0
    V <- pracma::gramSchmidt(Vproj)$Q[,1:(ncolV-1)]
  } else {
    Vproj <- V
    Vproj[item] <- 0
    Vproj <- Vproj / sqrt(sum(Vproj^2))
  }
  ncolV <- ncolV - 1
}
states
#Y <- sort(Y)
Y

# low fi movie
for(i in 1:length(Y)) {
  Sys.sleep(1)
  plot(x[Y[1:i],],xlim=c(1,N),ylim=c(1,N))
}
