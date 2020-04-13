# upperBound.R
# (c) Alex Bank 2020

# Load library for drawing circles
library(plotrix)

# Function to generate Jordan blocks
genJ <- function(N) {
  J <- matrix(rep(0,N*N),nrow = N,ncol = N)
  for (i in 1:(N-1)) J[i,i+1] <- 1
  return(J)
}

# Function for generating a complex number
rcomp <- function(x){
  return(complex(
    real = runif(1,-1,1),
    imaginary = runif(1,-1,1)))
}

# Function for generating an error matrix
# with norm <= eps
genE <- function(eps,N){
  errM <- matrix(
    sapply(rep(0,N*N),rcomp),
    nrow = N,ncol = N)
  m <- norm(errM,type='2')
  if (m>eps) {
    x <- eps/norm(errM,type = "2")}
  else {
    x <- 1}
  return(x*errM)
}

# Function for returning eigenvalues
# of A+E
genPseudo <- function(A,E){
  return(eigen(A+E)[[1]])
}

# Function for generating the plots
genPlot <- function(A,eps,ntimes,ax=TRUE,toplab=''){
  N <- dim(A)[1]
  eigList <- lapply(rep(0,ntimes),
                    function(x) 
                      {return(genPseudo(A,genE(eps,N)))})
  origEigs <- eigen(A)[[1]]
  lims <- 1+eps+0.1
  if (ax) {
    plot(1,type="n", xlab="", ylab="",
         xlim=c(-lims, lims), ylim=c(-lims, lims),main = toplab)}
  else {
    plot(1,type="n", xlab="", ylab="",
         xlim=c(-lims, lims), ylim=c(-lims, lims),
         xaxt='n',yaxt='n',main=toplab)
  }
  for (e in eigList) {
    points(Re(e),Im(e),pch=20,cex=.25)
  }
  points(Re(origEigs),Im(origEigs),col = 'gold',pch=20)
  draw.circle(0,0,1+eps,border='red',lty=2)
  draw.circle(0,0,(eps+eps^((N+1)/N))^(1/N),
              border='blue',lty=2)
  if (ax) {
  legend("topleft",legend = c("Old Upper Bound","New Upper Bound"),
         col=c('red','blue'),lty=2,cex=.7)}
  else {
    legend("topleft", legend = paste0("eps=",round(eps,3)),cex = .75,bty='n')
  }
}

# Figure 1
par(mai =rep(.5,4))
J3 <- genJ(3)
genPlot(J3,.1,2000)

# Figure 2
par(mfrow = c(4,4), mai = rep(.1,4))
genPlot(J3,.1,500,FALSE)
for (i in 1:15) genPlot(J3,runif(1,0,.5),500,FALSE)

# Figure 3
par(mfrow = c(2,2), mai = rep(.2,4))
J4 <- genJ(4); J5 <- genJ(5); J6 <- genJ(6); J7 <- genJ(7)
genPlot(J4,.1,1000,FALSE,"J4")
genPlot(J5,.1,1000,FALSE,"J5")
genPlot(J6,.1,1000,FALSE,"J6")
genPlot(J7,.1,1000,FALSE,"J7")

