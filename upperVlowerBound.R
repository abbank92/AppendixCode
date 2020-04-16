# upperVlowerBound.R
# (c) Alex Bank 2020

genPlot2 <- function(A,eps,ntimes,ax=TRUE,toplab=''){
  N <- dim(A)[1]
  eigList <- lapply(rep(0,ntimes),
                    function(x) 
                    {return(genPseudo(A,genE(eps,N)))})
  origEigs <- eigen(A)[[1]]
  lims <- (eps+eps^((N+1)/N))^(1/N)+0.1
  lowerB <- (eps+eps^((N+2)/N))^(1/N)
  if (ax) {
    plot(1,type="n", xlab="", ylab="",
         xlim=c(-lims, lims), ylim=c(-lims, lims),main = toplab)}
  else {
    plot(1,type="n", xlab="", ylab="",
         xlim=c(-lims, lims), ylim=c(-lims, lims),
         xaxt='n',yaxt='n',main=toplab)
  }
  for (e in eigList) {
    for(indEig in e) {
      if (abs(indEig) > lowerB){
        points(Re(indEig),Im(indEig),pch=20,cex=.4,col='red')
      }
      else {points(Re(indEig),Im(indEig),pch=20,cex=.25)}  
    }
  }
  points(Re(origEigs),Im(origEigs),col = 'gold',pch=20)
  draw.circle(0,0,(eps+eps^((N+2)/N))^(1/N),border='forestgreen',lty=2)
  draw.circle(0,0,(eps+eps^((N+1)/N))^(1/N),
              border='blue',lty=2)
  if (ax) {
    legend("topleft",legend = c("Leading two terms of Puiseux","Upper Bound"),
           col=c('forestgreen','blue'),lty=2,cex=.6)}
  else {
    legend("topleft", legend = paste0("eps=",round(eps,3)),cex = .75,bty='n')
  }
}

genPlot2(J3,.2,20000)
