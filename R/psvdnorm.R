###Description
##
###Computes the probability than an n-dimensional (n>1) multivariate
###normal distribution lies inside a rectangles of R^n. Suited for large
###probabilities. In particular, it can be used to compute p-values of
###tests.
##
###Usage
##
###msvdnorm(lower=Null,upper=Null,Sigma=Null,Nsample=10000,Ncov=1)
##
###Arguments
##
###lower    the vector of lower limits of length n.
###upper    the vector of upper limits of length n.
###Sigma    the covariance matrix of dimension n.
###Nsample  the sample size.
###Ncov     the number of control variates to be used (<=n).
##
##
###Details
##
###This program involves the computation of multivariate normal
###probabilities with arbitrary covariance matrices. The methodology is
###described in Phinikettos and Gandy (2010), Computational Statistics &
###Data Analysis.
##
###Value
##
###The evaluated probability is returned.
##
###Examples
##  
## # source("psvdnorm.R")
##  #lower <- rep(-5, 10)
##  #upper <- rep(5, 10)
##  #a  <-  matrix(rnorm(100),10,10)
##  #Sigma <-  a %*% t(a)
##  #psvdnorm(lower,upper,Sigma,Nsample=10000,Ncov=5)
##
###Disclaimer
##
###This is experimental code - use at your own risk. There is no
###guarantee that it is stable for all possible inputs. The code is
###designed to work for high-dimensional cases and for cases with 
###high results (close to 1).
##
##
##
###THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
###"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
###LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
###A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
###OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
###SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
###LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
###DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
###THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
###(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
###OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE





########################################################################################

########################################################################################
########### Importance Sampling Section 2.2 ############################################
########################################################################################

SVD_cond_IS_fast <- function(initialsteps,maxsteps,var,U,d,Lower,Upper) {
  
  lower <- Lower + (Upper-Lower)* (U[,1]<0)
  upper <- Upper + (Lower-Upper)* (U[,1]<0)
  nIS=length(d)-1
  Ucond <- U[,-1,drop=FALSE]
  dcond <- d[-1]
  p <- vector("numeric")
  S <- vector("numeric")
  slong <- vector("numeric")
  sumIS <- vector("numeric")
  S[1] <- var
  while(length(slong)<maxsteps){
    i <- length(S)
    n <- max(initialsteps,length(slong))
    nsamples <- min(n,maxsteps-length(slong))
    if (nsamples<0) break
    Zcond <- matrix(rnorm(nsamples*nIS,0,sqrt(S[i])),nIS,nsamples)
    h <- Ucond %*% (dcond*Zcond)
    a <-  (lower - h)/(d[1]*U[,1])
    b <-  (upper - h)/(d[1]*U[,1])
    sumISS <- colSums(Zcond^2)
    w <-exp(-1/2*sumISS*(1-1/S[i]) ) * (sqrt(S[i]))^nIS
    maxa <- apply(a,2,max)
    minb <- apply(b,2,min)
    phi <- ((pnorm(minb)-pnorm(maxa))*(maxa < minb)) 
    P <- phi*w
    slong <- c(slong,rep(1,nsamples))
    p <- c(p,P)
    sumIS <- c(sumIS,sumISS)
    if (initialsteps < maxsteps) {
      z <- function(s0){
        s <- s0[1]
        if (s<=0)  return(Inf);
        W <- exp( -1/2*sumIS*(1-1/s)  + nIS/2*log(s) )
        res <- mean(W*p)
        res
      }
      optimisation <- optimize(z,c(0,5),tol=0.00001)
      S[i+1] <- optimisation$min  
      if (   all( abs(p-p[1])<1e-15 )  ) S[i+1]=S[i]
    }   
  }
  list(p=mean(p),sd=sd(p)/sqrt(length(p)),S=S[length(S)])
}
###########################################################################
SVD_cond_IS_fast2 <- function(initialsteps,maxsteps,var,U,d,Lower,Upper) {
  
  lower <- Lower + (Upper-Lower)* (U[,1]<0)
  upper <- Upper + (Lower-Upper)* (U[,1]<0)
  nIS=length(d)-1
  Ucond <- U[,-1,drop=FALSE]
  dcond <- d[-1]
  p <- vector("numeric")
  S <- vector("numeric")
  slong <- vector("numeric")
  sumIS <- vector("numeric")
  S[1] <- var   
  while(length(slong)<maxsteps){
    i <- length(S)
    n <- max(initialsteps,length(slong))
    nsamples <- min(n,maxsteps-length(slong))
    if (nsamples<0) break
    Zcond <- matrix(rnorm(nsamples*nIS,0,sqrt(S[i])),nIS,nsamples)
    h <- Ucond %*% (dcond*Zcond)
    a <-  (lower - h)/(d[1]*U[,1])
    b <-  (upper - h)/(d[1]*U[,1])
    sumISS <- colSums(Zcond^2)
    w <-exp(-1/2*sumISS*(1-1/S[i]) ) * (sqrt(S[i]))^nIS
    maxa <- apply(a,2,max)
    minb <- apply(b,2,min)
    phi <- (1-(pnorm(minb)-pnorm(maxa))*(maxa < minb)) 
    P <- phi*w
    slong <- c(slong,rep(1,nsamples))
    p <- c(p,P)
    sumIS <- c(sumIS,sumISS)
    if (initialsteps < maxsteps) {
      z <- function(s0){
        s <- s0[1]
        if (s<=0)  return(Inf);
        W <- exp( -1/2*sumIS*(1-1/s)  + nIS/2*log(s) )
        res <- mean(W*p)
        res
      }
      optimisation <- optimize(z,c(0,5),tol=0.00001)
      S[i+1] <- optimisation$min  
      if (   all( abs(p-p[1])<1e-15 )  ) S[i+1]=S[i]
    }   
  }
  list(p=mean(p),sd=sd(p)/sqrt(length(p)),S=S[length(S)])
}


################################################################################
########### Optimal splitting parameters Section 2.3 ###########################
################################################################################

SVD_cond_split_optimal_few <- function(R,U,d,var,Lower,Upper){
  if (length(Lower) > 2 ){
    lower <- Lower + (Upper-Lower)* (U[,1]<0)
    upper <- Upper + (Lower-Upper)* (U[,1]<0)
    Ucond <- U[,-1,drop=FALSE]
    dcond <- d[-1]
    n <- length(d)-1
    cumk <- cumsum(dcond^2)
    sumk <- sum(dcond^2)
    nSplit <- which(cumk>0.85*sumk)[1] 
    if (nSplit > floor(n/2) & n >1 ) nSplit = floor(n/2)
    Zcond <- matrix(rnorm(R*n,0,sqrt(var)),n,R)
    Zcond <- matrix(Zcond,n,R*2)  
    hunsplit <- Ucond[,(nSplit+1):n,drop=FALSE] %*% (dcond[(nSplit+1):n]*Zcond[(nSplit+1):n,]) 
    Zcond[1:nSplit,] <- rnorm(2*R*nSplit,0,sqrt(var))
    hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zcond[1:nSplit,,drop=FALSE])
    h <- hsplit + hunsplit
    a <-  (lower -h)/(d[1]*U[,1])
    b <-  (upper -h)/(d[1]*U[,1])
    sumIS <- colSums(Zcond^2)
    w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
    maxa <- apply(a,2,max)
    minb <- apply(b,2,min)
    P <- ((pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
    A <- P[1:R]
    B <- P[(R+1):(2*R)]
    if (abs(sum(A-B))>0) rho <- abs(cor(A,B)) 
    else rho <- 0.01
    a1= n-nSplit
    b1= nSplit
    timeSplit = (a1)*(  sqrt(b1/(a1))  + sqrt(rho) )^2
    optS=sqrt(a1/(rho*b1))    
    list(timeSplit=nSplit,optimalS=optS)
  }
  else list(timeSplit=1,optimalS=1)
}

#######################################################################

SVD_cond_split_optimal_few2 <- function(R,U,d,var,Lower,Upper){
  if (length(Lower) > 2 ){
    lower <- Lower + (Upper-Lower)* (U[,1]<0)
    upper <- Upper + (Lower-Upper)* (U[,1]<0)
    Ucond <- U[,-1,drop=FALSE]
    dcond <- d[-1]
    n <- length(d)-1
    cumk <- cumsum(dcond^2)
    sumk <- sum(dcond^2)
    nSplit <- which(cumk>0.85*sumk)[1] 
    if (nSplit > floor(n/2)& n >1) nSplit = floor(n/2)
    Zcond <- matrix(rnorm(R*n,0,sqrt(var)),n,R)
    Zcond <- matrix(Zcond,n,R*2)
    hunsplit <- Ucond[,(nSplit+1):n,drop=FALSE] %*% (dcond[(nSplit+1):n]*Zcond[(nSplit+1):n,,drop=FALSE]) 
    Zcond[1:nSplit,] <- rnorm(2*R*nSplit,0,sqrt(var))
    hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zcond[1:nSplit,])
    h <- hsplit + hunsplit
    a <-  (lower -h)/(d[1]*U[,1])
    b <-  (upper -h)/(d[1]*U[,1])
    sumIS <- colSums(Zcond^2)
    w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
    maxa <- apply(a,2,max)
    minb <- apply(b,2,min)
    P <- (1-(pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
    A <- P[1:R]
    B <- P[(R+1):(2*R)]
    if (abs(sum(A-B))>0) rho <- abs(cor(A,B)) 
    else rho <- 0.01
    a1= n-nSplit
    b1= nSplit
    timeSplit = (a1)*(  sqrt(b1/(a1))  + sqrt(rho) )^2
    optS=sqrt(a1/(rho*b1))    
    list(timeSplit=nSplit,optimalS=optS)
  } 
  else list(timeSplit=1,optimalS=1)  
}

########################################################################
############ Splitting + control variates Sections 2.3-2.4 #############
########################################################################

SVD_cond_split_control_row <- function(ncov,nSplit,R,S,U,d,var,Lower,Upper){
  n <- length(d)-1 
  lower <- Lower + (Upper-Lower)* (U[,1]<0)
  upper <- Upper + (Lower-Upper)* (U[,1]<0)   
  sdcontrol <- vector("numeric")
  pcontrol <- matrix(0,R*S,ncov)
  pmeancontrol <- vector("numeric")   
  p <- vector("numeric")   
  Zsplit <- matrix(rnorm(R*nSplit,0,sqrt(var)),nSplit,R)
  if (n>1) Zunsplit <- matrix(rnorm(R*(n-nSplit),0,sqrt(var)),n-nSplit,R)
  Ucond <- U[,-1,drop=FALSE]
  dcond <- d[-1]  
  hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zsplit)
  hunsplit = 0   
  if (n>1) hunsplit <- Ucond[,(nSplit+1):n,drop=FALSE] %*% (dcond[(nSplit+1):n]*Zunsplit)
  h <- hsplit + hunsplit   
  a <-  (lower -h)/(d[1]*U[,1])
  b <-  (upper -h)/(d[1]*U[,1])
  sumISsplit <- colSums(Zsplit^2)
  sumISunsplit = 0
  if (n>1) sumISunsplit <- colSums(Zunsplit^2)
  sumIS <- sumISunsplit + sumISsplit  
  w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
  maxa <- apply(a,2,max)
  minb <- apply(b,2,min)   
  argbs <- apply(b,2,order)[1,]
  freq <- sapply(1:length(d),function(i) sum(argbs==i))
  argb <- order(freq,decreasing=T)[1:ncov]
  if (ncov>0) {
      for (i in 1:ncov) sdcontrol[i] <- sqrt(sum( (U[argb[i],] * d)^2 ))
      for (i in 1:ncov) pmeancontrol[i] <- (pnorm(Upper[argb[i]],0, sdcontrol[i])- pnorm(Lower[argb[i]],0, sdcontrol[i]))
  }
  P <- ((pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
  p <- c(p,P)   
  if (ncov>0) {
     for (i in 1:ncov) {
        maxca <- (lower[argb[i]] - h[argb[i],])/(d[1]*U[argb[i],1])
        mincb <- (upper[argb[i]] - h[argb[i],])/(d[1]*U[argb[i],1]) 
        pcontrol[1:R,i] <- ((pnorm(mincb)-pnorm(maxca))*(maxca <= mincb))*w -pmeancontrol[i]
     }  
  }
  if (S >1){
    for (i in 2:S){
      Zsplit[1:nSplit,] <- rnorm(nSplit*R,0,sqrt(var))
      hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zsplit)
      h <- hsplit + hunsplit    
      a <-  (lower -h)/(d[1]*U[,1])
      b <-  (upper -h)/(d[1]*U[,1])
      sumISsplit <- colSums(Zsplit^2)
      sumIS <- sumISunsplit + sumISsplit
      w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
      maxa <- apply(a,2,max)
      minb <- apply(b,2,min)
      P <- ((pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
      p <- c(p,P)
      if (ncov>0) {
         for (j in 1:ncov) {
           maxca <- (lower[argb[j]] - h[argb[j],])/(d[1]*U[argb[j],1])
           mincb <- (upper[argb[j]] - h[argb[j],])/(d[1]*U[argb[j],1]) 
          pcontrol[((i-1)*R+1):(i*R),j] <- ((pnorm(mincb)-pnorm(maxca))*(maxca <= mincb))*w -pmeancontrol[j]
         }      
     }  
     }
  } 
  pc = p
  if (ncov>0) {
    m <- lm(p~pcontrol)
    beta <- as.vector(m$coefficients)    
    pc <- beta[1]
  }
  list(p=p,pc=pc)
}
#################################################################

SVD_cond_split_control_row2 <- function(ncov,nSplit,R,S,U,d,var,Lower,Upper){
  n <- length(d)-1
  lower <- Lower + (Upper-Lower)* (U[,1]<0)
  upper <- Upper + (Lower-Upper)* (U[,1]<0)
  sdcontrol <- vector("numeric")
  pcontrol <- matrix(0,R*S,ncov)
  pmeancontrol <- vector("numeric")
  p <- vector("numeric")   
  Zsplit <- matrix(rnorm(R*nSplit,0,sqrt(var)),nSplit,R)
  if (n>1)   Zunsplit <- matrix(rnorm(R*(n-nSplit),0,sqrt(var)),n-nSplit,R)
  Ucond <- U[,-1,drop=FALSE]
  dcond <- d[-1]  
  hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zsplit)
  hunsplit = 0 
  if (n>1)   hunsplit <- Ucond[,(nSplit+1):n,drop=FALSE] %*% (dcond[(nSplit+1):n]*Zunsplit)
  h <- hsplit + hunsplit   
  a <-  (lower -h)/(d[1]*U[,1])
  b <-  (upper -h)/(d[1]*U[,1])
  sumISsplit <- colSums(Zsplit^2)
  sumISunsplit = 0 
  if (n>1)   sumISunsplit <- colSums(Zunsplit^2)
  sumIS <- sumISunsplit + sumISsplit  
  w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
  maxa <- apply(a,2,max)
  minb <- apply(b,2,min)   
  argbs <- apply(b,2,order)[1,]
  freq <- sapply(1:length(d),function(i) sum(argbs==i))
  argb <- order(freq,decreasing=T)[1:ncov]
  if (ncov>0) {
     for (i in 1:ncov) sdcontrol[i] <- sqrt(sum( (U[argb[i],] * d)^2 ))
     for (i in 1:ncov) pmeancontrol[i] <- 1-(pnorm(Upper[argb[i]],0, sdcontrol[i])- pnorm(Lower[argb[i]],0, sdcontrol[i]))
  }
  P <- (1-(pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
  p <- c(p,P)
  if (ncov>0) {
    for (i in 1:ncov) {
      maxca <- (lower[argb[i]] - h[argb[i],])/(d[1]*U[argb[i],1])
      mincb <- (upper[argb[i]] - h[argb[i],])/(d[1]*U[argb[i],1]) 
      pcontrol[1:R,i] <- (1-(pnorm(mincb)-pnorm(maxca))*(maxca <= mincb))*w -pmeancontrol[i]
    }  
  }
  if (S >1){
    for (i in 2:S){
      Zsplit[1:nSplit,] <- rnorm(nSplit*R,0,sqrt(var))
      hsplit <- Ucond[,1:nSplit,drop=FALSE] %*% (dcond[1:nSplit]*Zsplit)
      h <- hsplit + hunsplit    
      a <-  (lower -h)/(d[1]*U[,1])
      b <-  (upper -h)/(d[1]*U[,1])
      sumISsplit <- colSums(Zsplit^2)
      sumIS <- sumISunsplit + sumISsplit
      w <-exp(-1/2*sumIS*(1-1/var) ) * sqrt(var)^n
      maxa <- apply(a,2,max)
      minb <- apply(b,2,min)
      P <- (1-(pnorm(minb)-pnorm(maxa))*(maxa <= minb))*w
      p <- c(p,P)
      if (ncov>0) {
        for (j in 1:ncov) {
          maxca <- (lower[argb[j]] - h[argb[j],])/(d[1]*U[argb[j],1])
          mincb <- (upper[argb[j]] - h[argb[j],])/(d[1]*U[argb[j],1]) 
          pcontrol[((i-1)*R+1):(i*R),j] <- (1-(pnorm(mincb)-pnorm(maxca))*(maxca <= mincb))*w -pmeancontrol[j]
         }      
      }  
    }
  } 
  pc = p
  if (ncov>0) {
     m <- lm(p~pcontrol)
     beta <- as.vector(m$coefficients)    
     pc <- beta[1]
  }
  list(p=p,pc=pc)
}
########################################################################
######### coordinating the four variance reduction techniques ########## 
########################################################################

SVD_split_cond_IS_row2 <- function(ncov,sampleIS,R,U,d,var,lower,upper){
  n=length(d)-1
  checkp <- SVD_cond_IS_fast(50,50,1,U,d,lower,upper)
  checkp <- mean(checkp$p)
  if (checkp > 0.5){    
    timeIS <- system.time(optimalIS <- SVD_cond_IS_fast2(floor(sampleIS/(3*10)),floor(sampleIS/3),var,U,d,lower,upper))[1]
    s <- optimalIS$S
    optimalSplit <- SVD_cond_split_optimal_few2(R,U,d,s,lower,upper)
    nSplit= optimalSplit$timeSplit
    a=n-nSplit
    b=nSplit
    optS <- max(1,floor(optimalSplit$optimalS))
    ratio <- (optS * nSplit +(n-nSplit))/n
    optR <- floor(sampleIS/ratio)
    timeSplit <- system.time(Split <- SVD_cond_split_control_row2(ncov,nSplit,optR,optS,U,d,s,lower,upper))[1]
    psplit <- Split$p
    pcov   <- Split$pc
    pcov <- 1-pcov
  }
  if (checkp <= 0.5){    
    timeIS <- system.time(optimalIS <- SVD_cond_IS_fast(floor(sampleIS/(3*10)),floor(sampleIS/3),var,U,d,lower,upper))[1]
    s  <- optimalIS$S
    optimalSplit <- SVD_cond_split_optimal_few(R,U,d,s,lower,upper)
    nSplit= optimalSplit$timeSplit
    a=n-nSplit
    b=nSplit
    optS <- max(1,floor(optimalSplit$optimalS))
    ratio <- (optS * nSplit +(n-nSplit))/n
    optR <- floor(sampleIS/ratio)
    timeSplit <- system.time(Split <- SVD_cond_split_control_row(ncov,nSplit,optR,optS,U,d,s,lower,upper))[1]
    psplit <- Split$p
    pcov   <- Split$pc
    }
  mean(pcov)
}

####################################################################
############# main function, computing the svd #####################
####################################################################

psvdnorm <- function(lower,upper,Sigma=diag(rep(1,length(lower))),Nsample=10000,Ncov=1){
  svdS <- svd(Sigma)
  U <- svdS$u
  d <- sqrt(svdS$d) 
  SVD_split_cond_IS_row2(Ncov,Nsample,1000,U,d,1,lower,upper)
}




