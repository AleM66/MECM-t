###########################################################################################
# INFLU?NCIA LOCAL PARA O MEMC-t (Modelo com erro de medida censurado baseado na dist. t) #
###########################################################################################

 ## A. PERTURBA??O DE PONDERA??O DE CASOS ##

# 1ro rodar o ECM para usar as estimativas dos par?metros e as matrizes geradas no ECM

##############################################
### MATRIZ DELTA: 1ra derivada da fun??o Q ###
##############################################

deltaPondMEMC.t<-function(est){
  alpha <- est$alpha1
  beta <- est$beta1
  mux <- est$mux
  sigmax <- est$sigmax
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  uz <- est$tuzi
  uxz <- est$tuxzi
  uzz <- est$tuzzi
  n <- ncol(uz)
  delta <- matrix(0,k,n)
  u <- est$tui
  ux <- est$tuxi
  uxx <- est$tuxxi
  d_omega <- est$dd
  phi <- sqrt(d_omega)
  iomega22 <- solve(diag(d_omega[-1]))
  for (i in 1:n){
    uzi <- uz[((i-1)*p+1):(i*p),i]
    uz_s <- as.matrix(uzi[-1])
    uxzi <- uxz[i,((i-1)*p+1):(i*p)]
    uxz_s <- t(as.matrix(uxzi[-1]))
    uzzi <- uzz[((i-1)*p+1):(i*p),((i-1)*p+1):(i*p)]
    delta[1:r,i] <- iomega22%*%(uz_s-u[i]*alpha-ux[i]*beta)
    delta[(r+1):(2*r),i] <- iomega22%*%(t(uxz_s)-ux[i]*alpha-uxx[i]*beta)
    delta[2*r+1,i] <- (ux[i]-mux*u[i])/sigmax
    delta[2*r+2,i] <- -1/(2*sigmax)+(uxx[i]-2*mux*ux[i]+mux^2*u[i])/(2*sigmax^2)
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -1/phi[j]+(uzzi[j,j]-2*uxzi[j]+uxx[i])/phi[j]^3
      }
      else{
        delta[2*p+j,i] <- -1/phi[j]+(uzzi[j,j]-2*uzi[j]*alpha[j-1]-2*uxzi[j]*beta[j-1]+u[i]*alpha[j-1]^2+2*ux[i]*alpha[j-1]*beta[j-1]+uxx[i]*beta[j-1]^2)/phi[j]^3
      }
    }
  }
  return(delta)
}

#####################################################
### MATRIZ HESSIANA: Segunda derivada da fun??o Q ###
#####################################################

Q2MEMC.t<-function(est){
  
  alpha <- est$alpha1
  beta <- est$beta1
  mux <- est$mux
  sigmax <- est$sigmax
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  uz <- est$tuzi
  n <- ncol(uz)
  uxz <- est$tuxzi
  uzz <- est$tuzzi
  u <- est$tui
  ux <- est$tuxi
  uxx <- est$tuxxi
  d_omega <- est$dd
  phi <- sqrt(d_omega)
  iomega22 <- solve(diag(d_omega[-1]))
  
  QI<-matrix(0,k,k)
  
  QI[1:r,1:r]<- -sum(u)*iomega22
  QI[1:r,p:(2*r)]<- -sum(ux)*iomega22
  for (j in 1:r){
    QI[j,(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[j,(2*p+1+j)] <- QI[j,(2*p+1+j)]-2*(uz[((i-1)*p+1+j),i]-u[i]*alpha[j]-ux[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[p:(2*r),1:r]<- -sum(ux)*iomega22
  QI[p:(2*r),p:(2*r)]<- -sum(uxx)*iomega22
  for (j in 1:r){
    QI[(r+j),(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[(r+j),(2*p+1+j)] <- QI[(r+j),(2*p+1+j)]-2*(uxz[i,((i-1)*p+1+j)]-ux[i]*alpha[j]-uxx[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[(2*r+1),(2*r+1)] <- -sum(u)/sigmax
  QI[(2*r+1),(2*p)] <- -(sum(ux)-mux*sum(u))/sigmax^2
  QI[(2*p),(2*r+1)] <- -(sum(ux)-mux*sum(u))/sigmax^2
  QI[(2*p),(2*p)] <- n/(2*sigmax^2)-(sum(uxx)-2*mux*sum(ux)+mux^2*sum(u))/sigmax^3
  QI[(2*p+1):k,1:r] <- t(QI[1:r,(2*p+1):k])
  QI[(2*p+1):k,p:(2*r)] <- t(QI[p:(2*r),(2*p+1):k])
  for (j in 1:p){
    QI[(2*p+j),(2*p+j)] <- 0
    for (i in 1:n){
      if(j==1){
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)]+1/phi[j]^2-3*(uzz[((i-1)*p+j),((i-1)*p+j)]-2*uxz[i,((i-1)*p+1)]+uxx[i])/(phi[j]^4)
      }
      else{
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)]+1/phi[j]^2-3*(uzz[((i-1)*p+j),((i-1)*p+j)]-2*uz[((i-1)*p+j),i]*alpha[j-1]-2*uxz[i,((i-1)*p+j)]*beta[j-1]+u[i]*alpha[j-1]^2+2*ux[i]*alpha[j-1]*beta[j-1]+uxx[i]*beta[j-1]^2)/(phi[j]^4)
      }
    }
  }
  return(QI)
}

###########################################################
## Influ?ncia local (Perturba??o de pondera??o de casos) ##
###########################################################

Q2.t<-Q2MEMC.t(ECM)
DeltaPond.t<-deltaPondMEMC.t(ECM)
If.t=t(DeltaPond.t)%*%solve(-Q2.t)%*%DeltaPond.t
Mo.t<-diag(If.t)/sum(diag(If.t))
xrange.t <- c(0,length(ECM$tui))
yrange.t <- c(0,0.187)
plot(Mo.t,xlim=xrange.t,ylim=yrange.t,ylab="M(0)",xlab="?ndice",main="Pondera??o de Casos: MEMC-t")
abline(h=mean(Mo.t)+3*sd(Mo.t),lty=2)
identify(Mo.t, n=2)

##################################################
## B. PERTURBA??O NA COVARI?VEL MEDIDA COM ERRO ##
##################################################

# 1ro rodar o ECM para usar as estimativas dos par?metros e as matrizes geradas no ECM

######################################################################
# MATRIZ DELTA: 1ra derivada da fun??o Q (Perturba??o na covari?vel) #
######################################################################

deltaPCovMEMC.t<-function(est){
  
  alpha <- est$alpha1
  beta <- est$beta1
  mux <- est$mux
  sigmax <- est$sigmax
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  u <- est$tui
  ux <- est$tuxi
  d_omega <- est$dd
  phi <- sqrt(d_omega)
  uz <- est$tuzi
  n <- length(u)
  delta <- matrix(0,k,n)
  iomega22 <- solve(diag(d_omega[-1]))
  for (i in 1:n){
    uzi <- uz[((i-1)*p+1):(i*p),i]
    uz_s <- as.matrix(uzi[-1])
    delta[1:r,i] <- -iomega22%*%beta*u[i]
    delta[(r+1):(2*r),i] <- iomega22%*%(uz_s-u[i]*alpha-2*ux[i]*beta)
    delta[2*r+1,i] <- u[i]/sigmax
    delta[2*r+2,i] <- (ux[i]-mux*u[i])/sigmax^2
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -2*(uzi[j]-ux[i])/phi[j]^3
      }
      else{
        delta[2*p+j,i] <- -2*(uzi[j]*beta[j-1]-u[i]*alpha[j-1]*beta[j-1]-ux[i]*beta[j-1]^2)/phi[j]^3
      }
    }
  }
  return(delta)
}