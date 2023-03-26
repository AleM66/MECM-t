##################################################################################
## Modelo com erro de medida multivariado com censura baseado na distribui??o t ##
##################################################################################

## Algoritmo EM ##

require(mvtnorm)

ECM.t<-function(cc,Z,X,r,nu){
  t <- proc.time()
  GB = GenzBretz(maxpts = 50000, abseps = 1e-9, releps = 0)
  n <- length(X)
  nj <- rep(r+1,n)
  p <- r+1
  N <- p*n
  
  #valores iniciais
  alpha1 <- matrix(c(0.05,-0.25,0.10,1.50),r,1)
  beta1 <- matrix(c(0.90,1,1.15,1.10),r,1)
  a1 <- matrix(c(0,alpha1),p,1)
  b1 <- matrix(c(1,beta1),p,1)
  omega <- diag(p)
  iomega<-solve(omega)
  mux <- 4
  sigmax <- 25
  teta <- c(alpha1,beta1,mux,sigmax,diag(omega))
  
  criterio<-1
  cont<-0
  
  while(criterio > 0.00005){
    cont <- cont + 1
    print(cont)
    soma1 <- 0
    soma2 <- 0
    soma3 <- 0
    soma4 <- 0
    soma5 <- matrix(0,p-1,1)
    soma6 <- matrix(0,p-1,p-1)
    soma7 <- matrix(0,1,p-1)
    soma8 <- matrix(0,p,1)
    soma9 <- matrix(0,p,p)
    soma10 <- matrix(0,1,p)
    
    tui <- rep(0,n)
    tuxi <- rep(0,n)
    tuxxi <- rep(0,n)
    tuxzi <- matrix(0,n,N)
    tuzzi <- matrix(0,N,N)
    tuzi <- matrix(0,N,n)
    ver <- matrix(0,n,1)
    
    for (j in 1:n ){
      cc1 <- cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      Z1 <- Z[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      muZ <- a1+b1*mux
      sigmaZ <- sigmax*(b1 %*% t(b1))+omega
      dm <- t(Z1-muZ)%*%solve(sigmaZ)%*%(Z1-muZ)
      cdm <- as.numeric((nu+p)/(nu+dm))
      Lambda <- as.numeric(sigmax/(1+sigmax*t(b1) %*% iomega %*% b1))
      Phi <- Lambda*t(b1) %*% iomega
      
      if(sum(cc1)==0){
        tu <- cdm
        tuz <- matrix(Z1,nj[j],1)*cdm
        tuzz <- Z1%*%t(Z1)*cdm
        tux <- mux*tu+Phi%*%(tuz-muZ*tu)
        tuxx <- Lambda+mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi)
        tuxz <- mux*t(tuz)+Phi%*%(tuzz-muZ%*%t(tuz))
        ver[j,] <- dmvt(as.vector(Z1),as.vector(muZ),as.matrix(sigmaZ),df=nu,log=FALSE)
      }
      
      if(sum(cc1)>0){
        if(sum(cc1)==nj[j]){
          muUi <- muZ
          sigmaUi <- sigmaZ
          sigmaUiA <- sigmaUi*nu/(nu+2)
          auxupper <- Z1-muUi
          auxU1 <- pmvt(upper=c(auxupper), sigma=sigmaUiA, df=nu+2, algorithm=GB)[1]
          auxU2 <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu, algorithm=GB)[1]
          MoMT <- Mtmvt(muUi, sigmaUiA, nu+2, rep(-Inf, nj[j]), Z1)
          U0 <- as.numeric(auxU1/auxU2)
          U1 <- auxU1/auxU2*MoMT$Ey
          U2 <- auxU1/auxU2*MoMT$Eyy
          
          tu <- U0
          tuz <- U1
          tuzz <- U2
          tux <- mux*tu + Phi%*%(tuz-muZ*tu)
          tuxx <- Lambda + mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi)
          tuxz <- mux*t(tuz) + Phi%*%(tuzz-muZ%*%t(tuz))
          ver[j,] <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu, algorithm=GB)[1]
        }
        else {
          Psi <- sigmaZ
          PsiA <- Psi*nu/(nu+2)
          nu1 <- (nu+length(cc1[cc1==0]))
          
          muc <- muZ[cc1==1]+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0])
          Sc <- Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]
          ScA <- nu/(nu+2)*Sc
          
          Qz1 <- t(Z1[cc1==0]-muZ[cc1==0])%*%solve(Psi[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0])
          Qz2 <- t(Z1[cc1==0]-muZ[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0])
          
          auxcte <- as.numeric((nu+Qz1)/(nu+length(cc1[cc1==0])))
          auxcte1 <- as.numeric((nu+2+Qz2)/(nu+2+length(cc1[cc1==0])))
          
          Sc22 <- auxcte*Sc
          muUi <- muc
          sigmaUi <- Sc22
          sigmaUiA <- auxcte1*ScA
          auxupper <- Z1[cc1==1]-muUi
          
          auxU1 <- pmvt(upper=c(auxupper), sigma=sigmaUiA, df=nu1+2, algorithm=GB)[1]
          auxU2 <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu1, algorithm=GB)[1]
          MoMT <- Mtmvt(muUi, sigmaUiA, nu1+2, rep(-Inf,length(cc1[cc1==1])), Z1[cc1==1])
          
          U0 <- as.numeric(auxU1/auxU2)/auxcte
          U1 <- U0 * MoMT$Ey
          U2 <- U0 * MoMT$Eyy
          
          Auxtuz <- matrix(Z1,nj[j],1)
          
          tuz <- Auxtuz*U0
          tuz[cc1==1] <- U1
          tuzz <- Auxtuz%*%t(Auxtuz)
          
          AAx <- tuzz[cc1==0,cc1==0]*U0
          ABx <- Auxtuz[cc1==0]%*%t(U1)
          BAx <- t(ABx)
          BBx <- U2
          
          tuzz[cc1==0,cc1==0] <- AAx
          tuzz[cc1==0,cc1==1] <- ABx
          tuzz[cc1==1,cc1==0] <- BAx
          tuzz[cc1==1,cc1==1] <- BBx
          
          tu <- U0
          tux <- mux*tu + Phi%*%(tuz-muZ*tu)
          tuxx <- Lambda + mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi)
          tuxz <- mux*t(tuz)+Phi%*%(tuzz-muZ%*%t(tuz))
          
          sigmaUi <- (sigmaUi+t(sigmaUi))/2
          Psi[cc1==0,cc1==0] <- (Psi[cc1==0,cc1==0]+t(Psi[cc1==0,cc1==0]))/2
          ver[j,] <- dmvt(Z1[cc1==0],muZ[cc1==0],as.matrix(Psi[cc1==0,cc1==0]),df=nu,log=FALSE)*pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu1, algorithm=GB)[1]
        }
      }
      tuz_s <- as.matrix(tuz[2:length(cc1)])
      tuxz_s <- tuxz[2:length(cc1)]
      tuzz_s <- tuzz[2:length(cc1),2:length(cc1)]
      
      soma1 <- as.numeric(soma1 + tuzz[1,1]-2*tuxz[1]+tuxx)
      soma2 <- as.numeric(soma2 + tux)
      soma3 <- as.numeric(soma3 + tuxx)
      soma4 <- as.numeric(soma4 + tu)
      soma5 <- soma5 + tuz_s
      soma6 <- soma6 + tuzz_s
      soma7 <- soma7 + tuxz_s
      soma8 <- soma8 + tuz
      soma9 <- soma9 + tuzz
      soma10 <- soma10 + tuxz
      
      tui[j] <- tu
      tuxi[j] <- tux
      tuxxi[j] <- tuxx
      tuxzi[j, (sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- tuxz
      tuzzi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- tuzz
      tuzi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j] <- tuz
    }
    
    beta1 <- (t(soma7)*soma4 - soma5*soma2)/(soma3*soma4 - soma2*soma2)
    alpha1 <- soma5/soma4 - beta1*soma2/soma4
    omega1 <- soma1/n
    mux <- soma2/soma4
    sigmax <- (soma3 - mux*soma2)/n
    omega <- matrix(0,p,p)
    for(k in 1:p){
      if(k==1){
        omega[k,k] <- omega1
      }
      else{
        omega[k,k] <- (soma9[k,k]-2*soma8[k]*alpha1[k-1]-2*soma10[k]*beta1[k-1]+2*soma2*alpha1[k-1]*beta1[k-1]+soma4*(alpha1[k-1])^2+soma3*(beta1[k-1])^2)/n
      }
    }
    iomega<-solve(omega)
    
    teta1 <- c(alpha1,beta1,mux,sigmax,diag(omega))
    logver <- sum(log(ver))

    criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
        
    a1 <- matrix(c(0,alpha1),p,1)
    b1 <- matrix(c(1,beta1),p,1)
    
    teta<-teta1
    print(teta1)
    logver1<-logver
  }
  
  dd<-diag(omega)
  npar<-length(c(teta1))
  ni<-sum(nj)
  logver_t<-logver1
  
  AIC_t<- -2*logver_t +2*npar
  AICcorr_t<- AIC_t + 2*npar*(npar+1)/(ni-npar-1)
  BIC_t <- -2*logver_t +log(ni)*npar
  
  tempo <- proc.time()-t

  resul <- list(alpha1=alpha1, beta1=beta1, mux=mux, sigmax=sigmax, dd=dd, logver=logver_t, AIC=AIC_t, BIC=BIC_t, AICcorr=AICcorr_t, iter=cont,
                tui=tui, tuzzi=tuzzi, tuzi=tuzi, tuxi=tuxi, tuxxi=tuxxi, tuxzi=tuxzi, iomega=iomega, tempo=tempo)
  return(resul)
}