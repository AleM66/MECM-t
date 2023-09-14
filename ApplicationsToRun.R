#####################################################################################################
# INFLUÊNCIA LOCAL: MEMC-t e MEMC-N (Modelo com erro de medida censurado baseado na dist. t/Normal) #
#####################################################################################################

############################
## Algoritmo ECM (MEMC-t) ##
############################

#rm(list=ls(all=TRUE))
require(mvtnorm)

ECM.t<-function(cc,Z,X,r,nu){   # OK
  t <- proc.time() # Linha para medir o tempo de ejecução com a função proc.time(). Aqui inicia o cronómetro
  GB = GenzBretz(maxpts = 50000, abseps = 1e-9, releps = 0)   # função do pacote mvtnorm
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
  iomega<-solve(omega)      # inversa de omega
  mux <- 4	              	# média de X
  sigmax <- 25		          # var(X)
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
    
    tui <- rep(0,n)         # vetor com c/u dos ui's.
    tuxi <- rep(0,n)        # vetor com c/u dos uxi's.
    tuxxi <- rep(0,n)       # vetor com c/u dos uxxi's (ux^2).
    tuxzi <- matrix(0,n,N)  # matriz com c/u dos uxzi's. nxN (exemplo em documento base.tex)
    tuzzi <- matrix(0,N,N)  # matriz com c/u dos uzzi's. NxN (exemplo em documento base.tex)
    tuzi <- matrix(0,N,n)   # matriz com c/u dos uzi's. Nxn (exemplo em documento base.tex)
    ver <- matrix(0,n,1)    # vetor com c/u das verossimilhan?as. nx1
    
    for (j in 1:n ){
      cc1 <- cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] # j=1=>cc1=cc[1:5], j=2=>cc1=cc[6:10], etc
      Z1 <- Z[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]  # j=1=>Z1=Z[1:5], j=2=>Z1=Z[6:10], etc
      
      muZ <- a1+b1*mux                     # média de Z
      sigmaZ <- sigmax*b1%*%t(b1)+omega    # matriz de escala (dispersão) de Z. Cov(Z)=Matriz de escala*nu/(nu-2)
      dm <- t(Z1-muZ)%*%solve(sigmaZ)%*%(Z1-muZ) # Delta de Z (Dist. de Mahalanobis)
      cdm <- as.numeric((nu+p)/(nu+dm))
      Lambda <- as.numeric(sigmax/(1+sigmax*t(b1) %*% iomega %*% b1))
      Phi <- Lambda*t(b1) %*% iomega
      
      if(sum(cc1)==0){                     # Se não tem censura
        tu <- cdm                          # estimador de ui
        tuz <- matrix(Z1,nj[j],1)*cdm      # estimador de uzi
        tuzz <- Z1%*%t(Z1)*cdm             # estimador de uz^2i=E[UiZi.t(Zi)|Vi]
        tux <- mux*tu+Phi%*%(tuz-muZ*tu)   # estimador de uxi
        tuxx <- Lambda+mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi)
        tuxz <- mux*t(tuz)+Phi%*%(tuzz-muZ%*%t(tuz))  # estimador de uxzi=E[Uixi.t(Zi)|Vi]
        ver[j,] <- dmvt(as.vector(Z1),as.vector(muZ),as.matrix(sigmaZ),df=nu,log=FALSE) # Densidade da dist. t Multivariada dmvt(x, mu, sigma, nu, lambda, log=FALSE)
      }
      
      if(sum(cc1)>0){                  # se tem censura
        if(sum(cc1)==nj[j]){           # se todas as componentes são censuradas
          muUi <- muZ
          sigmaUi <- sigmaZ
          sigmaUiA <- sigmaUi*nu/(nu+2)
          auxupper <- Z1-muUi            # pmvt Calcula a FDA da t multivariada para limites, gl e matriz de correlação arbitraria baseado no algoritmo de Genz e Bretz
          auxU1 <- pmvt(upper=c(auxupper), sigma=sigmaUiA, df=nu+2, algorithm=GB)[1]   # FDA da t.
          auxU2 <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu, algorithm=GB)[1]   # FDA da t.
          MoMT <- Mtmvt(muUi, sigmaUiA, nu+2, rep(-Inf, nj[j]), Z1)
          U0 <- as.numeric(auxU1/auxU2)
          U1 <- auxU1/auxU2*MoMT$Ey
          U2 <- auxU1/auxU2*MoMT$Eyy
          
          tu <- U0               # estimador de ui.
          tuz <- U1              # estimador de uzi
          tuzz <- U2             # estimador de uzzi
          tux <- mux*tu + Phi%*%(tuz-muZ*tu)
          tuxx <- Lambda + mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi)
          tuxz <- mux*t(tuz) + Phi%*%(tuzz-muZ%*%t(tuz))
          ver[j,] <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu, algorithm=GB)[1]
        }
        else {
          Psi <- sigmaZ
          PsiA <- Psi*nu/(nu+2)
          nu1 <- (nu+length(cc1[cc1==0]))   # nu + p0.
          
          muc <- muZ[cc1==1]+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0]) # muZ^co
          Sc <- Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1] # sigmaZ^cc.o
          ScA <- nu/(nu+2)*Sc
          
          Qz1 <- t(Z1[cc1==0]-muZ[cc1==0])%*%solve(Psi[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0]) # dist. Mahalanobis1 t(Zi^o-muZ^o)*inv(sigmaZ^oo)*(Zi^o-muZ^o)
          Qz2 <- t(Z1[cc1==0]-muZ[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0])# dist. Mahalanobis2 t(Zi^o-muZ^o)*inv(nu*sigmaZ^oo/(nu+2))*(Zi^o-muZ^o)
          
          auxcte <- as.numeric((nu+Qz1)/(nu+length(cc1[cc1==0])))       # (nu+dist. Mahalanobis1)/(nu+p0)
          auxcte1 <- as.numeric((nu+2+Qz2)/(nu+2+length(cc1[cc1==0])))  # (nu+2+dist. Mahalanobis2)/(nu+p0+2)
          
          Sc22 <- auxcte*Sc         # Sz^co
          muUi <- muc               # muZ^co
          sigmaUi <- Sc22           # Sz^co
          sigmaUiA <- auxcte1*ScA   # S~z^co
          auxupper <- Z1[cc1==1]-muUi
          
          auxU1 <- pmvt(upper=c(auxupper), sigma=sigmaUiA, df=nu1+2, algorithm=GB)[1]
          auxU2 <- pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu1, algorithm=GB)[1]
          MoMT <- Mtmvt(muUi, sigmaUiA, nu1+2, rep(-Inf,length(cc1[cc1==1])), Z1[cc1==1])
          
          U0 <- as.numeric(auxU1/auxU2)/auxcte   # ui
          U1 <- U0 * MoMT$Ey                     # ui*E[Y]
          U2 <- U0 * MoMT$Eyy                    # ui*E[YY']
          
          Auxtuz <- matrix(Z1,nj[j],1)           # Zi
          
          tuz <- Auxtuz*U0                       # ui * Zi
          tuz[cc1==1] <- U1                      # uzi^c = ui*E[Y]
          tuzz <- Auxtuz%*%t(Auxtuz)             # ZZi = Zi*t(Zi)
          
          AAx <- tuzz[cc1==0,cc1==0]*U0          # ui * ZZi^o
          ABx <- Auxtuz[cc1==0]%*%t(U1)          # Zi^o * ui * E[Y']
          BAx <- t(ABx)                          # ui * E[Y] * t(Zi^o)
          BBx <- U2                              # ui*E[YY']
          
          tuzz[cc1==0,cc1==0] <- AAx             # ui * ZZi^o
          tuzz[cc1==0,cc1==1] <- ABx             # ui * Zi^o * E[Y']
          tuzz[cc1==1,cc1==0] <- BAx             # ui * E[Y] * t(Zi^o)
          tuzz[cc1==1,cc1==1] <- BBx             # ui*E[YY']
          
          tu <- U0                               # ui
          tux <- mux*tu + Phi%*%(tuz-muZ*tu)     # uxi
          tuxx <- Lambda + mux^2*tu + 2*mux*Phi%*%(tuz-muZ*tu) + Phi%*%(tuzz-tuz%*%t(muZ)-muZ%*%t(tuz)+muZ%*%t(muZ)*tu)%*%t(Phi) # uxi^2
          tuxz <- mux*t(tuz)+Phi%*%(tuzz-muZ%*%t(tuz))
          
          sigmaUi <- (sigmaUi+t(sigmaUi))/2
          Psi[cc1==0,cc1==0] <- (Psi[cc1==0,cc1==0]+t(Psi[cc1==0,cc1==0]))/2
          ver[j,] <- dmvt(Z1[cc1==0],muZ[cc1==0],as.matrix(Psi[cc1==0,cc1==0]),df=nu,log=FALSE)*pmvt(upper=c(auxupper), sigma=sigmaUi, df=nu1, algorithm=GB)[1]
        }
      }
      tuz_s <- as.matrix(tuz[2:length(cc1)])       # uzi*
      tuxz_s <- tuxz[2:length(cc1)]                # uxzi*
      tuzz_s <- tuzz[2:length(cc1),2:length(cc1)]  # uzzi*
      
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
    
    if(cont==250) criterio<-0.00000001
    
    a1 <- matrix(c(0,alpha1),p,1)
    b1 <- matrix(c(1,beta1),p,1)
    
    teta<-teta1
    print(teta1)
    logver1<-logver
  }
  MI <- matrix(0,3*p,3*p)
  
  for(i in 1:n){
    tu1 <- tui[i]
    tux1 <- tuxi[i]
    tuxx1 <- tuxxi[i]
    tuxz1 <- t(as.matrix(tuxzi[i,(sum(nj[1:i-1])+1) : (sum(nj[1:i]))]))
    tuzz1 <- tuzzi[(sum(nj[1:i-1])+1) : (sum(nj[1:i])), (sum(nj[1:i-1])+1) : (sum(nj[1:i]))]
    tuz1 <- tuzi[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),i]
    
    ai <- tuzz1-2*tuz1%*%t(a1)-2*t(tuxz1)%*%t(b1)+2*tux1*a1%*%t(b1)+tu1*a1%*%t(a1)+tuxx1*b1%*%t(b1) # eu
    dai <- as.matrix(diag(ai),p,1)
    I1 <- matrix(1,p,1)
    Ip <- matrix(0,p-1,p)
    Ip[,2:p] <- diag(p-1)   # ou matrix(c(rep(0,p-1),diag(p-1)),p-1,p) em vez das 2 ultimas linhas
    dalpha <- Ip%*%iomega%*%(tuz1-tu1*a1-tux1*b1)
    dbeta <- Ip%*%iomega%*%(t(tuxz1)-tux1*a1-tuxx1*b1)
    dmux <- (tux1-tu1*mux)/sigmax
    dsigmax <- -1/(2*sigmax) + (tuxx1 - 2*tux1*mux + tu1*mux^2)/(2*sigmax^2)
    dphi <- (-1/2)*iomega%*%I1 + (1/2)*iomega^2%*%(dai)
    
    si <- as.matrix(c(t(dalpha),t(dbeta),dmux,dsigmax,t(dphi)),3*p,1)
    MI <- MI + si%*%t(si)
  }
  dd<-diag(omega)
  npar<-length(c(teta1))
  ni<-sum(nj)
  loglik<-logver1
  
  AICc<- -2*loglik +2*npar
  AICcorr<- AICc + 2*npar*(npar+1)/(ni-npar-1)
  BICc <- -2*loglik +log(ni)*npar
  
  tempo <- proc.time()-t    # Detiene el cron?metro
  time <- tempo[3]
  horas <- ifelse(time<3600,time/60,time/3600)
  
  if (time < 3600){
    cat(time%/%60, "minutos e", time-(time%/%60)*60, "segundos\n")
  } else if (time < 86400) {
    cat(time%/%3600, "horas e", (time-(time%/%3600)*3600)/60, "minutos\n")
  } else {
    cat(time%/%86400, "dias,",(time-(time%/%86400)*86400)%/%3600, "horas e", (time-(time%/%3600)*3600)%/%60, "minutos\n")}
  
  resul <- list(alpha1=alpha1, beta1=beta1, mux=mux, sigmax=sigmax, dd=dd, loglik=loglik, AIC=AICc, BIC=BICc, AICcorr=AICcorr, iter=cont,
                tui=tui, tuzzi=tuzzi, tuzi=tuzi, tuxi=tuxi, tuxxi=tuxxi, tuxzi=tuxzi, MI=MI, EP = sqrt(diag(solve(MI))), tempo=horas)
  return(resul)
}

## FUNÇÕES AUXILIARES PARA O ALGORITMO ECM (MEMC-t) ##

## Momentos da Distribuição t Truncada ##

## Função de distribuição acumulada da Normal e da t ##

fdaNT<-function(x,mu,sigma2,nu,type="Normal"){
  resp<-matrix(0,length(x),1)
  if(type=="Normal"){
    resp<-pnorm(x,mu,sqrt(sigma2))	# fda da normal no ponto x.
  }
  if(type=="T"){
    z=(x-mu)/sqrt(sigma2)
    resp=pt(z,df=nu)		# fda da t no ponto z=(x-mu)/sqrt(sigma2).
  }
  return(resp)
}

# Calcula os dois primeiros momentos quando mu=0. Usa-se a matriz R (modificação da função TT.moment do R)
# A função original TT.moment do R calcula os 2 primeiros momentos da t multivariada truncada

TT.moment = function(a,b,R,nu){
  require(mvtnorm)
  GB = GenzBretz(maxpts = 50000, abseps = 1e-9, releps = 0)
  p = length(a)
  
  if(p==1){
    if(a== -Inf)  a <- -1e12
    if(b== Inf)   b <- 1e12
    
    G1<- 0.5*gamma((nu-1)/2)*nu^(nu/2)/((fdaNT(b,0,1,nu,"T")-fdaNT(a,0,1,nu,"T"))*gamma(nu/2)*gamma(1/2))
    EX<- G1*((nu+a^2)^(-(nu-1)/2)-(nu+b^2)^(-(nu-1)/2))
    EXX<- nu/(nu-2)+G1*(a*(nu+a^2)^(-(nu-1)/2)-b*(nu+b^2)^(-(nu-1)/2))
  }
  else{
    a = ifelse(a==-Inf, rep(-1e12,p),a)
    b = ifelse(b== Inf, rep( 1e12,p),b)
    al0 = pmvt(lower = a, upper = b, sigma = R, df = nu, algorithm = GB)[1]  # fda da t multivariada
    
    ### fdp & fda
    
    la1 = (nu-2)/nu; la2 = (nu-4)/nu   # lambda1 e lambda2 (paper do Ho et al (2012))
    da = (nu-1)/(nu+a^2); db = (nu-1)/(nu+b^2)  # delta(a) e delta(b) (paper do Ho et al (2012))
    f1a = sqrt(la1)*dt(sqrt(la1)*a,df=nu-2)   # primeira parcela de q*(a)
    f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)   # primeira parcela de q*(b)
    f2 = matrix(NA, p, p)
    G1a = G1b = rep(NA, p)
    G2 = matrix(NA, p, p)
    H = matrix(0,p,p)
    
    for(r in 1:(p-1)){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
      G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)
                      ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])
      mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)
                      ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])
      
      for(s in (r+1):p){
        rs = c(r,s)
        pdf.aa = dmvt(c(a[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
        pdf.ab = dmvt(c(a[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
        pdf.ba = dmvt(c(b[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
        pdf.bb = dmvt(c(b[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
        
        if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
        
        if(p>2){
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          mu.aa = c(tmp%*%c(a[r],a[s]))
          mu.ab = c(tmp%*%c(a[r],b[s]))
          mu.ba = c(tmp%*%c(b[r],a[s]))
          mu.bb = c(tmp%*%c(b[r],b[s]))
          daa = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
          dab = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
          dba = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
          dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          
          cdf.aa = ifelse(p==3,pt((b[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((a[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)
                          ,pmvt(lower = a[-rs]-mu.aa, upper = b[-rs]-mu.aa, sigma = R21/daa, df=nu-2, algorithm = GB)[1])
          
          cdf.ab = ifelse(p==3,pt((b[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((a[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)
                          ,pmvt(lower = a[-rs]-mu.ab, upper = b[-rs]-mu.ab, sigma = R21/dab, df=nu-2, algorithm = GB)[1])
          
          cdf.ba = ifelse(p==3,pt((b[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((a[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)
                          ,pmvt(lower = a[-rs]-mu.ba, upper = b[-rs]-mu.ba, sigma = R21/dba, df=nu-2, algorithm = GB)[1])
          
          cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((a[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)
                          ,pmvt(lower = a[-rs]-mu.bb, upper = b[-rs]-mu.bb, sigma = R21/dbb, df=nu-2, algorithm = GB)[1])
        }
        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
      }
    }
    ##Ultima parte do bucle (loop)
    r <- p
    temp = R[-r,r]
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
    G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)
                    ,pmvt(lower=low, upper=upp, sigma=S1/da[r], df=nu-1, algorithm=GB)[1])
    mub = temp*b[r]; low = a[-r]-mub; upp = b[-r]-mub
    G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)
                    ,pmvt(lower=low, upper=upp, sigma=S1/db[r], df=nu-1, algorithm=GB)[1])
    qa = f1a*G1a; qb = f1b*G1b
    EX = c(R %*% (qa-qb))/al0/la1
    H = H / la2
    D = matrix(0,p,p)
    diag(D) = a * qa - b * qb - diag(R%*%H)
    al1 = pmvt(lower=a, upper=b, sigma=R/la1, df=nu-2, algorithm=GB)[1]
    EXX = (al1 * R + R %*% (H + D) %*% R)/al0/la1
  }
  return(list(EX=EX,EXX=EXX))
}

## Calcula os dois primeiros momentos da t truncada quando mu não é 0 e sigma não é R ##

Mtmvt <- function(mu,sigma,nu,lower,upper){   # mu: vetor de médias
  p=length(lower)                             # sigma: matriz de dispersão
  if(p==1){
    if(lower== -Inf)  lower <- -1e12
    if(upper== Inf)  upper <- 1e12
    a1<-(lower-mu)/sqrt(sigma)
    b1<-(upper-mu)/sqrt(sigma)
    M <- TT.moment(a1, b1, 1, nu)
    Ey<- mu+sqrt(sigma)*M$EX
    Eyy<- mu^2+sigma*M$EXX+2*mu*sqrt(sigma)*M$EX
    Vary<- Eyy - Ey^2
  }
  else{
    Lambda <- diag(1/sqrt(diag(sigma)))   # inversa de Lambda (paper de Ho et al, 2012).
    
    if(length(which(upper == Inf)) != 0)  upper[which(upper == Inf)] <- 1e12	# substituir as entradas com Inf do vetor upper pelo valor 1e12
    b <- as.vector(Lambda %*% (upper - mu))
    
    if(length(which(lower == -Inf)) != 0) lower[which(lower == -Inf)] <- -1e12
    a <- as.vector(Lambda %*% (lower - mu))
    
    R <- Lambda %*% sigma %*% Lambda	# Cálculo da matriz R
    
    M <- TT.moment(a, b, R, nu)
    Ey <- mu + solve(Lambda) %*% M$EX
    Eyy <- mu %*% t(mu) + solve(Lambda) %*% M$EX %*% t(mu) + mu %*% t(M$EX) %*% solve(Lambda) + solve(Lambda) %*% M$EXX %*% solve(Lambda)
    Vary<- Eyy- Ey%*%t(Ey)
  }
  return(list(Ey=Ey,Eyy=Eyy,Vary=Vary))
}

################################
## ESQUEMAS DE PERTURBAÇÃO: t ##
################################

## A. PERTURBAÇÃO DE PONDERAÇÃO DE CASOS: t ##

## MATRIZ DELTA: t ##

deltaPondMEMC.t<-function(est){    # est: saida do EM para o MEMC_t     (OK)
  
  alpha <- est$alpha1                   # vetor alpha (rx1)
  beta <- est$beta1                     # vetor beta (rx1)
  mux <- est$mux                        # média de x (1x1)
  sigmax <- est$sigmax                  # variância de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  uz <- est$tuzi                        # matriz (npxn) formada pelos uzi (px1)
  uxz <- est$tuxzi                      # matriz (nxnp) formada pelos uxzi (1xp)
  uzz <- est$tuzzi                      # matriz (npxnp) formada pelos uzzi (pxp)
  n <- ncol(uz)
  delta <- matrix(0,k,n)
  u <- est$tui                          # vetor ui (nx1)
  ux <- est$tuxi                        # vetor uxi (nx1)
  uxx <- est$tuxxi                      # vetor uxxi (nx1)
  d_omega <- est$dd                     # diagonal da matriz omega (pxp)
  phi <- sqrt(d_omega)                  # vetor phi=(phi1, phi2, phi3, phi4, phi5)
  iomega22 <- solve(diag(d_omega[-1]))  # inversa da matriz omega22 (rxr)
  for (i in 1:n){
    uzi <- uz[((i-1)*p+1):(i*p),i]                   # matriz (px1)
    uz_s <- as.matrix(uzi[-1])                       # matriz uzi* (rx1) 
    uxzi <- uxz[i,((i-1)*p+1):(i*p)]                 # vetor linha (1xp)
    uxz_s <- t(as.matrix(uxzi[-1]))                  # vetor linha uxzi* (1xr)
    uzzi <- uzz[((i-1)*p+1):(i*p),((i-1)*p+1):(i*p)] # matriz (pxp)
    delta[1:r,i] <- iomega22%*%(uz_s-u[i]*alpha-ux[i]*beta)               # matriz delta.alpha (r*n)
    delta[(r+1):(2*r),i] <- iomega22%*%(t(uxz_s)-ux[i]*alpha-uxx[i]*beta) # matriz delta.beta (r*n)
    delta[2*r+1,i] <- (ux[i]-mux*u[i])/sigmax                             # matriz delta.mux (1*n)
    delta[2*r+2,i] <- -1/(2*sigmax)+(uxx[i]-2*mux*ux[i]+mux^2*u[i])/(2*sigmax^2) # matriz delta.sigmax (1*n)
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -1/phi[j]+(uzzi[j,j]-2*uxzi[j]+uxx[i])/phi[j]^3 # matriz delta.phi1 (1*n)
      }
      else{
        delta[2*p+j,i] <- -1/phi[j]+(uzzi[j,j]-2*uzi[j]*alpha[j-1]-2*uxzi[j]*beta[j-1]+u[i]*alpha[j-1]^2+2*ux[i]*alpha[j-1]*beta[j-1]+uxx[i]*beta[j-1]^2)/phi[j]^3 # matriz delta.phij (1*n)
      }
    }
  }
  return(delta)
}

## MATRIZ HESSIANA: t ##

Q2MEMC.t<-function(est){
  
  alpha <- est$alpha1                   # vetor alpha (rx1)
  beta <- est$beta1                     # vetor beta (rx1)
  mux <- est$mux                        # média de x (1x1)
  sigmax <- est$sigmax                  # variância de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  uz <- est$tuzi                        # matriz (npxn) formada pelos uzi (px1)
  n <- ncol(uz)
  uxz <- est$tuxzi                      # matriz (nxnp) formada pelos uxzi (1xp)
  uzz <- est$tuzzi                      # matriz (npxnp) formada pelos uzzi (pxp)
  u <- est$tui                          # vetor ui (nx1)
  ux <- est$tuxi                        # vetor uxi (nx1)
  uxx <- est$tuxxi                      # vetor uxxi (nx1)
  d_omega <- est$dd                     # diagonal da matriz omega (pxp)
  phi <- sqrt(d_omega)                  # vetor phi=(phi1, phi2, phi3, phi4, phi5)
  iomega22 <- solve(diag(d_omega[-1]))  # inversa da matriz omega22 (rxr)
  
  QI<-matrix(0,k,k)
  
  QI[1:r,1:r]<- -sum(u)*iomega22           #Q..alpha,alpha
  QI[1:r,p:(2*r)]<- -sum(ux)*iomega22      #Q..alpha,beta
  for (j in 1:r){                          #Q..alpha,phi
    QI[j,(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[j,(2*p+1+j)] <- QI[j,(2*p+1+j)]-2*(uz[((i-1)*p+1+j),i]-u[i]*alpha[j]-ux[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[p:(2*r),1:r]<- -sum(ux)*iomega22           #Q..beta,alpha
  QI[p:(2*r),p:(2*r)]<- -sum(uxx)*iomega22      #Q..beta,beta
  for (j in 1:r){                               #Q..beta,phi
    QI[(r+j),(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[(r+j),(2*p+1+j)] <- QI[(r+j),(2*p+1+j)]-2*(uxz[i,((i-1)*p+1+j)]-ux[i]*alpha[j]-uxx[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[(2*r+1),(2*r+1)] <- -sum(u)/sigmax                 #Q..mux,mux
  QI[(2*r+1),(2*p)] <- -(sum(ux)-mux*sum(u))/sigmax^2   #Q..mux,sigmax
  QI[(2*p),(2*r+1)] <- -(sum(ux)-mux*sum(u))/sigmax^2   #Q..sigmax,mux
  QI[(2*p),(2*p)] <- n/(2*sigmax^2)-(sum(uxx)-2*mux*sum(ux)+mux^2*sum(u))/sigmax^3   #Q..sigmax,sigmax
  QI[(2*p+1):k,1:r] <- t(QI[1:r,(2*p+1):k])             #Q..phi,alpha
  QI[(2*p+1):k,p:(2*r)] <- t(QI[p:(2*r),(2*p+1):k])     #Q..phi,beta
  for (j in 1:p){                                       #Q..phi,phi
    QI[(2*p+j),(2*p+j)] <- 0
    for (i in 1:n){
      if(j==1){     # 1ra entrada de Q..phi,phi   #n ou 1 (corresponde 1)
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)]+1/phi[j]^2-3*(uzz[((i-1)*p+j),((i-1)*p+j)]-2*uxz[i,((i-1)*p+1)]+uxx[i])/(phi[j]^4)
      }
      else{     # resto de entradas de Q..phi,phi #n ou 1 (corresponde 1)
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)]+1/phi[j]^2-3*(uzz[((i-1)*p+j),((i-1)*p+j)]-2*uz[((i-1)*p+j),i]*alpha[j-1]-2*uxz[i,((i-1)*p+j)]*beta[j-1]+u[i]*alpha[j-1]^2+2*ux[i]*alpha[j-1]*beta[j-1]+uxx[i]*beta[j-1]^2)/(phi[j]^4)
      }
    }
  }
  return(QI)
}

## B. PERTURBAÇÃO NA COVARIÁVEL MEDIDA COM ERRO: t ##

## MATRIZ DELTA: t ##

deltaPCovMEMC.t<-function(est){    # est: saida do EM para o MEMC_t     (OK)
  alpha <- est$alpha1                   # vetor alpha (rx1)
  beta <- est$beta1                     # vetor beta (rx1)
  mux <- est$mux                        # média de x (1x1)
  sigmax <- est$sigmax                  # variância de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  u <- est$tui                         # vetor ui (nx1)
  ux <- est$tuxi                       # vetor uxi (nx1)
  d_omega <- est$dd                    # diagonal da matriz omega (pxp)
  phi <- sqrt(d_omega)
  uz <- est$tuzi                       # matriz (npxn) formada pelos uzi (px1)
  n <- length(u)                       # total de sujeitos analisados
  delta <- matrix(0,k,n)
  iomega22 <- solve(diag(d_omega[-1])) # inversa da matriz omega22 (rxr)
  for (i in 1:n){
    uzi <- uz[((i-1)*p+1):(i*p),i]                   # vetor (px1)
    uz_s <- as.matrix(uzi[-1])                       # vetor uzi* (rx1) 
    delta[1:r,i] <- -iomega22%*%beta*u[i]             # matriz delta.alpha (r*n)
    delta[(r+1):(2*r),i] <- iomega22%*%(uz_s-u[i]*alpha-2*ux[i]*beta) # matriz delta.beta (r*n)
    delta[2*r+1,i] <- u[i]/sigmax                    # matriz delta.mux (1*n)
    delta[2*r+2,i] <- (ux[i]-mux*u[i])/sigmax^2      # matriz delta.sigmax (1*n)
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -2*(uzi[j]-ux[i])/phi[j]^3 # matriz delta.phi1 (1*n)
      }
      else{
        delta[2*p+j,i] <- -2*(uzi[j]*beta[j-1]-u[i]*alpha[j-1]*beta[j-1]-ux[i]*beta[j-1]^2)/phi[j]^3 # matriz delta.phi1 (1*n)
      }
    }
  }
  return(delta)
}


####################################################################################
## Algoritmo ECM (Modelo com EM multivariado com censura baseado na dist. NORMAL) ##
####################################################################################

library(mnormt)
library(mvtnorm)

ECM.N<-function(cc,Z,X,r){
  t <- proc.time() # Linha para medir o tempo de ejecu??o com a fun??o proc.time(). Aqui inicia o cron?metro
  n <- length(X)
  nj <- rep(r+1,n)
  p <- r+1
  N <- p*n
  #valores iniciais
  alpha1 <- matrix(c(0.05,-0.25,0.10,1.50),r,1)
  beta1 <- matrix(c(0.90,1,1.15,1.10),r,1)
  a1 <- matrix(c(0,alpha1),p,1)
  b1 <- matrix(c(1,beta1),p,1)
  omega <- diag(p)          # Dp = omega
  iomega<-solve(omega)      # iDp = iomega
  mux <- 4
  sigmax <- 25              # phix = sigmax
  teta <- c(alpha1, beta1, mux, sigmax, diag(omega))
  criterio<-1
  cont<-0
  while(criterio > 0.00005){
    cont <- cont + 1
    print(cont)
    soma1 <- 0
    soma2 <- 0
    soma3 <- 0
    soma5 <- matrix(0,p-1,1)
    soma6 <- matrix(0,p-1,p-1)
    soma7 <- matrix(0,1,p-1)
    soma8 <- matrix(0,p,1)
    soma9 <- matrix(0,p,p)
    soma10 <- matrix(0,1,p)
    #ui <- rep(0,m)
    xi <- rep(0,n)       # vetor com c/u dos xi's.
    xxi <- rep(0,n)      # vetor com c/u dos xxi's (x^2).
    xzi <- matrix(0,n,N) # matriz com c/u dos xzi's. nxN (exemplo em documento base.tex)
    zzi <- matrix(0,N,N) # matriz com c/u dos zzi's. NxN (exemplo em documento base.tex)
    zi <- matrix(0,N,n)  # matriz com c/u dos zi's. Nxn (exemplo em documento base.tex)
    ver <- matrix(0,n,1) # vetor com c/u das verossimilhan?as. nx1
    for (j in 1:n ){
      cc1 <- cc[(sum(nj[1:j-1])+1) : sum(nj[1:j])] # j=1 => cc1=cc[1:5], j=2 => cc1=cc[6:10], etc
      Z1 <- Z[(sum(nj[1:j-1])+1) : sum(nj[1:j])]   # j=1 => Z1=Z[1:5], j=2 => Z1=Z[6:10], etc
      
      muZ <- a1+b1*mux                                    # m?dia de Z
      SigmaZ <- sigmax*b1%*%t(b1)+omega                   # Matriz de covariancias de Z
      Lambda <- as.numeric(sigmax/(1 + sigmax * t(b1)%*%iomega%*%b1))
      Phi <- Lambda*t(b1) %*% iomega
      
      if(sum(cc1)==0){               # Se n?o tem censura
        z <- matrix(Z1,nj[j],1)      # estimador de zi
        zz <- Z1%*%t(Z1)             # estimador de z^2i=E[Zi*t(Zi)|Vi]
        x <- mux + Phi%*%(z-muZ)     # estimador de xi
        xx <- Lambda+mux^2 + 2*mux*Phi%*%(z-muZ) + Phi%*%(zz-z%*%t(muZ)-muZ%*%t(z)+muZ%*%t(muZ))%*%t(Phi)
        xz <- mux*t(z) + Phi%*%(zz-muZ%*%t(z))  # estimador de xzi=E[xi*t(Zi)|Vi]
        ver[j,]<- dmvnorm(Z1,muZ,SigmaZ) # Densidade da dist. normal Multivariada (pacote mvtnorm)
      }
      if(sum(cc1)>=1){                 # se tem censura
        if(sum(cc1)==nj[j]){           # se todas a suas componentes s?o censuradas
          muc <- muZ
          Sc <- SigmaZ
          aux <- MomemNT(muc,Sc,Z1)
          z <- aux$Ey
          zz <- aux$Eyy
          x <- mux + Phi%*%(z-muZ)     # estimador de xi
          xx <- Lambda+mux^2 + 2*mux*Phi%*%(z-muZ) + Phi%*%(zz-z%*%t(muZ)-muZ%*%t(z)+muZ%*%t(muZ))%*%t(Phi)
          xz <- mux*t(z) + Phi%*%(zz-muZ%*%t(z))  # estimador de xzi=E[xi*t(Zi)|Vi]
          ver[j,]<-pmnorm(Z1,muc,Sc) # FDA da dist. normal Multivariada (pacote mnormt)
        }
        else {
          Psi <- SigmaZ
          muc <- muZ[cc1==1]+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(Z1[cc1==0]-muZ[cc1==0])  # muZ^co (m?dia da normal truncada)
          Sc <- Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]  # S_z^co (matriz de covar. da normal truncada)
          aux <- MomemNT(muc,Sc,Z1[cc1==1])
          z <- matrix(Z1,nj[j],1)
          z[cc1==1] <- aux$Ey
          zz <- matrix(0,nj[j],nj[j])#??????
          zz[cc1==1,cc1==1] <- aux$Vary#??????
          zz <- zz+z%*%t(z)#??????
          x <- mux + Phi%*%(z-muZ)     # estimador de xi
          xx <- Lambda+mux^2 + 2*mux*Phi%*%(z-muZ) + Phi%*%(zz-z%*%t(muZ)-muZ%*%t(z)+muZ%*%t(muZ))%*%t(Phi)
          xz <- mux*t(z) + Phi%*%(zz-muZ%*%t(z))  # estimador de xzi=E[xi*t(Zi)|Vi]
          ver[j,]<-dmnorm(Z1[cc1==0],muZ[cc1==0],Psi[cc1==0,cc1==0])*pmnorm(Z1[cc1==1],t(muc),Sc)
        }
      }
      z_s <- as.matrix(z[2:length(cc1)])  
      xz_s <- xz[2:length(cc1)] 
      zz_s <- zz[2:length(cc1),2:length(cc1)]
      
      soma1 <- as.numeric(soma1 + zz[1,1]-2*xz[1]+xx)
      soma2 <- as.numeric(soma2 + x)
      soma3 <- as.numeric(soma3 + xx)
      soma5 <- soma5 + z_s
      soma6 <- soma6 + zz_s
      soma7 <- soma7 + xz_s
      soma8 <- soma8 + z
      soma9 <- soma9 + zz
      soma10 <- soma10 + xz
      
      xi[j] <- x
      xxi[j] <- xx
      xzi[j, (sum(nj[1:j-1])+1) : sum(nj[1:j])] <- xz
      zzi[(sum(nj[1:j-1])+1) : sum(nj[1:j]), (sum(nj[1:j-1])+1) : sum(nj[1:j])] <- zz
      zi[(sum(nj[1:j-1])+1) : sum(nj[1:j]), j] <- z
    }
    beta1 <- (t(soma7)*n - soma5*soma2)/(soma3*n - soma2*soma2)
    alpha1 <- soma5/n - beta1*soma2/n
    omega1 <- soma1/n
    mux <- soma2/n
    sigmax <- (soma3 - mux*soma2)/n
    omega <- matrix(0,p,p)
    for(k in 1:p){
      if(k==1){	
        omega[k,k] <- omega1
      }
      else{
        omega[k,k] <-(soma9[k,k]-2*soma8[k]*alpha1[k-1]-2*soma10[k]*beta1[k-1]+2*soma2*alpha1[k-1]*beta1[k-1]+n*(alpha1[k-1])^2+soma3*(beta1[k-1])^2)/n
      }
    }
    iomega<-solve(omega)
    teta1 <- c(alpha1,beta1,mux,sigmax,diag(omega))
    logver <- sum(log(ver))
    criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
    if(cont==250) criterio<-0.00000001
    a1 <- matrix(c(0,alpha1),p,1)
    b1 <- matrix(c(1,beta1),p,1)
    teta<-teta1
    print(teta1)
    logver1<-logver
  }
  MI <- matrix(0,3*p,3*p)  
  for(i in 1:n){
    x1 <- xi[i] 
    xx1 <- xxi[i] 
    xz1 <- t(as.matrix(xzi[i, (sum(nj[1:i-1])+1) : sum(nj[1:i])]))
    zz1 <- zzi[(sum(nj[1:i-1])+1) : sum(nj[1:i]), (sum(nj[1:i-1])+1) : sum(nj[1:i])] 
    z1 <- zi[(sum(nj[1:i-1])+1) : sum(nj[1:i]), i]
    ai <- zz1-2*z1%*%t(a1)-2*t(xz1)%*%t(b1)+2*x1*a1%*%t(b1)+a1%*%t(a1)+xx1*b1%*%t(b1)
    dai <- as.matrix(diag(ai),p,1)
    I1 <- matrix(1,p,1) 
    Ip <- matrix(0,p-1,p)
    Ip[,2:p] <- diag(p-1)
    dalpha <- Ip%*%iomega%*%(z1-a1-x1*b1)
    dbeta <- Ip%*%iomega%*%(t(xz1)-x1*a1-xx1*b1)
    dmux <- (x1-mux)/sigmax
    dsigmax <- -1/(2*sigmax) + (xx1-2*x1*mux+mux^2)/(2*sigmax^2)
    dphi <- -iomega%*%I1/2 + iomega^2%*%dai/2 
    si <- as.matrix(c(t(dalpha),t(dbeta),dmux,dsigmax,t(dphi)),3*p,1)
    MI <- MI + si%*%t(si)
  }
  dd<-diag(omega)
  npar<-length(c(teta1))
  ni<-sum(nj)
  loglik<-logver1
  AICc<- -2*loglik +2*npar
  AICcorr<- AICc + 2*npar*(npar+1)/(ni-npar-1)
  BICc <- -2*loglik +log(ni)*npar
  tempo <- proc.time()-t    # Detiene el cron?metro
  time <- tempo[3]
  horas <- ifelse(time<3600,time/60,time/3600)
  if (time < 3600){
    cat(time%/%60, "minutos e", time-(time%/%60)*60, "segundos\n")
  } else if (time < 86400) {
    cat(time%/%3600, "horas e", (time-(time%/%3600)*3600)/60, "minutos\n")
  } else {
    cat(time%/%86400, "dias,",(time-(time%/%86400)*86400)%/%3600, "horas e", (time-(time%/%3600)*3600)%/%60, "minutos\n")}
  resul <- list(alpha1=alpha1, beta1=beta1, mux=mux, sigmax=sigmax, dd=dd,  loglik=loglik, AIC=AICc, BIC=BICc, AICcorr=AICcorr, iter=cont,
                zzi=zzi, zi=zi, xi=xi, xxi=xxi, xzi=xzi, MI=MI, EP=sqrt(diag(solve(MI))), tempo=tempo)
  return(resul)
}

## FUNÇÃO AUXILIAR PARA O ALGORITMO ECM (MEMC-N) ##

## Momentos da Normal Truncada ##

MomemNT<-function(u=c(0,0),S=diag(2),qc=c(1,2)) {
  nic=length(u)
  if (nic==1) {
    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)
    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2
  }
  else {
    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)  # FDA da normal multivariada (pacote mvtnorm).
    dd <- rep(0, nic)   #derivative vector
    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]  #drop=F mantem o objeto como matriz, quando ele ? formado por uma linha/coluna
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V) #fdp da normal
    }
    H <- matrix(0, nrow=nic, ncol=nic)
    RH <- matrix(0, nic, nic)
    if(nic==2)     {
      H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }
    else {
      for( s in 1:(nic-1)){
        for (t in (s+1):nic){
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
        }
      }
    }
    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)
  }
  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))
}

##########################################################################################
# INFLUÊNCIA LOCAL: MEMC-N (Modelo com erro de medida censurado baseado na dist. Normal) #
##########################################################################################

## A. PERTURBAÇÃO DE PONDERAÇÃO DE CASOS: NORMAL ##

## MATRIZ DELTA: Normal ##

deltaPondMEMC.N<-function(est){    # est: saida do EM para o MEMC_N     (OK)???
  alpha <- est$alpha1                   # vetor alpha (rx1)
  beta <- est$beta1                     # vetor beta (rx1)
  mux <- est$mux                        # m?dia de x (1x1)
  sigmax <- est$sigmax                  # vari?ncia de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  z <- est$zi                           # matriz (npxn) formada pelos zi (px1)
  xz <- est$xzi                         # matriz (nxnp) formada pelos xzi (1xp)
  zz <- est$zzi                         # matriz (npxnp) formada pelos zzi (pxp)
  n <- ncol(z)
  delta <- matrix(0,k,n)
  x <- est$xi                           # vetor xi (nx1)
  xx <- est$xxi                         # vetor xxi (nx1)
  d_omega <- est$dd                     # diagonal da matriz omega (d_omega=(phi1^2,...,phip^2))
  iomega22 <- solve(diag(d_omega[-1]))  # inversa da matriz omega22 (rxr)
  for (i in 1:n){
    zi <- z[((i-1)*p+1):(i*p),i]                   # matriz (px1)
    z_s <- as.matrix(zi[-1])                       # matriz zi* (rx1) 
    xzi <- xz[i,((i-1)*p+1):(i*p)]                 # vetor linha (1xp)
    xz_s <- t(as.matrix(xzi[-1]))                  # vetor linha xzi* (1xr)
    zzi <- zz[((i-1)*p+1):(i*p),((i-1)*p+1):(i*p)] # matriz (pxp)
    delta[1:r,i] <- iomega22%*%(z_s-alpha-x[i]*beta)  # matriz delta.alpha (r*n)
    delta[(r+1):(2*r),i] <- iomega22%*%(t(xz_s)-x[i]*alpha-xx[i]*beta) # matriz delta.beta (r*n)
    delta[2*r+1,i] <- (x[i]-mux)/sigmax # matriz delta.mux (1*n)
    delta[2*r+2,i] <- -1/(2*sigmax)+(xx[i]-2*mux*x[i]+mux^2)/(2*sigmax^2) # matriz delta.sigmax (1*n)
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -1/sqrt(d_omega[j])+(zzi[j,j]-2*xzi[j]+xx[i])/sqrt(d_omega[j])^3 # matriz delta.phi1 (1*n)
      }
      else{
        delta[2*p+j,i] <- -1/sqrt(d_omega[j])+(zzi[j,j]-2*zi[j]*alpha[j-1]-2*xzi[j]*beta[j-1]+alpha[j-1]^2+2*x[i]*alpha[j-1]*beta[j-1]+xx[i]*beta[j-1]^2)/sqrt(d_omega[j])^3 # matriz delta.phij (1*n)
      }
    }
  }
  return(delta)
}

## MATRIZ HESSIANA: Normal ##

Q2MEMC.N<-function(est){   # Revisar as parcelas derivadas com phi
  alpha <- est$alpha1                   # vetor alpha (rx1)
  beta <- est$beta1                     # vetor beta (rx1)
  mux <- est$mux                        # m?dia de x (1x1)
  sigmax <- est$sigmax                  # vari?ncia de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  z <- est$zi                           # matriz (npxn) formada pelos zi (px1)
  n <- ncol(z)
  xz <- est$xzi                         # matriz (nxnp) formada pelos xzi (1xp)
  zz <- est$zzi                         # matriz (npxnp) formada pelos zzi (pxp)
  x <- est$xi                           # vetor xi (nx1)
  xx <- est$xxi                         # vetor xxi (nx1)
  d_omega <- est$dd                     # diagonal da matriz omega (pxp)
  phi <- sqrt(d_omega)                  # vetor phi=(phi1, phi2, phi3, phi4, phi5)
  iomega22 <- solve(diag(d_omega[-1]))  # inversa da matriz omega22 (rxr)
  
  QI<-matrix(0,k,k)
  
  QI[1:r,1:r]<- -n*iomega22               #Q..alpha,alpha
  QI[1:r,p:(2*r)]<- -sum(x)*iomega22      #Q..alpha,beta
  for (j in 1:r){                          #Q..alpha,phi
    QI[j,(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[j,(2*p+1+j)] <- QI[j,(2*p+1+j)]-2*(z[((i-1)*p+1+j),i]-alpha[j]-x[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[p:(2*r),1:r]<- -sum(x)*iomega22           #Q..beta,alpha
  QI[p:(2*r),p:(2*r)]<- -sum(xx)*iomega22      #Q..beta,beta
  for (j in 1:r){                               #Q..beta,phi
    QI[(r+j),(2*p+1+j)] <- 0
    for (i in 1:n){
      QI[(r+j),(2*p+1+j)] <- QI[(r+j),(2*p+1+j)]-2*(xz[i,((i-1)*p+1+j)]-x[i]*alpha[j]-xx[i]*beta[j])/(phi[j+1]^3)
    }
  }
  QI[(2*r+1),(2*r+1)] <- -n/sigmax                 #Q..mux,mux
  QI[(2*r+1),(2*p)] <- -(sum(x)-n*mux)/sigmax^2    #Q..mux,sigmax
  QI[(2*p),(2*r+1)] <- -(sum(x)-n*mux)/sigmax^2    #Q..sigmax,mux
  QI[(2*p),(2*p)] <- n/(2*sigmax^2)-(sum(xx)-2*mux*sum(x)+n*mux^2)/sigmax^3   #Q..sigmax,sigmax
  QI[(2*p+1):k,1:r] <- t(QI[1:r,(2*p+1):k])             #Q..phi,alpha
  QI[(2*p+1):k,p:(2*r)] <- t(QI[p:(2*r),(2*p+1):k])     #Q..phi,beta
  for (j in 1:p){                                       #Q..phi,phi
    QI[(2*p+j),(2*p+j)] <- 0
    for (i in 1:n){
      if(j==1){     # 1ra entrada de Q..phi,phi    #n ou 1??????
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)] + 1/phi[j]^2 - 3*(zz[((i-1)*p+j),((i-1)*p+j)] - 2*xz[i,((i-1)*p+1)] + xx[i])/(phi[j]^4)
      }
      else{     # resto de entradas de Q..phi,phi  #n ou 1??????
        QI[(2*p+j),(2*p+j)] <- QI[(2*p+j),(2*p+j)]+1/phi[j]^2-3*(zz[((i-1)*p+j),((i-1)*p+j)]-2*z[((i-1)*p+j),i]*alpha[j-1]-2*xz[i,((i-1)*p+j)]*beta[j-1]+alpha[j-1]^2+2*x[i]*alpha[j-1]*beta[j-1]+xx[i]*beta[j-1]^2)/(phi[j]^4)
      }
    }
  }
  return(QI)
}

## B. PERTURBAÇÃO NA COVARIÁVEL MEDIDA COM ERRO: NORMAL ##

## MATRIZ DELTA: Normal ##

deltaPCovMEMC.N<-function(est){    # est: saida do EM para o MEMC_N     (OK)
  alpha <- est$alpha1                  # vetor alpha (rx1)
  beta <- est$beta1                    # vetor beta (rx1)
  mux <- est$mux                       # m?dia de x (1x1)
  sigmax <- est$sigmax                 # vari?ncia de x (1x1)
  r <- length(alpha)
  p <- r+1
  k <- 3*p
  x <- est$xi                          # vetor xi (nx1)
  d_omega <- est$dd                    # diagonal da matriz omega (pxp)
  z <- est$zi                          # matriz (npxn) formada pelos zi (px1)
  n <- ncol(z)                         # total de sujeitos analisados
  delta <- matrix(0,k,n)
  iomega22 <- solve(diag(d_omega[-1])) # inversa da matriz omega22 (rxr)
  for (i in 1:n){
    zi <- z[((i-1)*p+1):(i*p),i]               # vetor (px1)
    z_s <- as.matrix(zi[-1])                   # vetor zi* (rx1) 
    delta[1:r,i] <- -iomega22%*%beta           # matriz delta.alpha (r*n)
    delta[(r+1):(2*r),i] <- iomega22%*%(z_s-alpha-2*x[i]*beta) # matriz delta.beta (r*n)
    delta[2*r+1,i] <- 1/sigmax                 # matriz delta.mux (1*n)
    delta[2*r+2,i] <- (x[i]-mux)/sigmax^2      # matriz delta.sigmax (1*n)
    for(j in 1:p){
      if(j==1){
        delta[2*p+j,i] <- -2*(zi[j]-x[i])/sqrt(d_omega[j])^3 # matriz delta.phi1 (1*n)
      }
      else{
        delta[2*p+j,i] <- -2*(zi[j]*beta[j-1]-alpha[j-1]*beta[j-1]-x[i]*beta[j-1]^2)/sqrt(d_omega[j])^3 # matriz delta.phi1 (1*n)
      }
    }
  }
  return(delta)
}

######################################################################
## Influência local (Perturbação de ponderação de casos: Normal, t) ##
######################################################################

# 1ro rodar o ECM para usar as estimativas dos parâmetros e as matrizes geradas no ECM

data<-read.table("https://raw.githubusercontent.com/AleM66/MECM-t/main/dataChipkevitchSC.txt", header=T)
cens <- 0.1 # Para outros porcentagens, modificar aqui: 0.1, 0.3, 0.5, 0.7
kij <- quantile(as.matrix(data), cens)
data[data<=kij] <- kij
X<-data[,1]
Z<-as.vector(t(data))
cc<-rep(0,length(as.matrix(data)))
cc[Z<=kij]<-1

est10 <- ECM.t(cc=cc,Z=Z,X=X,r=4,nu=6)	# Tem que ser rodado antes de rodar o resto.
estN10 <- ECM.N(cc=cc,Z=Z,X=X,r=4)	    # Tem que ser rodado antes de rodar o resto.

# Perturbação de ponderação de casos: Normal
Q2.N<-Q2MEMC.N(estN10)
DeltaPond.N<-deltaPondMEMC.N(estN10)
If.Pond.N=t(DeltaPond.N)%*%solve(-Q2.N)%*%DeltaPond.N
Mo.Pond.N<-diag(If.Pond.N)/sum(diag(If.Pond.N))
xrange.Pond.N <- c(0,ncol(estN10$zi))
yrange.Pond.N <- c(0,max(Mo.Pond.N)+0.01)
plot(Mo.Pond.N, xlim=xrange.Pond.N, ylim=yrange.Pond.N, ylab="M(0)", xlab="Index", main="N-MEC (10%)")
abline(h=mean(Mo.Pond.N)+3*sd(Mo.Pond.N), lty=2)
identify(Mo.Pond.N, n=2)

# Perturbação de ponderação de casos: t
Q2.t<-Q2MEMC.t(est10)
DeltaPond.t<-deltaPondMEMC.t(est10)
If.Pond.t=t(DeltaPond.t)%*%solve(-Q2.t)%*%DeltaPond.t
Mo.Pond.t<-diag(If.Pond.t)/sum(diag(If.Pond.t))
xrange.Pond.t <- c(0,length(est10$tui))
yrange.Pond.t <- c(0,max(Mo.Pond.N)+0.01)
plot(Mo.Pond.t, xlim=xrange.Pond.t, ylim=yrange.Pond.t, ylab="M(0)", xlab="Index", main="t-MEC (10%)")
abline(h=mean(Mo.Pond.t)+3*sd(Mo.Pond.t),lty=2)

#############################################################################
## Influência local (Perturbação da covariável medida com erro: Normal, t) ##
#############################################################################

# Normal
Q2.N<-Q2MEMC.N(estN10)
DeltaPCov.N<-deltaPCovMEMC.N(estN10)
If.Cov.N=t(DeltaPCov.N)%*%solve(-Q2.N)%*%DeltaPCov.N
Mo.Cov.N<-diag(If.Cov.N)/sum(diag(If.Cov.N))
xrange.Cov.N <- c(0,ncol(estN10$zi))
yrange.Cov.N <- c(0,max(Mo.Cov.N)+0.01)
plot(Mo.Cov.N, xlim=xrange.Cov.N, ylim=yrange.Cov.N, ylab="M(0)", xlab="Index", main="N-MEC (10%)")
abline(h=mean(Mo.Cov.N)+3*sd(Mo.Cov.N), lty=2)
identify(Mo.Cov.N, n=1)

# t
Q2.t <- Q2MEMC.t(est10)
DeltaPCov.t <- deltaPCovMEMC.t(est10)
If.Cov.t = t(DeltaPCov.t)%*%solve(-Q2.t)%*%DeltaPCov.t
Mo.Cov.t <- diag(If.Cov.t)/sum(diag(If.Cov.t))
xrange.Cov.t <- c(0,length(est10$tui))
yrange.Cov.t <- c(0,max(Mo.Cov.N)+0.01)
plot(Mo.Cov.t, xlim=xrange.Cov.t, ylim=yrange.Cov.t, ylab="M(0)", xlab="Index",main="t-MEC (10%)")
abline(h=mean(Mo.Cov.t)+3*sd(Mo.Cov.t), lty=2)
