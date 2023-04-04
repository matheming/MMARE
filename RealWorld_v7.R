## a version based on v6

setwd('E:/Learn/Paper_REMA/RealData')
source('ff2n.R')

wages = wooldridge::wage2
wage1 = na.omit(wages)
wage1 = subset(wage1, select = c('wage','educ','married','tenure','age','exper'))

wage1 = cbind(wage1, wage1$educ^2, wage1$exper^2,wage1$tenure^2,wage1$age^2, wage1$age*wage1$married,
              wage1$educ*wage1$married)

wage1$wage = log(wage1$wage)
p = 2

#Cor = numeric(15)
#for (i in 1:15) {
#  Cor[i] = cor(wage1$wage,wage1[,i+1])
#}
#Cor_ind = order(abs(Cor),decreasing = T)
#wage1 = wage1[,c(1,Cor_ind+1)]

Index = ff2n(ncol(wage1)-1-p)
Index = cbind(matrix(1,nrow(Index),1+p),Index)


realworld = function(n_train){
  MSPE=rep(0,7)
  #n_train = 500
  for (rt in 1:100) {
    #print(Sys.time())
    n = nrow(wage1)
    n_test = n - n_train
    #set.seed(seed)
    ind_train = sample(seq(n),n_train)
    wage1_train = wage1[ind_train,]
    wage1_test = wage1[-ind_train,]
    
    
    ## Lasso
    #library(lars)
    #Lars_b1<-lars(as.matrix(wage1_train[,-1]),wage1_train[,1],type = "lar",trace = F)
    
    #Index_b1<-matrix(0,nrow = nrow(Lars_b1$beta),ncol = ncol(Lars_b1$beta))
    #for (i in 1:nrow(Lars_b1$beta)) {
    #  Zeros<-rep(0,ncol(Lars_b1$beta))
    #  Zeros[which(Lars_b1$beta[i,]!=0,arr.ind = TRUE)]<-1
    #  Index_b1[i,]<-Zeros
    #}
    #Index = cbind(matrix(1,nrow(Index_b1),1),Index_b1)
    
    
    wage1_train_X = cbind(rep(1,n_train),wage1_train[,-1])
    wage1_train_Y = wage1_train[,1]
    wage1_test_X = cbind(rep(1,n_test),wage1_test[,-1])
    wage1_test_Y = wage1_test[,1]
    
    #Index = ff2n(ncol(wage1)-1-p)
    #Index = cbind(matrix(1,nrow(Index),1+p),Index)
    # Index = matrix(1,ncol(wage1),ncol(wage1))
    # Index[upper.tri(Index,diag = F)] = 0
    maxs = nrow(Index)
    
    # estimate
    sigmas = numeric(maxs)
    logyhat = matrix(NA,n_train,maxs)
    logy_fore = matrix(NA,n_test,maxs)
    hatmaxtrix= matrix(NA,n_train,maxs)
    for (k in 1:maxs) {
      Ind = which(Index[k,]==1)
      pn = sum(Index[k,])
      newX=as.matrix(wage1_train_X[,Ind],n_train,pn)
      newX_fore = as.matrix(wage1_test_X[,Ind],n_test,pn)
      betahat = solve(t(newX)%*%newX)%*%t(newX)%*%wage1_train_Y
      hatmaxtrix[,k]=diag(newX%*%solve(t(newX)%*%newX)%*%t(newX))
      logyhat[,k] = newX%*%betahat
      logy_fore[,k] = newX_fore%*%betahat
      sigmas[k]=sum((logyhat[,k]-wage1_train_Y)^2)/(n_train-pn)
    }
    # model averaging
    sigmar=sigmas[maxs]
    # print(sigmar)
    logyest=matrix(NA,n_test,7)
    f1=exp(logyhat-matrix(wage1_train_Y,n_train,maxs)+2*sigmar*hatmaxtrix-sigmar)
    f2=exp(logyhat-matrix(wage1_train_Y,n_train,maxs)+sigmar*hatmaxtrix-sigmar/2)
    cv<-function(w){
      mean((f1%*%w)^2)-2*mean(f2%*%w)+1
    }
    w0<-rep(1,maxs)/maxs
    eqfun<-function(w){
      sum(w)-1
    }
    library(Rsolnp)
    #w<-solnp(w0,cv,eqfun = eqfun,eqB = 0,
    #         LB=rep(0,maxs),UB=rep(1,maxs),control = list(trace=0))$pars
    
    library(quadprog)
    a1 = (t(f1)%*%f1)*2/n_train+diag(1e-10,maxs,maxs)
    a2 = (t(f2)%*%rep(1,n_train))*2/n_train
    a3 <- t(rbind(matrix(1,nrow=1,ncol=maxs),diag(maxs),-diag(maxs)))
    a4 <- rbind(1,matrix(0,nrow=maxs,ncol=1),matrix(-1,nrow=maxs,ncol=1))
    QP <- solve.QP(a1,a2,a3,a4,1)
    w <- QP$solution
    w <- w*(w>0)
    w <- w/sum(w0)
    #print(c(which.max(w),max(w)))
    
    logyest[,1]=log(exp(logy_fore)%*%w)+sigmar/2
    
    # AIC and BIC
    sigmamle=sigmas*(n_train-rowSums(Index))
    aic=n_train*log(sigmamle)+2*rowSums(Index)
    i_aic=which(aic==min(aic),arr.ind = T)
    i_aic=i_aic[1]
    logyest[,2]=logy_fore[,i_aic]+sigmas[i_aic]/2
    w_aic=exp(-(aic-min(aic))/2)/sum(exp(-(aic-min(aic))/2))
    logyest[,3]=log(exp(logy_fore)%*%w_aic)+sigmar/2
    bic=n_train*log(sigmamle)+log(n_train)*rowSums(Index)
    i_bic=which(bic==min(bic),arr.ind = T)
    i_bic=i_bic[1]
    logyest[,4]=logy_fore[,i_bic]+sigmas[i_bic]/2
    w_bic=exp(-(bic-min(bic))/2)/sum(exp(-(bic-min(bic))/2))
    logyest[,5]=log(exp(logy_fore)%*%w_bic)+sigmar/2
    # equal weight
    
    logyest[,6]=log(exp(logy_fore)%*%w0)+sigmar/2
    
    logyest[,7]=logy_fore[,maxs]+sigmar/2
    
    #print(c('aic',i_aic))
    #print(c('bic',i_bic))
    
    
    for (i in 1:7) {
      MSPE[i]=MSPE[i]+mean((exp(logyest[,i]-wage1_test_Y)-1)^2)
    }
  }
  MSPE = MSPE/rt
  MSPE
}

realworld(600)
