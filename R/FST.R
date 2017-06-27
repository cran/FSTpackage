
#########################
#   Hypothesis Testing
#########################

FST.prelim<-function(Y, X=NULL, id=NULL, out_type="C"){
  ##Preliminary
  Y<-as.matrix(Y);n<-nrow(Y)

  if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
  X0<-cbind(rep(1,n),X0)

  if(out_type=="C"){nullglm<-glm(Y~0+X0,family=gaussian)}
  if(out_type=="D"){nullglm<-glm(Y~0+X0,family=binomial)}

  if (length(id)==0){id<-1:n}

  #prepare the intermediate results
  result.prelim<-list(Y=Y,id=id,n=n,X0=X0,nullglm=nullglm,out_type=out_type)
  return(result.prelim)
}

#G is an n*p matrix, each row coresponds to one subject.
#Z is an p*q matrix, each row coresponds to one genetic variant.
FST.test<-function(result.prelim,G,Z,Gsub.id=NULL,weights=NULL,B=5000,impute.method='fixed'){
  ## Load preliminary data
  Y<-result.prelim$Y;X0<-result.prelim$X0
  n<-result.prelim$n;id<-result.prelim$id
  #X<-result.prelim$X #Z.A0=svd(cbind(rep(1,N),X, M.E))$u
  nullglm<-result.prelim$nullglm;out_type<-result.prelim$out_type
  mu<-nullglm$fitted.values;Y.res<-Y-mu

  ## Deal with the genotype
  SNP.list<-colnames(G)
  # match the phenotype and genotype subject ID
  if(length(Gsub.id)==0){Gsub.id<-id}
  G<-as.matrix(G[match(id,Gsub.id),])

  # missing genotype imputation
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,impute.method)
  }
  # genotype
  G<-as.matrix(G);center.G<-t(t(G)-colMeans(G))
  MAF<-colMeans(as.matrix(G[,colMeans(center.G^2)!=0]))/2;#MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  G<-as.matrix(G[,colMeans(center.G^2)!=0])
  SNP.name<-SNP.list[colMeans(center.G^2)!=0]
  # annotation
  Z<-matrix(Z[colMeans(center.G^2)!=0,],ncol(G),ncol(Z))#matrix(Z[colMeans(center.G^2)!=0,],ncol(G),ncol(Z))
  # make negative scores positive
  if(min(Z,na.rm=T)<0){
    msg<-sprintf('Functional scores contain negative values, x-min(x) is applied to that score')
    warning(msg,call.=F)
    positive<-function(x){if (min(x,na.rm=T)<0){x<-x-min(x,na.rm=T)};return(x)}
    Z<-apply(Z,2,positive)
  }
  # remove na in functional scores
  na.index<-apply(is.na(Z),1,sum)!=0
  G<-Matrix(G[,!na.index]);Z<-matrix(Z[!na.index,],ncol(G),ncol(Z));MAF<-as.matrix(MAF[!na.index])
  SNP.name<-SNP.name[!na.index]

  p<-ncol(G);q<-ncol(Z);colnames(G)<-SNP.name

  if(ncol(G) == 0){
    msg<-sprintf("G does not include any heterozygous variant")
    stop(msg)
  }

  #inference
  FST.fit<-FST.score(nullglm,Y,X0,G,out_type)

  if(length(weights)>1){weights<-weights[colMeans(center.G^2)!=0]}else{
    temp.MAF<-apply(G,2,mean,na.rm=T)/2;temp.MAF[temp.MAF>0.5]<-1-temp.MAF[temp.MAF>0.5]
    weights<-dbeta(temp.MAF,1,25)}

  score<-FST.fit$score;Sigma<-FST.fit$Sigma
  
  result<-FST.SummaryStat.test(score,Sigma,Z,weights,B)
  p.single<-result$p.single
  p.value<-result$p.value

  return(list(MAF=MAF,score=score,Sigma=Sigma,n.marker=p,p.single=p.single,p.value=p.value))
}


FST.SummaryStat.test<-function(score,Sigma,Z,weights,B=5000){
  p<-nrow(Sigma);q<-ncol(Z)
  re.score<-mvrnorm(B,rep(0,p),Sigma)
  
  if(min(Z,na.rm=T)<0){
    msg<-sprintf('Functional scores contain negative values, x-min(x) is applied to that score')
    warning(msg,call.=F)
    positive<-function(x){if (min(x,na.rm=T)<0){x<-x-min(x,na.rm=T)};return(x)}
    Z<-apply(Z,2,positive)
  }
  
  weighted.Z<-Z*weights

  ###
  #   Calculate p-values and resampled p-values
  ###

  #calculate scores for all functional scores
  temp.score<-weighted.Z*as.vector(score)
  #calculate test statistics for all functional scores
  Q1<-apply(temp.score^2,2,sum);Q2<-apply(temp.score,2,sum)^2

  #weight.Q1<-1/apply(apply(G^2,2,sum)*weighted.Z^2,2,sum)
  #weight.Q2<-1/apply((G%*%weighted.Z)^2,2,sum)

  #Qsum1<-sum(Q1*weight.Q1)
  #Qsum2<-sum(Q2*weight.Q2)
  #Qsum<-Qsum1+Qsum2

  #resampling method to generate statistics
  re.Q1<-(re.score^2)%*%(weighted.Z^2);re.Q2<-(re.score%*%weighted.Z)^2

  #re.Qsum1<-apply(t(re.Q1)*weight.Q1,2,sum)
  #re.Qsum2<-apply(t(re.Q2)*weight.Q2,2,sum)
  #re.Qsum<-re.Qsum1+re.Qsum2

  #p.sum<-Get.p(c(Qsum1,Qsum2,Qsum),cbind(re.Qsum1,re.Qsum2,re.Qsum))

  #cauculated resampled p-values using resampled test statistics
  #rho.class<-c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
  rho.class<-c(0,1)

  all.p<-array(0,c(B+1,length(rho.class),q))
  for (rho in rho.class){
    Q.rho<-(1-rho)*Q1+rho*Q2#(1-rho)*Q1+rho*Q2
    re.Q.rho<-(1-rho)*re.Q1+rho*re.Q2
    all.p[,which(rho.class==rho),]<-Get.p(rbind(Q.rho,re.Q.rho),re.Q.rho)
  }
  re.p<-all.p[-1,,] # resampled p-values
  temp.p<-all.p[1,,] # p-values

  ###
  #   Combine p-values using Fisher's method
  ###
  p.Fisher<-matrix(0,3,q+1)
  p.Fisher[1:2,1:q]<-temp.p[c(1,length(rho.class)),]
  for (j in 1:q){p.Fisher[3,j]<-FCombine.p(temp.p[,j],re.p[,,j])}
  for (k in 1:2){p.Fisher[k,q+1]<-FCombine.p(temp.p[1+(k-1)*(length(rho.class)-1),],re.p[,1+(k-1)*(length(rho.class)-1),])}
  p.Fisher[3,q+1]<-FCombine.p(temp.p,re.p)
  ###
  #   Combine p-values using MinP method
  ###
  p.MinP<-matrix(0,3,q+1)
  p.MinP[1:2,1:q]<-temp.p[c(1,length(rho.class)),]
  for (j in 1:q){p.MinP[3,j]<-MCombine.p(temp.p[,j],re.p[,,j])}
  for (k in 1:2){p.MinP[k,q+1]<-MCombine.p(temp.p[1+(k-1)*(length(rho.class)-1),],re.p[,1+(k-1)*(length(rho.class)-1),])}
  p.MinP[3,q+1]<-MCombine.p(as.vector(temp.p),t(apply(re.p,1,as.vector)))

  p.value<-cbind(p.MinP,p.Fisher[,q+1])

  ###
  #   Summerize results
  ###
  #p.value<-rbind(p.dispersion,p.burden,p.DB);p.value<-cbind(p.value,p.unified)
  colnames(p.value)<-c(paste0('Score',1:q),'MinP','Fisher')
  rownames(p.value)<-c('Dispersion','Burden','D+B')
  
  p.single<-pchisq(as.matrix(score^2/diag(Sigma)),1,lower.tail=F)

  return(list(p.single=p.single,p.value=p.value))
}

#################################################
#   Calculate Score Statistics and their Variance
#################################################

FST.score<-function(nullglm,Y,X0,G,out_type){
  n<-nrow(Y);p<-ncol(G);X0<-svd(X0)$u
  mu<-nullglm$fitted.values;Y.res<-Y-mu
  if(out_type=="C"){sigma<-summary(nullglm)$dispersion;v<-rep(sigma,n)}
  if(out_type=="D"){v<-mu*(1-mu)}
  #g<-Matrix(G)
  score<-t(G)%*%Y.res
  Sigma<-t(v*G)%*%G-(t(v*G)%*%X0)%*%solve(t(v*X0)%*%X0)%*%(t(v*X0)%*%G)
  #score<-t(g)%*%Y.res
  #Sigma<-t(v*g)%*%g-(t(v*g)%*%X0)%*%solve(t(v*X0)%*%X0)%*%(t(v*X0)%*%g)

  return(list(score=score,Sigma=Sigma))
}

#cauculated p-values using resampled test statistics
Get.p<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean<-apply(re.Q,2,mean)
  re.variance<-apply(re.Q,2,var)
  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
  return(re.p)
}

#Resample.p<-function(re.Q){ #re.Q is a B*q matrix, each row is one replicate
#  re.mean<-apply(re.Q,2,mean)
#  re.variance<-apply(re.Q,2,var)
#  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
#  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
#  re.p<-t(1-pchisq((t(re.Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df))
#  return(re.p)
#}

# combine p-values
FCombine.p<-function(p,re.p){ #re.p: a B*b*c array; function to combine multipla p-values
  Fisher.stat<--2*sum(log(p))
  re.Fisher<--2*apply(log(re.p),1,sum);re.Fisher[re.Fisher==Inf]<-NA
  Fisher.mean<-mean(re.Fisher,na.rm=T)
  Fisher.variance<-var(re.Fisher,na.rm=T)
  Fisher.kurtosis<-mean((re.Fisher-Fisher.mean)^4,na.rm=T)/Fisher.variance^2-3
  if (Fisher.kurtosis>0){df<-12/Fisher.kurtosis}else{
    df<-100000
  }
  p.combined<-pchisq((Fisher.stat-Fisher.mean)*sqrt(2*df)/sqrt(Fisher.variance)+df,df,lower.tail=F)
  return(p.combined)
}

MCombine.p<-function(p,re.p){ #re.p: a b*c matrix; function to combine multipla p-values for a minP test
  MinP.stat<-min(p)#;re.p[re.p==1]<-0.99
  re.Normal<-qnorm(re.p);re.Normal[re.Normal==Inf]<-NA
  D<-cor(re.Normal,use='complete.obs')
  #diag(D)<-1
  p.combined<-as.numeric(1-pmvnorm(lower=rep(qnorm(MinP.stat),length(p)),sigma=D))
  if (p.combined==0){
    msg<-sprintf('Extreme p-value < 1e-15 is not supported by MinP method, p-value is recorded as 1e-15')
    warning(msg,call.=F)
    p.combined<-1e-15
  }
  return(p.combined)
}


# calculate number of independent tests
effect.n<-function(x,MinP.adjust){
  temp<-0;sum.EV<-sum(x) #summation of eigen values
  for (i in 1:length(x)){
    temp<-temp+x[i];if (temp>sum.EV*MinP.adjust){break}
  }
  return(i)
}

#Data management
Minor.allele<-function(x){if(mean(x)>1){x<-2-x};return(x)}
Standardize<-function(x){return((x-mean(x))/sd(x))}
Center<-function(x){return(x-mean(x))}
Variation<-function(x){x<-apply(x,2,Minor.allele);x<-x[,which(apply(x,2,sd)>0)];return(x)}
Common<-function(x){x<-apply(x,2,Minor.allele);freq<-apply(x,2,mean)/2;x<-x[,which(freq>0.05)];return(x)}
##Matrix calculation
#Get sqrt.root of a matrix
Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Get inverse of a matrix; numerically robust
Get.inverse<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix inverse!")}
  a.inverse <- a.eig$vectors[,ID1] %*% diag(a.eig$values[ID1]^-1) %*% t(a.eig$vectors[,ID1])
  return(a.inverse)
}
#Get -1/2 of a matrix; numerically robust
Get.inv.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])^-1,length(ID1)) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Time similarity
Check.same<-function(x,y){return(as.numeric(x==y))}
Check.near<-function(x,y){return(as.numeric(abs(x-y)==1))}
D.matrix<-function(X){
  n<-length(X[,1]);temp<-matrix(0,n,n)
  #checking time
  temp<-outer(X[,1],X[,1],Check.same)*outer(X[,2],X[,2],Check.near)
  diag(temp)<-0
  return(temp)
}
# Simple Imputation (from SKAT package)
# Z : an m x p genotype matrix with m samples and p SNPs

Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}



