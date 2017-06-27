
#########################
#   Hypothesis Testing
#########################

#G is an n*p matrix, each row coresponds to one subject.
#Z is an p*q matrix, each row coresponds to one genetic variant.
FST.GeneSet.test<-function(result.prelim,G,Z,GeneSetID,Gsub.id=NULL,weights=NULL,B=5000,impute.method='fixed'){
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
  MAF<-colMeans(as.matrix(G[,colMeans(center.G^2)!=0]))/2;MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  G<-as.matrix(G[,colMeans(center.G^2)!=0])
  SNP.name<-SNP.list[colMeans(center.G^2)!=0]
  # annotation
  Z<-as.matrix(Z[colMeans(center.G^2)!=0,])#matrix(Z[colMeans(center.G^2)!=0,],ncol(G),ncol(Z))
  # make negative scores positive
  if(min(Z,na.rm=T)<0){
    msg<-sprintf('Functional scores contain negative values, x-min(x) is applied to that score')
    warning(msg,call.=F)
    positive<-function(x){if (min(x,na.rm=T)<0){x<-x-min(x,na.rm=T)};return(x)}
    Z<-apply(Z,2,positive)
  }
  # remove na in functional scores
  na.index<-apply(is.na(Z),1,sum)!=0
  G<-Matrix(G[,!na.index]);Z<-as.matrix(Z[!na.index,])

  p<-ncol(G);q<-ncol(Z)

  if(ncol(G) == 0){
    msg<-sprintf("G does not include any heterozygous variant")
    stop(msg)
  }

  #
  GeneSetID<-GeneSetID[match(colnames(G),GeneSetID[,2]),]

  #inference
  FST.fit<-FST.GeneSet.score(nullglm,Y,X0,G,out_type,GeneSetID,B)
  re.score<-FST.fit$re.score
  #B.norm<-matrix(rnorm(B*p,0,1),B,p)
  #re.score<-mvrnorm(B,rep(0,p),FST.fit$Sigma)

  if(length(weights)>1){weighted.Z<-Z*weights[colMeans(center.G^2)!=0]}else{
    weights<-dbeta(apply(G,2,mean,na.rm=T)/2,1,25);weighted.Z<-Z*weights}

  ###
  #   Calculate p-values and resampled p-values
  ###

  #calculate scores for all functional scores
  temp.score<-weighted.Z*as.vector(FST.fit$score)
  #calculate test statistics for all functional scores
  Q1<-apply(temp.score^2,2,sum)/n;Q2<-apply(temp.score,2,sum)^2/n

  #weight.Q1<-1/apply(apply(G^2,2,sum)*weighted.Z^2,2,sum)
  #weight.Q2<-1/apply((G%*%weighted.Z)^2,2,sum)

  #Qsum1<-sum(Q1*weight.Q1)
  #Qsum2<-sum(Q2*weight.Q2)
  #Qsum<-Qsum1+Qsum2

  #resampling method to generate statistics
  re.Q1<-(re.score^2)%*%(weighted.Z^2)/n;re.Q2<-(re.score%*%weighted.Z)^2/n

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

  return(list(n.marker=p,p.value=p.value))
}

#################################################
#   Calculate Score Statistics and their Variance
#################################################

FST.GeneSet.score<-function(nullglm,Y,X0,G,out_type,GeneSetID,B){
  n<-nrow(Y);p<-ncol(G);X0<-svd(X0)$u
  mu<-nullglm$fitted.values;Y.res<-Y-mu
  if(out_type=="C"){sigma<-summary(nullglm)$dispersion;v<-rep(sigma,n)}
  if(out_type=="D"){v<-mu*(1-mu)}
  #g<-Matrix(G)
  score<-t(G)%*%Y.res

  re.score<-matrix(NA,B,p)
  n.total<-1;n.rep<-as.numeric(table(as.character(GeneSetID[,1])))
  for (i in 1:length(n.rep)){
    print(i)
    ni<-n.rep[i];index<-n.total:(n.total+ni-1)
    n.total<-n.total + ni
    temp.G<-Matrix(G[,index])
    temp.sigma<-t(v*temp.G)%*%temp.G-(t(v*temp.G)%*%X0)%*%solve(t(v*X0)%*%X0)%*%(t(v*X0)%*%temp.G)
    re.score[,index]<-mvrnorm(B,rep(0,ni),temp.sigma)
  }
  return(list(score=score,re.score=re.score))
}

