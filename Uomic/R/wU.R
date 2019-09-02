#'  A U-statistics method for association test using multi-layer omics data.
#'
#' @param Y The continuous or discrete outcomes of subjects.
#' @param Ks A list with each element being a pair-wise similarity matrix derived from a layer of omcis data. 
#' @param X A matrix of demographic variables (e.g., age and gender).
#' @param scheme Specify the weighting schemes (0: average weighting, 1: consensus weighting, 2: permutation-based weighting) 
#' @param B The number of permutations (this is ignored if \code{scheme} is not equal 2). It can be set as null, and the algorithm will determine B.
#' @return The P-value of our proposed method.
#' @export
#'
#' @examples
#' data(Omics)
#' data(outcome)
#' data(Cov)
#' Input<-ComputeKs(Omics,Y=outcome$all,Xcov=Cov,Kernel=c("Lin0","Lin0","IBS"))
#' p_ave<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=0)
#' p_con<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=1)
#' p_perm<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=2,B=1000)
#' Input<-ComputeKs(Omics,Y=outcome$expr,Xcov=Cov,Kernel=c("Lin0","Lin0","IBS"))
#' p_ave<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=0)
#' p_con<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=1)
#' p_perm<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=2,B=1000)
#' Input<-ComputeKs(Omics,Y=outcome$null,Xcov=Cov,Kernel=c("Lin0","Lin0","IBS"))
#' p_ave<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=0)
#' p_con<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=1)
#' p_perm<-wU(Y=Input$Y,Ks=Input$Ks,X=Input$X,scheme=2)
wU<-function(Y,Ks,X=NULL,scheme=0,B=NULL)
{
  n=length(Y);
  if(!is.null(X)) X=cbind(1,X)
  if(is.null(X)) X=matrix(1,nrow=n,ncol=1);
  In=diag(n);
  Pn=X %*% solve(t(X) %*% X) %*% t(X)
  PY=(In-Pn) %*% Y
  R=rank(PY);
  R=scale(R);
  if(length(Ks)<1 | class(Ks)!="list") stop("Ks must be a list with at least one component")
  for(i in 1:length(Ks)){Ks[[i]]=Ks[[i]]/mean(Ks[[i]]^2)  }

  if(scheme==0 | scheme==1)
  {
    if(scheme==0)
    {
      W=Ks[[1]];
      if(length(Ks)>1) {
        for(i in 2:length(Ks)) W=W+Ks[[i]]
      }
      diag(W)=0;
    }
    if(scheme==1)
    {
      X.scaled <- lapply(Ks, function(x) {
        x.cosinus <- sweep(sweep(x, 2, sqrt(diag(x)), "/"), 1, sqrt(diag(x)), "/")
        t(t(x.cosinus - colSums(x.cosinus)/nrow(x.cosinus)) - rowSums(x.cosinus)/nrow(x.cosinus)) + sum(x.cosinus)/nrow(x.cosinus)^2
      })
      similarities <- outer(1:length(X.scaled), 1:length(X.scaled), FUN = Vectorize(function(i, j) {
        tr(X.scaled[[i]] %*% X.scaled[[j]])/(norm(X.scaled[[i]], type = "F") * norm(X.scaled[[j]], type = "F"))}))
      weights <- eigen(similarities, symmetric = TRUE)$vectors[, 1]
      weights <- weights/sum(weights)
      W=Ks[[1]]*weights[1];
      if(length(Ks)>1) {
        for(i in 2:length(Ks)) W=W+Ks[[i]]*weights[i]
      }
      diag(W)=0;
    }
    Q=t(R) %*% W %*% R;
    Jn=matrix(1/n,ncol=n,nrow=n)
    WJ=matrix(rep(apply(W,1,sum)/n,n),nrow=n);
    WJW=matrix(rep(apply(WJ,2,sum)/n,each=n),nrow=n)
    V=W-WJ-t(WJ)+WJW
    lambda=eigen(V)$values
    pvalue=CompQuadForm::davies(Q,lambda)$Qq
  }
  
  if(scheme==2)
  {
    pvalue<-permWeightedU(Y,Ks,X=X[,-1],B=B,Bmax=1e6)
  }
  return(pvalue)
}






permWeightedU<-function(Y,Ks,X=NULL,B=NULL,Bmax=1e6)
{
  n=length(Y);
  if(!is.null(X)) X=cbind(1,X)
  if(is.null(X)) X=matrix(1,nrow=n,ncol=1);
  In=diag(n);
  Pn=X %*% solve(t(X) %*% X) %*% t(X)
  PY=(In-Pn) %*% Y
  R=rank(PY);
  R=scale(R);
  if(length(Ks)<1 | class(Ks)!="list") stop("Ks must be a list with at least one component")
  for(i in 1:length(Ks)){
    Ks[[i]]=Ks[[i]]/mean(Ks[[i]]^2)
    diag(Ks[[i]])=0;
  }
  Con=matrix(c(0,1),ncol=1);
  for(i in 2:length(Ks))
  {
    weight1=matrix(c(rep(c(0,1),each=nrow(Con))),ncol=1)
    Con=rbind(Con,Con);
    Con=cbind(Con,weight1)
  }
  Con=Con[-1,];
  
  # Get pvalues for each combinations #
  Eigen=list(); porg=NULL; Wlist=list();
  for(i in 1:nrow(Con))
  {
    W=Reduce("+",Map("*",Ks,Con[i,]))
    Q=t(R) %*% W %*% R;
    Jn=matrix(1/n,ncol=n,nrow=n)
    WJ=matrix(rep(apply(W,1,sum)/n,n),nrow=n);
    WJW=matrix(rep(apply(WJ,2,sum)/n,each=n),nrow=n)
    V=W-WJ-t(WJ)+WJW  
    # Vlist[[i]]=V;
    Wlist[[i]]=W;
    Eigen[[i]]=eigen(V)$values;
    porg=c(porg,CompQuadForm::davies(Q,Eigen[[i]])$Qq)
  }
  PvaluesB=matrix(porg,nrow=1);
  if(!is.null(B)) It <- Perm_U(B = B, R, Wlist,Eigen,PvaluesB,Con)
  if (is.null(B)) {
    It <- Perm_U(B = 1000, R, Wlist,Eigen,PvaluesB,Con)
    B1 = round(1/min(abs(It[-1]))) * 50; 
    B1=min(B1,Bmax)
    #print(B1)
    if (B1 > 1000 )  It <- Perm_U(B = B1, R, Wlist,Eigen,PvaluesB,Con)
  }
  return(It)
}
tr<-function(X)
{
  sum(diag(X))
}

Perm_U<-function(B = 200, R, Wlist,Eigen,PvaluesB,Con)
{
  for(i in 1:B)
  {
    if(i %% 100000==0) print(i);
    RR=matrix(sample(1:length(R),length(R)),ncol=1);
    RR=scale(RR);
    Q=sapply(1:length(Wlist),function(jj){t(RR) %*% Wlist[[jj]] %*% RR})
    ptmp=sapply(1:length(Eigen),function(jj)CompQuadForm::davies(Q[jj],Eigen[[jj]])$Qq)
    PvaluesB=rbind(PvaluesB,ptmp)
  }
  pvalueb=apply(abs(PvaluesB),1,min);
  per=mean(pvalueb[-1]<=pvalueb[1]);
  pvalue=c(per,PvaluesB[1,])
  names(pvalue)=c("per",apply(Con,1,paste,collapse="_"));
  return(pvalue);
}