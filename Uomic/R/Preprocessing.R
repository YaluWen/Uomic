#'  A U-statistics method for association test using multi-layer omics data.
#'
#' @param Y The continuous or discrete outcomes for $n$ subjects.
#' @param Omics A list with each element being a n by p matrix. If there is missing, then should be NA for the subject.
#' @param Xcov A matrix of demographic variables (e.g., age and gender). The default is NULL. 
#' @param Kernel The method to calculate the pair-wise sample similarity matrix. It takes the method \code{AM}, \code{IBS}, \code{Lin0}, 
#' \code{Quad1},\code{Minkowski}, \code{Polyk}.Details can be found in \code{varComp}package.
#' @param impute Whether to impute the missing values with the mean. Default is \code{TRUE}.
#' @return A list with the outcome vector \code{Y}, the demographic variable \code{X}, and the pair-wise similarities \code{Ks}
#' @export
#'
#' @examples
#' data(Omics)
#' data(outcome)
#' data(Cov)
#' Input<-ComputeKs(Omics,Y=outcome$all,Xcov=NULL,Kernel=c("Lin0","Lin0","IBS"))
#' Input<-ComputeKs(Omics,Y=outcome$all,Xcov=Cov,Kernel=c("Lin0","Lin0","IBS"))
#'
ComputeKs<-function(Omics,Y,Xcov=NULL,Kernel,c=0,d=2,impute=TRUE)
{
  if(length(Omics)!=length(Kernel)) stop("The number of kernel functions should be the same as the number of omics data")
  n=length(Y);
  if(!is.null(Xcov))
  {
    if(!is.null(dim(Xcov))) {if(nrow(Xcov)!=n) stop("The number of subjects for demographic variables and that of the outcomes are not the same")}
    if(is.null(dim(Xcov))) {if(length(Xcov)!=n) stop("The number of subjects for demographic variables and that of the outcomes are not the same")}
  }
  for(i in 1:length(Omics))
  {
    X=Omics[[i]]
    if(!is.null(dim(X))) {if(nrow(X)!=n) stop("The number of subjects for omics data and that of the outcomes are not the same")}
    if(is.null(dim(X))){if(length(X)!=n) stop("The number of subjects for omics data and that of the outcomes are not the same")}
  }
  ### Use the mean to compute ###
  if(impute) Omics<-lapply(Omics,function(X){
    if(is.null(dim(X))) X[is.na(X)]=mean(X,na.rm = T)
    if(!is.null(dim(X))){est=apply(X,2,mean,na.rm=T);for(i in 1:ncol(X)) X[is.na(X[,i]),i]=est[i]}
    X
  })
  if(!is.null(Xcov))
  {
    if(!is.null(dim(Xcov)))  {exc=(apply(is.na(Xcov),1,sum)!=0) | (is.na(Y)); Y=Y[!exc];Xcov=Xcov[!exc,]}
    if(is.null(dim(Xcov))) {exc=(is.na(Xcov)) | (is.na(Y)); Y=Y[!exc]; Xcov=Xcov[!exc]}
    for(i in 1:length(Omics)) {
      tmp=Omics[[i]];
      if(is.null(dim(tmp))) tmp=matrix(tmp,ncol=1)
      Omics[[i]]=tmp[!exc,]
    }
  }
  exc=is.na(Y);
  if(sum(exc)>0)
  {
    Y=Y[!exc];
    for(i in 1:length(Omics)) {
      tmp=Omics[[i]];
      if(is.null(dim(tmp))) tmp=matrix(tmp,ncol=1)
      Omics[[i]]=tmp[!exc,]
    }  
  }
  Ks=list();
  for(i in 1:length(Omics))
  {
    Ks[[i]]=calKernel(Omics[[i]],Kernel[i],c=c,d=d)
  }

  Input=list();
  Input$Y=Y;
  Input$Ks=Ks;
  Input$X=Xcov;
  Input
}

calKernel<-function(X,method,c=0,d=2)
{
  K=NULL;
  if(tolower(method)=="am") K=varComp::AM(X)
  if(tolower(method)=="ibs") K=varComp::IBS(X)
  if(tolower(method)=="lin0") K=varComp::Lin0(X)
  if(tolower(method)=="quad1") K=varComp::Quad1(X)
  if(tolower(method)=="minkowski") K=varComp::Minkowski(X)
  if(tolower(method)=="polyk") K=varComp::Polyk(X,c=c,d=d)
  K
}

