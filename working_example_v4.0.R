Rprof("Rprof.out",memory.profiling=TRUE)

library(sommer)
library(BGLR)

################################################################################
################################################################################
#Functionsforthecalculationofthegenomicvarianceindifferentset-ups
GenVarCur<-function(GRM,varg,eBLUP,eBLUPcov){
n<-dim(GRM)[[1]]
  V<-varg*sum(diag(GRM))/(n-1)#formula(8)
  W<-V+sum(eBLUP^2)/(n-1)-sum(diag(eBLUPcov))/(n-1)
  return(c(V,W))
}
GenVarBase_Ped<-function(PedMat,GRM,varg,eBLUP,eBLUPcov){
  n<-dim(GRM)[[1]]
  P<-diag(1,n,n)-matrix(1,n,n)/n #Equation(2)
 # Q<-solve(chol(P))
  
   E<-eigen(PedMat)
   Q<-E$vectors
   D<-E$values
  Ped.inv.halve<-Q%*%diag(1/sqrt(D),nrow(Q),nrow(Q))%*%t(Q)
  B_b<-Ped.inv.halve%*%P%*%Ped.inv.halve #equation(16)
  V<-varg*sum(diag(B_b%*%GRM))/(n-1)#formula(17)
  W<-V+t(eBLUP)%*%B_b%*%eBLUP/
    (n-1)-sum(diag(B_b%*%eBLUPcov))/(n-1)#formula(18)
  return(c(V,W))
}
GenVarBase_GRM<-function(GRM,phenos,mu,varg,eBLUP,eBLUPcov,vare){
  n<-dim(GRM)[[1]]
  P<-diag(1,n,n)-matrix(1,n,n)/n#Equation(2)
  E_G<-eigen(GRM)
  Q_G<-E_G$vectors
  D_G<-E_G$values
  D_G[D_G<0]<-0#LastEVmightbe-0
  G_halve<-Q_G%*%diag(sqrt(D_G),nrow(Q_G),nrow(Q_G))%*%t(Q_G)
  Cov.y.minus<-chol2inv(GRM*varg+diag(vare,n,n))
  eBLUP.s<-varg*G_halve%*%Cov.y.minus%*%(phenos-mu)#formula(20)
  eBLUPcov.s<-varg^2*G_halve%*%Cov.y.minus%*%G_halve-(varg^2/sum(Cov.y.minus))*G_halve%*%Cov.y.minus%*%matrix(1,n,n)%*%Cov.y.minus%*%G_halve#formula(21)
  V<-varg#formula(19)
  W<-V+t(eBLUP.s)%*%P%*%eBLUP.s/(n-1)-sum(diag(P%*%eBLUPcov.s))/(n-1)#formula(22)
  return(c(V,W))
}
data(mice)




y<-mice.pheno$Obesity.BodyLength#fixphenotype
y.scaled<-y/sd(y)#yscaledtovarianceof1,seeequation(5)
n<-length(y)
ones<-matrix(1,nrow=n,ncol=1)
P<-diag(1,n,n)-matrix(1,n,n)/n#Equation(2)
X<-mice.X#fixmarker-genotypes
p<-dim(X)[2]
c<-sum(colMeans(X)*(1-colMeans(X)/2))#vanRadenc,seeequation(3)
G<-P%*%tcrossprod(X)%*%P/c#equation(1)
rm(X)
############################################################################
#FitgBLUPmodel
lm.equi<-mmer(Y~X,random=~G,data=data.frame(Y=y.scaled,X=ones),iters=50,tolpar=1e-7,tolparinv=1e-9,date.warning=FALSE)
var.g.hat<-as.numeric(lm.equi$sigma$G)#variancecomponentg
var.e.hat<-as.numeric(lm.equi$sigma$units)#variancecomp.eps
mu.hat<-rep(lm.equi$Beta$Estimate,n)#EstimateofIntercept
g.hat<-as.vector(lm.equi$U$G[[1]])#eBLUP,equation(6)
cov.g.hat<-as.matrix(lm.equi$VarU$G$Y)#CovofeBLUP,equation(7)
rm(lm.equi)
#Genomicvariancesinthecurrentpopulation
gv_cur<-GenVarCur(GRM=G,varg=var.g.hat,eBLUP=g.hat,eBLUPcov=cov.g.hat)
V.hat<-gv_cur[1]
h2.V<-V.hat#formula(10)
sum.V<-V.hat+var.e.hat#formula(12)
h2.V.sum<-V.hat/sum.V#formula(14)
W.hat<-gv_cur[2]
h2.W<-W.hat#(11)
sum.W<-W.hat+var.e.hat#formula(13)
h2.W.sum<-W.hat/sum.W#formula(15)
#Genomicvariancesinthebasepopulation(viaPedigree)
Ped<-mice.A
if(!any(is.na(Ped))){
gv_base_ped<-GenVarBase_Ped(PedMat=Ped,GRM=G,varg=var.g.hat,eBLUP=g.hat,eBLUPcov=cov.g.hat)
V.hat.star<-gv_base_ped[1]
W.hat.star<-gv_base_ped[2]
}else{
V.hat.star<-W.hat.star<-NA
}
#Genomicvariancesinthebasepopulation(viaGRM)
gv_base_grm<-GenVarBase_GRM(GRM=G,phenos=y.scaled,mu=mu.hat,
varg=var.g.hat,eBLUP=g.hat,eBLUPcov=cov.g.hat,
vare=var.e.hat)
V.hat.star.s<-gv_base_grm[1]
W.hat.star.s<-gv_base_grm[2]
##############################################################################
result<-c(var.e.hat,V.hat,h2.V,sum.V,h2.V.sum,W.hat,h2.W,
sum.W,h2.W.sum,V.hat.star,W.hat.star,V.hat.star.s,
W.hat.star.s)
rm(var.e.hat,V.hat,h2.V,sum.V,h2.V.sum,W.hat,h2.W,sum.W,h2.W.sum,
V.hat.star,W.hat.star,V.hat.star.s,W.hat.star.s,y,y.scaled,var.g.hat,
g.hat,cov.g.hat)

Rprof(NULL)
