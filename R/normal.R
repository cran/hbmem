.packageName='hbmem'
.First.lib=function(lib,pkg) library.dynam('hbmem',pkg,lib)

#Truncated Normal
rtnorm=function(N,mu,sigma,a,b)
  {
    y=1:N*0
    .C("rtruncnorm",as.double(y),as.integer(N),as.double(mu),as.double(sigma),as.double(a),as.double(b),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }

#Extract predicted from block
getPred=function(block,sub,item,lag,I,J,R)
  {
    pred=1:R
    tmp=.C("getPred",as.double(pred),as.double(block),as.integer(sub),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),PACKAGE=.packageName)
    return(tmp[[1]])
  }

#Simulate Hierarchical Normal Data
normalSim=function(I=30,J=300,mu=-.5,s2a=.2,s2b=.2,muS2=0,s2aS2=.2,s2bS2=.2)
{
  R=I*J
  alpha=rnorm(I,0,sqrt(s2a))
  beta=rnorm(J,0,sqrt(s2b))
  alphaS2=rnorm(I,0,sqrt(s2aS2))
  betaS2=rnorm(J,0,sqrt(s2bS2))
  subj=rep(0:(I-1),each=J)
  item=rep(0:(J-1),I)
  lag=rep(0,R)
  resp=1:R
    for(r in 1:R)
      {
        mean=mu+alpha[subj[r]+1]+beta[item[r]+1]
        sd=sqrt(exp(muS2+alphaS2[subj[r]+1]+betaS2[item[r]+1]))
        resp[r]=rnorm(1,mean,sd)
      }
return(as.data.frame(cbind(subj,item,lag,resp)))
  
 } 
#BLOCK MEANS WITH ONE VARIANCE
sampleNorm = function(sample,y,subj,item,lag,I,J,R,nsub,nitem,s2mu,s2a,s2b,meta,metb,sigma2,sampLag)
{
b0=c(0,0)
tmp=.C("sampleNormal", as.double(sample),as.double(y),as.integer(subj),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(meta),as.double(metb),as.integer(b0),as.double(sigma2),as.integer(sampLag),PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[16]]
return(list(samp,b0))
}


#BLOCK MEANS WITH BLOCK VARIANCES
sampleNormb = function(sample,y,subj,item,lag,I,J,R,nsub,nitem,s2mu,s2a,s2b,meta,metb,blockSigma2,sampLag)
{
b0=c(0,0)
tmp=.C("sampleNormalb", as.double(sample),as.double(y),as.integer(subj),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(meta),as.double(metb),as.integer(b0),as.double(blockSigma2),as.integer(sampLag),NAOK=TRUE,PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[16]]
return(list(samp,b0))
}

#BLOCK VARIANCES
sampleSig2b = function(sample,y,sub,item,lag,I,J,R,nsub,nitem,s2mu,s2a,s2b,met,blockMean,sampLag)
{
b0=rep(0,I+J+3)
tmp=.C("sampleSigma2b", as.double(sample),as.double(y),as.integer(sub),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(met),as.integer(b0),as.double(blockMean),as.integer(sampLag),PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[15]]
return(list(samp,b0))
}


#ONE SIGMA
sampleSig2=function(sig2,block,y,sub,item,lag,I, J, R,a,b) return(.C("sampleSigma2",as.double(sig2),as.double(block),as.double(y),as.integer(sub),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),as.double(a),as.double(b),PACKAGE=.packageName)[[1]])

