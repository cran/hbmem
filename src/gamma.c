#include <R.h>
#include <Rmath.h>
#include "normal.h"

void rtruncgamma(double *y,int *N,double *shape,double *scale,double *a,double *b)
{
  double lower,upper,u,s;
  int n;
  GetRNGstate();
  for(n=0;n<N[0];n++)
    {
      lower=pgamma(a[n],shape[0],scale[n],1,0);
      if(b[n]==INFINITY) upper=1;
      else upper=pgamma(b[n],shape[0],scale[n],1,0);
      u=runif(lower,upper);
      y[n]=qgamma(u,shape[0],scale[n],1,0);
    }
  PutRNGstate();
}

void logLikeGamma(double *like,int *R,int *NN, int *NS,int *I,int *J,int 
*K,int *dat,int *cond,int *Scond, int *subj, int *item,double *lag,double 
*blockN,double *blockS, double *critVec,double *shape)
{
  int r;
  double scale,crit[7];
  crit[0]=0;
  crit[3]=1;
  crit[K[0]]=INFINITY;

  for(r=0;r<R[0];r++)
    {
      crit[1]=critVec[1+subj[r]*7];
      crit[2]=critVec[2+subj[r]*7];
      crit[4]=critVec[4+subj[r]*7];
      crit[5]=critVec[5+subj[r]*7];
      
      if(Scond[r]==0) 
	scale=exp(blockN[cond[r]]+blockN[NN[0]+subj[r]]+blockN[NN[0]+I[0]+item[r]]+lag[r]*blockN[NN[0]+I[0]+J[0]+2]);
      if(Scond[r]==1)
	scale=exp(blockS[cond[r]]+blockS[NS[0]+subj[r]]+blockS[NS[0]+I[0]+item[r]]+lag[r]*blockS[NS[0]+I[0]+J[0]+2]);
      like[r]=log(pgamma(crit[dat[r]+1],shape[0],scale,1,0)-pgamma(crit[dat[r]],shape[0],scale,1,0));
    }
}

void sampleGamma(double *block,double *dat,int *cond,int *subj, int *item,double *lag,int *N,int *I,int *J,int *R,int *ncond,int *nsub, int *nitem,double *sig2Mu,double *sig2A,double *sig2B,double *met,int *b0,double *shape,int *sampLag)
{
  int n,i,j,r;
  double accept;

  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=block[n];
  for(i=0;i<I[0];i++) alpha[i]=block[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=block[N[0]+I[0]+j];
  double sig2Alpha=block[N[0]+I[0]+J[0]];
  double sig2Beta=block[N[0]+I[0]+J[0]+1];
  double theta=block[N[0]+I[0]+J[0]+2];
  GetRNGstate(); 
  

  //Sample Mu
  double muProp[N[0]];
  for(n=0;n<N[0];n++) muProp[n]=mu[n]+rnorm(0,met[n]);    
  double nLike[N[0]],nLikeProp[N[0]];
  for(n=0;n<N[0];n++)
    {
      nLike[n]=0;
      nLikeProp[n]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      nLike[cond[r]]-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      nLikeProp[cond[r]]-=dat[r]/exp(muProp[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
    }
  
  for(n=0;n<N[0];n++)
    {	  
      nLike[n]-=shape[0]*ncond[n]*mu[n] + (mu[n]*mu[n])/(2*sig2Mu[0]);
      nLikeProp[n]-=shape[0]*ncond[n]*muProp[n] + (muProp[n]*muProp[n])/(2*sig2Mu[0]);
      if(rbinom(1,fmin2(1,exp(nLikeProp[n]-nLike[n]))))
	{
	  mu[n]=muProp[n];
	  b0[n]++;
	}
    }

  //Sample Alpha
  double alphaProp[I[0]];
  if(I[0]>1){
  for(i=0;i<I[0];i++) 
    alphaProp[i]=alpha[i]+rnorm(0,met[N[0]+i]);    
  double aLike[I[0]],aLikeProp[I[0]];
  for(i=0;i<I[0];i++)
    {
      aLike[i]=0;
      aLikeProp[i]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      aLike[subj[r]]-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      aLikeProp[subj[r]]-=dat[r]/exp(mu[cond[r]]+alphaProp[subj[r]]+beta[item[r]]+theta*lag[r]);
    }
  
  for(i=0;i<I[0];i++)
    {	  
      aLike[i]-=shape[0]*nsub[i]*alpha[i] + (alpha[i]*alpha[i])/(2*sig2Alpha);
      aLikeProp[i]-=shape[0]*nsub[i]*alphaProp[i] + (alphaProp[i]*alphaProp[i])/(2*sig2Alpha);
      if(rbinom(1,fmin2(1,exp(aLikeProp[i]-aLike[i]))))
	{
	  alpha[i]=alphaProp[i];
	  b0[N[0]+i]++;
	}
    }
  }

   //Sample Beta
  double betaProp[J[0]];
  if(J[0]>1){
  for(j=0;j<J[0];j++) 
    betaProp[j]=beta[j]+rnorm(0,met[N[0]+I[0]+j]);    
  double bLike[J[0]],bLikeProp[J[0]];
  for(j=0;j<J[0];j++)
    {
      bLike[j]=0;
      bLikeProp[j]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      bLike[item[r]]-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      bLikeProp[item[r]]-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+betaProp[item[r]]+theta*lag[r]);
    }
  
  for(j=0;j<J[0];j++)
    {	  
      bLike[j]-=shape[0]*nitem[j]*beta[j] + (beta[j]*beta[j])/(2*sig2Beta);
      bLikeProp[j]-=shape[0]*nitem[j]*betaProp[j] + (betaProp[j]*betaProp[j])/(2*sig2Beta);
      if(rbinom(1,fmin2(1,exp(bLikeProp[j]-bLike[j]))))
	{
	  beta[j]=betaProp[j];
	  b0[N[0]+I[0]+j]++;
	}
    }
  }

  //Sample s2alpha
  double rate=0;
  for(i=0;i<I[0];i++) rate+=alpha[i]*alpha[i];
  rate=rate/2 + sig2B[0];
  sig2Alpha=1.0/rgamma(I[0]/2+sig2A[0],1/rate);
  //Sample s2beta
  rate=0;
  for(j=0;j<J[0];j++) rate+=beta[j]*beta[j];
  rate=rate/2 + sig2B[0];
  sig2Beta=1/rgamma(J[0]/2+sig2A[0],1/rate);
 
  //Sample theta
  if(sampLag[0]){
  double thetaProp;
    thetaProp=theta+rnorm(0,met[N[0]+I[0]+J[0]+2]);    
    double thetaLike=0,thetaLikeProp=0;
    for(r=0;r<R[0];r++)
    {  
      thetaLike-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      thetaLikeProp-=dat[r]/exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+thetaProp*lag[r]);
    }
    thetaLike-=(theta*theta)/(2*sig2Mu[0]);
    thetaLikeProp-= (thetaProp*thetaProp)/(2*sig2Mu[0]);
    //note: element in shape goes away since sum(lags)=0
    if(rbinom(1,fmin2(1,exp(thetaLikeProp-thetaLike))))
	{
	  theta=thetaProp;
	  b0[N[0]+I[0]+J[0]+2]++;
	} 
  }

  double p,shift;
  //DECORR Mu and Alpha
  if(I[0]>1){
    shift=rnorm(0,met[N[0]+I[0]+J[0]]);
    for(n=0;n<N[0];n++) muProp[n]=mu[n]+shift;
    for(i=0;i<I[0];i++) alphaProp[i]=alpha[i]-shift;
    p=(sumsqr(alpha,I[0])-sumsqr(alphaProp,I[0]))/(2*sig2Alpha) + (sumsqr(mu,N[0])-sumsqr(muProp,N[0]))/(2*sig2Mu[0]);
    if (rbinom(1,fmin2(1,exp(p))))
      {
	for(n=0;n<N[0];n++) mu[n]=muProp[n];
	for(i=0;i<I[0];i++) alpha[i]=alphaProp[i];
	b0[N[0]+I[0]+J[0]]+=1;
      }
  }
  
  //Decorrelate Alpha and Beta
  if((I[0]>1) && (J[0]>1)){
    shift=rnorm(0,met[N[0]+I[0]+J[0]+1]);
    for(i=0;i<I[0];i++) alphaProp[i]=alpha[i]+shift;
    for(j=0;j<J[0];j++) betaProp[j]=beta[j]-shift;
    
    p=(sumsqr(alpha,I[0])-sumsqr(alphaProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(betaProp,J[0]))/(2*sig2Beta);
    
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=alphaProp[i];
      for(j=0;j<J[0];j++) beta[j]=betaProp[j];
      b0[N[0]+I[0]+J[0]+1]+=1;
    }
  }

  //write block
  for(n=0;n<N[0];n++) block[n]=mu[n];
  for(i=0;i<I[0];i++) block[i+N[0]]=alpha[i];
  for(j=0;j<J[0];j++) block[N[0]+I[0]+j]=beta[j];
  block[N[0]+I[0]+J[0]]=sig2Alpha;
  block[N[0]+I[0]+J[0]+1]=sig2Beta;
  block[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();

}


//Likelihood for likelihood model
void logLikeGammaLike(double *like,int *R,int *NN, int *NS,int *I,int *J,int 
*K,int *dat,int *cond,int *Scond, int *subj, int *item,double *lag,double 
*blockN,double *blockS, double *critVec,double *shape)
{
  int r;
  double scale,scaleN=1,scaleS,crit[7];
  crit[0]=0;
  crit[K[0]]=INFINITY;

  for(r=0;r<R[0];r++)
    {
      scaleS=1+exp(blockS[cond[r]]+blockS[NS[0]+subj[r]]+blockS[NS[0]+I[0]+item[r]]+lag[r]*blockS[NS[0]+I[0]+J[0]+2]);

      crit[1]=critVec[1+subj[r]*7];
      crit[2]=critVec[2+subj[r]*7];
      crit[3]=critVec[3+subj[r]*7];
      crit[4]=critVec[4+subj[r]*7];
      crit[5]=critVec[5+subj[r]*7];
      
      if(Scond[r]==0) scale=1-1/scaleS;
      else scale=scaleS-1;

   like[r]=log(pgamma(crit[dat[r]+1],shape[0],scale,1,0)-pgamma(crit[dat[r]],shape[0],scale,1,0));
    }
}

void samplePosGamma(double *block,double *dat,int *cond,int *subj, int *item,double *lag,int *N,int *I,int *J,int *R,int *ncond,int *nsub, int *nitem,double *sig2Mu,double *sig2A,double *sig2B,double *met,int *b0,double *shape,int *sampLag)
{
  int n,i,j,r;
  double accept;

  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=block[n];
  for(i=0;i<I[0];i++) alpha[i]=block[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=block[N[0]+I[0]+j];
  double sig2Alpha=block[N[0]+I[0]+J[0]];
  double sig2Beta=block[N[0]+I[0]+J[0]+1];
  double theta=block[N[0]+I[0]+J[0]+2];
  GetRNGstate(); 
  

  //Sample Mu
  double muProp[N[0]],scale;
  for(n=0;n<N[0];n++) muProp[n]=mu[n]+rnorm(0,met[n]);    
  double nLike[N[0]],nLikeProp[N[0]];
  for(n=0;n<N[0];n++)
    {
      nLike[n]=0;
      nLikeProp[n]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      nLike[cond[r]]-=dat[r]/scale + shape[0]*log(scale);
      scale=1+exp(muProp[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      nLikeProp[cond[r]]-=dat[r]/scale + shape[0]*log(scale);
    }
  
  for(n=0;n<N[0];n++)
    {	  
      nLike[n]-=(mu[n]*mu[n])/(2*sig2Mu[0]);
      nLikeProp[n]-=(muProp[n]*muProp[n])/(2*sig2Mu[0]);
      if(rbinom(1,fmin2(1,exp(nLikeProp[n]-nLike[n]))))
	{
	  mu[n]=muProp[n];
	  b0[n]++;
	}
    }

  //Sample Alpha
  double alphaProp[I[0]];
  if(I[0]>1){
  for(i=0;i<I[0];i++) 
    alphaProp[i]=alpha[i]+rnorm(0,met[N[0]+i]);    
  double aLike[I[0]],aLikeProp[I[0]];
  for(i=0;i<I[0];i++)
    {
      aLike[i]=0;
      aLikeProp[i]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      aLike[subj[r]]-=dat[r]/scale + shape[0]*log(scale);
      scale=1+exp(mu[cond[r]]+alphaProp[subj[r]]+beta[item[r]]+theta*lag[r]);
      aLikeProp[subj[r]]-=dat[r]/scale + shape[0]*log(scale);
    }
  
  for(i=0;i<I[0];i++)
    {	  
      aLike[i]-=(alpha[i]*alpha[i])/(2*sig2Alpha);
      aLikeProp[i]-=(alphaProp[i]*alphaProp[i])/(2*sig2Alpha);
      if(rbinom(1,fmin2(1,exp(aLikeProp[i]-aLike[i]))))
	{
	  alpha[i]=alphaProp[i];
	  b0[N[0]+i]++;
	}
    }
  }

   //Sample Beta
  double betaProp[J[0]];
  if(J[0]>1){
  for(j=0;j<J[0];j++) 
    betaProp[j]=beta[j]+rnorm(0,met[N[0]+I[0]+j]);    
  double bLike[J[0]],bLikeProp[J[0]];
  for(j=0;j<J[0];j++)
    {
      bLike[j]=0;
      bLikeProp[j]=0;
    }
  for(r=0;r<R[0];r++)
    {  
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      bLike[item[r]]-=dat[r]/scale + shape[0]*log(scale);
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+betaProp[item[r]]+theta*lag[r]);
      bLikeProp[item[r]]-=dat[r]/scale + shape[0]*log(scale); 
    }
  
  for(j=0;j<J[0];j++)
    {	  
      bLike[j]-=(beta[j]*beta[j])/(2*sig2Beta);
      bLikeProp[j]-=(betaProp[j]*betaProp[j])/(2*sig2Beta);
      if(rbinom(1,fmin2(1,exp(bLikeProp[j]-bLike[j]))))
	{
	  beta[j]=betaProp[j];
	  b0[N[0]+I[0]+j]++;
	}
    }
  }

  //Sample s2alpha
  double rate=0;
  for(i=0;i<I[0];i++) rate+=alpha[i]*alpha[i];
  rate=rate/2 + sig2B[0];
  sig2Alpha=1.0/rgamma(I[0]/2+sig2A[0],1/rate);
  //Sample s2beta
  rate=0;
  for(j=0;j<J[0];j++) rate+=beta[j]*beta[j];
  rate=rate/2 + sig2B[0];
  sig2Beta=1/rgamma(J[0]/2+sig2A[0],1/rate);
 
  //Sample theta
  if(sampLag[0]){
  double thetaProp;
    thetaProp=theta+rnorm(0,met[N[0]+I[0]+J[0]+2]);    
    double thetaLike=0,thetaLikeProp=0;
    for(r=0;r<R[0];r++)
    {  
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+theta*lag[r]);
      thetaLike-=dat[r]/scale + shape[0]*log(scale);
      scale=1+exp(mu[cond[r]]+alpha[subj[r]]+beta[item[r]]+thetaProp*lag[r]);
	thetaLikeProp-=dat[r]/scale + shape[0]*log(scale);
    }
    thetaLike-=(theta*theta)/(2*sig2Mu[0]);
    thetaLikeProp-= (thetaProp*thetaProp)/(2*sig2Mu[0]);
    if(rbinom(1,fmin2(1,exp(thetaLikeProp-thetaLike))))
	{
	  theta=thetaProp;
	  b0[N[0]+I[0]+J[0]+2]++;
	} 
  }

  double p,shift;
  //DECORR Mu and Alpha
  if(I[0]>1){
    shift=rnorm(0,met[N[0]+I[0]+J[0]]);
    for(n=0;n<N[0];n++) muProp[n]=mu[n]+shift;
    for(i=0;i<I[0];i++) alphaProp[i]=alpha[i]-shift;
    p=(sumsqr(alpha,I[0])-sumsqr(alphaProp,I[0]))/(2*sig2Alpha) + (sumsqr(mu,N[0])-sumsqr(muProp,N[0]))/(2*sig2Mu[0]);
    if (rbinom(1,fmin2(1,exp(p))))
      {
	for(n=0;n<N[0];n++) mu[n]=muProp[n];
	for(i=0;i<I[0];i++) alpha[i]=alphaProp[i];
	b0[N[0]+I[0]+J[0]]+=1;
      }
  }
  
  //Decorrelate Alpha and Beta
  if((I[0]>1) && (J[0]>1)){
    shift=rnorm(0,met[N[0]+I[0]+J[0]+1]);
    for(i=0;i<I[0];i++) alphaProp[i]=alpha[i]+shift;
    for(j=0;j<J[0];j++) betaProp[j]=beta[j]-shift;
    
    p=(sumsqr(alpha,I[0])-sumsqr(alphaProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(betaProp,J[0]))/(2*sig2Beta);
    
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=alphaProp[i];
      for(j=0;j<J[0];j++) beta[j]=betaProp[j];
      b0[N[0]+I[0]+J[0]+1]+=1;
    }
  }

  //write block
  for(n=0;n<N[0];n++) block[n]=mu[n];
  for(i=0;i<I[0];i++) block[i+N[0]]=alpha[i];
  for(j=0;j<J[0];j++) block[N[0]+I[0]+j]=beta[j];
  block[N[0]+I[0]+J[0]]=sig2Alpha;
  block[N[0]+I[0]+J[0]+1]=sig2Beta;
  block[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();

}

