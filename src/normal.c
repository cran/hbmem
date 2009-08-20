#include "normal.h"

//Calculate the sum of squares of vector
double sumsqr(double *invec, int size)
{
  int i;
  double sum=0;
  
  for (i=0;i<size;i++)
    sum+=(invec[i]*invec[i]);
  return (sum);
}

double sum(double *invec, int size)
{
  int i;
  double sum=0;
  
  for (i=0;i<size;i++)
    sum+=(invec[i]);
  return (sum);
}

//Sample Random Truncated Normal.  Needs better sampling
//in low density tails
void rtruncnorm(double *y, int *N,double *mu, double *sigma, double *a,double *b)
{
  long double lower,upper,sample,unif;
  int n;
  GetRNGstate(); 
  for(n=0;n<N[0];n++)
    {
      lower=pnorm(a[n],mu[n],sigma[n],1,0);
      upper=pnorm(b[n],mu[n],sigma[n],1,0);
      unif=runif(lower, upper);
      y[n]=qnorm(unif,mu[n],sigma[n],1,0);

               
      //sloppy fix
      if(y[n]==INFINITY | y[n]==-INFINITY)
	{
	  printf("Warning! RTNORM = INF; %lf %lf %lf %lf\n",mu[n],sigma[n],a[n],b[n]);
	  y[n]=runif(a[n],b[n]);
	}       
      
    }
  PutRNGstate();
}

//Take a block of paramters and spit back the predicted mean for each trial

void getPred(double *pred,double *block,int *cond,int *subj, int *item,double *cov, int *N, int *I, int *J, int *R)
{
  int r;
  for(r=0;r<R[0];r++) pred[r]=(block[cond[r]]+block[N[0]+subj[r]]+block[N[0]+I[0]+item[r]]+cov[r]*block[N[0]+I[0]+J[0]+2]);
}

//Sample block means with ONE sigma
void sampleNormal(double *block,double *dat,int *cond, int *subj, int *item,double *cov,int *N,int *I,int *J,int *R,int *Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *sig2,int *sampCov,int *Hier)
{
  int r,n,i,j;
  double sumY=0,sumYN[N[0]],sumYI[I[0]],sumYJ[J[0]];
  double a,postB;
  double p,shift,muProp[N[0]],aProp[I[0]],bProp[J[0]];

  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;

  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=block[n];
  for(i=0;i<I[0];i++) alpha[i]=block[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=block[N[0]+I[0]+j];
  double sig2Alpha=block[N[0]+I[0]+J[0]];
  double sig2Beta=block[N[0]+I[0]+J[0]+1];
  double theta=block[N[0]+I[0]+J[0]+2];
  GetRNGstate(); 

  //sample mu
  for(n=0;n<N[0];n++) sumYN[n]=0;
  for(r=0;r<R[0];r++) sumYN[cond[r]]+=dat[r]-alpha[subj[r]]-beta[item[r]]-theta*cov[r];
  
  for(n=0;n<N[0];n++) 
    {
      a=1/(Ncond[n]/sig2[0]+1/sig2Mu[0]);
      mu[n]=rnorm(sumYN[n]/sig2[0]*a,sqrt(a));
    }
 
  //sample alpha
  if(I[0]>1){
  for(i=0;i<I[0];i++) sumYI[i]=0;
  for(r=0;r<R[0];r++) sumYI[subj[r]]+=dat[r]-mu[cond[r]]-beta[item[r]]-theta*cov[r];
  
  for(i=0;i<I[0];i++) 
    {
      a=1/(Nsub[i]/sig2[0]+1/sig2Alpha);
      alpha[i]=rnorm(sumYI[i]/sig2[0]*a,sqrt(a));
    }
  if(Hier[0]==1)
    {
      postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
      sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
    }
  }
  if(J[0]>1){
  //sample beta
  for(j=0;j<J[0];j++) sumYJ[j]=0;
  for(r=0;r<R[0];r++) sumYJ[item[r]]+=dat[r]-mu[cond[r]]-alpha[subj[r]]-theta*cov[r];
  
  for(j=0;j<J[0];j++)
    {
      a=1/(Nitem[j]/sig2[0]+1/sig2Beta);
      beta[j]=rnorm(sumYJ[j]/sig2[0]*a,sqrt(a));
    }
  if(Hier[0]==1)
    {
      postB=sumsqr(beta,J[0])/2.0+sig2B[0];
      sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
    }
  }
  //sample theta
  if(sampCov[0])
    {
      sumY=0;
      for(r=0;r<R[0];r++) sumY+=cov[r]*(dat[r]-mu[cond[r]]-alpha[subj[r]]-beta[item[r]]);
      a=1/(sumsqr(cov,R[0])/sig2[0]+1/sig2Mu[0]);
      theta=rnorm(sumY/sig2[0]*a,sqrt(a)); 
    }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
  shift=rnorm(0,sqrt(s2decor1[0]));
  for(n=0;n<N[0];n++) muProp[n]=mu[n]+shift;
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(mu,N[0])-sumsqr(muProp,N[0]))/(2*sig2Mu[0]);
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(n=0;n<N[0];n++) mu[n]=muProp[n];
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      b0[0]+=1;
    }
  }

  //Decorrelate Alpha and Beta
  if(I[0]>1 & J[0]>1){
  shift=rnorm(0,sqrt(s2decor2[0]));
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]+shift;
  for(j=0;j<J[0];j++) bProp[j]=beta[j]-shift;
  
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(bProp,J[0]))/(2*sig2Beta);
  
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      for(j=0;j<J[0];j++) beta[j]=bProp[j];
      b0[1]+=1;
    }
  }

  for(n=0;n<N[0];n++) block[n]=mu[n];
  for(i=0;i<I[0];i++) block[i+N[0]]=alpha[i];
  for(j=0;j<J[0];j++) block[N[0]+I[0]+j]=beta[j];
  block[N[0]+I[0]+J[0]]=sig2Alpha;
  block[N[0]+I[0]+J[0]+1]=sig2Beta;
  block[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();
}

//sample one sigma
void sampleSigma2(double *sigma2,double *block, double *dat, int *cond, int *subj,int *item,double *cov, int *N,int *ncond,int *I, int *J,double *a, double *b)
{
double tmp,postA,postB,resid[N[0]];
int n,r,R=0;
for(n=0;n<N[0];n++)
{
 resid[n]=0;
 R+=ncond[n];
}

for(r=0;r<R;r++)
{
tmp=dat[r]-(block[cond[r]]+block[N[0]+subj[r]]+block[N[0]+I[0]+item[r]]+cov[r]*block[N[0]+I[0]+J[0]+2]); 
resid[cond[r]]+=tmp*tmp;
}
GetRNGstate();
for(n=0;n<N[0];n++){
postA=a[0] + ncond[n]/2.0;
postB=resid[n]/2.0 + b[0];
sigma2[n]=1.0/rgamma(postA,1.0/postB);
}
PutRNGstate();
}


void sampleNormalb(double *block,double *dat,int *cond,int *subj, int *item,double *cov,int *N,int *I,int *J,int *R,int *Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *blocks2,int *sampCov,int *Hier)
{
  int n,r,i,j,v;
  double a,postB;
  double p,shift,muProp[N[0]],aProp[I[0]],bProp[J[0]];

  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;

  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=block[n];
  for(i=0;i<I[0];i++) alpha[i]=block[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=block[N[0]+I[0]+j];
  double sig2Alpha=block[N[0]+I[0]+J[0]];
  double sig2Beta=block[N[0]+I[0]+J[0]+1];
  double theta=block[N[0]+I[0]+J[0]+2];

  double sigma2[R[0]];
  for(r=0;r<R[0];r++) 
    {
      sigma2[r]=exp(blocks2[cond[r]]+blocks2[N[0]+subj[r]]+blocks2[N[0]+I[0]+item[r]]+cov[r]*blocks2[N[0]+I[0]+J[0]+2]);     
    }

  GetRNGstate(); 

  //sample Mu
  double sumYN[N[0]],sumInvS2N[N[0]];
  for(n=0;n<N[0];n++) 
    {
      sumYN[n]=0;
      sumInvS2N[n]=0;
    }

  for(r=0;r<R[0];r++) 
    {
      sumYN[cond[r]]+=(dat[r]-alpha[subj[r]]-beta[item[r]]-theta*cov[r])/sigma2[r];  
      sumInvS2N[cond[r]]+=1/sigma2[r];
    }
  for(n=0;n<N[0];n++) 
    {
      a=1/(sumInvS2N[n] + 1/sig2Mu[0]);
      mu[n]=rnorm(a*sumYN[n],sqrt(a));
    }
  
  //sample alpha
  if(I[0]>1){
  double sumYI[I[0]],sumInvS2I[I[0]];
  for(i=0;i<I[0];i++)
    {
      sumYI[i]=0;
      sumInvS2I[i]=0;
    }
  for(r=0;r<R[0];r++) 
    {
      sumYI[subj[r]]+=(dat[r]-mu[cond[r]]-beta[item[r]]-theta*cov[r])/sigma2[r];
      sumInvS2I[subj[r]]+=1/sigma2[r];
    }  
  for(i=0;i<I[0];i++) 
    {
      a=1/(sumInvS2I[i] + 1/sig2Alpha);
      alpha[i]=rnorm(a*sumYI[i],sqrt(a));
    }
  if(Hier[0]==1)
    {
      postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
      sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
    }
  }
  //sample beta
  if(J[0]>1){
  double sumYJ[J[0]],sumInvS2J[J[0]];
  for(j=0;j<J[0];j++) 
    {
      sumYJ[j]=0;
      sumInvS2J[j]=0;
    }

  for(r=0;r<R[0];r++) 
    {
      sumYJ[item[r]]+=(dat[r]-mu[cond[r]]-alpha[subj[r]]-theta*cov[r])/sigma2[r];
      sumInvS2J[item[r]]+=1/sigma2[r];
    }

  for(j=0;j<J[0];j++)
    {
      a=1/(sumInvS2J[j] + 1/sig2Beta);
      beta[j]=rnorm(a*sumYJ[j],sqrt(a));
    }
  if(Hier[0]==1)
    {
      postB=sumsqr(beta,J[0])/2.0+sig2B[0];
      sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
    }
  }
  //sample theta
  if(sampCov[0])
    {
      double sumY=0,ycov=0;
      a=0;
      for(r=0;r<R[0];r++) 
	{
	  sumY+=cov[r]*(dat[r]-mu[cond[r]]-alpha[subj[r]]-beta[item[r]])/sigma2[r];
	  a+=(cov[r]*cov[r])/sigma2[r];
	}

      a=1/(a+1/sig2Mu[0]);
      theta=rnorm(a*sumY,sqrt(a)); 
    }

  //Decorrelate Mu and Alpha
  if(I[0]>1){
  shift=rnorm(0,sqrt(s2decor1[0]));
  for(n=0;n<N[0];n++) muProp[n]=mu[n]+shift;
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(mu,N[0])-sumsqr(muProp,N[0]))/(2*sig2Mu[0]);
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(n=0;n<N[0];n++) mu[n]=muProp[n];
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      b0[0]+=1;
    }
  }

  //Decorrelate Alpha and Beta
  if(I[0]>1 & J[0]>1){
  shift=rnorm(0,sqrt(s2decor2[0]));
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]+shift;
  for(j=0;j<J[0];j++) bProp[j]=beta[j]-shift;
  
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(bProp,J[0]))/(2*sig2Beta);
  
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      for(j=0;j<J[0];j++) beta[j]=bProp[j];
      b0[1]+=1;
    }
  }

  //Write Output
  for(n=0;n<N[0];n++) block[n]=mu[n];
  for(i=0;i<I[0];i++) block[N[0]+i]=alpha[i];
  for(j=0;j<J[0];j++) block[N[0]+I[0]+j]=beta[j];
  block[N[0]+I[0]+J[0]]=sig2Alpha;
  block[N[0]+I[0]+J[0]+1]=sig2Beta;
  block[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();
}



////////////////////////////////////////////////
void sampleSigma2b(double *blockS2,double *dat,int *cond,int *subj, int *item,double *cov,int *N,int *I,int *J,int *R,int *Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *met,int *b0,double *blockMu,int *sampCov,int *Hier)
{
  int n,r,i,j;
  double likePropN[N[0]],likeOldN[N[0]],likePropI[I[0]],likeOldI[I[0]],likePropJ[J[0]],likeOldJ[J[0]];
  double a,S=0,postB;
  double p,shift,muProp[N[0]],aProp[I[0]],bProp[J[0]];
  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;
  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=blockS2[n];
  for(i=0;i<I[0];i++) alpha[i]=blockS2[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=blockS2[N[0]+I[0]+j];
  double sig2Alpha=blockS2[N[0]+I[0]+J[0]];
  double sig2Beta=blockS2[N[0]+I[0]+J[0]+1];
  double theta=blockS2[N[0]+I[0]+J[0]+2];
  double mean[R[0]],sigma2[R[0]];
 
  for(r=0;r<R[0];r++) mean[r]=blockMu[cond[r]]+blockMu[N[0]+subj[r]]+blockMu[N[0]+I[0]+item[r]]+cov[r]*blockMu[N[0]+I[0]+J[0]+2];
      
  GetRNGstate(); 

  //sample Mu
  for(n=0;n<N[0];n++)
    {
      likeOldN[n]=0;
      likePropN[n]=0;
      muProp[n]=mu[n]+rnorm(0,met[n]);
    }

  for(r=0;r<R[0];r++) 
    {
      likeOldN[cond[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + beta[item[r]]+theta*cov[r]);
      likePropN[cond[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(muProp[cond[r]] + alpha[subj[r]] + beta[item[r]]+theta*cov[r]);
    }

  for(n=0;n<N[0];n++)
    {
      likeOldN[n]=-.5*(Ncond[n]*mu[n] + (mu[n]*mu[n])/sig2Mu[0] + likeOldN[n]);
      likePropN[n]=-.5*(Ncond[n]*muProp[n] + (muProp[n]*muProp[n])/sig2Mu[0] + likePropN[n]);
      if (rbinom(1,fmin2(1,exp(likePropN[n]-likeOldN[n]))))
	{
	  mu[n]=muProp[n];
	  b0[n]+=1;
	}
    }
      
  //sample alpha
  if(I[0]>1){
    for(i=0;i<I[0];i++) 
      {
	likeOldI[i]=0;	
	likePropI[i]=0;
	aProp[i]=alpha[i]+rnorm(0,met[N[0]+i]);
      }
    for(r=0;r<R[0];r++) 
      {
	likeOldI[subj[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + beta[item[r]]+theta*cov[r]);
	likePropI[subj[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + aProp[subj[r]] + beta[item[r]]+theta*cov[r]);
      }
    for(i=0;i<I[0];i++) 
      {
	likeOldI[i]=-.5*(Nsub[i]*alpha[i] + (alpha[i]*alpha[i])/sig2Alpha + likeOldI[i]);
	likePropI[i]=-.5*(Nsub[i]*aProp[i] + (aProp[i]*aProp[i])/sig2Alpha + likePropI[i]);
	if (rbinom(1,fmin2(1,exp(likePropI[i]-likeOldI[i]))))
	  {
	    alpha[i]=aProp[i];
	    b0[N[0]+i]+=1;
	  }
      }
    if(Hier[0]==1)
      {
	postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
	sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
      }
  }
    
  //Sample Beta
  if(J[0]>1){
    for(j=0;j<J[0];j++) 
      {
	likeOldJ[j]=0;	
	likePropJ[j]=0;
	bProp[j]=beta[j]+rnorm(0,met[N[0]+I[0]+j]);
      }
    for(r=0;r<R[0];r++) 
      {
	likeOldJ[item[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + beta[item[r]]+theta*cov[r]);
      likePropJ[item[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + bProp[item[r]]+theta*cov[r]);
    }
    for(j=0;j<J[0];j++) 
      {
	likeOldJ[j]=-.5*(Nitem[j]*beta[j] + (beta[j]*beta[j])/sig2Beta + likeOldJ[j]);
	likePropJ[j]=-.5*(Nitem[j]*bProp[j] + (bProp[j]*bProp[j])/sig2Beta + likePropJ[j]);
	if (rbinom(1,fmin2(1,exp(likePropJ[j]-likeOldJ[j]))))
	  {
	    beta[j]=bProp[j];
	    b0[N[0]+I[0]+j]+=1;
	  }
      }
    if(Hier[0]==1)
      {
	postB=sumsqr(beta,J[0])/2.0+sig2B[0];
	sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
      }
  }
  //sample theta
  if(sampCov[0])
    {
      double thetaProp;
      likeOldN[0]=0;
      likePropN[0]=0;
      thetaProp=theta+rnorm(0,met[N[0]+I[0]+J[0]+2]);
     
      for(r=0;r<R[0];r++) 
	{
	  likeOldN[0]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + beta[item[r]]+theta*cov[r]);
	  likePropN[0]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu[cond[r]] + alpha[subj[r]] + beta[item[r]]+thetaProp*cov[r]);
	}
      
      likeOldN[0]=-.5*((theta*theta)/sig2Mu[0] + likeOldN[0]);
      likePropN[0]=-.5*((thetaProp*thetaProp)/sig2Mu[0] + likePropN[0]);
      if (rbinom(1,fmin2(1,exp(likePropN[0]-likeOldN[0]))))
	{
	  theta=thetaProp;
	  b0[N[0]+I[0]+J[0]+2]+=1;
	}
    }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
    shift=rnorm(0,met[N[0]+I[0]+J[0]]);
    for(n=0;n<N[0];n++) muProp[n]=mu[n]+shift;
    for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
    p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(mu,N[0])-sumsqr(muProp,N[0]))/(2*sig2Mu[0]);
    if (rbinom(1,fmin2(1,exp(p))))
      {
	for(n=0;n<N[0];n++) mu[n]=muProp[n];
	for(i=0;i<I[0];i++) alpha[i]=aProp[i];
	b0[N[0]+I[0]+J[0]]+=1;
      }
  }
  
  //Decorrelate Alpha and Beta
  if(I[0]>1 & J[0]>1){
    shift=rnorm(0,met[N[0]+I[0]+J[0]+1]);
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]+shift;
  for(j=0;j<J[0];j++) bProp[j]=beta[j]-shift;
  
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(bProp,J[0]))/(2*sig2Beta);
  
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      for(j=0;j<J[0];j++) beta[j]=bProp[j];
      b0[N[0]+I[0]+J[0]+1]+=1;
    }
  }

  for(n=0;n<N[0];n++) blockS2[n]=mu[n];
  for(i=0;i<I[0];i++) blockS2[N[0]+i]=alpha[i];
  for(j=0;j<J[0];j++) blockS2[N[0]+I[0]+j]=beta[j];
  blockS2[N[0]+I[0]+J[0]]=sig2Alpha;
  blockS2[N[0]+I[0]+J[0]+1]=sig2Beta;
  blockS2[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();
}

/////////////////////////////////////////////////
void sampleNormalR(double *block,double *phi,double *blockD,double *dat,int *subj, int *item,double *cov,int *I,int *J,int *R,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *sig2,int *sampCov)
{
  int r,i,j;
  double sumY=0,sumYI[I[0]],sumYJ[J[0]],sumYprop;
  double a,postB;
  double p,shift,muProp,aProp[I[0]],bProp[J[0]];

  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;
  double phiA=phi[0],phiB=phi[1];
  double alpha[I[0]],beta[J[0]];
  double mu=block[0];
  for(i=0;i<I[0];i++) alpha[i]=block[i+1];
  for(j=0;j<J[0];j++) beta[j]=block[1+I[0]+j];
  double sig2Alpha=block[I[0]+J[0]+1];
  double sig2Beta=block[I[0]+J[0]+2];
  double theta=blockD[I[0]+J[0]+3]; 
  GetRNGstate(); 

  for(r=0;r<R[0];r++) sumY+=dat[r]-alpha[subj[r]]-beta[item[r]]-theta*cov[r];  
  a=1/(R[0]/sig2[0]+1/sig2Mu[0]);
  mu=rnorm(sumY/sig2[0]*a,sqrt(a));
  
  //sample alpha
  if(I[0]>1){
  for(i=0;i<I[0];i++) sumYI[i]=0;
  for(r=0;r<R[0];r++) sumYI[subj[r]]+=dat[r]-mu-beta[item[r]]-theta*cov[r];
  
  for(i=0;i<I[0];i++) 
    {
      a=1/(Nsub[i]/sig2[0]+1/sig2Alpha);
      alpha[i]=rnorm((sumYI[i]/sig2[0]+(phiA*blockD[i+1])/sig2Alpha)*a,sqrt(a));
    }
  sumY=0;
  for(i=0;i<I[0];i++) sumY+=(alpha[i]-(phiA*blockD[i+1]))*(alpha[i]-(phiA*blockD[i+1]));
  postB=sumY/2.0+sig2B[0];
  sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
  }
 
  if(J[0]>1){
  //sample beta
  for(j=0;j<J[0];j++) sumYJ[j]=0;
  for(r=0;r<R[0];r++) sumYJ[item[r]]+=dat[r]-mu-alpha[subj[r]]-theta*cov[r];
  
  for(j=0;j<J[0];j++)
    {
      a=1/(Nitem[j]/sig2[0]+1/sig2Beta);
      beta[j]=rnorm((sumYJ[j]/sig2[0]+(phiB*blockD[j+I[0]+1])/sig2Beta)*a,sqrt(a));
    }
  sumY=0;
  for(j=0;j<J[0];j++) sumY+=(beta[j]-phiB*blockD[j+I[0]+1])*(beta[j]-phiB*blockD[j+I[0]+1]);
  postB=sumY/2.0+sig2B[0];
  sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
  }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
    sumY=0;
    sumYprop=0;
    shift=rnorm(0,sqrt(s2decor1[0]));
    muProp=mu+shift;
    for(i=0;i<I[0];i++) 
      {
	aProp[i]=alpha[i]-shift;
	sumY+=(alpha[i]-phiA*blockD[i+1])*(alpha[i]-phiA*blockD[i+1]);
	sumYprop+=(aProp[i]-phiA*blockD[i+1])*(aProp[i]-phiA*blockD[i+1]);
      }  
	  
    p=(sumY-sumYprop)/(2*sig2Alpha) + (mu*mu-muProp*muProp)/(2*sig2Mu[0]);
    if (rbinom(1,fmin2(1,exp(p))))
      {
	mu=muProp;
	for(i=0;i<I[0];i++) alpha[i]=aProp[i];
	b0[0]+=1;
      }
  }
  
  //Decorrelate Alpha and Beta
  if(1==0){
  if(I[0]>1 & J[0]>1){
    sumY=0;
    sumYprop=0;
    shift=rnorm(0,sqrt(s2decor2[0]));
  for(i=0;i<I[0];i++) 
    {
      aProp[i]=alpha[i]+shift;
      sumY+=(alpha[i]-phiA*blockD[i+1])*(alpha[i]-phiA*blockD[i+1]);
      sumYprop+=(aProp[i]-phiA*blockD[i+1])*(aProp[i]-phiA*blockD[i+1]);
    }
  p=(sumY-sumYprop)/(2*sig2Alpha);
  
  sumY=0;
  sumYprop=0;
  for(j=0;j<J[0];j++)
    {
      bProp[j]=beta[j]-shift;
      sumY+=(beta[j]-phiB*blockD[j+I[0]+1])*(beta[j]-phiB*blockD[j+I[0]+1]);
      sumYprop+=(bProp[j]-phiB*blockD[j+I[0]+1])*(bProp[j]-phiB*blockD[j+I[0]+1]);
    }
  p+=(sumY-sumYprop)/(2*sig2Beta);
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      for(j=0;j<J[0];j++) beta[j]=bProp[j];
      b0[1]+=1;
    }
  }
  }


  block[0]=mu;
  for(i=0;i<I[0];i++) block[i+1]=alpha[i];
  for(j=0;j<J[0];j++) block[1+I[0]+j]=beta[j];
  block[I[0]+J[0]+1]=sig2Alpha;
  block[I[0]+J[0]+2]=sig2Beta;
  block[I[0]+J[0]+3]=theta;

  PutRNGstate();
}

