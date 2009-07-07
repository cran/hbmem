#include <R.h>
#include <Rmath.h>

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

void getPred(double *pred,double *block,int *subj, int *item,double *cov, int *I, int *J, int *R)
{
  int r;
  for(r=0;r<R[0];r++) pred[r]=(block[0]+block[1+subj[r]]+block[1+I[0]+item[r]]+cov[r]*block[I[0]+J[0]+3]);
}

void logLikeUvsd(double *like,int *R,int *I,int *J,int *dat,int *subj, int *item,double *cov,int *cond,double *blockN,double *blockS,double *blockS2, double *critVec)
{
  int r;
  double mu,sigma,crit[7];
  crit[0]=-INFINITY;
  crit[3]=0;
  crit[6]=INFINITY;

  for(r=0;r<R[0];r++)
    {
      crit[1]=critVec[1+subj[r]*7];
      crit[2]=critVec[2+subj[r]*7];
      crit[4]=critVec[4+subj[r]*7];
      crit[5]=critVec[5+subj[r]*7];
	
      if(cond[r]==0)
	{
	  mu=blockN[0]+blockN[1+subj[r]]+blockN[1+I[0]+item[r]]+cov[r]*blockN[I[0]+J[0]+3];
	  sigma=1;
	  like[r]=log(pnorm(crit[dat[r]+1],mu,sigma,1,0)-pnorm(crit[dat[r]],mu,sigma,1,0));
	}
      
      if(cond[r]==1)
	{	
	  mu=blockS[0]+blockS[1+subj[r]]+blockS[1+I[0]+item[r]]+cov[r]*blockS[I[0]+J[0]+3];
	  sigma=sqrt(exp(blockS2[0]+blockS2[1+subj[r]]+blockS2[1+I[0]+item[r]]));
	  like[r]=log(pnorm(crit[dat[r]+1],mu,sigma,1,0)-pnorm(crit[dat[r]],mu,sigma,1,0));
	}  

      if(like[r]==0) printf("Warning: Log Likelihood = 0");
    }
}

void logLikeDpsd(double *like,int *R,int *I,int *J,int *K,int *dat,int *cond,int *subj, int *item,double *lag,double *blockN,double *blockS,double *blockR, double *critVec)
{
  int r;
  double mu,pR,crit[7];
  crit[0]=-INFINITY;
  crit[3]=0;
  crit[K[0]]=INFINITY;

  for(r=0;r<R[0];r++)
    {
      crit[1]=critVec[1+subj[r]*7];
      crit[2]=critVec[2+subj[r]*7];
      crit[4]=critVec[4+subj[r]*7];
      crit[5]=critVec[5+subj[r]*7];
	
      if(cond[r]==0)
	{
	  mu=blockN[0]+blockN[1+subj[r]]+blockN[1+I[0]+item[r]]+lag[r]*blockN[I[0]+J[0]+3];
	  like[r]=log(pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0));
	}
      
      if(cond[r]==1)
	{	
	  mu= blockS[0]+blockS[1+subj[r]]+blockS[1+I[0]+item[r]]+lag[r]*blockS[I[0]+J[0]+3];
	  pR=pnorm(blockR[0]+blockR[1+subj[r]]+blockR[1+I[0]+item[r]]+lag[r]*blockR[I[0]+J[0]+3],0,1,1,0);
	  
	  if(dat[r]<(K[0]-1)) like[r]=log((1-pR) * (pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0)));
	      
	  if(dat[r]==(K[0]-1)) like[r]=log(pR + (1-pR)*pnorm(crit[K[0]-1],mu,1,0,0));

	}
    }
}



//Sample Block with ONE sigma
void sampleNormal(double *block,double *dat,int *subj, int *item,double *cov,int *I,int *J,int *R,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *sig2,int *sampCov)
{
  int r,i,j;
  double sumY=0,sumYI[I[0]],sumYJ[J[0]];
  double a,postB;
  double p,shift,muProp,aProp[I[0]],bProp[J[0]];

  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;

  double alpha[I[0]],beta[J[0]];
  double mu=block[0];
  for(i=0;i<I[0];i++) alpha[i]=block[i+1];
  for(j=0;j<J[0];j++) beta[j]=block[1+I[0]+j];
  double sig2Alpha=block[I[0]+J[0]+1];
  double sig2Beta=block[I[0]+J[0]+2];
  double theta=block[I[0]+J[0]+3];
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
      alpha[i]=rnorm(sumYI[i]/sig2[0]*a,sqrt(a));
    }
  postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
  sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
  }
  
  if(J[0]>1){
  //sample beta
  for(j=0;j<J[0];j++) sumYJ[j]=0;
  for(r=0;r<R[0];r++) sumYJ[item[r]]+=dat[r]-mu-alpha[subj[r]]-theta*cov[r];
  
  for(j=0;j<J[0];j++)
    {
      a=1/(Nitem[j]/sig2[0]+1/sig2Beta);
      beta[j]=rnorm(sumYJ[j]/sig2[0]*a,sqrt(a));
    }
  postB=sumsqr(beta,J[0])/2.0+sig2B[0];
  sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
  }

 
  //sample theta
  if(sampCov[0])
    {
      sumY=0;
      for(r=0;r<R[0];r++) sumY+=cov[r]*(dat[r]-mu-alpha[subj[r]]-beta[item[r]]);
      a=1/(sumsqr(cov,R[0])/sig2[0]+1/sig2Mu[0]);
      theta=rnorm(sumY/sig2[0]*a,sqrt(a)); 
    }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
  shift=rnorm(0,sqrt(s2decor1[0]));
  muProp=mu+shift;
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (mu*mu-muProp*muProp)/(2*sig2Mu[0]);
  if (rbinom(1,fmin2(1,exp(p))))
    {
      mu=muProp;
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

  block[0]=mu;
  for(i=0;i<I[0];i++) block[i+1]=alpha[i];
  for(j=0;j<J[0];j++) block[1+I[0]+j]=beta[j];
  block[I[0]+J[0]+1]=sig2Alpha;
  block[I[0]+J[0]+2]=sig2Beta;
  block[I[0]+J[0]+3]=theta;

  PutRNGstate();
}



//sample one sigma
void sampleSigma2(double *sigma2,double *block, double *dat, int *subj, int *item,double *cov, int *I, int *J, int *R,double *a, double *b)
{
  double postA=a[0] + R[0]/2.0;
  double resid[R[0]];
  int r;

  for(r=0;r<R[0];r++) resid[r]=dat[r]-(block[0]+block[1+subj[r]]+block[1+I[0]+item[r]]+cov[r]*block[I[0]+J[0]+3]);
  double postB=sumsqr(resid,R[0])/2.0 + b[0];
  GetRNGstate(); 
  sigma2[0]=1.0/rgamma(postA,1.0/postB);
  PutRNGstate();
}


void sampleNormalb(double *block,double *dat,int *subj, int *item,double *cov,int *I,int *J,int *R,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *blocks2,int *sampCov)
{
  int r,i,j;
  double sumY=0,sumYI[I[0]],sumYJ[J[0]],v;
  double a,postB;
  double p,shift,muProp,aProp[I[0]],bProp[J[0]];

  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;

  double alpha[I[0]],beta[J[0]];
  double mu=block[0];
  for(i=0;i<I[0];i++) alpha[i]=block[i+1];
  for(j=0;j<J[0];j++) beta[j]=block[1+I[0]+j];
  double sig2Alpha=block[I[0]+J[0]+1];
  double sig2Beta=block[I[0]+J[0]+2];
  double theta=block[I[0]+J[0]+3];

  double sigma2[R[0]],sumInvS2=0;
  for(r=0;r<R[0];r++) 
    {
      sigma2[r]=exp(blocks2[0]+blocks2[1+subj[r]]+blocks2[1+I[0]+item[r]]+cov[r]*blocks2[I[0]+J[0]+3]);
      sumInvS2+=1/sigma2[r];
    }

  GetRNGstate(); 

  //sample Mu
  a=1/(sumInvS2 + 1/sig2Mu[0]);
  for(r=0;r<R[0];r++) sumY+=(dat[r]-alpha[subj[r]]-beta[item[r]]-theta*cov[r])/sigma2[r];  
  mu=rnorm(a*sumY,sqrt(a));
 
  //sample alpha
  if(I[0]>1){
  for(i=0;i<I[0];i++) sumYI[i]=0;
  for(r=0;r<R[0];r++) sumYI[subj[r]]+=(dat[r]-mu-beta[item[r]]-theta*cov[r])/sigma2[r];
  
  for(i=0;i<I[0];i++) 
    {
      a=1/(Nsub[i]/exp(blocks2[0]+blocks2[1+i])+1/sig2Alpha);
      alpha[i]=rnorm(a*sumYI[i],sqrt(a));
    }
  postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
  sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
  }
  
  if(J[0]>1){
  //sample beta
  for(j=0;j<J[0];j++) sumYJ[j]=0;
  for(r=0;r<R[0];r++) sumYJ[item[r]]+=(dat[r]-mu-alpha[subj[r]]-theta*cov[r])/sigma2[r];
  
  for(j=0;j<J[0];j++)
    {
      a=1/(Nitem[j]/exp(blocks2[0]+blocks2[1+I[0]+j])+1/sig2Beta);
      beta[j]=rnorm(a*sumYJ[j],sqrt(a));
    }
  postB=sumsqr(beta,J[0])/2.0+sig2B[0];
  sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
  }

 
  //sample theta
  if(sampCov[0])
    {
      double ycov=0;
      sumY=0;
      a=0;
      for(r=0;r<R[0];r++) 
	{
	  sumY+=cov[r]*(dat[r]-mu-alpha[subj[r]]-beta[item[r]])/sigma2[r];
	  a+=(cov[r]*cov[r])/sigma2[r];
	}

      a=1/(a+1/sig2Mu[0]);
      theta=rnorm(a*sumY,sqrt(a)); 
    }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
  shift=rnorm(0,sqrt(s2decor1[0]));
  muProp=mu+shift;
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (mu*mu-muProp*muProp)/(2*sig2Mu[0]);
  if (rbinom(1,fmin2(1,exp(p))))
    {
      mu=muProp;
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

  block[0]=mu;
  for(i=0;i<I[0];i++) block[i+1]=alpha[i];
  for(j=0;j<J[0];j++) block[1+I[0]+j]=beta[j];
  block[I[0]+J[0]+1]=sig2Alpha;
  block[I[0]+J[0]+2]=sig2Beta;
  block[I[0]+J[0]+3]=theta;

  PutRNGstate();
}



////////////////////////////////////////////////
void sampleSigma2b(double *blockS2,double *dat,int *subj, int *item,double *cov,int *I,int *J,int *R,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *met,int *b0,double *blockMu,int *sampCov)
{
  int r,i,j;
  double S=0,likeProp=0,likeOld=0,likePropI[I[0]],likeOldI[I[0]],likePropJ[J[0]],likeOldJ[J[0]];
  double a,postB;
  double p,shift,muProp,aProp[I[0]],bProp[J[0]];
  double postAlphaA=sig2A[0]+(1.0*I[0])/2.0;
  double postBetaA=sig2A[0]+(1.0*J[0])/2.0;
  double alpha[I[0]],beta[J[0]];
  double mu=blockS2[0];
  for(i=0;i<I[0];i++) alpha[i]=blockS2[i+1];
  for(j=0;j<J[0];j++) beta[j]=blockS2[1+I[0]+j];
  double sig2Alpha=blockS2[I[0]+J[0]+1];
  double sig2Beta=blockS2[I[0]+J[0]+2];
  double theta=blockS2[I[0]+J[0]+3];

  double mean[R[0]],sigma2[R[0]];
 
  for(r=0;r<R[0];r++) mean[r]=blockMu[0]+blockMu[1+subj[r]]+blockMu[1+I[0]+item[r]]+cov[r]*blockMu[I[0]+J[0]+3];
      
  GetRNGstate(); 

  //sample Mu
  muProp=mu+rnorm(0,met[0]);
  for(r=0;r<R[0];r++) 
    {
      likeOld+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu + alpha[subj[r]] + beta[item[r]]-theta*cov[r]);
      likeProp+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(muProp + alpha[subj[r]] + beta[item[r]]-theta*cov[r]);
    }
  likeOld=-.5*(R[0]*mu + (mu*mu)/sig2Mu[0] + likeOld);
  likeProp=-.5*(R[0]*muProp + (muProp*muProp)/sig2Mu[0] + likeProp);
  if (rbinom(1,fmin2(1,exp(likeProp-likeOld))))
    {
      mu=muProp;
      b0[0]+=1;
    }


  //sample alpha
  if(I[0]>1){
    for(i=0;i<I[0];i++) 
      {
	likeOldI[i]=0;	
	likePropI[i]=0;
	aProp[i]=alpha[i]+rnorm(0,met[i+1]);
      }
    
  for(r=0;r<R[0];r++) 
    {
      likeOldI[subj[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu + alpha[subj[r]] + beta[item[r]]-theta*cov[r]);
      likePropI[subj[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu + aProp[subj[r]] + beta[item[r]]-theta*cov[r]);
    }

    for(i=0;i<I[0];i++) 
      {
	likeOldI[i]=-.5*(Nsub[i]*alpha[i] + (alpha[i]*alpha[i])/sig2Alpha + likeOldI[i]);
	likePropI[i]=-.5*(Nsub[i]*aProp[i] + (aProp[i]*aProp[i])/sig2Alpha + likePropI[i]);
	if (rbinom(1,fmin2(1,exp(likePropI[i]-likeOldI[i]))))
	  {
	    alpha[i]=aProp[i];
	    b0[i+1]+=1;
	  }
      }

    postB=sumsqr(alpha,I[0])/2.0+sig2B[0];
    sig2Alpha=1.0/rgamma(postAlphaA,1.0/postB);
  }


  //Sample Beta
  if(J[0]>1){
    for(j=0;j<J[0];j++) 
      {
	likeOldJ[j]=0;	
	likePropJ[j]=0;
	bProp[j]=beta[j]+rnorm(0,met[j+I[0]+1]);
      }
    

    for(r=0;r<R[0];r++) 
      {
	likeOldJ[item[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu + alpha[subj[r]] + beta[item[r]]-theta*cov[r]);
      likePropJ[item[r]]+=(dat[r]-mean[r])*(dat[r]-mean[r])/exp(mu + alpha[subj[r]] + bProp[item[r]]-theta*cov[r]);
    }


    for(j=0;j<J[0];j++) 
      {
	likeOldJ[j]=-.5*(Nitem[j]*beta[j] + (beta[j]*beta[j])/sig2Beta + likeOldJ[j]);
	likePropJ[j]=-.5*(Nitem[j]*bProp[j] + (bProp[j]*bProp[j])/sig2Beta + likePropJ[j]);
	if (rbinom(1,fmin2(1,exp(likePropJ[j]-likeOldJ[j]))))
	  {
	    beta[j]=bProp[j];
	    b0[j+I[0]+1]+=1;
	  }
      }

    postB=sumsqr(beta,J[0])/2.0+sig2B[0];
    sig2Beta=1.0/rgamma(postBetaA,1.0/postB);
  }



  ///NEED TO FIX THIS
  //sample theta
  if(sampCov[0])
    {
      double ycov=0,sumY=0;
      a=0;
      for(r=0;r<R[0];r++) 
	{
	  sumY+=cov[r]*(dat[r]-mu-alpha[subj[r]]-beta[item[r]])/sigma2[r];
	  a+=(cov[r]*cov[r])/sigma2[r];
	}

      a=1/(a+1/sig2Mu[0]);
      theta=rnorm(a*sumY,sqrt(a)); 
      printf("WARNING!!! Not setup to sample theta in Sigma!");
    }

  //DECORRELATE
  //Decorrelate Mu and Alpha
  if(I[0]>1){
    shift=rnorm(0,met[1+I[0]+J[0]]);
    muProp=mu+shift;
    for(i=0;i<I[0];i++) aProp[i]=alpha[i]-shift;
    p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (mu*mu-muProp*muProp)/(2*sig2Mu[0]);
    if (rbinom(1,fmin2(1,exp(p))))
      {
	mu=muProp;
	for(i=0;i<I[0];i++) alpha[i]=aProp[i];
	b0[1+I[0]+J[0]]+=1;
      }
  }
  
  //Decorrelate Alpha and Beta
  if(I[0]>1 & J[0]>1){
    shift=rnorm(0,met[I[0]+J[0]+2]);
  for(i=0;i<I[0];i++) aProp[i]=alpha[i]+shift;
  for(j=0;j<J[0];j++) bProp[j]=beta[j]-shift;
  
  p=(sumsqr(alpha,I[0])-sumsqr(aProp,I[0]))/(2*sig2Alpha) + (sumsqr(beta,J[0])-sumsqr(bProp,J[0]))/(2*sig2Beta);
  
  if (rbinom(1,fmin2(1,exp(p))))
    {
      for(i=0;i<I[0];i++) alpha[i]=aProp[i];
      for(j=0;j<J[0];j++) beta[j]=bProp[j];
      b0[I[0]+J[0]+2]+=1;
    }
  }

  blockS2[0]=mu;
  for(i=0;i<I[0];i++) blockS2[i+1]=alpha[i];
  for(j=0;j<J[0];j++) blockS2[1+I[0]+j]=beta[j];
  blockS2[I[0]+J[0]+1]=sig2Alpha;
  blockS2[I[0]+J[0]+2]=sig2Beta;
  blockS2[I[0]+J[0]+3]=theta;

  PutRNGstate();
}

/////////////////////////////////////////////////
