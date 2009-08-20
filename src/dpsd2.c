#include "normal.h"

//Sample mu in alphas, for non-hierarchical & no items
void sampleNormal2(double *block,double *dat, int *subj,int *N,int *I,int *J,int *R,int *Nsub)
{
  int r,n,i,j;
  double sumYI[I[0]];
  double a;
  double mu[N[0]],alpha[I[0]],beta[J[0]];
  for(n=0;n<N[0];n++) mu[n]=block[n];
  for(i=0;i<I[0];i++) alpha[i]=block[N[0]+i];
  for(j=0;j<J[0];j++) beta[j]=block[N[0]+I[0]+j];
  double sig2Alpha=block[N[0]+I[0]+J[0]];
  double sig2Beta=block[N[0]+I[0]+J[0]+1];
  double theta=block[N[0]+I[0]+J[0]+2];
  GetRNGstate(); 

  //sample alpha
  for(i=0;i<I[0];i++) sumYI[i]=0;
  for(r=0;r<R[0];r++) sumYI[subj[r]]+=dat[r];
  
  for(i=0;i<I[0];i++) 
    {
      a=1/(Nsub[i]+1/sig2Alpha);
      alpha[i]=rnorm(sumYI[i]*a,sqrt(a));
    }
    

  for(n=0;n<N[0];n++) block[n]=mu[n];
  for(i=0;i<I[0];i++) block[i+N[0]]=alpha[i];
  for(j=0;j<J[0];j++) block[N[0]+I[0]+j]=beta[j];
  block[N[0]+I[0]+J[0]]=sig2Alpha;
  block[N[0]+I[0]+J[0]+1]=sig2Beta;
  block[N[0]+I[0]+J[0]+2]=theta;

  PutRNGstate();
}

