#include <R.h>
#include <Rmath.h>


void logLikeDpsd(double *like,int *R,int *NN, int *NS,int *I,int *JN,int *JS,int *K,int *dat,int *cond,int *Scond, int *subj, int *item,double *lag,double *blockN,double *blockS,double *blockR, double *critVec)
{
  int r,k;
  double mu,pR,crit[K[0]+1];
  crit[0]=-INFINITY;
  crit[K[0]]=INFINITY;
  for(r=0;r<R[0];r++)
    {
      for(k=1;k<K[0];k++)
	{
	  crit[k]=critVec[k+subj[r]*(K[0]+1)];
	}
      
      if(Scond[r]==0)
	{
	  mu=blockN[cond[r]]+blockN[NN[0]+subj[r]]+blockN[NN[0]+I[0]+item[r]]+lag[r]*blockN[NN[0]+I[0]+JN[0]+2];
	  like[r]=log(pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0));
	}
      
      if(Scond[r]==1)
	{	
	  mu=blockS[cond[r]]+blockS[NS[0]+subj[r]]+blockS[NS[0]+I[0]+item[r]]+lag[r]*blockS[NS[0]+I[0]+JS[0]+2];
	  pR=pnorm(blockR[cond[r]]+blockR[NS[0]+subj[r]]+blockR[NS[0]+I[0]+item[r]]+lag[r]*blockR[NS[0]+I[0]+JS[0]+2],0,1,1,0);
	  
	  if(dat[r]<(K[0]-1)) like[r]=log((1-pR) * (pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0)));	      
	  if(dat[r]==(K[0]-1)) like[r]=log(pR + (1-pR)*pnorm(crit[K[0]-1],mu,1,0,0));
	}
    }
}




void logLikeDpsdPos(double *like,int *R,int *NN, int *NS,int *I,int *JN,int *JS,int *K,int *dat,int *cond,int *Scond, int *subj, int *item,double *lag,double *blockN,double *blockD,double *blockR, double *critVec)
{
  int r,k;
  double muN,muS,pR,crit[K[0]+1];
  crit[0]=-INFINITY;
  crit[K[0]]=INFINITY;
  for(r=0;r<R[0];r++)
    {
      for(k=1;k<K[0];k++)
	{
	  crit[k]=critVec[k+subj[r]*(K[0]+1)];
	}

      muN=blockN[cond[r]]+blockN[NN[0]+subj[r]]+blockN[NN[0]+I[0]+item[r]]+lag[r]*blockN[NN[0]+I[0]+JN[0]+2];

      if(Scond[r]==0)
	{
	  like[r]=log(pnorm(crit[dat[r]+1],muN,1,1,0)-pnorm(crit[dat[r]],muN,1,1,0));
	}
      
      if(Scond[r]==1)
	{	
	  muS=muN+exp(blockD[cond[r]]+blockD[NS[0]+subj[r]]+blockD[NS[0]+I[0]+item[r]]+lag[r]*blockD[NS[0]+I[0]+JS[0]+2]);
	  pR=pnorm(blockR[cond[r]]+blockR[NS[0]+subj[r]]+blockR[NS[0]+I[0]+item[r]]+lag[r]*blockR[NS[0]+I[0]+JS[0]+2],0,1,1,0);
	  
	  if(dat[r]<(K[0]-1)) like[r]=log((1-pR) * (pnorm(crit[dat[r]+1],muS,1,1,0)-pnorm(crit[dat[r]],muS,1,1,0)));	      
	  if(dat[r]==(K[0]-1)) like[r]=log(pR + (1-pR)*pnorm(crit[K[0]-1],muS,1,0,0));
	}
    }
}



