#include <R.h>
#include <Rmath.h>


void logLikeDpsd(double *like,int *R,int *NN, int *NS,int *I,int *J,int *K,int *dat,int *cond,int *Scond, int *subj, int *item,double *lag,double *blockN,double *blockS,double *blockR, double *critVec)
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
	
      if(Scond[r]==0)
	{
	  mu=blockN[cond[r]]+blockN[NN[0]+subj[r]]+blockN[NN[0]+I[0]+item[r]]+lag[r]*blockN[NN[0]+I[0]+J[0]+2];
	  like[r]=log(pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0));
	}
      
      if(Scond[r]==1)
	{	
	  mu= blockS[cond[r]]+blockS[NS[0]+subj[r]]+blockS[NS[0]+I[0]+item[r]]+lag[r]*blockS[NS[0]+I[0]+J[0]+2];
	  pR=pnorm(blockR[cond[r]]+blockR[NS[0]+subj[r]]+blockR[NS[0]+I[0]+item[r]]+lag[r]*blockR[NS[0]+I[0]+J[0]+2],0,1,1,0);
	  
	  if(dat[r]<(K[0]-1)) like[r]=log((1-pR) * (pnorm(crit[dat[r]+1],mu,1,1,0)-pnorm(crit[dat[r]],mu,1,1,0)));	      
	  if(dat[r]==(K[0]-1)) like[r]=log(pR + (1-pR)*pnorm(crit[K[0]-1],mu,1,0,0));

	}
    }
}

void logLikeDpsdPos(double *like,int *R,int *I,int *J,int *K,int *dat,int *cond,int *subj, int *item,double *lag,double *blockN,double *blockD,double *blockR, double *critVec)
{
  int r;
  double muN,muS,pR,crit[7];
  crit[0]=-INFINITY;
  crit[3]=0;
  crit[K[0]]=INFINITY;

  for(r=0;r<R[0];r++)
    {
      crit[1]=critVec[1+subj[r]*7];
      crit[2]=critVec[2+subj[r]*7];
      crit[4]=critVec[4+subj[r]*7];
      crit[5]=critVec[5+subj[r]*7];
	  
      muN=blockN[0]+blockN[1+subj[r]]+blockN[1+I[0]+item[r]]+lag[r]*blockN[I[0]+J[0]+3];
	
      if(cond[r]==0)
	{
	  like[r]=log(pnorm(crit[dat[r]+1],muN,1,1,0)-pnorm(crit[dat[r]],muN,1,1,0));
	}
      
      if(cond[r]==1)
	{	
	  muS=muN + exp(blockD[0]+blockD[1+subj[r]]+blockD[1+I[0]+item[r]]+lag[r]*blockD[I[0]+J[0]+3]);
	  pR=pnorm(blockR[0]+blockR[1+subj[r]]+blockR[1+I[0]+item[r]]+lag[r]*blockR[I[0]+J[0]+3],0,1,1,0);
	  
	  if(dat[r]<(K[0]-1)) like[r]=log((1-pR) * (pnorm(crit[dat[r]+1],muS,1,1,0)-pnorm(crit[dat[r]],muS,1,1,0)));
	      
	  if(dat[r]==(K[0]-1)) like[r]=log(pR + (1-pR)*pnorm(crit[K[0]-1],muS,1,0,0));

	}
    }
}

