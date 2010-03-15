#include <R.h>
#include <Rmath.h>

void logLikeUvsd(double *like,int *R,int *NN, int *NS,int *I,int *J,int *dat,int *subj, int *item,double *cov,int *cond,int *Scond,double *blockN,double *blockS,double *blockS2, double *critVec)
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
	
      if(Scond[r]==0)
	{
	  mu=blockN[cond[r]]+blockN[NN[0]+subj[r]]+blockN[NN[0]+I[0]+item[r]]+cov[r]*blockN[NN[0]+I[0]+J[0]+2];
	  sigma=1;
	  like[r]=log(pnorm(crit[dat[r]+1],mu,sigma,1,0)-pnorm(crit[dat[r]],mu,sigma,1,0));
	}
      
      if(Scond[r]==1)
	{	
mu=blockS[cond[r]]+blockS[NS[0]+subj[r]]+blockS[NS[0]+I[0]+item[r]]+cov[r]*blockS[NS[0]+I[0]+J[0]+2];	  
sigma=sqrt(exp(blockS2[cond[r]]+blockS2[NS[0]+subj[r]]+blockS2[NS[0]+I[0]+item[r]]+cov[r]*blockS2[NS[0]+I[0]+J[0]+2]));
	  like[r]=log(pnorm(crit[dat[r]+1],mu,sigma,1,0)-pnorm(crit[dat[r]],mu,sigma,1,0));
	}  

      //if(like[r]==-INFINITY) printf("LL = 0, DP= %lf, S2= %lf, c1= %lf, c2=%lf\n",mu,sigma,crit[dat[r]+1],crit[dat[r]]);
    }
}


