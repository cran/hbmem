#include <R.h>
#include <Rmath.h>

double sumsqr(double *invec, int size);

double sum(double *invec, int size);

void rtruncnorm(double *y, int *N,double *mu, double *sigma, double 
*a,double *b);

void getPred(double *pred,double *block,int *cond,int *subj, int *item,double *cov, int *N, int *I, int *J, int *R);

void sampleNormal(double *block,double *dat,int *cond, int *subj, int *item,double *cov,int *N,int *I,int *J,int *R,int 
*Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double 
*s2decor1,double *s2decor2,int *b0,double *sig2,int *sampCov,int *Hier);

void sampleSigma2(double *sigma2,double *block, double *dat, int *cond, int *subj,int *item,double *cov,int *N,int 
*ncond,int *I, int *J,double *a, double *b);

void sampleNormalb(double *block,double *dat,int *cond,int *subj, int *item,double *cov,int *N,int *I,int *J,int *R,int 
*Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double 
*s2decor1,double *s2decor2,int *b0,double *blocks2,int *sampCov,int *Hier);

void sampleSigma2b(double *blockS2,double *dat,int *cond,int *subj, int *item,double *cov,int *N,int *I,int *J,int 
*R,int *Ncond,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double 
*met,int *b0,double *blockMu,int *sampCov,int *Hier);

void sampleNormalR(double *block,double *phi,double *blockD,double *dat,int *subj, int *item,double *cov,int *I,int 
*J,int *R,int *Nsub, int *Nitem, double *sig2Mu,double *sig2A,double *sig2B,double *s2decor1,double *s2decor2,int *b0,double *sig2,int *sampCov);








