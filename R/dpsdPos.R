.packageName='hbmem'
.First.lib=function(lib,pkg) library.dynam('hbmem',pkg,lib)

dpsdPosSim=function(I=30,J=200,K=6,muN=-.7,s2aN=.2,s2bN=.2,muD=0,s2aD=.2,s2bD=.2,muR=qnorm(.25),s2aR=.2,s2bR=.2,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
  {
    R=I*J
    alphaN=rnorm(I,0,sqrt(s2aN))
    betaN=rnorm(J,0,sqrt(s2bN))
    alphaD=rnorm(I,0,sqrt(s2aD))
    betaD=rnorm(J,0,sqrt(s2bD))
    alphaR=rnorm(I,0,sqrt(s2aR))
    betaR=rnorm(J,0,sqrt(s2bR))
    subj=rep(0:(I-1),each=J)
    item=rep(0:(J-1),I)
    lag=rep(0,R)
    cond.sub.A=rep(0:1,J/2)
    cond.sub.B=rep(1:0,J/2)
    cond=rep(c(cond.sub.A,cond.sub.B),I/2)

    dat=1:R
    for(r in 1:R)
      {
        meanN=muN+alphaN[subj[r]+1]+betaN[item[r]+1]
        if(cond[r]==0) p=dpsdProbs(0,meanN,crit[subj[r]+1,])
        if(cond[r]==1)
          {
            delta=exp(muD + alphaD[subj[r]+1]+betaD[item[r]+1])
            meanS=meanN+delta
            pR=pnorm(muR+alphaR[subj[r]+1]+betaR[item[r]+1])
            p=dpsdProbs(pR,meanS,crit[subj[r]+1,])
          }
           dat[r]=which.max(rmultinom(1,1,p))-1
      }
     
    ret=new("dpsdSim")
    ret@cond=cond
    ret@subj=subj
    ret@item=item
    ret@lag=lag
    ret@resp=dat
    ret@muN=muN
    ret@muS=muD
    ret@muR=muR
    ret@alphaN=alphaN
    ret@alphaS=alphaD
    ret@alphaR=alphaR
    ret@betaN=betaN
    ret@betaS=betaD
    ret@betaR=betaR
    return(ret)
  }


dpsdPosLogLike=function(R,I,J,K,dat,cond,sub,item,lag,blockN,blockD,blockR,crit)
  {
    l=1:R*0
    .C("logLikeDpsdPos",as.double(l),as.integer(R),as.integer(I),as.integer(J),as.integer(K),as.integer(dat),as.integer(cond),as.integer(sub),as.integer(item),as.double(lag),as.double(blockN),as.double(blockD),as.double(blockR),as.double(as.vector(t(crit))),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }


dpsdPosSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.01)
{
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag
resp=dat$resp


#DEFINE CONSTANTS
I=length(levels(as.factor(sub)))
J=length(levels(as.factor(item)))
K=length(levels(as.factor(resp)))
nsubT=table(sub)
nitemT=table(item)
nsubS=table(sub[cond==1])
nitemS=table(item[cond==1])
R=dim(dat)[1]
RS=sum(cond)
RN=R-sum(cond)
B=I+J+4

#INDEXING
mu=1
alpha=2:(I+1)
beta=(max(alpha)+1):(I+J+1)
s2alpha=max(beta)+1
s2beta=s2alpha+1
theta=s2beta+1

#SPACE AND STARTING VALUES
blockN=blockD=blockR=matrix(0,nrow=M,ncol=B)
wS=wR=rnorm(RS)
allwN=rnorm(R)
s.crit=array(dim=c(I,7,M))
s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)
b0=matrix(0,2,2)
met=matrix(.01,2,2)
b0D=rep(0,B)
metD=c(.01,rep(.2,(B-4)),.05,.05,.01)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]

notk=cond==1 & resp<(K-1)
Rnotk=sum(notk)
isk= cond==1 & resp==(K-1) 
Risk=sum(isk)

pb=txtProgressBar(min=1,max=M,style=3,width=10)
print("Starting MCMC")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
meanN=getPred(blockN[m-1,],sub,item,lag,I,J,R)
meanD=exp(getPred(blockD[m-1,],sub,item,lag,I,J,R))
meanS=meanN+meanD
meanR=getPred(blockR[m-1,],sub,item,lag,I,J,R)

wN=rtnorm(RN,meanN[cond==0],rep(1,RN),s.crit[cbind(sub[cond==0]+1,resp[cond==0]+1,m-1)],s.crit[cbind(sub[cond==0]+1,resp[cond==0]+2,m-1)])

wS[notk[cond==1]]=rtnorm(Rnotk,meanS[notk],rep(1,RN),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[cond==1]]=rtnorm(Rnotk,meanR[notk],rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[cond==1] & wR>0
wS[ind]=rnorm(sum(ind),meanS[cond==1][ind],1)
ind=isk[cond==1] & wR<0
wS[ind]=rtnorm(sum(ind),meanS[cond==1][ind],rep(1,sum(ind)),s.crit[cbind(sub[cond==1][ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[cond==1] & wS>s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),meanR[cond==1][ind],1)
ind=isk[cond==1] & wS<s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),meanR[cond==1][ind],rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))


#Sample Blocks
wD=wS-meanN[cond==1]
tmp=samplePosNorm(blockD[m-1,],wD,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,10,.01,.01,metD,1,1)
blockD[m,]=tmp[[1]]
b0D=b0D+tmp[[2]]

meanD=exp(getPred(blockD[m,],sub,item,lag,I,J,R))
allwN[cond==0]=wN
allwN[cond==1]=wS-meanD[cond==1]
tmp=sampleNorm(blockN[m-1,],allwN,sub,item,lag,I,J,R,nsubT,nitemT,10,.01,.01,met[1,1],met[1,2],1,0)
blockN[m,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

tmp=sampleNorm(blockR[m-1,],wR,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,nsubS,nitemS,10,.01,.01,met[2,1],met[2,2],1,1)
blockR[m,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[2]]

#Sample Criteria
PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rnorm(I*4,0,rep(met.crit,4))
violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,m-1]

likeOld=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
likeProp=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],PropCrit),sub,sum)
accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
b0.crit=b0.crit+accept
s.crit[accept==1,,m]=PropCrit[accept==1,]
s.crit[accept==0,,m]=s.crit[accept==0,,m-1]


#AUTOTUNING
if(m>20 & m<keep[1] & m%%10==0)
  {
    met=met+(b0/m<.3)*matrix(-jump,2,2) +(b0/m>.5)*matrix(jump,2,2)
    metD=metD+(b0D/m<.3)*-jump + (b0D/m>.5)*jump
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
    met[met<jump]=jump
    metD[metD<jump]=jump
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    print("")
    print("Continuning with following acceptance probabilities:")
    print(b0/m)
    print(b0D/m)
    print(b0.crit/m)
    print("If they are <.2 or >.6, adjust jump or keep")
  }

setTxtProgressBar(pb, m)
}
close(pb)

estN=colMeans(blockN[keep,])
estD=colMeans(blockD[keep,])
estR=colMeans(blockR[keep,])
estCrit=apply(s.crit[,,keep],c(1,2),mean)

#GET DIC
DIC=pD=NA
if(getDIC)
{
print("Getting DIC")
pb=txtProgressBar(min=keep[1],max=tail(keep,1),style=3,width=10)

D0=0
for(m in keep)
  {
    D0=D0+sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,estN,estD,estR,estCrit))
pD=D0-Dhat
DIC=pD+D0
}


u=new("dpsd")
u@mu=mu
u@alpha=alpha
u@beta=beta
u@s2alpha=s2alpha
u@s2beta=s2beta
u@theta=theta
u@estN=estN
u@estS=estD
u@estR=estR
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockD[keep,]
u@blockR=blockR[keep,]
u@s.crit=s.crit[,,keep]
u@pD=pD
u@DIC=DIC
u@M=M
u@keep=keep
u@b0=b0/M
u@b0Crit=b0.crit/M
print("")
print("Finished")
return(u)
}
