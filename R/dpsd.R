.packageName='hbmem'
.First.lib=function(lib,pkg) library.dynam('hbmem',pkg,lib)

setClass("dpsd",representation(
                               mu="numeric",
                               alpha="numeric",
                               beta="numeric",
                               s2alpha="numeric",
                               s2beta="numeric",
                               theta="numeric",
                               estN="numeric",
                               estS="numeric",
                               estR="numeric",
                               estCrit="matrix",
                               blockN="matrix",
                               blockS="matrix",
                               blockR="matrix",
                               s.crit="array",
                               pD="numeric",
                               DIC="numeric",
                               M="numeric",
                               keep="numeric",
                               b0="matrix",
                               b0Crit="numeric"))

setClass("dpsdSim",representation(
                                  cond="numeric",
                                  subj="numeric",
                                  item="numeric",
                                  lag="numeric",
                                  resp="numeric",
                                  muN="numeric",
                                  muS="numeric",
                                  muR="numeric",
                                  alphaN="numeric",
                                  betaN="numeric",
                                  alphaS="numeric",
                                  betaS="numeric",
                                  alphaR="numeric",
                                  betaR="numeric"))

dpsdProbs=function(r,d,crit)
{
K=length(crit)+1
theta=diff(c(0,pnorm(crit,d,1),1))
ind=c(rep(0,K-1),1)
p=r*ind+(1-r)*theta
return(p)
}

dpsdSim=function(I=30,J=200,K=6,muN=-.7,s2aN=.2,s2bN=.2,muS=0,s2aS=.2,s2bS=.2,muR=qnorm(.25),s2aR=.2,s2bR=.2,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
  {
    R=I*J
    alphaN=rnorm(I,0,sqrt(s2aN))
    betaN=rnorm(J,0,sqrt(s2bN))
    alphaS=rnorm(I,0,sqrt(s2aS))
    betaS=rnorm(J,0,sqrt(s2bS))
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
        if(cond[r]==0) p=dpsdProbs(0,muN+alphaN[subj[r]+1]+betaN[item[r]+1],crit[subj[r]+1,])
        if(cond[r]==1) p=dpsdProbs(pnorm(muR+alphaR[subj[r]+1]+betaR[item[r]+1]),muS+alphaS[subj[r]+1]+betaS[item[r]+1],crit[subj[r]+1,])
        dat[r]=which.max(rmultinom(1,1,p))-1
      }

     
    ret=new("dpsdSim")
    ret@cond=cond
    ret@subj=subj
    ret@item=item
    ret@lag=lag
    ret@resp=dat
    ret@muN=muN
    ret@muS=muS
    ret@muR=muR
    ret@alphaN=alphaN
    ret@alphaS=alphaS
    ret@alphaR=alphaR
    ret@betaN=betaN
    ret@betaS=betaS
    ret@betaR=betaR
    return(ret)
  }


dpsdLogLike=function(R,I,J,K,dat,cond,sub,item,lag,blockN,blockS,blockR,crit)
  {
    l=1:R*0
    .C("logLikeDpsd",as.double(l),as.integer(R),as.integer(I),as.integer(J),as.integer(K),as.integer(dat),as.integer(cond),as.integer(sub),as.integer(item),as.double(lag),as.double(blockN),as.double(blockS),as.double(blockR),as.double(as.vector(t(crit))),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }


dpsdSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.01)
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
nsubN=table(sub[cond==0])
nitemN=table(item[cond==0])
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
blockN=blockS=blockR=matrix(0,nrow=M,ncol=B)
wS=wR=rnorm(RS)
s.crit=array(dim=c(I,7,M))
s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)
b0=matrix(0,3,2)
met=matrix(rep(.01,6),3,2)
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
wN=rtnorm(RN,getPred(blockN[m-1,],sub[cond==0],item[cond==0],lag[cond==0],I,J,RN),rep(1,RN),s.crit[cbind(sub[cond==0]+1,resp[cond==0]+1,m-1)],s.crit[cbind(sub[cond==0]+1,resp[cond==0]+2,m-1)])

wS[notk[cond==1]]=rtnorm(Rnotk,getPred(blockS[m-1,],sub[notk],item[notk],lag[notk],I,J,Rnotk),rep(1,Rnotk),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[cond==1]]=rtnorm(Rnotk,getPred(blockR[m-1,],sub[notk],item[notk],lag[notk],I,J,Rnotk),rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[cond==1] & wR>0
wS[ind]=rnorm(sum(ind),getPred(blockS[m-1,],sub[cond==1][ind],item[cond==1][ind],lag[cond==1][ind],I,J,sum(ind)),1)
ind=isk[cond==1] & wR<0
wS[ind]=rtnorm(sum(ind),getPred(blockS[m-1,],sub[cond==1][ind],item[cond==1][ind],lag[cond==1][ind],I,J,sum(ind)),rep(1,sum(ind)),s.crit[cbind(sub[cond==1][ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[cond==1] & wS>s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),getPred(blockR[m-1,],sub[cond==1][ind],item[cond==1][ind],lag[cond==1][ind],I,J,sum(ind)),1)
ind=isk[cond==1] & wS<s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),getPred(blockR[m-1,],sub[cond==1][ind],item[cond==1][ind],lag[cond==1][ind],I,J,sum(ind)),rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))


#Sample Blocks
tmp=sampleNorm(blockN[m-1,],wN,sub[cond==0],item[cond==0],lag[cond==0],I,J,RN,nsubN,nitemN,10,.01,.01,met[1,1],met[1,2],1,1)
blockN[m,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

tmp=sampleNorm(blockS[m-1,],wS,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,nsubS,nitemS,10,.01,.01,met[2,1],met[2,2],1,1)
blockS[m,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[2]]

tmp=sampleNorm(blockR[m-1,],wR,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,nsubS,nitemS,10,.01,.01,met[3,1],met[3,2],1,1)
blockR[m,]=tmp[[1]]
b0[3,]=b0[3,]+tmp[[2]]

#Sample Criteria
PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rnorm(I*4,0,rep(met.crit,4))
violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,m-1]

likeOld=tapply(dpsdLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
likeProp=tapply(dpsdLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],PropCrit),sub,sum)
accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
b0.crit=b0.crit+accept
s.crit[accept==1,,m]=PropCrit[accept==1,]
s.crit[accept==0,,m]=s.crit[accept==0,,m-1]


#AUTOTUNING
if(m>20 & m<tail(keep,1) & m%%10==0)
  {
    met=met+(b0/m<.3)*matrix(-jump,3,2) +(b0/m>.5)*matrix(jump,3,2)
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
    met[met<jump]=jump
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    print("")
    print("Continuning with following acceptance probabilities:")
    print(b0/m)
    print(b0.crit/m)
    print("If they are <.2 or >.6, adjust jump or keep")
  }

setTxtProgressBar(pb, m)
}
close(pb)

estN=colMeans(blockN[keep,])
estS=colMeans(blockS[keep,])
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
    D0=D0+sum(-2*dpsdLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockS
      [m,],blockR[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*dpsdLogLike(R,I,J,K,resp,cond,sub,item,lag,estN,estS,estR,estCrit))
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
u@estS=estS
u@estR=estR
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockS[keep,]
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
