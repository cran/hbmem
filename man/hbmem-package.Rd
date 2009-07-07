\name{hbmem-package}
\alias{hbmem-package}
\alias{hbmem}
\docType{package}
\title{Hierarchical Models of Recognition Memory}
\description{
Contains functions for fitting hierarchial versions of EVSD, UVSD, and DPSD
to recognition memory confidence-ratings data.
}
\details{
\tabular{ll}{
Package: \tab hbmem\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-07-07\cr
License: \tab LGPL\cr
LazyLoad: \tab yes\cr
}
}
\author{Michael S. Pratte <prattems@gmail.com>}

\references{Morey, Pratte, and Rouder (2008); Pratte, Rouder, and Morey (2009)}
\keyword{ package }
\seealso{'uvsdSample' to fit hierarchical UVSD model, 'uvsdSim' to
  simulate data from the hierarchical UVSD model, dpsdSample to fit the
  hierarchial DPSD model, 'dpsdSim' to simulate data from the
  hierarchial DPSD model, dataset 'prm09' for some real data to fit.}
\examples{
#In this example data are simulated from EVSD
#They are then fit by both UVSD and DPSD

library(hbmem)
sim=uvsdSim(s2aS2=0,s2bS2=0) #Simulate data from hierarchical EVSD
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","lag","resp")

M=500 #Set low for speed
evsd=uvsdSample(dat,M=M,equalVar=TRUE) #Fit EVSD
uvsd=uvsdSample(dat,M=M,freeSig2=FALSE) #Fit UVSD w/1 Sigma2
dpsd=dpsdSample(dat,M=M) #Fit DPSD 

#Look at available information
slotNames(uvsd)
slotNames(dpsd)

#Compare DIC; smaller is better
evsd@DIC
uvsd@DIC
dpsd@DIC

#PLOTS FROM UVSD FIT
par(mfrow=c(3,2),pch=19,pty='s')
#Make sure chains look OK
plot(uvsd@blockN[,uvsd@mu],t='l',xlab="Iteration",ylab="Mu-N")
abline(h=sim@muN,col="blue")
plot(uvsd@blockS[,uvsd@mu],t='l',xlab="Iteration",ylab="Mu-S")
abline(h=sim@muS,col="blue")

#Estimates of Alpha as function of true values
plot(uvsd@estN[uvsd@alpha]~sim@alphaN,xlab="True
Alpha-N",ylab="Est. Alpha-N");abline(0,1,col="blue")
plot(uvsd@estS[uvsd@alpha]~sim@alphaS,xlab="True
Alpha-S",ylab="Est. Alpha-S");abline(0,1,col="blue")
#Estimates of Beta as function of true values
plot(uvsd@estN[uvsd@beta]~sim@betaN,xlab="True
Beta-N",ylab="Est. Beta-N");abline(0,1,col="blue")
plot(uvsd@estS[uvsd@beta]~sim@betaS,xlab="True
Beta-S",ylab="Est. Beta-S");abline(0,1,col="blue")

###Look at Sigma2 and Recollection from UVSD and DPSD###
par(mfrow=c(2,3),pch=19,pty='s')
plot(sqrt(exp(uvsd@blockS2[,uvsd@mu])),
t='l',ylab="Sigma")
abline(h=1,col="blue")
hist(uvsd@blockS2[,uvsd@s2alpha],main="Participant Effect")
hist(uvsd@blockS2[,uvsd@s2beta],main="Item Effect")

plot(pnorm(dpsd@blockR[,dpsd@mu]),
t='l',ylim=c(0,.2),ylab="P(Recollection)")
abline(h=0,col="blue")
hist(dpsd@blockR[,dpsd@s2alpha],main="Participant Effect")
hist(dpsd@blockR[,dpsd@s2beta],main="Item Effect")


#See what DPSD does with EVSD effects
par(mfrow=c(2,2))
plot(dpsd@estN[dpsd@alpha]~sim@alphaS,xlab="True
Alpha-S",ylab="DPSD Alpha-S");abline(0,1,col="blue")
plot(dpsd@estR[dpsd@alpha]~sim@alphaS,xlab="True
Alpha-S",ylab="DPSD. Alpha-R");abline(0,1,col="blue")
plot(dpsd@estN[dpsd@beta]~sim@betaS,xlab="True
Beta-S",ylab="DPSD Beta-S");abline(0,1,col="blue")
plot(dpsd@estR[dpsd@beta]~sim@betaS,xlab="True
Beta-S",ylab="DPSD. Beta-R");abline(0,1,col="blue")
}