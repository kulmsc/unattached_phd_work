import numpy as np
import gzip
import os
import sys
import pdb
from scipy.stats import norm
from scipy.optimize import minimize

method=sys.argv[1]
inputLambda=float(sys.argv[2])
dd=sys.argv[3]
pLim=float(sys.argv[4])

ssName=dd+"/specific_ss"
#pLim=0.05

def normRead(fileName,delimiter):
    totalData=[]
    with open(fileName,"r") as f:
        for line in f.read().splitlines():
            totalData.append(line.split())
    totalData=np.array(totalData)
    return(totalData)


def justRead(fileName,delimiter):
    rList=[]
    with open(fileName,"r") as f:
        for line in f.read().splitlines():
            rList.append(line.split(delimiter))
    return(np.array(rList))


def likelihood(beta,betaHat,se,alpha):
    density=norm.pdf((betaHat-beta)/se)
    lambdaVal=norm.ppf(1-alpha/2)*se
    cumu1=norm.cdf(beta/se - lambdaVal/se)
    cumu2=norm.cdf(-beta/se - lambdaVal/se)
    if abs(beta)>lambdaVal:
        indi=1
    else:
        indi=0
    like=((density/se)/(cumu1+cumu2))*indi
    return(like)


def lasso(beta,lambdaVal):
    newBeta=np.sign(beta)*abs(abs(beta)-lambdaVal)*(1 if abs(beta)>lambdaVal else 0)
    return(newBeta)


ss=normRead(ssName,'\t')
ss=ss[ss[:,5] != "Inf", :]
ss=ss[ss[:,5] != "0", :]
#pdb.set_trace()


ssText=ss[:,2:5]
ssNums=ss[:,(0,1,5,6,7,8)]
ssNums=ssNums.astype(float)
goodBetaHolder=np.zeros(ss.shape[0])


for i in range(ss.shape[0]):
    ssTextRow=ssText[i,:]
    ssNumsRow=ssNums[i,:]
    zVal=abs(norm.ppf(ssNumsRow[4]))
    betaSE=ssNumsRow[3]/zVal

    if method=="like":
        y=minimize(likelihood,1,args=(ssNumsRow[3],betaSE,pLim),method='nelder-mead')
        betaAns=y.x
    elif method=="lasso":
        betaAns=lasso(ssNumsRow[3],inputLambda)
    goodBetaHolder[i]=betaAns


goodBetaHolder.astype("str")
ss[:,6]=goodBetaHolder

np.savetxt(dd+'/summStat', ss, fmt='%s', delimiter='\t')
