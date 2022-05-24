import numpy as np

def Random(Min,Max,N_par,ens):
    X=np.random.uniform(low=Min, high=Max, size=(N_par,ens))
    return X
def PdfGamma(aMin,aMax,kMin,kMax,nMax,nMin,bMin,bMax,x,N_par,ens):
    from scipy.stats import gamma
    a=np.random.uniform(low=aMin, high=aMax, size=(ens))  #first value
    k=np.random.uniform(low=kMin, high=kMax, size=(ens))  #scale parameter
    n=np.random.uniform(low=nMin, high=nMax, size=(ens))  #shape parameter
    b=np.random.uniform(low=bMin, high=bMax, size=(ens))  #volume
    X=np.zeros((N_par,ens))
    for i in range(0,ens):
        X[:,i]=gamma.pdf(x, n[i], loc=0, scale=k[i])*b[i]+a[i]
    return X
def PdfGammaNpeaks(aMin,aMax,kMin,kMax,nMax,nMin,
                   bMin,bMax,x,N_par,ens,N_peaks):
    from scipy.stats import gamma
    a=np.zeros((ens,N_peaks))
    k=np.zeros((ens,N_peaks))
    n=np.zeros((ens,N_peaks))
    b=np.zeros((ens,N_peaks))
    for i in range(0,N_peaks):
        a[:,i]=np.random.uniform(low=aMin, high=aMax, size=(ens))  #first value
        k[:,i]=np.random.uniform(low=kMin, high=kMax, size=(ens))  #scale parameter
        n[:,i]=np.random.uniform(low=nMin, high=nMax, size=(ens))  #shape parameter
        b[:,i]=np.random.uniform(low=bMin, high=bMax, size=(ens))  #volume
        #nMin=nMax/1.3
        #nMax=nMax+nMax/1.90
    X=np.zeros((N_par,ens))
    for j in range(0,ens):
        for i in range(0,N_peaks):
            X[:,j]+=gamma.pdf(x, n[j,i], loc=0, scale=k[j,i])*\
                    b[j,i]+a[j,i]
    return X
def PdfNormal(aMin,aMax,muMin,muMax,sigmaMin,sigmaMax,
              bMin,bMax,x,N_par,ens):
    from scipy.stats import norm
    a=np.random.uniform(low=aMin, high=aMax, size=(ens))  #first value
    mu=np.random.uniform(low=muMin, high=muMax, size=(ens))  #scale parameter
    sigma=np.random.uniform(low=sigmaMin, high=sigmaMax, size=(ens))  #shape parameter
    b=np.random.uniform(low=bMin, high=bMax, size=(ens))  #volume
    X=np.zeros((N_par,ens))
    for i in range(0,ens):
        X[:,i]=norm.pdf(x, mu[i], sigma[i])*b[i]+a[i]
    return X
def ConstantRandom(Min,Max,N_par,ens):       
    values=np.random.uniform(low=Min, high=Max, size=(ens))
    values=np.atleast_2d(values)
    vector=np.ones([N_par])
    vector=np.atleast_2d(vector).T
    X=vector@values
    return X
def ConstantNormal(mu,sigma,N_par,ens):       
    values=np.random.normal(mu, sigma, size=(ens))
    values=np.atleast_2d(values)
    vector=np.ones([N_par])
    vector=np.atleast_2d(vector).T
    X=vector@values
    return X          
#X=CreateEnsemble.Gamma(0.000001,0.05,15,60,5,20,600,900,np.arange(5, 2500, 35),72,5)
#X=RandomEnsemble.GammaNpeaks(0.000001,0.05,20,50,5,20,600,900,np.arange(5, 2500, 35),72,5,2)
#x=np.arange(5, 2500, 35)
#X=PdfNormal(0.000001,0.05,np.quantile(x,0.2),np.quantile(x,0.55),np.quantile(x,0.03),np.quantile(x,0.15),600,900,x,N_par=72,ens=5)

