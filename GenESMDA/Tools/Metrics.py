import math
import numpy as np

def RMSE(actual, predicted):
    # root mean square error
    N=actual.shape[0]
    rmse=math.sqrt(np.sum((predicted-actual)**2)/N)
    return rmse

def AES(ensemble):
    # average ensemble spread
    N=ensemble.shape[0]
    aes=math.sqrt(np.sum(np.var(ensemble,axis=1))/N)
    return aes

def NSE(actual, predicted):
    # Nash-Sutcliffe efficiency criterion
    nse=(1-(np.sum((predicted-actual)**2))/\
         (np.sum((actual-np.mean(actual))**2)))*100
    return nse

def RSS(actual, predicted):
    # residual sum of squares
    rss=sum((predicted-actual)**2)
    return rss

def spatial_distance(actual, prediction):
     dd=math.sqrt(np.sum(prediction-actual)**2)
     return dd
 