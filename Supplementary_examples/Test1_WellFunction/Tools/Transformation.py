import numpy as np

def Log_forward(xx):
    X_t=np.log(xx)
    return X_t
def Log_backward(xx):
    X_bt=np.exp(xx)
    return X_bt    

def LogLim_forward(xx,Xmin,Xmax):
    X_t=np.log((xx-Xmin)/(Xmax-xx))
    return X_t
def LogLim_backward(xx,Xmin,Xmax):
    X_bt=(np.exp(xx)*Xmax+Xmin)/(1+np.exp(xx))
    return X_bt
    
def SquareRoot_forward(xx):   
    X_t=xx**(1/2)
    return X_t
def SquareRoot_backward(xx):   
    X_bt=xx**(2)
    return X_bt
      
def SquareRootLim_forward(xx,Xmin,Xmax):   
    X_t=((xx-Xmin)/(Xmax-xx))**(1/2)
    return X_t
def SquareRootLim_backward(xx,Xmin,Xmax):  
    X_t=(Xmax-Xmin)*xx**2/(1+xx**2)+Xmin
    return X_t
    
#A=np.array([1,2,3])
#X=Transform.Log(A)
#XX=BackTransform.Log(X)
#print(XX)