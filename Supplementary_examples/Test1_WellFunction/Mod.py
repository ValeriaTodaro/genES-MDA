import numpy as np
import os

def write_input(par):
    with open('par.txt','wb') as file:
        np.savetxt(file, par, fmt='%.4e', newline='\n') 
      
def run():
    #Read parameters
    par=np.loadtxt('par.txt', dtype=float)
    T=par[0] # Transmissivity in m²/s
    S=par[1] # Storativity (dimensionless)
    
    #Read observation times
    os.chdir('..')
    obs=np.loadtxt('Obs.txt')
    os.chdir('Model')
    t=obs[:,3]
    
    #Run forward model
    from Forward_model import well_function
    s=well_function(T, S, t)  
    
    #Write response variable at obsevation locations
    with open('pred.txt','wb') as file:
        data = np.column_stack((t, s))
        np.savetxt(file, data, fmt=['%.0f','%.3e'], newline='\n') 
   
def read_output():
    _, pred =np.loadtxt('pred.txt', dtype=float, unpack=True)
    return pred