import numpy as np
import os

def write_input(par):
    with open('par.txt','wb') as file:
        np.savetxt(file, par, fmt='%.4e', newline='\n') 
      
def run():
    par=np.loadtxt('par.txt', dtype=float)
    T=par[0] # Transmissivity in m²/s
    S=par[1] # Storativity (dimensionless)
    
    from Forward_model import well_function
    well_function(T, S)     
   
def read_output():
    os.chdir('..')
    Obs_file=np.loadtxt('Obs.txt')
    os.chdir('Model')
    t_obs=Obs_file[:,3]
    t,s=np.loadtxt('pred.txt', dtype=float, unpack=True)
    idx=np.where(np.isin(t, t_obs))[0]
    pred=s[idx]
    return pred