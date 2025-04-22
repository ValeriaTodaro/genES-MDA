import numpy as np

def NormalError(var_err,N_obs,ens):
    eps=np.zeros((N_obs,ens))
    R=np.identity(N_obs)*var_err
    for i in range(N_obs):
        eps[i,:]=np.random.normal(0, var_err**(1/2), size=(ens))
    return (eps,R)

def PercNormalError(obs,perc,var_err_lim,
                    N_obs,ens):
    var_err=(perc/100*obs/3)**2
    var_err[var_err < var_err_lim]=var_err_lim 
    R=np.diag(var_err[:,0])
    eps=np.zeros((N_obs,ens))
    for i in range(N_obs):
        eps[i,:]=np.random.normal(0, (var_err[i])**(1/2), size=(ens))
    return (eps,R)

