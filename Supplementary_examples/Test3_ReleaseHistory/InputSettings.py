import numpy as np

new_ens='n'            #Do you want to generate a new ensemble? (y/n)
if new_ens=='y':
    ens=20             #Ensemble size
new_err='n'            #Do you want to generate new random errors? (y/n)
N_iter=5
alpha_geo=1
w=0                    #Relaxation coefficient
inflation='y'          #Do you want to apply the inflation? (y/n)
if inflation=='y':
    rr=1.01            #Inflation coefficient
localize='y'
if localize=='y':
    loc_space='y'      #Localization in space? (y/n)
    loc_time='y'       #Localization in time? (y/n)
    iter_loc='n'       #Do iterative localization? (y/n)
space_transform='y'    #Do you want to work in transformed space? (y/n)
ref_solution='y'       #Do you have the reference solution? (y/n)
last_iter_pred=0       #Get predictions for the final parameter estimation: 0(no), 1(ensemble mean), 2(all realizations)

def Func_ens(par, ens):
#function for the generation of the initial ensemble
    from Tools import EnsembleGenerator
    N_par=par.shape[0]

    # Ensemble=EnsembleGenerator.PdfGammaNpeaks(0, 0.05, 2, 10, 2, 20, 5, 15, par[:,3]-50, N_par, ens, 3)

    # Ensemble=EnsembleGenerator.ConstantRandom(0.2,1,N_par,ens)+ \
    #               EnsembleGenerator.Random(0.001,0.03,N_par,ens)
     
    Ensemble=EnsembleGenerator.PdfNormal(0, 0.05, 80, 180, 10, 30, 10, 35, par[:,3], N_par, ens)                          
    
    return Ensemble
    
def Func_err(Obs, ens):
    from Tools import ErrorGenerator
    N_obs=Obs.shape[0]
    perc=1
    var_err_lim=1e-10
    eps,R=ErrorGenerator.PercNormalError(Obs, perc, var_err_lim, N_obs, ens)

    return (eps,R)
    
def forward_transf(xx):
    from Tools import Transformation as T
    (Xmin,Xmax)=(1e-8,30.0)
    xx[xx==Xmin]=Xmin+1e-8
    xx[xx==Xmax]=Xmax-1e-8
    xx=T.LogLim_forward(xx,Xmin,Xmax)
    return xx

def backward_transf(xx):
    from Tools import Transformation as T
    (Xmin,Xmax)=(1e-8,30.0)
    xx[xx>5]=5
    xx=T.LogLim_backward(xx,Xmin,Xmax)
    return xx

def localization(ens_par,par,obs,loc_space,loc_time,iter_loc):
    a_space=300   #correlation length in space
    a_time=75    #correlation length in time
    time_par=par[:,3]
    pos_obs=obs[:,0:3]
    time_obs=obs[:,3]
    if iter_loc=='n':
        pos_par=par[:,0:3]
    else:
        ens_m=ens_par.mean(1)
        pos_par=np.tile(ens_m[0:3],(time_par.shape[0],1))
    from Tools import Localization as Loc
    [rho_yy_sp,rho_xy_sp,rho_xx_sp]=Loc.SpaceLocal(a_space,pos_par,pos_obs)
    [rho_yy_tm,rho_xy_tm,rho_xx_tm]=Loc.TimeLocal(a_time,time_par,time_obs)
    if loc_space=='y' and loc_time=='n':
        rho_yy=rho_yy_sp
        rho_xy=rho_xy_sp
        rho_xx=rho_xx_sp
    elif loc_space=='n' and loc_time=='y':
        rho_yy=rho_yy_tm
        rho_xy=rho_xy_tm
        rho_xx=rho_xx_tm        
    elif loc_space=='y' and loc_time=='y':
        rho_yy=rho_yy_sp*rho_yy_tm
        rho_xy=rho_xy_sp*rho_xy_tm
        rho_xx=rho_xx_sp*rho_xx_tm     
    rho_yy[np.isnan(rho_yy)]=1
    rho_xy[np.isnan(rho_xy)]=1
    rho_xx[np.isnan(rho_xx)]=1
    return (rho_yy,rho_xy,rho_xx)

def Metrics_obs(Xprev,pred,True_par,obs):
    from Tools import Metrics as m
    par=Xprev
    RMSE_obs=m.RMSE(obs, np.atleast_2d(pred.mean(1)).T)
    AES=m.AES(par)
    metrics_dict ={}
    for variable in ['RMSE_obs','AES']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict
    metrics_dict['rmse_obs']=eval('rmse_obs')
    metrics_dict['rmse_par']=eval('rmse_par')
    return metrics_dict

def Metrics_obs_par(Xprev,pred,True_par,obs):
    from Tools import Metrics as m
    RMSE_obs=[m.RMSE(obs.flatten(), pred.mean(1))]
    RMSE_par=[m.RMSE(True_par, Xprev.mean(1))]
    AES=[m.AES(Xprev)]
    NSE_par=[m.NSE(True_par, Xprev.mean(1))]
    metrics_dict ={}
    for variable in ['RMSE_obs', 'RMSE_par', 'AES', 'NSE_par']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict
