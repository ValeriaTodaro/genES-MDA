#InputSettings is the module that includes all the input information required to perform the inverse procedure.
#It requires to be modified by the user to set up the best algorithm configuration for the problem at hand.

#Import here the required Python packages and libraries

new_ens='n'            #Do you want to generate a new ensemble? (y/n)
if new_ens=='y':
    ens=35             #Ensemble size
new_err='n'            #Do you want to generate new random errors? (y/n)
N_iter=6               #Number of iterations
alpha_geo=3            #Alpha_geo
w=0                    #Relaxation coefficient
inflation='y'          #Do you want to apply the inflation? (y/n)
if inflation=='y':
    rr=1.01            #Inflation coefficient
localize='y'           #Do you want to apply the localization? (y/n)
if localize=='y':
    loc_space='n'      #Localization in space? (y/n)
    loc_time='y'       #Localization in time? (y/n)
    iter_loc='n'       #Do iterative localization? (y/n)
space_transform='y'    #Do you want to work in transformed space? (y/n)
ref_solution='n'       #Do you have the reference solution? (y/n)
last_iter_pred=0       #Get predictions for the final parameter estimation: 0(no), 1(ensemble mean), 2(all realizations)

'''
def Func_ens(par, ens):
    #Function used to set up the initial ensemble of parameter
    from Tools import EnsembleGenerator  
    return Ensemble
    
def Func_err(Obs, ens):
    #Function used to set up the ensemble of observation errors and its covariance matrix.
    from Tools import ErrorGenerator
    return (eps,R)
    
def forward_transf(xx):
    #Function used to set up the parameter space transformation
    from Tools import Transformation as T
    return xx

def backward_transf(xx):
    #Function used to back-transform the parameters in their physical space
    from Tools import Transformation as T
    return xx

def localization(ens_par,par,obs,loc_space,loc_time,iter_loc):
    #Function used to set up the covariance localization
    from Tools import Localization as Loc    
    a_space=#      #correlation length in space
    a_time=#       #correlation length in time
    time_par=par[:,3]
    pos_obs=obs[:,0:3]
    time_obs=obs[:,3]
    if iter_loc=='n':
        pos_par=par[:,0:3]
    else:
        ens_m=ens_par.mean(1)
        #if localization is performed in the itarive form, define the ensemble rows pertaining the coordinates
        pos_par=#
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

def Metrics_obs(Xprev,pred,obs):
    #Function used to select the metrics based on actual observations and predicted values
    from Tools import Metrics as m
    metric1=#
    metric2=#
    metric3=#
    metrics_dict ={}
    for variable in ['metric1', 'metric2', 'metric3']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict

def Metrics_obs_par(Xprev,pred,True_par,obs):
    #Function used to select the metrics when a reference solution is available
    from Tools import Metrics as m
    metric1=#
    metric2=#
    metric3=#
    metrics_dict ={}
    for variable in ['metric1', 'metric2', 'metric3']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict
'''
