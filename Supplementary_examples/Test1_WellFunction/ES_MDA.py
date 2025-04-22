import numpy as np
import Mod
import os
#input
Default= (input('Do you want to use default settings? (y/n) '))
if Default=='n':
    new_ens=(input('Do you want to generate a new ensemble? (y/n) '))
    if new_ens=='y':
       ens=int(input('Ensemble size: '))
    new_err=(input('Do you want to generate new random errors? (y/n) '))
    N_iter=int(input('Number of iterations: '))
    alpha_geo=float(input('Alpha_geo= '))
    w=float(input('Relaxation coefficient: '))
    inflation=(input('Do you want to apply the inflation? (y/n) '))
    if inflation=='y':
       rr=float(input('Inflation coefficient: ')) 
    localize=(input('Do you want to apply the localization? (y/n) '))
    if localize=='y':
        loc_space=(input('Localization in space? (y/n) '))
        loc_time=(input('Localization in time? (y/n) '))
        iter_loc=(input('Do iterative localization? (y/n) '))
    
    space_transform=(input('Do you want to work in transformed space? (y/n) '))
    ref_solution=(input('Do you have the reference solution? (y/n) '))
    last_iter_pred=int(input('Get predictions for the final parameter estimation: 0(no), 1(ensemble mean), 2(all realizations) '))
else:
    from InputSettings import(new_ens,new_err,N_iter,alpha_geo,w,inflation,localize,space_transform,ref_solution,last_iter_pred)
    if new_ens=='y':
        from InputSettings import(ens)
    if inflation=='y':
        from InputSettings import(rr)    
    if localize=='y':
        from InputSettings import(loc_space,loc_time,iter_loc)
   
Obs_file=np.atleast_2d(np.loadtxt("Obs.txt", dtype=float))
x_obs=Obs_file[:,0]
y_obs=Obs_file[:,1]
z_obs=Obs_file[:,2]
time_obs=Obs_file[:,3]
Obs=np.atleast_2d(Obs_file[:,4]).T
N_obs=Obs.shape[0]

True_par_file=np.atleast_2d(np.loadtxt('Par.txt', dtype=float))
x_par=True_par_file[:,0]
y_par=True_par_file[:,1]
z_par=True_par_file[:,2]
time_par=True_par_file[:,3]
True_par=True_par_file[:,4]

if new_ens=='n':
    X=np.loadtxt('Ens.txt')
    ens=X.shape[1]
elif new_ens=='y':
    from InputSettings import Func_ens
    X=Func_ens(True_par_file, ens)
else:
    print('Error, invalid input')
N_par=X.shape[0]

if new_err=='n':
    R=np.loadtxt('R.txt')
    eps=np.loadtxt('Errors.txt')
elif new_err=='y':
    from InputSettings import Func_err
    eps,R=Func_err(Obs, ens)
else:
    print('Error, invalid input')

al_i = np.ones((N_iter, 1),float)
for i in range(1,N_iter):
    al_i[i]=al_i[i-1]/alpha_geo
sum_al_i=sum(1./al_i)
alpha=al_i*sum_al_i
sum_alpha=sum(1./alpha)

if localize=='y':
    from InputSettings import localization
if localize=='y' and iter_loc=='n':
    (rho_yy,rho_xy,rho_xx)=localization(X,True_par_file[:,0:4],Obs_file[:,0:4],loc_space, loc_time, iter_loc)

r=[]
pred=np.zeros((N_obs,ens))

Xprev=np.copy(X)
for i in range(0,N_iter):
    print('Iteration ' + str(i+1))
    R_corr=alpha[i]*R
    r.append(R_corr[i,i])
    if localize=='y' and iter_loc=='y':
        (rho_yy,rho_xy,rho_xx)=localization(Xprev,True_par_file[:,0:4],Obs_file[:,0:4],loc_space, loc_time, iter_loc)
        
    os.chdir('Model')
    for j in range(0,ens):
        print('Iteration ' + str(i+1) + '- Realization ' +str(j+1))
        Mod.write_input(Xprev[:,j])
        Mod.run()  #run forward model
        pred[:,j]=Mod.read_output()
    os.chdir('..')
          
    if space_transform=='y':
        from InputSettings import forward_transf
        Xprev=forward_transf(Xprev)

    xm=np.atleast_2d(Xprev.mean(1)).T
    ym=np.atleast_2d(pred.mean(1)).T
    Qx=Xprev-xm*np.ones((1,ens))
    Qy=pred-ym*np.ones((1,ens))
    Qxy=Qx@Qy.T/(ens-1)
    Qyy=Qy@Qy.T/(ens-1)
    Qxx=Qx@Qx.T/(ens-1)
    
    if localize=='y':
       Qxy=rho_xy*Qxy
       Qyy=rho_yy*Qyy
       Qxx=rho_xx*Qxx

    Gain=Qxy @ np.linalg.inv(Qyy+R_corr)

    Xnew=Xprev+Gain@(Obs@np.ones((1,ens))+(alpha[i])**(1/2)*eps-pred)
    Xnew=(1-w)*Xnew+w*Xprev
    
    if inflation=='y':
        Xnew_mean=np.atleast_2d(np.mean(Xnew,axis=1)).T
        Xnew=Xnew_mean*np.ones((1,ens))+rr*(Xnew-Xnew_mean*np.ones((1,ens)));
        
    if space_transform=='y':
       from InputSettings import backward_transf
       Xnew=backward_transf(Xnew)
       Xprev=backward_transf(Xprev)     
  
    Xprev=np.copy(Xnew)
    Xp=np.mean(Xprev,axis=1)
    pred_mean=np.mean(pred,axis=1)
    
    if ref_solution=='y':
        from InputSettings import Metrics_obs_par
        metrics_iter=Metrics_obs_par(Xprev,pred,True_par,Obs)
    else:
        from InputSettings import Metrics_obs
        metrics_iter=Metrics_obs(Xprev,pred,Obs)
        
    if i==0:
        metrics_name=list(metrics_iter.keys())
        metrics_dict=metrics_iter.copy()
    else:
        for m in metrics_name:
            metrics_dict[m]=metrics_dict[m]+metrics_iter[m]
    
    if i==0:
        with open('_Xprev_iter0.txt','wb') as f:
            np.savetxt(f, X, fmt='%.6e')
    with open('_Xprev_iter'+str(i+1)+'.txt','wb') as f:
        np.savetxt(f, Xprev, fmt='%.6e', newline='\r\n') 
    with open('_pred_iter'+str(i+1)+'.txt','wb') as f:
        np.savetxt(f, pred, fmt='%.6e', newline='\n')

if  last_iter_pred==1:
    os.chdir('Model')
    Mod.write_input(Xp)
    Mod.run()  #run forward model
    pred_ParMean=Mod.read_output()
    os.chdir('..')
    with open('_pred_last_iter_mean.txt','wb') as f:
        np.savetxt(f, pred_ParMean, fmt='%.6e', newline='\n')
elif last_iter_pred==2:
    os.chdir('Model')
    for j in range(0,ens):
        print('Last_iter_pred ' + '- Realization ' +str(j+1))
        Mod.write_input(Xprev[:,j])
        Mod.run()  #run forward model
        pred[:,j]=Mod.read_output()
    os.chdir('..')   
    with open('_pred_last_iter.txt','wb') as f:
        np.savetxt(f, pred, fmt='%.6e', newline='\n')