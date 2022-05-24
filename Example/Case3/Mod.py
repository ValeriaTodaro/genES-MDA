import numpy as np
import math, os
import flopy
import flopy.modflow as mf
import flopy.mt3d as mt
import flopy.utils as fu

Obs_file=np.loadtxt('Obs.txt', dtype=float)
True_par_file=np.loadtxt('Par.txt', dtype=float)

os.chdir('Model')

modelname='Ayvaz_flow'
namemt3d='Ayvaz_trans'

workspace='ModelFiles'

mfMod=mf.Modflow.load(f'{modelname}.nam', model_ws=workspace,exe_name='mf2005dbl.exe')
mtMod=mt.Mt3dms.load(f'{namemt3d}.nam',model_ws=workspace,modflowmodel=mfMod,exe_name='mt3d-usgs_1.1.0_64.exe')
ssm=mfMod.get_package('ssm')
dis=mfMod.get_package('dis')

delr=dis.delr.array[0]
delc=dis.delc.array[0]
top=dis.top.array[0,0]
botm=dis.botm.array[-1,-1,-1]
nlay=dis.nlay
nrow=dis.nrow
delv=(top-botm)/nlay

os.chdir('..')
def write_input(par):
    itype = 15
    S1_loc=[4,4,1]
    S2_loc=[4,7,1]
    S1_flux=np.concatenate((par[0:4],np.array([0.,0.,0.,0.,0.,0.])))
    S2_flux=np.concatenate((par[4:],np.array([0.,0.,0.,0.,0.,0.])))
    dtype = np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('css', '<f4'), ('itype', '<i8')])
    ssm_data={}
    ssm_par=[]
    nper=len(S1_flux)
    for n in range(0,nper,1):
        ssm_par=[]
        ssm_par.append([S1_loc[2]-1]+[S1_loc[0]-1]+[S1_loc[1]-1]+[S1_flux[n],itype])
        ssm_par.append([S2_loc[2]-1]+[S2_loc[0]-1]+[S2_loc[1]-1]+[S2_flux[n],itype])
        ssm_data[n]=ssm_par
    mtMod.remove_package("ssm")
    ssm_par = flopy.mt3d.Mt3dSsm(mtMod, mxss=15, stress_period_data=ssm_data,dtype=dtype)
    ssm_par.write_file()
    
def run():
    try:
        os.remove(os.path.join(workspace,'MT3D001.UCN'))
    except:
        pass
    #mfMod.run_model()
    mtMod.run_model() 
      
def read_output():
    N_obs=Obs_file.shape[0]
    x_obs=Obs_file[:,0]
    y_obs=Obs_file[:,1]
    z_obs=Obs_file[:,2]
    time_obs=Obs_file[:,3]
    obs_col=np.zeros((N_obs),int)
    obs_lay=np.zeros((N_obs),int)
    obs_row=np.zeros((N_obs),int)
    for ll in range(0,N_obs):
        obs_col[ll]=math.ceil(x_obs[ll]/delr)
        obs_row[ll]=nrow-math.floor(y_obs[ll]/delc)
        obs_lay[ll]=math.ceil(z_obs[ll]/delv)
    concobj=fu.UcnFile(os.path.join(workspace,'MT3D001.UCN'))
    times_mod = np.array(concobj.get_times())
    conc = concobj.get_alldata()
    concobj.close()
    
    C=np.zeros((N_obs))
    tt=[]
    for ii in range(0,N_obs):
        tt+=np.where(times_mod==time_obs[ii])
        C[ii]=conc[tt[ii],obs_lay[ii]-1,obs_row[ii]-1,obs_col[ii]-1]
    return C