import numpy as np
import math, os
import flopy.modflow as mf
import flopy.utils as fu

Obs_file=np.loadtxt('Obs.txt', dtype=float)
True_par_file=np.loadtxt('Par.txt', dtype=float)

os.chdir('Model') 

modelname='3Dexample'
workspace='ModelFiles'

mfMod=mf.Modflow.load(f'{modelname}.nam', model_ws=workspace,exe_name='mf2005dbl.exe')

lpf=mfMod.get_package('lpf')
dis=mfMod.get_package('dis')
wel=mfMod.get_package('wel')

delr=dis.delr.array[0]
delc=dis.delc.array[0]

top=dis.top.array[0,0]
botm=dis.botm.array[-1,-1,-1]
nlay=dis.nlay
nrow=dis.nrow
ncol=dis.ncol
delv=(top-botm)/nlay

X=True_par_file[:,0]
Y=True_par_file[:,1]
Z=True_par_file[:,2]
var_val=True_par_file[:,4]

x_vals=np.sort(np.unique(X))
y_vals=np.sort(np.unique(Y))
z_vals=np.sort(np.unique(Z))

z_idx = np.searchsorted(z_vals, Z)[::-1]
x_idx = np.searchsorted(x_vals, X)
y_idx = np.searchsorted(y_vals, Y)[::-1]

os.chdir('..')

def write_input(par):

    kfield = np.full((nlay, nrow, ncol), np.nan)
    kfield[z_idx, y_idx, x_idx,] = par
    mfMod.remove_package("lpf")
    lpf_par=mf.ModflowLpf(mfMod, laytyp=0, chani=-1, layvka=1, hk=kfield,vka=1, hdry=-888.0)
    lpf_par.write_file()  #write file
    
def run():
    #mfMod.run_model()
    mfMod.run_model(silent=True)
   
def read_output():
    N_obs=Obs_file.shape[0]
    x_obs=Obs_file[:,0]
    y_obs=Obs_file[:,1]
    z_obs=Obs_file[:,2]
    obs_col=np.zeros((N_obs),int)
    obs_raw=np.zeros((N_obs),int)
    obs_lay=np.zeros((N_obs),int)
    for ll in range(0,N_obs):
        obs_col[ll]=math.ceil(x_obs[ll]/delr)
        obs_raw[ll]=nrow-math.floor(y_obs[ll]/delc)
        obs_lay[ll]=nlay-math.floor(z_obs[ll]/delv)
    headobj=fu.binaryfile.HeadFile(os.path.join(workspace,modelname+'.hds'))
    head = headobj.get_alldata()
    headobj.close()
    H=np.zeros((N_obs))
    for i in range(0,N_obs):
        H[i]=head[0,obs_lay[i]-1,obs_raw[i]-1,obs_col[i]-1]
    return H