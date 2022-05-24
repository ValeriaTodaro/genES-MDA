import numpy as np
import math, os
import flopy.modflow as mf
import flopy.utils as fu

Obs_file=np.loadtxt('Obs.txt', dtype=float)
True_par_file=np.loadtxt('Par.txt', dtype=float)

os.chdir('Model') 

modelname='Ayvaz_flow'
workspace='ModelFiles'

mfMod=mf.Modflow.load(f'{modelname}.nam', model_ws=workspace,exe_name='mf2005dbl.exe')

lpf=mfMod.get_package('lpf')
dis=mfMod.get_package('dis')
wel=mfMod.get_package('wel')
# AA=mfMod.wel.stress_period_data.to_array(kper=0)
delr=dis.delr.array[0]
delc=dis.delc.array[0]
top=dis.top.array[0,0]
botm=dis.botm.array[-1,-1,-1]
nlay=dis.nlay
nrow=dis.nrow
delv=(top-botm)/nlay

os.chdir('..')

def write_input(par):
    kfield_ind=np.loadtxt('Kfield_indices.txt')
    Kpar=np.copy(kfield_ind)
    
    for ind in np.arange(1,len(par)+1,1):
        Kpar[Kpar==ind]=par[ind-1]
        mfMod.remove_package("lpf")
        lpf_par=mf.ModflowLpf(mfMod, laytyp=0, chani=-1, layvka=1, hk=Kpar,vka=1, hdry=-888.0)
        lpf_par.write_file()  #write file
    
def run():
    try:
        os.remove(os.path.join(workspace,'MT3D001.UCN'))
    except:
        pass
    mfMod.run_model()
   
def read_output():
    N_obs=Obs_file.shape[0]
    x_obs=Obs_file[:,0]
    y_obs=Obs_file[:,1]
    z_obs=Obs_file[:,2]
    obs_col=np.zeros((N_obs),int)
    obs_raw=np.zeros((N_obs),int)
    for ll in range(0,N_obs):
        obs_col[ll]=math.ceil(x_obs[ll]/delr)
        obs_raw[ll]=nrow-math.floor(y_obs[ll]/delc)
    headobj=fu.binaryfile.HeadFile(os.path.join(workspace,modelname+'.hds'))
    #times_mod = np.array(headobj.get_times())
    # conc = concobj.get_data(totim=times[150])
    head = headobj.get_alldata()
    headobj.close()
    H=np.zeros((N_obs))
    for i in range(0,N_obs):
        H[i]=head[0,0,obs_raw[i]-1,obs_col[i]-1]
    return H