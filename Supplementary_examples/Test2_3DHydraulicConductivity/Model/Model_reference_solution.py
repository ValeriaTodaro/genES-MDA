import numpy as np
import matplotlib.pyplot as plt
import flopy, math, os
import flopy.modflow as mf
import flopy.mt3d as mt
import flopy.utils as fu
import os

modelname='3Dexample'
workspace='ModelFiles'
mfMod=mf.Modflow(modelname=modelname,model_ws=workspace,exe_name='mf2005dbl.exe')

#DIS file
Lx=100
Ly=100
ztop=30
zbot=0.
nlay=3
nrow=10
ncol=10
delr=Lx/ncol
delc=Ly/nrow
delv=(ztop-zbot)/nlay
lay_dis=np.linspace(ztop, zbot, nlay+1)
#TemporalDiscretization
nper=1
perlen=np.array(np.ones(nper-1)*10)
perlen=np.insert(perlen,0,0)
nstp=1
steady=np.zeros(nper,dtype=bool)
steady[0]=True
#UnitySystem
itmuni=4   #time units   (4=year, 1=seconds)
lenuni=2   #lenght units (2=m,  3=centimeters)
dis=flopy.modflow.ModflowDis(mfMod,nlay,nrow,ncol,delr=delr,delc=delc,
                              top=ztop,botm=lay_dis[1:],
                              nper=nper, perlen=perlen, nstp=nstp, steady=steady,
                              itmuni=itmuni, lenuni=lenuni)

#LPF file
os.chdir('..')
par=np.loadtxt('Par.txt')
os.chdir('Model')
X=par[:,0]
Y=par[:,1]
Z=par[:,2]
var_val=par[:,4]

x_vals=np.sort(np.unique(X))
y_vals=np.sort(np.unique(Y))
z_vals=np.sort(np.unique(Z))

kfield_true = np.full((nlay, nrow, ncol), np.nan)

z_idx = np.searchsorted(z_vals, Z)[::-1]
x_idx = np.searchsorted(x_vals, X)
y_idx = np.searchsorted(y_vals, Y)[::-1]

kfield_true[z_idx, y_idx, x_idx,] = var_val

lpf=mf.ModflowLpf(mfMod, laytyp=1, chani=-1, layvka=1, hk=kfield_true,vka=1, hdry=-888.0)

#BAS file
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)  # Active cells (1).
ibound[:,:,0] = -1  # Set inactive cells (0).
ibound[:,:,-1] = -1  # Set inactive cells (0).
# Initial heads equal to top elevations.
strt = np.ones((nlay, nrow, ncol))*ztop
strt[:,:,0] = 29.0
strt[:,:,-1] = 28.5

bas = flopy.modflow.ModflowBas(mfMod, ibound=ibound, strt=strt)

#Well file
pumping_rate=-900.0
x_pos1=35
y_pos1=75
z_pos1=25
source_col=math.ceil(x_pos1/delr)
source_row=nrow-math.floor(y_pos1/delc)
source_lay=nlay-math.floor(z_pos1/delv)
well1 = [source_lay-1,source_row-1,source_col-1,pumping_rate]

#well 2
pumping_rate=-650.0
x_pos2=75
y_pos2=35
z_pos2=5
source_col=math.ceil(x_pos2/delr)
source_row=nrow-math.floor(y_pos2/delc)
source_lay=nlay-math.floor(z_pos2/delv)
well2 = [source_lay-1,source_row-1,source_col-1,pumping_rate]

sp_data={0: [well1, well2]}

wel=mf.ModflowWel(mfMod,stress_period_data=sp_data)

#PCG file
pcg=mf.mfpcg.ModflowPcg(mfMod, mxiter=20, iter1=50, hclose=1e-03, rclose=1e-03, relax=1.0)

#OC file
stress_period_data={}
for kper in range(nper):
    for kstp in range(nstp):
        stress_period_data[(kper, kstp)] = ['save head']
oc=mf.ModflowOc(mfMod, stress_period_data=stress_period_data, compact=True)

#Write input
mfMod.write_input()

#Try to delete the output files, to prevent accidental use of older files
try:
    os.remove(os.path.join(workspace,modelname+'.hds'))
except:
    pass

#Run the model
success,buff=mfMod.run_model()

#write observations
os.chdir('..')
Obs_file=np.loadtxt('Obs.txt', dtype=float)
os.chdir('Model')

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
head_all = headobj.get_alldata()
headobj.close()
H=np.zeros((N_obs))
for i in range(0,N_obs):
    H[i]=head_all[0,obs_lay[i]-1,obs_raw[i]-1,obs_col[i]-1]
Obs_file[:,-1]=H 

os.chdir('..')
with open('Obs.txt','wb') as f:
    np.savetxt(f, Obs_file, fmt='%.2f')
os.chdir('Model')

#Results
import flopy.utils.binaryfile as bf
headobj=bf.HeadFile(os.path.join(workspace,modelname+'.hds'))
times = headobj.get_times()
head = headobj.get_data(totim=times[-1])
 
plt.rcParams.update({'font.size': 20}) 
fig, axes = plt.subplots(nlay, 1, figsize=(18, 8 * nlay))
kfield_true_plot=kfield_true/(365.25) #convert from m/year to m/day
vmin, vmax = kfield_true_plot.min(), kfield_true_plot.max()
for i in range(nlay):
    ax = axes[i] if nlay > 1 else axes  # Adjust for single layer
    extent = (0, Lx, 0, Ly)  # Spatial extent of the plot
    im = ax.imshow(kfield_true_plot[i, :, :], extent=extent, cmap='plasma_r' ,vmin=vmin, vmax=vmax)
    ax.set_title(f'Layer {i + 1}')
    ax.set_title(f'Layer {i + 1}', fontsize=20)
    ax.set_xlabel('x (m)', fontsize=20)
    ax.set_ylabel('y (m)', fontsize=20)
    cbar = fig.colorbar(im, ax=ax, label='HK (m/d)', fraction=0.05, pad=0.02)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((-2, 2))  
    cbar.update_ticks()
plt.tight_layout()
plt.show()  
    
fig, axes = plt.subplots(nlay, 1, figsize=(18, 8 * nlay))
vmin, vmax = head.min(), head.max()
for i in range(nlay):
    ax = axes[i] if nlay > 1 else axes  # Handle single layer case
    extent = (0, Lx, 0, Ly)  # Spatial extent of the plot
    im = ax.imshow(head[i, :, :], extent=extent, cmap='GnBu', vmin=vmin, vmax=vmax)

    for x, y in zip(x_obs, y_obs):
        ax.scatter(x, y, c='red', marker='.', s=100, label='Obs Location')
    ax.set_title(f'Layer {i + 1}', fontsize=20)
    ax.set_xlabel('x (m)', fontsize=20)
    ax.set_ylabel('y (m)', fontsize=20)
    cbar = fig.colorbar(im, ax=ax, label='Hydraulic head (m)', fraction=0.05, pad=0.02)
axes[0].scatter(x_pos1, y_pos1, edgecolors='black', facecolors='none', marker='^', s=300, label='Well1')
axes[2].scatter(x_pos2, y_pos2, edgecolors='black', facecolors='none', marker='^', s=300, label='Well2')
plt.tight_layout()
plt.show()