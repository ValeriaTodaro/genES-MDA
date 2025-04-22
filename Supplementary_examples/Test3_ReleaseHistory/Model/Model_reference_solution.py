import numpy as np
import matplotlib.pyplot as plt

import flopy, math, os
import flopy.modflow as mf
import flopy.mt3d as mt
import flopy.utils as fu

modelname='ReleaseHistory_flow'
workspace='ModelFiles'
mfMod=mf.Modflow(modelname=modelname,model_ws=workspace,exe_name='mf2005dbl.exe')

os.chdir('..')
Obs_file=np.loadtxt('Obs.txt', dtype=float)
Par_file=np.loadtxt('Par.txt', dtype=float)
os.chdir('Model')

#DIS file
Lx=300
Ly=100
ztop=10
zbot=0.
nlay=1
nrow=50
ncol=150
delr=Lx/ncol
delc=Ly/nrow
delv=(ztop-zbot)/nlay
lay_dis=np.linspace(ztop, zbot, nlay+1)

#TemporalDiscretization
nper=101
perlen=np.array(np.ones(nper-1)*3)
perlen=np.insert(perlen,0,0)
nstp=2
steady=np.zeros(nper,dtype=bool)
steady[0]=True
#UnitySystem
itmuni=4   #time units (1=seconds, 4=days)
lenuni=2   #lenght units (2=meters, 3=centimeters)
dis=flopy.modflow.ModflowDis(mfMod,nlay,nrow,ncol,delr=delr,delc=delc,
                              top=ztop,botm=lay_dis[1:],
                              nper=nper, perlen=perlen, nstp=nstp, steady=steady,
                              itmuni=itmuni, lenuni=lenuni)

#LPF file
kfield=np.loadtxt('HK_field_md.txt')
kfield=kfield

lpf=mf.ModflowLpf(mfMod, laytyp=0, chani=-1, layvka=1, hk=kfield, vka=1, hdry=-888.0)

#BAS file
ibound=np.ones((nlay,nrow,ncol),dtype=np.int32)
ibound[:,:,0]=-1
ibound[:,:,-1]=-1
strt=ztop
bas=flopy.modflow.ModflowBas(mfMod, ibound=ibound,strt=strt)

#CHD file
chd_dataU=[]
chd_dataD=[]
up_cond=24
down_cond=21

rows = np.arange(nrow)
chd_dataU = [(0, row, 0, up_cond, up_cond) for row in rows]
chd_dataD = [(0, row, ncol - 1, down_cond, down_cond) for row in rows]

chd_data=chd_dataU+chd_dataD
stress_period_data = {0:chd_data}
chd=mf.ModflowChd(mfMod, stress_period_data=stress_period_data)

#Well file
#Par_file=np.loadtxt('Par.txt')
source_col=math.ceil(Par_file[0,0]/delc)
source_row=nrow-math.floor(Par_file[0,1]/delr)
source_lay=nlay-math.floor(Par_file[0,2]/delv)

t=Par_file[:,3]
rel=1/1000 #convert l/d to m3/d
sp_data={}
sp_data[0]=[source_lay-1,source_row-1,source_col-1,0]
for i in range(1,t.shape[0]+1):
    sp_data[i]=[source_lay-1,source_row-1,source_col-1,rel]
wel1=mf.ModflowWel(mfMod,stress_period_data=sp_data)

#LMT Linkage with MT3DMS for multi-species mass transport modeling
lmt = mf.ModflowLmt(mfMod, output_file_name='mt3d_link.ftl')

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

#Results
import flopy.utils.binaryfile as bf
headobj=bf.HeadFile(os.path.join(workspace,modelname+'.hds'))
times = headobj.get_times()
head = headobj.get_data(totim=times[-1])

# plot
levels = np.arange(0, 100, 0.1)
extent = (0, Lx, 0, Ly)
# Make the plots
fig = plt.figure(figsize=(18,8))
plt.subplot(1, 1, 1, aspect='equal')
plt.title('Head distribution (m)')
plt.imshow(head[0, :, :], extent=extent, cmap='YlGnBu')
plt.colorbar()
contours = plt.contour(np.flipud(head[0, :, :]), levels=levels, extent=extent, zorder=0.1)
plt.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
plt.show()

plt.rcParams.update({'font.size': 20}) 
fig = plt.figure(figsize=(18,5))
ax=plt.imshow(kfield, extent=extent,cmap='plasma_r')
plt.scatter(Par_file[0,0],Par_file[0,1],marker='^',c='black',s=120, label='Well')
plt.scatter(Obs_file[:,0],Obs_file[:,1],c='white', edgecolors='black', s=100, label='Obs. location')
cbar = fig.colorbar(ax, label='HK', fraction=0.05, pad=0.02)
plt.legend(loc='upper left', borderpad=0.25)
plt.show()

#MT3D-USGS
namemt3d='ReleaseHistory_trans'
mtMod=mt.Mt3dms(modelname=namemt3d, model_ws=workspace,version='mt3d-usgs',
                exe_name='mt3d-usgs_1.1.0_64.exe', modflowmodel=mfMod)
#BTN file
nprs=101
timprs=np.arange(0, 300+3, 3)
icbund=np.ones((nlay,nrow,ncol),dtype=np.int32)
btn=flopy.mt3d.Mt3dBtn(mtMod, sconc=0.0, prsity=0.3, thkmin=0.01, munit='MG', nprs=nprs, timprs=timprs, icbund=icbund,mxstrn=1000)

#ADV file
mixelm = -1 #Third-order TVD scheme (ULTIMATE)
percel = 1 #Courant number PERCEL is also a stability constraint
adv = flopy.mt3d.Mt3dAdv(mtMod, mixelm=mixelm, percel=percel)

#GCG file
mxiter = 1 #Maximum number of outer iterations
iter1 = 50 #Maximum number of inner iterations
isolve = 1 #Preconditioner = 1 Jacobi
gcg = mt.Mt3dGcg(mtMod, mxiter=mxiter, iter1=iter1, isolve=isolve)

#DSP file
al = 1 #longitudinal dispersivity
dmcoef = 0 #effective molecular diffusion coefficient
trpt = 1/10 #ratio of the horizontal transverse dispersivity to the longitudinal dispersivity
trpv = 1/10 #ratio of the vertical transverse dispersivity to the longitudinal dispersivity
dsp = mt.Mt3dDsp(mtMod, al=al, dmcoef=dmcoef, trpt=trpt, trpv=trpv)

#SSM file
itype = -1 #20

S_loc=[source_lay,source_row,source_col]

S_flux=Par_file[:,-1] #unit in par file is mg/l
S_flux=S_flux #convert mg/l to mg/m3

dtype = np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('cssm', '<f4'), ('itype', '<i8')])
ssm_data={}
ssm=[]
for i in range(0,nper,1):
    ssm=[]    
    ssm.append([S_loc[0]-1]+[S_loc[1]-1]+[S_loc[2]-1]+[S_flux[i],itype])
    ssm_data[i]=ssm
ssm = flopy.mt3d.Mt3dSsm(mtMod, stress_period_data=ssm_data,dtype=dtype)

#Write model input   sp_data[0]=[source_lay-1,source_row-1,source_col-1,0]
mtMod.write_input()

#Try to delete the output files, to prevent use of older files
try:
    os.remove(os.path.join(workspace,'MT3D001.UCN'))
except:
    pass

#Run the model
mtMod.run_model()

#Results
concobj=fu.UcnFile(os.path.join(workspace,'MT3D001.UCN'))
times = concobj.get_times()
conc_t = concobj.get_data(totim=times[5])

# observations
#Obs_file=np.loadtxt('Obs.txt', dtype=float)
N_obs=Obs_file.shape[0]
x_obs=Obs_file[:,0]
y_obs=Obs_file[:,1]
z_obs=Obs_file[:,2]
time_obs=Obs_file[:,3]
obs_col=np.zeros((N_obs),int)
obs_row=np.zeros((N_obs),int)
obs_lay=np.zeros((N_obs),int)

for ll in range(0,N_obs):
    obs_col[ll]=math.ceil(x_obs[ll]/delr)
    obs_row[ll]=nrow-math.floor(y_obs[ll]/delc)
    obs_lay[ll]=nlay-math.floor(z_obs[ll]/delv)
times_mod = np.array(concobj.get_times())
conc = concobj.get_alldata()
concobj.close()
# t=conc.shape[0]
C=np.zeros((N_obs))
tt=[]
for i in range(0,N_obs):
    tt+=np.where(times_mod==time_obs[i])
    C[i]=conc[tt[i][0],obs_lay[i]-1,obs_row[i]-1,obs_col[i]-1]
concobj.close()
Obs=np.copy(Obs_file)
Obs[:,4]=C

os.chdir('..')
with open('Obs.txt','wb') as file:
    np.savetxt(file, Obs, fmt='%.1f %.1f %.1f %.1f %.6f', newline='\n') 

    
#Plot
from matplotlib.colors import LinearSegmentedColormap
YlGnBu = plt.cm.YlGnBu
colors = [(1, 1, 1), *YlGnBu(np.linspace(0, 1, 700))]  # Adding white at the start
cmap = LinearSegmentedColormap.from_list("white_YlGnBu", colors, N=256)

plt.rcParams.update({'font.size': 20}) 
fig = plt.figure(figsize=(18,5))
ax=plt.imshow(conc[68, 0,:, :], extent=extent,cmap=cmap,vmin=0,vmax=0.12)
plt.scatter(Par_file[0,0],Par_file[0,1],marker='^',c='black',s=120, label='Well')
plt.scatter(Obs_file[:,0],Obs_file[:,1],c='white', edgecolors='black', s=100, label='Obs. location')
cbar = fig.colorbar(ax, label='Concentration', fraction=0.05, pad=0.02)
plt.legend(loc='upper left', borderpad=0.25)
plt.show()

x, y, z, time = Obs[:, 0], Obs[:, 1], Obs[:, 2], Obs[:, 3]
unique_positions = np.unique(np.column_stack((x, y, z)), axis=0)
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
axes = axes.flatten()  
num_plots = len(unique_positions)
for i in range(num_plots):
    ux, uy, uz = unique_positions[i]
    mask = (x == ux) & (y == uy) & (z == uz)   
    axes[i].plot(time[mask], Obs[mask,4], '-o', markersize=3, label=f'({ux}, {uy})')
    axes[i].set_ylabel('Concentration (mg/l)')
    axes[i].legend()
    axes[i].grid()
for ax in axes:
    ax.set_xlabel('Time')
plt.tight_layout()
plt.show()


