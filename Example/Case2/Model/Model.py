import numpy as np
import matplotlib.pyplot as plt
import flopy, math, os
import flopy.modflow as mf
import flopy.mt3d as mt
import flopy.utils as fu

modelname='Ayvaz_flow'
workspace='ModelFiles'
mfMod=mf.Modflow(modelname=modelname,model_ws=workspace,exe_name='mf2005dbl.exe')

#DIS file
Lx=2500
Ly=1500
ztop=30
zbot=0.
nlay=1
nrow=15
ncol=25
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
kfield_true={'ind1':34.56,'ind2':17.28,'ind3':8.64,'ind4':25.92,'ind5':60.48}
kfield_ind=np.loadtxt('Kfield_indices.txt')

Kpar=np.copy(kfield_ind)
for ind in np.arange(1,len(kfield_true)+1,1):
    Kpar[Kpar==ind]=kfield_true['ind'+str(ind)]

lpf=mf.ModflowLpf(mfMod, laytyp=0, chani=-1, layvka=1, hk=Kpar,vka=1, hdry=-888.0)

#BAS file
ibound=np.ones((nlay,nrow,ncol),dtype=np.int32)
ibound[0,np.where(kfield_ind == 0)[0], np.where(kfield_ind == 0)[1]]=0
strt=ztop
bas=flopy.modflow.ModflowBas(mfMod, ibound=ibound,strt=strt)

#CHD file
chd_dataU=[]
chd_dataD=[]
up_cond=100
down_cond=80
pos_up_cond=np.array([[1,1,4], [1, 2, 3], [1, 3, 3], [1, 4, 2], [1, 5, 1]])
pos_down_cond=np.array([[1,11,25], [1, 12, 24], [1, 13, 23], [1, 14, 22], [1, 15, 21]])
for c in range(0,pos_up_cond.shape[0]):
    lChd = np.append(pos_up_cond[c,:]-1, (up_cond, up_cond))
    chd_dataU.append(lChd)
for c in range(0,pos_down_cond.shape[0]):
    rChd = np.append(pos_down_cond[c,:]-1, (down_cond, down_cond))
    chd_dataD.append(rChd)
chd_data=chd_dataU+chd_dataD
stress_period_data = {0:chd_data}
chd=mf.ModflowChd(mfMod, stress_period_data=stress_period_data)

#Well file
#rel_input=np.loadtxt('True_input_model.txt', dtype=float)
x_pos=1150
y_pos=1050
z_pos=0
source_col=math.ceil(x_pos/delr)
source_row=nrow-math.floor(y_pos/delc)
source_lay=nlay-math.floor(z_pos/delv)
#t=rel_input[:,2]
#rel=rel_input[:,3]
sp_data={}
sp_data[0]=[source_lay-1,source_row-1,source_col-1,-3000.]
#for i in range(1,t.shape[0]+1):
#    sp_data[i]=[source_lay-1,source_row-1,source_col-1,rel[i-1]]
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

# #plot
levels = np.arange(0, 500, 1)
extent = (0, Lx, 0, Ly)
 # Make the plots
fig = plt.figure(figsize=(18,8))
plt.subplot(1, 1, 1, aspect='equal')
plt.title('Head distribution (cm)')
plt.imshow(head[0, :, :], extent=extent, cmap='YlGnBu', vmin=50, vmax=100)
plt.colorbar()
contours = plt.contour(np.flipud(head[0, :, :]), levels=levels, extent=extent, zorder=10)
plt.clabel(contours, inline=1, fontsize=10, fmt='%d')
plt.show()

fig = plt.figure(figsize=(18,8))
plt.imshow(Kpar)