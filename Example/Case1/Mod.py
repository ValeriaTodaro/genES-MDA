import numpy as np
import os
import win32com.client
from pydsstools.heclib.dss import HecDss

projectname='3ReachUnsteady'

def write_input(par):
    col=10
    row=int(np.floor(par.size/10))
    par_sq=par[0:row*col]
    par_sq=par_sq.reshape(row,col)

    with open(projectname+".u01","r") as f:
        text=f.readlines()
    keyword="Butte Cr."
    i=0
    while i < len(text):
        if keyword in text[i]:
            #Htext[i]=keyword+"= "+str(par.size)+"\n"
            i=i+3
            for j in range(0,row):
                text[i]=''.join(map(' {:7.2f}'.format, par_sq[j,:]))+"\n"
                i=i+1
            text[i]=''.join(map(' {:7.2f}'.format, par[row*col:]))+"\n"
        i=i+1    
        
    with open(projectname+".u01","w") as f:
        for item in text:
            f.write(item)

def run():
    HEC = win32com.client.Dispatch('RAS507.HECRASCONTROLLER')
    RASProject = os.path.join(os.getcwd(),projectname+'.prj')
    HEC.Project_Open(RASProject)
    HEC.Compute_HideComputationWindow()
    HEC.Compute_CurrentPlan(None, None, True)
    HEC.QuitRas()
    
def read_output():
    dss_file = projectname+".dss"
    pathname = "/FALL RIVER LOWER REACH/9.5/STAGE/01SEP2014/1HOUR/EXISTING/"
    startDate = "02SEP1998 24:00:00"
    endDate = "05SEP1998 24:00:00"
    
    with HecDss.Open(dss_file) as fid:
        ts = fid.read_ts(pathname,window=(startDate,endDate),trim_missing=True)
        values = ts.values
    return values
    
