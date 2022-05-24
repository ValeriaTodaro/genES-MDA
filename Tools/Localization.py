import numpy as np
# import matplotlib.pyplot as plt

def TimeLocal(a,time_par,time_obs):
    N_par=time_par.shape[0]
    N_obs=time_obs.shape[0]
    d_yy=np.zeros((N_obs,N_obs))
    d_xy=np.zeros((N_par,N_obs))
    d_xx=np.zeros((N_par,N_par))
    rho_yy=d_yy;
    rho_xy=d_xy;
    rho_xx=d_xx;              
    for row in range(0,d_yy.shape[0]):
        for col in range(0,d_yy.shape[1]):
           d_yy[row,col]=abs(time_obs[row]-time_obs[col])
           if d_yy[row,col]<=a:
              rho_yy[row,col]=-1/4*(abs(d_yy[row,col])/a)**5+\
                               1/2*(abs(d_yy[row,col])/a)**4+\
                               5/8*(abs(d_yy[row,col])/a)**3-\
                               5/3*(abs(d_yy[row,col])/a)**2+\
                               1
           elif d_yy[row,col]>a  and d_yy[row,col]<=2*a:
              rho_yy[row,col]=1/12*(abs(d_yy[row,col])/a)**5-\
                              1/2*(abs(d_yy[row,col])/a)**4+\
                              5/8*(abs(d_yy[row,col])/a)**3+\
                              5/3*(abs(d_yy[row,col])/a)**2-\
                              5*(abs(d_yy[row,col])/a)+\
                              4-\
                              2/3*a/abs(d_yy[row,col])
           elif d_yy[row,col]>2*a:
              rho_yy[row,col]=0
    for row in range(0,d_xy.shape[0]):
        for col in range(0,d_xy.shape[1]):
           d_xy[row,col]=abs(time_par[row]-time_obs[col])
           if d_xy[row,col]<=a:
              rho_xy[row,col]=-1/4*(abs(d_xy[row,col])/a)**5+\
                               1/2*(abs(d_xy[row,col])/a)**4+\
                               5/8*(abs(d_xy[row,col])/a)**3-\
                               5/3*(abs(d_xy[row,col])/a)**2+\
                               1
           elif d_xy[row,col]>a  and d_xy[row,col]<=2*a:
              rho_xy[row,col]=1/12*(abs(d_xy[row,col])/a)**5-\
                              1/2*(abs(d_xy[row,col])/a)**4+\
                              5/8*(abs(d_xy[row,col])/a)**3+\
                              5/3*(abs(d_xy[row,col])/a)**2-\
                              5*(abs(d_xy[row,col])/a)+\
                              4-\
                              2/3*a/abs(d_xy[row,col])
           elif d_xy[row,col]>2*a:
              rho_xy[row,col]=0  
    for row in range(0,d_xx.shape[0]):
        for col in range(0,d_xx.shape[1]):
           d_xx[row,col]=abs(time_par[row]-time_par[col])
           if d_xx[row,col]<=a:
              rho_xx[row,col]=-1/4*(abs(d_xx[row,col])/a)**5+\
                               1/2*(abs(d_xx[row,col])/a)**4+\
                               5/8*(abs(d_xx[row,col])/a)**3-\
                               5/3*(abs(d_xx[row,col])/a)**2+\
                               1
           elif d_xx[row,col]>a  and d_xx[row,col]<=2*a:
              rho_xx[row,col]=1/12*(abs(d_xx[row,col])/a)**5-\
                              1/2*(abs(d_xx[row,col])/a)**4+\
                              5/8*(abs(d_xx[row,col])/a)**3+\
                              5/3*(abs(d_xx[row,col])/a)**2-\
                              5*(abs(d_xx[row,col])/a)+\
                              4-\
                              2/3*a/abs(d_xx[row,col])
           elif d_xx[row,col]>2*a:
              rho_xx[row,col]=0
    return (rho_yy,rho_xy,rho_xx)

def SpaceLocal(a,pos_par,pos_obs):
    x_par=pos_par[:,0]
    y_par=pos_par[:,1]
    z_par=pos_par[:,2]
    x_obs=pos_obs[:,0]
    y_obs=pos_obs[:,1]
    z_obs=pos_obs[:,2]
    N_par=pos_par.shape[0]
    N_obs=pos_obs.shape[0]
    d_yy=np.zeros((N_obs,N_obs))
    d_xy=np.zeros((N_par,N_obs))
    d_xx=np.zeros((N_par,N_par))
    rho_yy=d_yy;
    rho_xy=d_xy;
    rho_xx=d_xx;       
    for row in range(0,d_yy.shape[0]):
        for col in range(0,d_yy.shape[1]):
           d_yy[row,col]=np.sqrt((x_obs[row]-x_obs[col])**2+\
                                 (y_obs[row]-y_obs[col])**2+\
                                 (z_obs[row]-z_obs[col])**2)
           if d_yy[row,col]<=a:
              rho_yy[row,col]=-1/4*(abs(d_yy[row,col])/a)**5+\
                               1/2*(abs(d_yy[row,col])/a)**4+\
                               5/8*(abs(d_yy[row,col])/a)**3-\
                               5/3*(abs(d_yy[row,col])/a)**2+\
                               1
           elif d_yy[row,col]>a  and d_yy[row,col]<=2*a:
              rho_yy[row,col]=1/12*(abs(d_yy[row,col])/a)**5-\
                              1/2*(abs(d_yy[row,col])/a)**4+\
                              5/8*(abs(d_yy[row,col])/a)**3+\
                              5/3*(abs(d_yy[row,col])/a)**2-\
                              5*(abs(d_yy[row,col])/a)+4-\
                              2/3*a/abs(d_yy[row,col])
           elif d_yy[row,col]>2*a:
              rho_yy[row,col]=0
    for row in range(0,d_xy.shape[0]):
        for col in range(0,d_xy.shape[1]):
           d_xy[row,col]=np.sqrt((x_par[row]-x_obs[col])**2+\
                                 (y_par[row]-y_obs[col])**2+\
                                 (z_par[row]-z_obs[col])**2)
           if d_xy[row,col]<=a:
              rho_xy[row,col]=-1/4*(abs(d_xy[row,col])/a)**5+\
                               1/2*(abs(d_xy[row,col])/a)**4+\
                               5/8*(abs(d_xy[row,col])/a)**3-\
                               5/3*(abs(d_xy[row,col])/a)**2+\
                               1
           elif d_xy[row,col]>a  and d_xy[row,col]<=2*a:
              rho_xy[row,col]=1/12*(abs(d_xy[row,col])/a)**5-\
                              1/2*(abs(d_xy[row,col])/a)**4+\
                              5/8*(abs(d_xy[row,col])/a)**3+\
                              5/3*(abs(d_xy[row,col])/a)**2-\
                              5*(abs(d_xy[row,col])/a)+4-\
                              2/3*a/abs(d_xy[row,col])
           elif d_xy[row,col]>2*a:
              rho_xy[row,col]=0    
    for row in range(0,d_xx.shape[0]):
        for col in range(0,d_xx.shape[1]):
           np.sqrt((x_par[row]-x_par[col])**2+\
                   (y_par[row]-y_par[col])**2+\
                   (z_par[row]-z_par[col])**2)
           if d_xx[row,col]<=a:
              rho_xx[row,col]=-1/4*(abs(d_xx[row,col])/a)**5+\
                              1/2*(abs(d_xx[row,col])/a)**4+\
                              5/8*(abs(d_xx[row,col])/a)**3-\
                              5/3*(abs(d_xx[row,col])/a)**2+\
                              1
           elif d_xx[row,col]>a  and d_xx[row,col]<=2*a:
              rho_xx[row,col]=1/12*(abs(d_xx[row,col])/a)**5-\
                              1/2*(abs(d_xx[row,col])/a)**4+\
                              5/8*(abs(d_xx[row,col])/a)**3+\
                              5/3*(abs(d_xx[row,col])/a)**2-\
                              5*(abs(d_xx[row,col])/a)+4-\
                              2/3*a/abs(d_xx[row,col])
           elif d_xx[row,col]>2*a:
              rho_xx[row,col]=0
    return (rho_yy,rho_xy,rho_xx)