import numpy as np
from scipy.special import expi

def well_function(T, S):
    """
    Computes the drawdown at a given radial distance and time due to a pumping well using Theis well function
    
    Parameters:
    Q  - Pumping rate (m³/s)
    T  - Transmissivity (m²/s)
    S  - Storativity (dimensionless)
    r  - Radial distance from the well (m)
    t  - Time since pumping started (s)
    
    Returns:
    s : float - Drawdown at distance r and time t
    """
    Q = 5e-4  # Pumping rate in (m³/s)
    r = 15    # Distance from well (m)
    t = np.arange(1, 2*60+1, 1)*60 #Time (s)
        
    u = (r**2 * S) / (4 * T * t)
    
    # Theis well function W(u) using the exponential integral
    W_u = -expi(-u)
    
    s = (Q / (4 * np.pi * T)) * W_u
    
    with open('pred.txt','wb') as file:
        data = np.column_stack((t, s))
        np.savetxt(file, data, fmt=['%.0f','%.3e'], newline='\n')