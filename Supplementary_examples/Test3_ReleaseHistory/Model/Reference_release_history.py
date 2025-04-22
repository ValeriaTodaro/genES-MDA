import numpy as np
import os

def F(t):
    term1 = np.exp(-((t - 130) ** 2) / 50)
    term2 = 0.3 * np.exp(-((t - 150) ** 2) / 200)
    term3 = 0.5 * np.exp(-((t - 190) ** 2) / 98)
    return term1 + term2 + term3

t_values = np.arange(0, 300+1, 3)

F_values = F(t_values)

import matplotlib.pyplot as plt

plt.plot(t_values, F_values)
plt.xlabel("t")
plt.ylabel("C")
plt.show()

X=50
Y=45
Z=5

os.chdir('..')
par_file=np.vstack([np.full(F_values.size,X),np.full(F_values.size,Y), np.full(F_values.size,Z),t_values,F_values]).T
with open('Par.txt','wb') as file:
    np.savetxt(file, par_file, fmt=('%.1f','%.1f','%s','%.1f','%.5f'), newline='\n')