import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import TMM_Sims_Starshot as tm

lamb = np.arange(2000, 14005, 5)
Unit = ['PDMS (IR)']
thick = [3*180.8]
Sub = ['Air']
SubThick = [10]
ff = 1

angles = np.arange(90)
Abs = np.zeros((angles.shape[0], lamb.shape[0]))

for i, phi in enumerate(angles*np.pi/180):
    Abs[i, :] = tm.Stack(lamb, Unit, thick, Sub, SubThick, TM=0, phi=phi, output='Absorptance')
    Abs[i, :] += tm.Stack(lamb, Unit, thick, Sub, SubThick, TM=1, phi=phi, output='Absorptance')

Abs = Abs/2

data = np.zeros((angles.shape[0]+1, lamb.shape[0]+1)).astype('str')
data[0, 0] = 'Angle (o)/Wavelength (nm)'
data[1:, 0] = angles
data[0, 1:] = lamb
data[1:, 1:] = Abs
pd.DataFrame(data).to_excel('Starshot/Whittam Angle Dependent Abs.xlsx', index=False, header=None)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.pcolor(lamb, angles, Abs, shading='auto', cmap=plt.get_cmap('jet'))
plt.colorbar()
plt.show()
