import matplotlib.pyplot as plt
import numpy as np
#name = "~/OpenFOAM/OpenFOAM-v1912/run/postProcessing/forceCoeffs/0/coefficient.dat"
#data = np.loadtxt(name, skiprows=0)
#np.savetxt(name[:-4]+'.csv', data[:,:-1], delimiter=",", header='Time,Cd,Cl', comments='')
#t = data[0,:,-1]
#Cd = data[1,:,-1]
#Cl = data[2,:,-1]
name = "/home/zetison/hugeFiles/IFEM/cylinder/Re100/1/chorin/Cyl2D_force.dat"
data = np.loadtxt(name, skiprows=0)
t = data[:,0]
Cd = data[:,1]/0.5
Cl = data[:,2]/0.5
drag = plt.plot(t,Cd,label='Drag')
lift = plt.plot(t,Cl,label=('Lift'))
plt.legend()
plt.xlabel('Time')
plt.ylabel('Coefficients')
plt.ylim(-2,1)
plt.show()
