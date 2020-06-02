import numpy as np
name = "postProcessing/forceCoeffs/0/coefficient.dat"
data = np.loadtxt(name, skiprows=0)
np.savetxt(name[:-4]+'.csv', data, delimiter=",", header='Time,Cd,Cs,Cl,CmRoll,CmPitch,CmYaw,Cd(f),Cd(r),Cs(f),Cs(r),Cl(f),Cl(r)', comments='')

