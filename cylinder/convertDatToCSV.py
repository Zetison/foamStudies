import pandas as pd
import itertools
import os

dir = os.getcwd()
name = "postProcessing/forceCoeffs/0/coefficient.dat"
data = np.loadtxt(name, skiprows=0)
np.savetxt(name[:-4]+'.csv', data[:,:-1], delimiter=",", header='Time,Cd,Cl', comments='')

