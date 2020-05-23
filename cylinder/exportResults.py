import numpy as np
import os
import scipy.signal as signal
import matplotlib.pyplot as plt

dir = os.getcwd()
name = "postProcessing/forceCoeffs/0/coefficient.dat"

data = np.loadtxt(name, skiprows=0)

t_start = 130
t_start = 0
u_inf = 1.0
D = 1.0

t  = data[:,0]
Cd = data[:,2]
Cl = data[:,3]

N  = len(t)
dt = t[1] - t[0]

nmax=512
freq, Cl_amp   = signal.welch(Cl, 1./dt, nperseg=nmax)
Cl_max_fft_idx = np.argmax(abs(Cl_amp))  
freq_shed      = freq[Cl_max_fft_idx]
St             = freq_shed*D/u_inf

idx = np.where(t > t_start)
Cd = Cd[idx]
Cl = Cl[idx]
Clzeros = np.where(np.multiply(np.sign(Cl[1:]),np.sign(Cl[:-1])) == -1)
idx = range(int(Clzeros[0][0]),int(Clzeros[0][-1]))

avgDragCoeff = (Cd[idx].max() + Cd[idx].min())/2.0
avgLiftCoeff = (Cl[idx].max() + Cl[idx].min())/2.0
peakToPeakLift = Cl[idx].max()-Cl[idx].min()
rmsLiftCoeff = np.sqrt(np.mean(Cl[idx]**2))
peakToPeakDrag = Cd[idx].max()-Cd[idx].min()

file1 = open('results.txt', 'w') 

file1.write("Average drag coefficient: %.5f (1.3353 from Rajani2008nso, 1.411 from NSCM_Cylinder)\n" % avgDragCoeff)
file1.write("Rms of lift coefficient: %.5f (0.1792 from Rajani2008nso)\n" % rmsLiftCoeff)
file1.write("Peak to peak lift: %.5f (0.727 from NSCM_Cylinder)\n" % peakToPeakLift)
file1.write("Peak to peak drag: %.5f (0.0203 from NSCM_Cylinder)\n" % peakToPeakDrag)
file1.write("Strouhals number: %.5f (0.1569 from Rajani2008nso, 0.173 from NSCM_Cylinder)\n" % St)
#file1.write("Average lift coefficient: %.15f\n" % avgLiftCoeff)

file1.close()




