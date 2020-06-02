import numpy as np
import os
import scipy.signal as signal
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

dir = os.getcwd()

u_inf = 1.0
D = 1.0
def printResults(inputName,outputName,scaling,idx):
    # Rajani2008nso: When the statistically stationary state is indicated in the computation, the calculation is continued for the next 50 vortex shedding cycles and the time-averaged values are obtained by averaging the instantaneous field data on the flow variables for 50 shedding cycles. 
    data = np.loadtxt(inputName, skiprows=0)
    
    t  = data[:,idx[0]]
    Cd = data[:,idx[1]]*scaling
    Cl = data[:,idx[2]]*scaling
    
    # Compute Strouhals number
    dt = t[1] - t[0]
    freq, Cl_amp   = signal.welch(Cl, 1/dt)
    Cl_max_fft_idx = np.argmax(abs(Cl_amp))  
    f              = freq[Cl_max_fft_idx]    # approximated shedding frequency
    T              = 1/f                     # approximated shedding period
    
    # Find the final 50 shedding cycles
    t_start = t[-1] - 51*T # start at the time at which 50+1 periods remains (an extra period added for later cut off)
    if t_start < 0:
        print('Data set does not contain enough periods: t_end = '+str(t[-1])+', T = '+str(T))

    idx = np.where(t > t_start)
    t = t[idx]
    Cd = Cd[idx]
    Cl = Cl[idx]

    # Cut off data set to achieve full periods
    Clzeros = np.where(np.multiply(np.sign(Cl[1:]),np.sign(Cl[:-1])) == -1)
    idx = range(int(Clzeros[0][0]),int(Clzeros[0][-1]))
    t = t[idx]
    Cd = Cd[idx]
    Cl = Cl[idx]
    
    # Interpolate data set to obtain acurate extremas
    def calcExtremas(x,y):
        spl = interpolate.UnivariateSpline(x, y, k=4, s=0) # NB: sproot only works for order 3 splines, hence the 4th order
        extremas = spl.derivative().roots()
        T = extremas[-1] - extremas[-3]
        return spl(extremas).min(), spl(extremas).max(), T

    Cd_min, Cd_max, Td = calcExtremas(t,Cd)
    Cl_min, Cl_max, T = calcExtremas(t,Cl)
    
    f              = 1/T                     # shedding frequency
    St             = f*D/u_inf               # Strouhals number
    avgDragCoeff   = (Cd_max + Cd_min)/2.0
    peakToPeakLift = Cl_max - Cl_min
    peakToPeakDrag = Cd_max - Cd_min
    rmsLiftCoeff = np.sqrt(np.mean(Cl**2))
    
    file1 = open(outputName, 'w') 
    
    file1.write("Average drag coefficient: %.5f (1.3353 from Rajani2008nso, 1.411 from NSCM_Cylinder)\n" % avgDragCoeff)
    file1.write("Rms of lift coefficient: %.5f (0.1792 from Rajani2008nso)\n" % rmsLiftCoeff)
    file1.write("Peak to peak lift: %.5f (0.727 from NSCM_Cylinder)\n" % peakToPeakLift)
    file1.write("Peak to peak drag: %.5f (0.0203 from NSCM_Cylinder)\n" % peakToPeakDrag)
    file1.write("Strouhals number: %.5f (0.1569 from Rajani2008nso, 0.173 from NSCM_Cylinder)\n" % St)
    
    file1.close()
    
    
printResults("postProcessing/forceCoeffs/0/coefficient.dat","results_OF.txt",1.0,[0,1,3])
#printResults("/home/zetison/hugeFiles/IFEM/cylinder/Re100/1/chorin/Cyl2D_force.dat","results_IFEM.txt",1/0.5,[0,1,2])
    
    
