import numpy as np
import os
import scipy.signal as signal
import matplotlib.pyplot as plt

dir = os.getcwd()

u_inf = 1.0
D = 1.0
def printResults(inputName,outputName,scaling,idx):
    
    data = np.loadtxt(inputName, skiprows=0)
    
    t_start = 100
    #t_start = 0
    
    t  = data[:,idx[0]]
    Cd = data[:,idx[1]]*scaling
    Cl = data[:,idx[2]]*scaling
    
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
    
    file1 = open(outputName, 'w') 
    
    file1.write("Average drag coefficient: %.5f (1.3353 from Rajani2008nso, 1.411 from NSCM_Cylinder)\n" % avgDragCoeff)
    file1.write("Rms of lift coefficient: %.5f (0.1792 from Rajani2008nso)\n" % rmsLiftCoeff)
    file1.write("Peak to peak lift: %.5f (0.727 from NSCM_Cylinder)\n" % peakToPeakLift)
    file1.write("Peak to peak drag: %.5f (0.0203 from NSCM_Cylinder)\n" % peakToPeakDrag)
    file1.write("Strouhals number: %.5f (0.1569 from Rajani2008nso, 0.173 from NSCM_Cylinder)\n" % St)
    #file1.write("Average lift coefficient: %.15f\n" % avgLiftCoeff)
    
    file1.close()
    
    
#printResults("postProcessing/forceCoeffs/0/coefficient.dat","results_OF.txt",1.0,[0,1,3])
printResults("/home/zetison/hugeFiles/IFEM/cylinder/Re100/1/chorin/Cyl2D_force.dat","results_IFEM.txt",1/0.5,[0,1,2])
    
    
