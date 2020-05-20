import pandas as pd
import numpy as np
import itertools
import os
import numpy.fft as fft
import matplotlib.pyplot as plt

dir = os.getcwd()
name = "postProcessing/forces/0/coefficient.dat"

with open(name) as handle:
    *_comments, names = itertools.takewhile(lambda line: line.startswith('#'), handle)
    # This is not the most robust way, adjust for your needs :)
    names = names[1:].split()

headers = pd.read_table(name, header=0, names=names, sep='\t', comment='#')
headers.to_csv(name[:-4]+'.csv', index=False)

t_start = 130
t_start = 130
u_inf = 1.0
D = 1.0

t = headers.Time.to_numpy()
Cd = headers.Cd.to_numpy()
Cl = headers.Cl.to_numpy()
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

#spectrum = fft.fft(Cl[idx] - avgLiftCoeff)
p = np.abs(np.fft.rfft(Cl[idx]))
#freq = fft.fftfreq(len(spectrum))
rate = 1/(t[1]-t[0])
freq = np.linspace(0, rate/2, len(p))
#threshold = 0.5 * max(abs(spectrum))
threshold = 0.99 * max(p)
#mask = abs(spectrum) > threshold
mask = p > threshold
print(freq[mask])
f = freq[mask]
StrouhalNumber = f[-1]*D/u_inf

#plt.plot(freq,p)
#plt.plot(t[idx],(Cd[idx] - avgDragCoeff)/peakToPeakDrag*2,t[idx],Cl[idx]/peakToPeakLift*2)
t = t[idx]
#print(np.any(abs(t[1:]-t[:-1]-0.005) > 1e-10))
#print(rate)
plt.show()

file1 = open('results.txt', 'w') 

file1.write("Average drag coefficient: %.5f (1.3353 from Rajani2008nso, 1.411 from NSCM_Cylinder)\n" % avgDragCoeff)
file1.write("Rms of lift coefficient: %.5f (0.1792 from Rajani2008nso)\n" % rmsLiftCoeff)
file1.write("Peak to peak lift: %.5f (0.727 from NSCM_Cylinder)\n" % peakToPeakLift)
file1.write("Peak to peak drag: %.5f (0.0203 from NSCM_Cylinder)\n" % peakToPeakDrag)
file1.write("Strouhals number: %.5f (0.1569 from Rajani2008nso, 0.173 from NSCM_Cylinder)\n" % StrouhalNumber)
#file1.write("Average lift coefficient: %.15f\n" % avgLiftCoeff)

file1.close()




