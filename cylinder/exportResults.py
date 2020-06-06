import numpy as np
import click
import scipy.signal as signal
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

@click.command()
@click.option('--ifem/--openfoam', default=False)
@click.option('--no-periods', default=50)
@click.option('--plots', default=False)
@click.option('--u_inf', default=1.0)
@click.option('--diam', default=1.0)
@click.option('--rho', default=1.0)
@click.option('--length', default=1.0)
@click.option('--inputname', default="postProcessing/forceCoeffs/0/coefficient.dat")
@click.option('--outputname', default="results_OF.txt")
def main(ifem,no_periods,plots,u_inf,diam,rho,length,inputname,outputname):
	D = diam	
	L = length
	if ifem:
	    scaling = 1/(0.5*rho*u_inf**2*D*L) # ifem results must be scaled as the inputs are forces and not coefficients
	    idx = [0,1,2]
	else:
	    scaling = 1.0
	    idx = [0,1,3]
	
	# Rajani2008nso: When the statistically stationary state is indicated in the computation, the calculation is continued for the next 50 vortex shedding cycles and the time-averaged values are obtained by averaging the instantaneous field data on the flow variables for 50 shedding cycles. 
	data = np.loadtxt(inputname, skiprows=0)
	
	t  = data[:,idx[0]]
	Cd = data[:,idx[1]]*scaling
	Cl = data[:,idx[2]]*scaling
	
	# Compute Strouhals number
	dt = t[1] - t[0]
	freq, Cl_amp   = signal.welch(Cl, 1/dt, nperseg=512)
	Cl_max_fft_idx = np.argmax(abs(Cl_amp))  
	f              = freq[Cl_max_fft_idx]    # approximated shedding frequency
	T              = 1/f                     # approximated shedding period
	
	# Find the final no_periods shedding cycles
	t_start = t[-1] - (no_periods+2)*T # start at the time at which no_periods+1 periods remains (two extra periods added for later cut off)
	if t_start < 0:
		print('Data set does not contain enough shedding periods: t_end = '+str(t[-1])+', T = '+str(T))
		return -1
	
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
	    spl = interpolate.UnivariateSpline(x, y, k=4, s=0.001) # NB: sproot only works for order 3 splines, hence the 4th order
	    extremas = spl.derivative().roots()
	    T = extremas[-1] - extremas[-3]
	    return spl, spl(extremas).min(), spl(extremas).max(), T
	
	spl_d, Cd_min, Cd_max, Td = calcExtremas(t,Cd)
	spl_l, Cl_min, Cl_max, T = calcExtremas(t,Cl)
	
	
	f              = 1/T                     # shedding frequency
	St             = f*D/u_inf               # Strouhals number
	avgDragCoeff   = (Cd_max + Cd_min)/2.0
	peakToPeakLift = Cl_max - Cl_min
	peakToPeakDrag = Cd_max - Cd_min
	rmsLiftCoeff = np.sqrt(np.mean(Cl**2))
	
	if plots:
	    t_fine = np.linspace(t[0],t[-1],100000)
	    plt.plot(t_fine,spl_l(t_fine),t,Cl,t,Cd, \
	             t,avgDragCoeff*np.ones(t.shape), \
	             t,Cd_max*np.ones(t.shape), \
	             t,Cd_min*np.ones(t.shape), \
	             t,Cl_max*np.ones(t.shape), \
	             t,Cl_min*np.ones(t.shape))
	    plt.legend(('Spline approximation','Cl','Cd','Average drag','Cd_max','Cd_min','Cl_max','Cl_min'))
	    plt.show()
	
	file1 = open(outputname, 'w') 
	
	file1.write("Average drag coefficient: %.5f (1.3353 from Rajani2008nso, 1.411 from NSCM_Cylinder)\n" % avgDragCoeff)
	file1.write("Rms of lift coefficient: %.5f (0.1792 from Rajani2008nso)\n" % rmsLiftCoeff)
	file1.write("Peak to peak lift: %.5f (0.727 from NSCM_Cylinder)\n" % peakToPeakLift)
	file1.write("Peak to peak drag: %.5f (0.0203 from NSCM_Cylinder)\n" % peakToPeakDrag)
	file1.write("Strouhals number: %.5f (0.1569 from Rajani2008nso, 0.173 from NSCM_Cylinder)\n" % St)
	
	file1.close()
	return 0

if __name__ == '__main__':
	main()
