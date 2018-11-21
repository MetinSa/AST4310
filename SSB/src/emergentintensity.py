import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Matplotlib estetics
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath1 = "../figures/2.2/"
savepath2 = "../figures/2.3/"
savepath3 = "../figures/2.4/"
filename1 = "falc.dat"
filename2 = "solspect.dat"

# Reading data
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt(path + filename1,usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, Fsmoothed, F, Ismoothed, I = np.loadtxt(path + filename2, usecols=(0,1,2,3,4), unpack=True)


"""
Constants/Variables/Parameters

"""

c = const.c.cgs.value
# h = const.h.cgs.value
k_B = const.k_B.cgs.value
wav_to_freq = (wav**2/(c))*1e6
sigmaT= 6.648e-25


"""
Functions

"""
def Planck(wav, T):

	"""
	Returns the planck curves for a given temperature and wavelength

	"""

	return (2*const.h.cgs.value*c**2*wav**(-5))*(1./(np.exp(const.h.cgs.value*c/(wav*k_B*T))-1.))

def Emergent_Intensity():
	contfuncCalc=np.zeros(len(wav))
	for j in range(len(wav)):
		ext = np.zeros(len(tau5))
		tau = np.zeros(len(tau5))
		integrand = np.zeros(len(tau5))
		contfunc = np.zeros(len(tau5))
		intt = 0.0
		hint = 0.0
		for ih in range(1, len(tau5)):
			ext[ih] = (exthmin(wav[j]*1e4, temp[ih], nel[ih])*(nhyd[ih]-nprot[ih])+ sigmaT*nel[ih])
			tau[ih] = tau[ih-1] + 0.5*(ext[ih] + ext[ih-1])*(h[ih-1]-h[ih])*1e5
			integrand[ih] = Planck(wav[j]*1e-4, temp[ih])*np.exp(-tau[ih])
			intt += 0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
			hint += h[ih]*0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
			contfunc[ih] = integrand[ih]*ext[ih]
		contfuncCalc[j]=intt
		# note : exthmin has wavelength in [Angstrom], planck in [cm]
		mean = hint / intt
		#tau5[69]=1
	return contfuncCalc, mean


def exthmin(wav,T,eldens):
	# other parameters
	c = const.c.cgs.value
	h = const.h.cgs.value
	k = const.k_B.cgs.value
	# wav_to_freq = (wav**2/(c))*1e6
	theta=5040./T
	elpress=eldens*k*T

	# evaluate H-min bound-free per H-min ion ? Gray (8.11)
	# his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
	sigmabf = (1.99654 -1.18267e-5*wav +2.64243e-6*wav**2-4.40524e-10*wav**3+3.23992e-14*wav**4-1.39568e-18*wav**5 +2.78701e-23*wav**6)
	sigmabf *= 1e-18 # cm^2 per H-min ion
	if np.size(wav) > 1:
		sigmabf[np.argwhere(wav>16444)] = 0 # H-min ionization limit at lambda = 1.6444 micron
	elif ( np.size(wav) == 1):
		if wav> 16444:
			sigmabf = 0

	# convert into bound-free per neutral H atom assuming Saha = Gray p135
	# units: cm2 per neutral H atom in whatever level (whole stage)
	graysaha=4.158e-10*elpress*theta**2.5*10.**(0.754*theta)# Gray (8.12)
	kappabf=sigmabf*graysaha # per neutral H atom
	kappabf=kappabf*(1.-np.exp(-h*c/(wav*1e-8*k*T)))# correct stimulated

	# check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
	# logratio=-0.1761-np.log10(elpress)+np.log10(2.)+2.5*np.log10(T)-theta*0.754
	# print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB

	# evaluate H-min free-free including stimulated emission = Gray p136
	lwav = np.log10(wav)
	f0 = - 2.2763 - 1.6850*lwav + 0.76661*lwav**2 - 0.0533464*lwav**3
	f1 =   15.2827 - 9.2846*lwav + 1.99381*lwav**2 - 0.142631*lwav**3
	f2 = - 197.789 + 190.266*lwav - 67.9775*lwav**2 + 10.6913*lwav**3 - 0.625151*lwav**4
	ltheta = np.log10(theta)
	kappaff = 1e-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2) #Gray(8.13)
	return kappabf+kappaff



def Emergent_IntensityPlot():

	contfuncCalc , mean = Emergent_Intensity()

	plt.title("Observed and computed continuum intensity")
	plt.plot(wav, contfuncCalc*1e-14, color = "royalblue", label = "FALC")
	plt.plot(wav,I, color = "crimson", label = "Observed")
	plt.legend()
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$ [$\mu$ m]")
	plt.ylabel(r"Intensity [$10^{14}$ erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]")
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath3 + "Observed_computed.pdf")
	plt.show()

def Contribution_single():
	wl = 0.5
	ext = np.zeros(np.size(tau5)) 
	tau = np.zeros(np.size(tau5))
	integrand = np.zeros(np.size(tau5)) 
	contfunc = np.zeros(np.size(tau5))
	intt = 0.0 
	hint = 0.0
	for ih in range(1, len(tau)): 
		ext[ih] = exthmin(wl*1e4, temp[ih], nel[ih])*(nhyd[ih]-nprot[ih])+sigmaT*nel[ih]
		tau[ih] = tau[ih-1] + 0.5*(ext[ih]+ext[ih-1])*(h[ih-1]-h[ih])*1e5
		integrand[ih] = Planck(wl*1e-4, temp[ih])*np.exp(-tau[ih])
		intt += 0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
		hint += h[ih]*0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
		contfunc[ih] = integrand[ih]*ext[ih]
	hmean = hint / intt

	plt.title("Peak-normalized contribution function")
	plt.plot(h, contfunc/np.amax(contfunc),color = "royalblue")
	plt.axvline(x=hmean, color = "black",linestyle = "--", alpha = 0.7, label = r"h$_\mathrm{mean} =$ %.1f"%hmean)
	plt.grid(linestyle = "--")
	plt.xlabel("Height[km]")
	plt.ylabel("Contribution function")
	plt.xlim(-100,500)
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath3 + "contfunc_single.pdf")
	plt.show()



def Contribution_multi():
	colors = ["royalblue", "darkorange", "mediumseagreen", "crimson"]
	idum = 0
	wl_list = [0.5, 1, 1.6, 5]

	for wl in wl_list:
		ext = np.zeros(np.size(tau5)) 
		tau = np.zeros(np.size(tau5))
		integrand = np.zeros(np.size(tau5)) 
		contfunc = np.zeros(np.size(tau5))
		intt = 0.0 
		hint = 0.0
		for ih in range(1, len(tau)): 
			ext[ih] = exthmin(wl*1e4, temp[ih], nel[ih])*(nhyd[ih]-nprot[ih])+sigmaT*nel[ih]
			tau[ih] = tau[ih-1] + 0.5*(ext[ih]+ext[ih-1])*(h[ih-1]-h[ih])*1e5
			integrand[ih] = Planck(wl*1e-4, temp[ih])*np.exp(-tau[ih])
			intt += 0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
			hint += h[ih]*0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])
			contfunc[ih] = integrand[ih]*ext[ih]
		hmean = hint / intt

		plt.plot(h, contfunc/np.amax(contfunc), color = colors[idum], label = r"$\lambda$ = %g, h$_\mathrm{mean}$ = %.1f" %(wl,hmean))
		idum += 1

	plt.title("Peak-normalized contribution function")
	plt.grid(linestyle = "--")
	plt.xlabel("Height[km]")
	plt.ylabel("Contribution function")
	plt.xlim(-100,500)
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath3 + "contfunc.pdf")
	plt.show()




# Emergent_IntensityPlot()
Contribution_single()
# Contribution_multi()