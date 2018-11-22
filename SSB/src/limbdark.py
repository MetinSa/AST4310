import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Matplotlib estetics
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath = "../figures/2.6/"
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
	mu = np.linspace(0.1,1.,10)
	contfuncMu = np.zeros((len(mu),len(wav)))
	contfuncCalc=np.zeros(len(wav))
	for j in range(len(wav)):
		
		for w in range(len(mu)):
			ext = np.zeros(len(tau5))
			tau = np.zeros(len(tau5))
			integrand = np.zeros(len(tau5))
			contfunc = np.zeros(len(tau5))
			intt = 0.0
			hint = 0.0
			for i in range(1, len(tau5)):
				ext[i] = (exthmin(wav[j]*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+ sigmaT*nel[i])
				tau[i] = tau[i-1] + 0.5*(ext[i] + ext[i-1])*(h[i-1]-h[i])*1e5
				integrand[i] = Planck(wav[j]*1e-4,temp[i])*np.exp(-tau[i]/mu[w])
				intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu[w]
				hint += h[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu[w]
				contfunc[i] = integrand[i]*ext[i]
			contfuncMu[w,j]=intt
	return contfuncMu, mu


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


def limb():

	contfuncMu, mu = Emergent_Intensity()


	plt.title("Limb darkening")
	for i in range(len(mu)):
		plt.plot(wav, contfuncMu[i,:] , label = r"$\mu =$ %.1f"%mu[i])
	plt.legend()
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$ [$\mu$ m]")
	plt.ylabel(r"Intensity [erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]")
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "limbdark_wave.pdf")
	plt.show()

def limbangle():

	contfuncMu, mu = Emergent_Intensity()

	ratio1 = np.zeros(len(mu))
	ratio2 = np.zeros(len(mu))
	ratio3 = np.zeros(len(mu))
	ratio4 = np.zeros(len(mu))	
	ratio5 = np.zeros(len(mu))	

	for i in range(len(mu)):
	    ratio1[i] = contfuncMu[i, np.argwhere(wav == 0.2)[0][0]]/contfuncMu[-1, np.argwhere(wav == 0.2)[0][0]]
	    ratio2[i] = contfuncMu[i, np.argwhere(wav == 0.5)[0][0]]/contfuncMu[-1, np.argwhere(wav == 0.5)[0][0]]
	    ratio3[i] = contfuncMu[i, np.argwhere(wav == 1.6)[0][0]]/contfuncMu[-1, np.argwhere(wav == 1.6)[0][0]]
	    ratio4[i] = contfuncMu[i, np.argwhere(wav == 2.5)[0][0]]/contfuncMu[-1, np.argwhere(wav == 2.5)[0][0]]
	    ratio5[i] = contfuncMu[i, np.argwhere(wav == 5.0)[0][0]]/contfuncMu[-1, np.argwhere(wav == 5.0)[0][0]]

	rRsun = np.sin(np.arccos(mu))

	plt.title("Limb darkening")
	plt.plot(rRsun, ratio1, label = r"$\lambda = 0.2 \mu$ m")
	plt.plot(rRsun, ratio2, label = r"$\lambda = 0.5 \mu$ m")
	plt.plot(rRsun, ratio3, label = r"$\lambda = 1.6 \mu$ m")
	plt.plot(rRsun, ratio4, label = r"$\lambda = 2.5 \mu$ m")
	plt.plot(rRsun, ratio5, label = r"$\lambda = 5.0 \mu$ m")
	plt.legend()
	plt.grid(linestyle = "--")
	plt.xlabel(r"sin $\theta = r/ R_\odot$")
	plt.ylabel(r"Normalized intensity")
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "limbdark_rRsun.pdf")
	plt.show()

def flux():

	xgauss=[-0.7745966692,0.0000000000,0.7745966692]
	wgauss=[ 0.5555555555,0.8888888888,0.5555555555]
	fluxspec = np.zeros(len(wav),dtype=float)
	intmu = np.zeros((3,len(wav)), dtype=float)
	for imu in range(3):
		mu=0.5+xgauss[imu]/2.
		# rescale xrange [-1,+1] to [0,1]
		wg=wgauss[imu]/2.
		# weights add up to 2 on [-1,+1]
		for iw in range(0,len(wav)):
			wl=wav[iw]
			ext = np.zeros(len(tau5))
			tau = np.zeros(len(tau5))
			integrand = np.zeros(len(tau5))
			intt = 0.0
			for i in range(1, len(tau5)):
				ext[i] = (exthmin(wl*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+ sigmaT*nel[i])
				tau[i] = (tau[i-1] + 0.5*(ext[i] + ext[i-1])*(h[i-1]-h[i])*1e5)
				integrand[i] = Planck(wl*1e-4, temp[i])*np.exp(-tau[i]/mu)
				intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu
			intmu[imu,iw]=intt
			fluxspec[iw]=fluxspec[iw] + wg*intmu[imu,iw]*mu
	fluxspec*= 2

	plt.title("Observed and computed continuum flux")
	plt.plot(wav, fluxspec*1e-14, label = "FALC", color = "royalblue")
	plt.plot(wav, F, label = "Observed", color = "crimson")
	plt.legend()
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$ [$\mu$ m]")
	plt.ylabel(r"Astrophysical flux [$10^{14}$ erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]")
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "flux_comparison.pdf")
	plt.show()


# limb()
# limbangle()
# flux()