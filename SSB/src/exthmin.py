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
filename1 = "falc.dat"
filename2 = "solspect.dat"

# Reading data
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt(path + filename1,usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, Fsmoothed, F, Ismoothed, I = np.loadtxt(path + filename2, usecols=(0,1,2,3,4), unpack=True)


"""
Constants/Variables/Parameters

"""
c = const.c.cgs.value
h_planck = const.h.cgs.value
k_B = const.k_B.cgs.value
wav_to_freq = (wav**2/(c))*1e6
sigmaT= 6.648e-25

def exthmin(wav,T,eldens):
	"""
	Extinction function copied from the SSB paper

	"""
	# Parameters
	theta=5040./T
	elpress=eldens*k_B*T

	# Info about the function:
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
	kappabf=kappabf*(1.-np.exp(-h_planck*c/(wav*1e-8*k_B*T)))# correct stimulated

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


def Extinction_H0():

	"""
	Plotting the extinction at the solar surface

	"""
	T = temp[np.argwhere(h ==0)[0][0]]
	n_el = nel[np.argwhere(h ==0)[0][0]]
	wavelength = np.linspace(wav[0],2,80)*1e4
	EXmin = exthmin(wavelength,T,n_el) 

	plt.figure()
	plt.title("H$^{-}$ extinction at h = 0")
	plt.plot(wavelength*1e-4,EXmin, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"H$^{-}$ extinction [cm$^2$ per H-atom]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath1 + "extinction_H0.pdf")

	plt.figure()
	plt.title("H$^{-}$ transmittance at h = 0")
	plt.plot(wavelength*1e-4,1/EXmin, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"H$^{-}$ extinction$^{-1}$ [cm$^2$ per H-atom]$^{-1}$")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath1 + "extinction_transmittance.pdf")
	plt.show()


def Extinction_comparison():

	"""
	Plotting the comparison of extinction profiles

	"""
	sigmaT= 6.648e-25
	nneutH=nhyd-nprot
	wavelength = 0.5*1e4
	Thomson = nel*sigmaT

	EXmin = exthmin(wavelength,temp,nel) *nneutH

	plt.title(r"H$^{-}$ extinction for $\lambda$ = 5000 $\mathrm{\AA}$")
	plt.semilogy(h,EXmin, color = "royalblue", label = r"H$^{-}$")
	plt.semilogy(h,Thomson, color = "mediumseagreen", label = "Thomson")
	plt.semilogy(h,EXmin+Thomson, color = "crimson", label = "Total ")
	plt.grid(linestyle = "--")
	plt.legend()
	plt.xlabel(r"Height [km]")
	plt.ylabel(r"$\alpha_\lambda$ [cm$^{-1}$]")
	plt.xlim(-200,2150)
	plt.ylim(1e-18,1e-5)
	plt.subplots_adjust(bottom = 0.12, left = 0.15)
	plt.savefig(savepath1 + "extinction_comparison.pdf")
	plt.show()


def Optical_Depth():

	"""
	Plotting the optical depth as a function of height

	"""
	tau = np.zeros(len(tau5), dtype=float) # initializing tau array 
	nneutH=nhyd-nprot
	sigmaT= 6.648e-25
	EXmin = exthmin(5000,temp,nel)*nneutH
	Thomson = nel*sigmaT
	ext = EXmin+Thomson

	for i in range(1,len(tau)):
		tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5

	plt.semilogy(h,tau5, color = "royalblue", label = r"$\tau_{5000}$ from FALC")
	plt.semilogy(h,tau, color = "crimson", label = r"$\tau_{5000}$ from H$^{-}$ and Thomson extinction")

	plt.grid(linestyle = "--")
	plt.legend()
	plt.xlabel(r"Height [km]")
	plt.ylabel(r"Optical depth $\tau_\lambda$")
	# plt.xlim(-200,2150)
	# plt.ylim(1e-18,1e-5)
	plt.subplots_adjust(bottom = 0.12, left = 0.15)
	plt.savefig(savepath2 + "extinction_5000.pdf")
	plt.show()


def PrintInfo():

	"""
	Printing height for tau = 1 for a set of wavelengths

	"""
	tau = np.zeros(len(tau5), dtype=float) # initializing tau array 
	nneutH=nhyd-nprot
	sigmaT= 6.648e-25
	wl_list = [0.5, 1, 1.6, 5]
	Thomson = nel*sigmaT
	for wl in wl_list:
		EXmin = exthmin(wl*1e4,temp,nel)*nneutH
		ext = EXmin+Thomson

		for i in range(1,len(tau)):
			tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5

		index = (np.abs(1-tau)).argmin()
		# print(index)

		H_tau1 = h[index]
		print(H_tau1)

# Extinction_H0()
# Extinction_comparison()
# Optical_Depth()
# PrintInfo()
