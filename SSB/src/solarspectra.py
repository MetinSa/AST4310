import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Matplotlib estetics
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath = "../figures/2.1/"
filename = "solspect.dat"

# Reading data
wav, Fsmoothed, F, Ismoothed, I = np.loadtxt(path + filename,usecols=(0,1,2,3,4), unpack=True)


"""
Constants/Variables/Parameters

"""
c = const.c.cgs.value
h = const.h.cgs.value
k_B = const.k_B.cgs.value
wav_to_freq = (wav**2/(c))*1e6
T_sun = 5770


"""
Functions

"""
def Planck(wav, temp):

	"""
	Returns the planck curves for a given temperature and wavelength

	"""
	return (2*h*c**2*wav**(-5))*(1./(np.exp(h*c/(wav*k_B*temp))-1.))


def BrightnessTemperature(wav, I):

	"""
	Inverts the Planck function and returns the brightness temperature

	"""
	return ((h*c)/(k_B*wav))/np.log(1 + ((2*h*c**2)/(I*wav**5))*1e-4)



"""
Plotting functions

"""
def SpectralDistributions_Wav():

	"""
	Plotting the four spectral distributions in wavelength

	"""
	plt.title("Solar continuum")
	plt.plot(wav,I, color = "royalblue", label = r"$I_\lambda^\mathrm{cont}$")
	plt.plot(wav,Ismoothed, color = "darkorange", label = r"$I_\lambda^\mathrm{cont}$ smoothed")
	plt.plot(wav,F, color = "crimson", label = r"$F_\lambda^\mathrm{cont}$")
	plt.plot(wav,Fsmoothed, color = "mediumseagreen", label = r"$F_\lambda^\mathrm{cont}$ smoothed")
	plt.grid(linestyle = "--")
	plt.legend()
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"Intensity and flux [$10^{10}$ erg cm$^{-2}$ s$^{-1}$ $\mu$m$^{-1}$ ster$^{-1}$]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "solarcontinuum.pdf")
	plt.show()

	print("max(Ic) = ",np.amax(I),"at",wav[np.argwhere(I == np.max(I))][0][0])


def SpectralDistributions_Freq():

	"""
	Plotting the four spectral distributions in frequeny

	"""
	plt.title("Solar continuum")
	plt.plot(wav,I*wav_to_freq, color = "royalblue", label = r"$I_\nu^\mathrm{cont}$")
	plt.plot(wav,Ismoothed*wav_to_freq, color = "darkorange", label = r"$I_\nu^\mathrm{cont}$ smoothed")
	plt.plot(wav,F*wav_to_freq, color = "crimson", label = r"$F_\nu^\mathrm{cont}$")
	plt.plot(wav,Fsmoothed*wav_to_freq, color = "mediumseagreen", label = r"$F_\nu^\mathrm{cont}$ smoothed")
	plt.grid(linestyle = "--")
	plt.legend()
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"Intensity and flux [$10^{10}$ erg cm$^{-2}$ s$^{-1}$ ster$^{-1}$]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "solarcontinuum_freq.pdf")
	plt.show()

	print("max(Ic) = ",np.amax(I*wav_to_freq),"at",wav[np.argwhere(I*wav_to_freq == np.max(I*wav_to_freq))][0][0])


def PlanckPlot():

	"""
	Plotting the planck fit

	"""
	T_test = 6300

	plt.title("Planck fit to observed continuum")
	plt.plot(wav,I, color = "royalblue", label = r"$I_\nu^\mathrm{cont}$")
	plt.plot(wav,Planck(wav/1e4,T_test)*1e-14, color = "mediumseagreen", label = r"$B_\lambda(T =$ %g)" %T2)
	plt.plot(wav,Planck(wav/1e4,T1)*1e-14, color = "crimson", label = r"$B_\lambda(T =$ %g)" %T1)
	plt.grid(linestyle = "--")
	plt.legend()
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"Intensity [$10^{10}$ erg cm$^{-2}$ s$^{-1}$ ster$^{-1}$]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "solarcontinuum_planckfit.pdf")
	plt.show()


def BrightnessTempPlot():

	"""
	Plotting the brightness temperature

	"""
	B_temp = BrightnessTemperature(wav*1e-4,I*1e10)
	plt.plot(wav,B_temp, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$[$\mu$m]")
	plt.ylabel(r"Brightness temperature [K]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "solarcontinuum_brighttemp.pdf")
	plt.show()


"""
Activating functions

"""
# SpectralDistributions_Wav()
# SpectralDistributions_Freq()
# PlanckPlot()
BrightnessTempPlot()