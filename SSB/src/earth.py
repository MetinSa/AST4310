import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Matplotlib estetics
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath = "../figures/1.3/"
filename = "earth.dat"

# Reading data
h, logP, temp, logdens, logN= np.loadtxt(path + filename, usecols=(0,1,2,3,4), unpack=True)

"""
Constants/Variables/Parameters

"""
D_sun = const.au.cgs.value
R_sun = const.R_sun.cgs.value
m_H = 1.67352e-24
N_phot = 1.351e+17


"""
Functions

"""
def ScaleHeight(h, logrho):

	"""
	Returns the scale height for at a given height and density

	"""
	return -h / (np.log((10**logrho)/(10**logdens[0])))


def SunShine():

	"""
	Returns the photon density

	"""
	return np.pi*((R_sun**2)/(D_sun**2))*N_phot


def PrintInfo():

	"""
	Prints usefull information

	"""
	print("H_P = %.4f" %ScaleHeight(100,logdens[np.argwhere(h == 8)[0][0]]))
	print("N_phoy = %.4E"%SunShine())


"""
Plotting functions

"""
def Temperature_Height():

	"""
	Plotting the temperature vs height

	"""
	plt.title("Temperature vs height (Earth)")
	plt.plot(h, temp, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel("Temperautre [K]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "Earth_temperature.pdf")
	plt.show()


def Pressure_Height():

	"""
	Plotting the pressure vs height

	"""
	plt.title("Pressure vs height (Earth)")
	plt.semilogy(h, 10**(logP), color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel("Pressure [dyne cm$^{-2}$]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "Earth_pressure.pdf")
	plt.show()


def GasDensity_Height():

	"""
	Plotting the particle density vs height

	"""
	plt.title("Gas density vs height (Earth)")
	plt.semilogy(h, 10**(logdens), color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel("Gas Density [g cm$^{-3}$]")
	plt.subplots_adjust(bottom = 0.12, left = 0.14)
	plt.savefig(savepath + "Earth_gasdensity.pdf")
	plt.show()


def ParticleDensity_Height():

	"""
	Plotting the particle density vs height

	"""
	plt.title("Particle density vs height (Earth)")
	plt.semilogy(h, 10**(logN), color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel("Particle Density [cm$^{-3}$]")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "Earth_particledensity.pdf")
	plt.show()


def NormalizedDensityPressure_Height():

	"""
	Co-Plotting the density and pressure in terms of normalized units

	"""
	plt.title("Normalized pressure and density (Earth)")
	plt.semilogy(h, 10**(logP)/10**(np.amax(logP)), color = "royalblue", label ="$P/P_0$")
	plt.semilogy(h, 10**(logdens)/10**(np.amax(logdens)), color = "crimson", label =r"$\rho/ \rho_0$")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.legend()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "Earth_Normalized_den_pres.pdf")
	plt.show()


def mu_Height():

	"""
	Plotting the mean molecular weight vs height 

	"""
	plt.title("Mean molecular weight (Earth)")
	plt.plot(h, 10**(logdens)/(10**(logN)*m_H), color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel(r"$\mu_E$")
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "Earth_mu.pdf")
	plt.show()


"""
Activating functions

"""
# Temperature_Height()
# Pressure_Height()
# GasDensity_Height()
# ParticleDensity_Height()
# NormalizedDensityPressure_Height()
# mu_Height()

# PrintInfo()