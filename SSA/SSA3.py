import numpy as np 
import matplotlib.pyplot as plt 
import astropy.constants as cc
from cycler import cycler


default_cycler = (cycler(color=['royalblue', 'crimson', 'mediumseagreen', 'darkorange', "darkorchid", "sienna"]))
plt.rc('axes', prop_cycle=default_cycler)

"""
Constants.

"""

k_eV = 8.61734e-5 							#boltzmann constant/eV
k_erg = 1.380658e-16						#boltzmann*erg
m_e = 9.10939e-28							#electron mass
h = 6.62607e-27								#plancks constant in erg
c = 2.99792e10								#speed of light

"""
Functions.

"""

def planck(T,wav):

	"""
	Plancks Law.
	"""

	return (2*h*c**2*wav**-5)*(1./(np.exp(h*c/(wav*k_erg*T))-1.))





"""
Plots.

"""
def plot_planck():

	"""
	Plotting Plancks law.

	"""

	wav = np.linspace(1000,20800,1000)
	B = np.zeros(len(wav))

	for T in range(5000,8000,200):
		B[:] = planck(T, wav[:]*1e-8)
		plt.plot(wav,B,color = "royalblue")

	plt.title(r"Planck's Law for T $\in$ [5000,8000] ")
	plt.xlabel(r'$\lambda$ [$\AA$]')
	plt.ylabel(r"$B_\lambda$")
	plt.xlim(0,20800)
	plt.grid(ls ="--")
	plt.savefig("planck.pdf")
	plt.show()


plot_planck()










