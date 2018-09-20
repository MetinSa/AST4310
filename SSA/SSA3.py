import numpy as np 
import matplotlib.pyplot as plt 
import astropy.constants as cc
from cycler import cycler
from scipy import special


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

	return (2*h*c**2*wav**(-5))*(1./(np.exp(h*c/(wav*k_erg*T))-1.))


def wien(T,wav):

	"""
	Wien approximation to Planck.

	"""

	return (2*h*c**2*wav**(-5))*np.exp(-h*c/(wav*k_erg*T))


def rayleighjeans(T,wav):

	"""
	Rayleigh-Jeans approximation to Planck.

	"""

	return 2*c*k_erg*T*wav**(-4)


def voigt(gamma,x):

	"""
	Voigt profil.

	"""

	return special.wofz(x+1j*gamma).real 




"""
Plots.

"""
def plot_planck():

	"""
	Plotting Plancks law.

	"""

	wav = np.linspace(1000,20800,1000)
	B = np.zeros(len(wav))
	dummy = 1
	for T in range(5000,8000,200):

		B[:] = planck(T, wav[:]*1e-8)

		if dummy == 1:
			plt.plot(wav,B,color = "royalblue", label = r"T $\in$ [5000,8000] K")
			dummy += 1

		else:
			plt.plot(wav,B,color = "royalblue")



	# plt.yscale('log')
	# plt.xscale('log')
	plt.title(r"Planck's Law ")
	plt.xlabel(r'$\lambda$ [$\AA$]')
	plt.ylabel(r"$B_\lambda$")
	plt.xlim(0,20800)
	plt.legend()
	plt.grid(ls ="--")
	# plt.savefig("planck.pdf")
	plt.show()


# plot_planck()


def plot_planck_approx():

	"""
	Plotting plancks law with the wien and rayleigh-Jeans approximations.

	"""

	wav = np.linspace(400,120080,5000)	
	B = np.zeros(len(wav))
	W = np.zeros(len(wav))
	R = np.zeros(len(wav))
	T = 10000

	B[:] = planck(T, wav[:]*1e-8)
	W[:] = wien(T, wav[:]*1e-8)
	R[:] = rayleighjeans(T, wav[:]*1e-8)


	plt.plot(wav,B, label = "Planck's law")
	plt.plot(wav,W, label = "Wien's law", ls ="--")
	plt.plot(wav,R, label = "Rayleigh-Jeans law", ls ="--")

	plt.ylim([8e10, 1e16])
	plt.yscale("log")
	plt.xscale("log")
	plt.title("Wien and Rayleigh-Jeans approximations (loglog plot)")
	plt.xlabel(r'$\lambda$ [$\AA$]')
	plt.ylabel(r"$B_\lambda$")
	plt.xlim([0,1e5])
	plt.legend()
	plt.grid(ls ="--")
	# plt.savefig("planckapproxlog.pdf")
	plt.show()


# plot_planck_approx()


def plot_optical_thickness():

	B = 2
	tau = np.linspace(0.01,10,100)
	intensity = np.zeros(len(tau))

	for I0 in range(4,-1,-1):
		intensity[:] = I0 * np.exp(-tau[:]) + B*(1-np.exp(-tau[:]))
		plt.plot(tau, intensity, label = 'I(0) = ' + str(I0))

	plt.yscale('log')
	plt.xscale('log')
	plt.title(r"Emergent Intensity for $B_\lambda = 2$")
	plt.xlabel(r'Optical depth $\tau$')
	plt.ylabel('Intensity')
	plt.grid(ls = "--")
	plt.legend()
	# plt.savefig("opticaldepthlog.pdf")
	plt.show()	
	
	
# plot_optical_thickness()


def plot_voigt():

	# u = np.arange(-10,10.1,0.1)
	u = np.linspace(-10,10,1000)
	a = np.array([0.001,0.01,0.1,1])
	vau = np.zeros((a.shape[0],u.shape[0]))

	for i in range(4):
		vau[i,:] = voigt(a[i],u[:])
		plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))

	plt.title("Voigt Profile for varying a")
	plt.ylim(0,1)
	plt.xlim(-10,10)
	plt.grid(ls = "--")
	plt.legend()
	plt.ylabel('Voigt profile')
	plt.xlabel("u")
	plt.savefig("voigt.pdf")
	plt.show()

	for i in range(4):
		vau[i,:] = voigt(a[i],u[:])
		plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))
	plt.title("Voigt Profile for varying a")
	plt.yscale('log')
	plt.grid(ls = "--")
	plt.legend()
	plt.xlabel('u')
	plt.ylabel('Logarithmic oigt profile')
	# plt.savefig("voigtlog.pdf")
	plt.show()


# plot_voigt()

def plot_spectral_lines(mode):	

	T_s = 5700 					#solar surface temperature
	T_l = 4200 					#solar T-min temperature - 'reversing layer'
	a = 0.1 						#damping parameter
	wav = 5000e-8 					#wavelength in cm
	tau0 = 1 						#reversing layer thickness at line center
	u  = np.linspace(-10,10,1000)
	intensity = np.zeros(len(u))	
	logtau0 = np.linspace(-2,2,9)

	if mode == "single":
		for i in range(len(u)):	

			tau = tau0 * voigt(a, u[i])
			intensity[i] = planck(T_s,wav) * np.exp(-tau) + planck(T_l,wav)*(1.-np.exp(-tau))	

		
		plt.plot(u,intensity, label = r"$\tau_0 = 1$")
		plt.grid(ls="--")
		plt.title(r"$T_\mathrm{surface} = 5700$ K, $T_\mathrm{layer} = 4200$ K, $a = 0.1$ , $\lambda = 5000$ Å", fontsize = 11)
		plt.xlabel("u")
		plt.legend()
		plt.ylabel(r"I$_\lambda$")
		# plt.savefig("ssline.pdf")
		plt.show()
	
	
	elif mode == "multi":

		for itau in range(len(logtau0)):
			for i in range(len(u)):
				tau = 10**(logtau0[itau]) * voigt(a, u[i])
				intensity[i] = planck(T_s,wav) * np.exp(-tau) + planck(T_l,wav)*(1-np.exp(-tau))
			plt.plot(u,intensity, label = r'$\tau_0 = %.2f$' %10**(logtau0[itau]))		

		plt.grid(ls="--")
		plt.title(r"$T_\mathrm{surface} = 5700$ K, $T_\mathrm{layer} = 4200$ K, $a = 0.1$ , $\lambda = 5000$ Å", fontsize = 11)
		plt.xlabel("u")
		plt.ylabel(r"I$_\lambda$")
		plt.legend()
		plt.savefig("sslines.pdf")
		plt.show()	
	
	elif mode == "var_wave":
		savefileindex = 0
		wavelist = [2000e-8, 5000e-8, 10000e-8 ]
		for wave in wavelist:

			for itau in range(len(logtau0)):
				for i in range(len(u)):	

					tau = 10**(logtau0[itau]) * voigt(a,u[i])
					intensity[i] = planck(T_s,wave) * np.exp(-tau) + planck(T_l,wave)*(1.-np.exp(-tau))	

				plt.plot(u,intensity, label = r'$\tau_0 = %.2f$' %10**(logtau0[itau]))
			plt.grid(ls="--")
			plt.title(r"$T_\mathrm{surface} = 5700$ K, $T_\mathrm{layer} = 4200$ K, $a = 0.1$ , $\lambda = %g$ Å" %(wave*1e8), fontsize = 11)
			plt.xlabel("u")
			plt.ylabel(r"I$_\lambda$")
			plt.legend()
			plt.savefig("ss"+str(savefileindex)+ ".pdf")
			plt.show()	
			savefileindex += 1
	
# plot_spectral_lines("single")
plot_spectral_lines("multi")
# plot_spectral_lines("var_wave")

	

