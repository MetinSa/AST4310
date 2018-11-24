import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Matplotlib estetics
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath = "../figures/3.2/"
filename1 = "int_nad.dat"
filename2 = "falc.dat"


# Reading data
wavenumber, E_spec, S_spec, S_spec_corr = np.loadtxt(path + filename1, usecols=(0,1,2,3), unpack=True)
height, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt(path + filename2,usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)

wave_vac = 1/wavenumber*1e8

# Constants
c = const.c.cgs.value
k_erg = const.k_B.cgs.value
keV = 8.61734E-5
erg2eV = 6.242e11
h = const.h.cgs.value
m_e = 9.109390E-28


def VacToAir(wavelength):

	"""
	Shifting wavelengths from vacuum to air

	"""
	return 0.99972683*wavelength + 0.0107 - 196.25/wavelength


def partfunc_Na(temp):
	# partition functions Na
	# input: temp (K)
	# output: float array(3) = partition functions U1,U2,U3
	u=np.zeros(3)
	# partition function Na I: follow Appendix D of Gray 1992
	# log(U1(T)) = c0 + c1 * log(theta) + c2 * log(theta)^2 +
	# c3 *log(theta)^3 + c4 log(theta)^4
	# with theta=5040./T
	theta=5040./temp
	# partition function Na I : Appendix D of Gray (1992)
	c0=0.30955
	c1=-0.17778
	c2=1.10594
	c3=-2.42847
	c4=1.70721
	logU1 = (c0 + c1 * np.log10(theta) + c2 * np.log10(theta)**2 +
	c3 * np.log10(theta)**3 + c4 * np.log10(theta)**4)
	u[0]=10**logU1
	# partition function Na II and Na III: approximate by the
	# statistical weights of the ion ground states
	u[1]=1 # from Allen 1976
	u[2]=6 # from Allen 1976
	return u

def boltzmann_na(temp, r, s):

	"""
	Computing the Boltzmann distribution for Na

	"""
	e_n1 = h*c / 5895.94e-8 * erg2eV
	e_n2 = h*c / 5889.97e-8 * erg2eV
	u = partfunc_Na(temp)
	chi = [0, e_n1, e_n2]
	g = [2, 2, 4]
	relnrs = g[s]/u[r]*np.exp(-(chi[s])/(keV*temp))
	return relnrs


def saha_na(temp, el_dens, ion_stage, e_ionization):

	"""
	Computing the Saha distribution for Na

	"""
	kevt = keV*temp
	kergt = k_erg * temp
	u = partfunc_Na(temp)
	u = np.append(u, 2)  # append element to array
	saha_const = (2.*np.pi*m_e*kergt/(h*h)) ** (3./2) * 2. / el_dens
	n_stage = np.zeros(4)
	n_stage[0] = 1.
	for r in range(3):
	    n_stage[r+1] = n_stage[r]*saha_const*u[r+1]/u[r] * np.exp(-e_ionization[r] / kevt)
	n_total = np.sum(n_stage)
	n_stage_rel = n_stage/n_total
	return n_stage_rel[ion_stage - 1]


def PlotSpectrum():

	"""
	Plotting the solar spectra

	"""
	plt.title(r"Solar NA I D lines ")
	plt.plot(wave_vac,S_spec, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$ [$\mathrm{\AA}$]")
	# plt.subplots_adjust(bottom = 0.12, left = 0.15)
	# plt.savefig(savepath + "solarspectran.pdf")
	plt.show()


def PlotSpectrumAir():

	"""
	Plotting the solar spectra seen through air

	"""

	wave_air = VacToAir(wave_vac)


	plt.title(r"Solar NA I D lines observed in air ")
	plt.plot(wave_air,S_spec, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel(r"Wavelength $\lambda$ [$\mathrm{\AA}$]")
	# plt.subplots_adjust(bottom = 0.12, left = 0.15)
	# plt.savefig(savepath + "solarspectranair.pdf")
	plt.show()


def PlotBoltz():

	"""
	Plotting the Boltzmann distribution

	"""

	b_l = 1.
	b_u = 1.
	A_Na = 1.8*1e-6
	f_lu = [0.318, 0.631]
	E_ionization = np.array([5.139, 47.29, 71.64])
	boltz = np.zeros((3, len(temp)))	

	for i in range(len(temp)):
	    boltz[0, i] = boltzmann_na(temp[i], 0, 0)
	    boltz[1, i] = boltzmann_na(temp[i], 0, 1)
	    boltz[2, i] = boltzmann_na(temp[i], 0, 2)	

	plt.plot(height, boltz[0], color = "royalblue", label='Na I s = 1: ground state')
	plt.plot(height, boltz[1], color = "crimson", label='Na I s = 2: level D1')
	plt.plot(height, boltz[2], color = "mediumseagreen", label='Na I s = 3: level D2')
	plt.title(r'Boltzmann distribution of Na I in Falc')
	plt.xlabel(r'Height [km]')
	plt.grid(linestyle="--")
	plt.ylabel(r'Population fraction $n_{1,s}/N_1$')
	plt.legend(loc='center')
	# plt.savefig(savepath + "boltz.pdf")
	plt.show()


def PlotSaha():

	"""
	Plotting the Saha distribution

	"""
	E_ionization = np.array([5.139, 47.29, 71.64])

	saha = np.zeros((2, len(temp)))
	for i in range(len(temp)):
	    saha[0, i] = saha_na(temp[i], nel[i], 1, E_ionization)
	    saha[1, i] = saha_na(temp[i], nel[i], 2, E_ionization)	

	plt.figure(37)
	plt.plot(height, saha[0], color = "royalblue", label='Na I')
	plt.plot(height, saha[1], color = "crimson", label='Na II')
	plt.xlim([np.min(height), 2000])
	plt.ylim([1e-4, 10])

	plt.yscale('log')
	plt.title(r'Saha distribution of Na in FALC')
	plt.xlabel(r'Height [km]')
	plt.ylabel(r'Ionization state fraction $N_r/N_\mathrm{total}$')
	plt.legend(loc='best')
	# plt.savefig(savepath + "saha.pdf")

	plt.show()



"""
Activating functions

"""
# PlotSpectrum()
# PlotSpectrumAir()
# PlotBoltz()
# PlotSaha()