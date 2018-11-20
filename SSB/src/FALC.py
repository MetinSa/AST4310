import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy.constants as const

# Chaning plot text
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size= 15)

# System specifics
path = "../data/"
savepath = "../figures/"
filename = "falc.dat"

# Reading file
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt(path + filename,usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)


"""
Constants/Parameters

"""

k_B = const.k_B.cgs.value
m_e = const.m_e.cgs.value
m_p = const.m_p.cgs.value
m_H = 1.67352e-24
m_He = 3.97*m_H
n_He = 0.1*nhyd

#Densities
dens_H = nhyd*m_H
dens_He = n_He*m_He
dens_Metals = dens - (dens_H + dens_He)
dens0 = dens[np.argwhere(h ==0)[0][0]]

#Pressure
p_idealgas = (nhyd+nel)*k_B*temp
p_idealgas_he = (nhyd+nel+n_He)*k_B*temp

# pgas = ptot-dens*vturb**2/2
# n_phot = 20*temp**3

"""
Functions

"""

def ScaleHeight(h, rho):

	"""
	Returns the scale height for at a given height and density

	"""

	return -h / (np.log(rho/dens0))


"""
Plotting Functions

"""

def Temperature_Height():

	"""
	Plotting atmosphere temperature for FALC model. 

	"""
	plt.title("FALC model of the solar atmosphere")
	plt.plot(h,temp, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlim(-500,2.5e3)
	plt.ylim(3e3,1e4)
	plt.xlabel("Height [km]")
	plt.ylabel("Temperature [K]")
	plt.tight_layout()
	plt.savefig(savepath + "FALCsolar_atmosphere.pdf")
	plt.show()


def TotalPressure_ColumnMass():

	"""
	Plotting pressure vs column mass both linearely and logarithmically

	"""
	plt.subplot(2,1,1)
	plt.title("Total pressure vs column mass")
	plt.plot(colm,ptot, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylabel(r"$P_\mathrm{tot}$ [dyn / cm$^2$]")

	plt.subplot(2,1,2)
	plt.loglog(colm,ptot, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Column mass [g / cm$^2$]")
	plt.ylabel(r"$P_\mathrm{tot}$ [dyn / cm$^2$]")
	plt.tight_layout()
	plt.savefig(savepath + "FALCpressure.pdf")
	plt.show()


def ColumnMass_Height():

	"""
	Plotting pressure vs column mass both linearely and logarithmically

	"""
	plt.subplot(2,1,1)
	plt.title("Column mass vs height")
	plt.plot(h,colm, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.ylabel("Column mass [g / cm$^2$]")

	plt.subplot(2,1,2)
	plt.semilogy(h,colm, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.ylabel("Column mass [g / cm$^2$]")
	plt.xlabel("Height [km]")
	# plt.tight_layout()
	plt.subplots_adjust(bottom = 0.12)
	plt.savefig(savepath + "FALCcolumnmass_height.pdf")
	plt.show()


def NumberDensity_Height():

	"""
	Plotting number density ratios of hydrogen, helium and metals

	"""
	plt.title("Comparing number densities in FALC")
	plt.semilogy(h,nel-nprot, color = "black", alpha = 0.7, linestyle = "--", label = r"$N_\mathrm{e} - N_\mathrm{p}$")
	plt.semilogy(h,nhyd, color = "royalblue", label = r"$N_\mathrm{H}$")
	plt.semilogy(h,nel, color = "mediumseagreen", label = r"$N_\mathrm{e}$")
	plt.semilogy(h,nprot, color = "crimson", label = r"$N_\mathrm{p}$")
	plt.grid(linestyle = "--")
	plt.legend()
	plt.xlabel("Height [km]")
	plt.ylabel("Number density [cm$^3$]")
	plt.xlim(-100, 2000)
	# plt.tight_layout()
	plt.savefig(savepath + "FALCnumberDensity.pdf")
	plt.show()


def gasDensity_Height():

	"""
	Plotting the gas density vs height along with the density heigh scale in the deep photosphere

	"""
	plt.title("Gas density vs height")
	plt.plot(h,dens, color = "royalblue")
	plt.grid(linestyle = "--")
	plt.xlabel("Height [km]")
	plt.ylabel(r" $\rho_\mathrm{gas}$ [g cm$^{-3}$]")
	plt.axvline(x=100, color = "black",linestyle = "--", alpha = 0.7, label = r"$H_P(h = 100$km, $\rho = $%.2E) $= $%.2f" %(dens[np.argwhere(h ==100)[0][0]], ScaleHeight(100,dens[np.argwhere(h ==100)[0][0]])))
	plt.legend()
	plt.savefig(savepath + "FALCgasDensity.pdf")
	plt.show()


def gasPressure_Height():

	"""
	Plotting the gas pressure vs height along with the product (n_H+n_el)k_B*T

	"""
	plt.figure()
	plt.title("Pressure vs height")
	plt.semilogx(h,ptot, color = "royalblue", label = r"$P_\mathrm{gas}$")
	plt.semilogx(h,p_idealgas, color = "crimson", label = r"$(n_\mathrm{H} + n_\mathrm{e})k_B T$")
	plt.semilogx(h,p_idealgas_he, color = "mediumseagreen", label = r"$(n_\mathrm{H} + n_\mathrm{He} + n_\mathrm{e})k_B T$")
	plt.grid(linestyle = "--")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel("Height [km]")
	plt.ylabel("Pressure [dyne cm$^{-2}$]")
	plt.legend()
	plt.savefig(savepath + "FALCgasPressure.pdf")
	# plt.show()

	plt.figure()
	plt.title("Deviations in gas pressure")
	plt.semilogx(h,p_idealgas/ptot, linestyle = "-",color = "crimson", label = r"$([n_\mathrm{H} + n_\mathrm{e}]k_B T)/P_\mathrm{gas}$")
	plt.semilogx(h,p_idealgas_he/ptot,linestyle = "-", color = "mediumseagreen", label = r"$([n_\mathrm{H} + n_\mathrm{He} n_\mathrm{e}]k_B T)/P_\mathrm{gas}$")
	plt.grid(linestyle = "--")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel("Height [km]")
	plt.ylabel("Pressure Ratios")
	plt.legend()
	plt.savefig(savepath + "FALCgasPressureDeviations.pdf")
	plt.show()


"""
Activating plots

"""
# Temperature_Height()
# TotalPressure_ColumnMass()
# ColumnMass_Height()
# NumberDensity_Height()
# gasDensity_Height()
gasPressure_Height()






# print ("Photon density at deepest model location: %g" % n_phot[np.argmin(h)])
# print ("Hydrogen density at deepest model location: %g" % nhyd[np.argmin(h)])
# print ("Photon density at highest model location: %g" % (20*5770**3/(2*np.pi)))
# print ("Hydrogen density at highest model location: %g" % nhyd[np.argmax(h)])

