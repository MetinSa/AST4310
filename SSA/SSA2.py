import numpy as np 
import matplotlib.pyplot as plt 
import astropy.constants as cc
from cycler import cycler

# chi_ion = np.array([7, 16, 31, 51]) 		#schadee ionization energies into numpy array
# chi_ion = np.array([6, 12, 51, 67]) 		#Calcium ionization
k_eV = 8.61734e-5 							#boltzmann constant/eV
k_erg = 1.380658e-16						#boltzmann*erg
m_e = 9.10939e-28							#electron mass
h = 6.62607e-27								#plancks constant in erg


default_cycler = (cycler(color=['royalblue', 'crimson', 'mediumseagreen', 'darkorange', "darkorchid", "sienna"]))
plt.rc('axes', prop_cycle=default_cycler)

def partfunc_E(T):

	"""
	Partition function of the Shchadee "E" element.

	"""

	u = np.zeros(4)
	for r in range(4):
		for s in range(chi_ion[r]):

			u[r] = u[r] + np.exp(-s/(k_eV*T))

	return u 


def boltz_E(T, r, s):

	"""
	Boltzman function.

	"""

	u = partfunc_E(T)

	return  (1/(u[r-1])) * np.exp(-(s-1)/(k_eV*T))



def saha_E(T, e_Press, ionstage):

	"""
	Saha routine.

	"""
	
	u = partfunc_E(T)
	u = np.append(u, 2) 
	sahaconst = ((2. * np.pi * m_e * k_erg*T)/(h**2))**(3/2) * (2./(e_Press/(k_erg*T)))
	nstage = np.zeros(5)
	nstage[0] = 1.

	for r in range(4):
		nstage[r+1] = nstage[r] * sahaconst * (u[r+1]/u[r]) * np.exp(-chi_ion[r]/(k_eV*T))

	ntotal = np.sum(nstage)
	nstagerel = (nstage/ntotal)
	return nstagerel[ionstage-1]



def sahabolt_E(T, e_Press, r, s):

	"""
	Saha-Boltzmann routine.

	"""

	return saha_E(T, e_Press, r) * boltz_E(T, r, s)



def sahabolt_H(T, e_Press, level):

	"""
	Saha-Boltzmann for Hydrogen.

	"""

	nrlevels = 100
	g = np.zeros((2,nrlevels))
	chi = np.zeros((2,nrlevels))

	for s in range(nrlevels):
		g[0,s] = 2*(s+1)**2
		chi[0,s] = 13.598*(1- (1/(s+1)**2))	

	g[1,0] = 1.
	chi[1,0] = 0. 

	# partition functions
	u = np.zeros([2])
	for s in range(nrlevels):
		u[0] = u[0] + g[0,s]*np.exp(-chi[0,s]/(k_eV*T))
	u[1] = g[1,0]
	
	# Saha
	sahaconst = ((2. * np.pi * m_e * k_erg*T)/(h**2))**(3/2) * (2./(e_Press/(k_erg*T)))
	nstage = np.zeros(2)
	nstage[0] = 1.
	nstage[1] = nstage[0] * sahaconst * (u[1]/u[0]) * np.exp(-13.598/(k_eV*T))
	ntotal = np.sum(nstage)

	# Boltzmann
	nlevel = (nstage[0]*g[0,level-1]/u[0])*np.exp(-chi[0,level-1]/(k_eV*T))
	nlevelrel = nlevel/ntotal
	return nlevelrel



#plots

def part_func_plot():

	temp = np.arange(1,30001,1000)
	U = np.zeros((4,len(temp)))
	for i in range(len(temp)):
		# print(temp[i])
		U[:,i] = partfunc_E(temp[i])


	plt.plot(temp,U[0,:], label = "E I: $(r = 1)$")
	plt.plot(temp,U[1,:], label = "E II: $(r = 2)$")
	plt.plot(temp,U[2,:], label = "E III: $(r = 3)$")
	plt.plot(temp,U[3,:], label = "E IV: $(r = 4)$")
	plt.grid(linestyle = "--")
	plt.title("Partition function $U_r$ for Schadeenium E")
	plt.xlabel("Temperature [K]")
	plt.ylabel("$U_r$")
	plt.legend()
	# plt.savefig("part_E.pdf")
	plt.show()

	# print(partfunc_E(5000))
	# print(partfunc_E(10000))
	# print(partfunc_E(20000))

# part_func_plot()


def boltz_e_plot():

	temp = np.linspace(1,30000,100)
	smax = 6
	boltz = np.zeros((smax,100))
	lab = ["s = 1", "s = 2", "s = 3", "s = 4", "s = 5", "s = 6"]

	for s in range(0,smax):
		for i in range(len(temp)):
			boltz[s,i] = boltz_E(temp[i],1,(s+1))

	# print(boltz)

	# ground-state plot
	for i in range(0,smax):
		plt.semilogy(temp,boltz[i,:],label = lab[i])

	plt.title("Boltzmann distribution of Schadeenium E")
	plt.grid(linestyle = "--")
	plt.xlabel('Temperature [k]')
	plt.ylabel('$n_{r,s}/N_r$')
	plt.ylim([1e-5, 2])
	plt.legend()
	# plt.savefig("boltz_E.pdf")
	plt.show()


# boltz_e_plot()

def saha_E_plot():

	temp = np.linspace(1,30000,1000)
	rmax = 5
	saha = np.zeros((rmax,1000))
	lab = ["E I: $(r = 1)$", "E II: $(r = 2)$", "E III: $(r = 3)$", "E IV: $(r = 4)$", "E V: $(r = 5)$"]



	for r in range(0,rmax):
		for i in range(len(temp)):
			saha[r,i] = saha_E(temp[i],1000,(r+1))

	for i in range(0,rmax):
		plt.plot(temp,saha[i,:],label = lab[i])

	plt.title("Saha distribution of element E, $P_e = 1000$ dyn cm$^{-2}$")
	plt.grid(linestyle = "--")
	plt.legend(loc = "left")
	plt.xlabel('Temperature [k]')
	plt.ylabel('$N_{r+1}/N_r$')
	# plt.savefig("saha_E.pdf")
	plt.show()

# saha_E_plot()


def saha_E_plot_comp():

	temp = np.linspace(1,30000,1000)
	rmax = 5
	saha1 = np.zeros((rmax,1000))
	saha2 = np.zeros((rmax,1000))
	lab = ["E I: $(r = 1)$", "E II: $(r = 2)$", "E III: $(r = 3)$", "E IV: $(r = 4)$", "E V: $(r = 5)$"]



	for r in range(0,rmax):
		for i in range(len(temp)):
			saha1[r,i] = saha_E(temp[i],1000,(r+1))
			saha2[r,i] = saha_E(temp[i],10,(r+1))

	# print(boltz)

	plt.subplot(2,1,1)
	# ground-state plot
	for i in range(0,rmax):
		plt.plot(temp,saha1[i,:],label = lab[i])

	plt.title("Saha distribution of element E")
	plt.grid(linestyle = "--")
	plt.text(-1200, 0.45, '$P_e = 1000$', fontsize=12)
	plt.legend(loc = "left")
	plt.ylabel('$N_{r+1}/N_r$')


	plt.subplot(2,1,2)
	# ground-state plot
	for i in range(0,rmax):
		plt.plot(temp,saha2[i,:],label = lab[i])


	plt.grid(linestyle = "--")
	plt.xlabel('Temperature [k]')
	plt.ylabel('$N_{r+1}/N_r$')
	plt.text(-1200, 0.45, '$P_e = 10$', fontsize=12)
	# plt.ylim([1e-5, 2])
	# plt.savefig("saha_sub.pdf")
	plt.show()

# saha_E_plot_comp()



def payne_curves_E():

	temp = np.arange(0,30001,1000)
	#print T
	pop = np.zeros((5,31))
	for T in np.arange(1,31):
		for r in np.arange(1,5):
			pop[r,T] = sahabolt_E(temp[T],131.,r,1)
			labellst = ['E I: (r=1)', 'E II: (r=2)', 'E III: (r=3)', 'E IV: (r=4)']

	# ground-state plot
	for i in range(1,5):
		plt.plot(temp,pop[i,:], label=labellst[i-1])

	plt.title("Grund state populations of element E, $P_e = 131$ dyn cm$^{-2}$")
	plt.xlabel('Temperature [k]', size=14)
	plt.ylabel('Population', size=14)
	plt.yscale('log')
	plt.ylim([1e-3, 1.1])
	plt.grid(linestyle = "--")
	plt.legend(loc = "right")
	# plt.savefig("payne_E.pdf")
	plt.show()


# payne_curves_E()

def payne_curves_E_multi():

	temp = np.linspace(1,30000,200)
	pop = np.zeros((4,4,200))

	for T in np.arange(1,200):
		for r in range(4):
			for s in range(4):
				pop[r,s,T] = sahabolt_E(temp[T],131.,r+1,s+1)

	col = ['royalblue', 'crimson', 'mediumseagreen', 'darkorange']
	labellst = ['S = 1', 'S = 2', 'S = 3', 'S = 4']
	# # ground-state plot
	for i in range(4):
		for k in range(4):
			if i == 0:
				plt.plot(temp,pop[i,k,:], color = col[k], label= labellst[k])

			else:
				plt.plot(temp,pop[i,k,:], color = col[k])

	plt.title("Element E: Ionization stages r and excitation states s")
	plt.xlabel('Temperature [k]', size=14)
	plt.ylabel('Population', size=14)
	plt.yscale('log')
	plt.ylim([1e-4, 3])
	plt.grid(linestyle = "--")
	plt.legend(loc = "bottom right")
	plt.text(350, 1.2, 'E I: r=1', fontsize=12)
	plt.text(5350, 1.2, 'E II: r=2', fontsize=12)
	plt.text(12350, 1.2, 'E III: r=3', fontsize=12)
	plt.text(21350, 1.2, 'E IV: r=4', fontsize=12)
	# plt.savefig("payne_E_multi.pdf")
	plt.show()

# payne_curves_E_multi()


def sahabolt_H_plot():

	T = np.arange(1,20000,100)

	pop_H = np.zeros((4,len(T)))
	pop_Ca = np.zeros(len(T))
	for i in range(len(T)):
		for s in range(4):
			pop_H[s,i] = sahabolt_H(T[i],131.,s+1)
		pop_Ca[i] = sahabolt_E(T[i],1e2,2,1)

	labellst = ['H I, (s=1)', 'H I, (s=2)', 'H I, (s=3)', 'H I, (s=4)']

	# ground-state plot
	for i in range(4):
		plt.plot(T,pop_H[i,:],label = labellst[i] )
	plt.plot(T,pop_Ca,label = "Ca II, (s=1)", linestyle="--" , color = "gray")
	plt.title("Saha-Boltzmann distribution for H I and Ca$^+$ K")
	plt.xlabel('Temperature [k]')
	plt.ylabel('Population')
	plt.yscale('log')
	plt.ylim([1e-10, 2])
	plt.grid(linestyle = "--")
	plt.legend()
	# plt.savefig("sahabolt_HCA.pdf")
	plt.show()


# sahabolt_H_plot()

def plot_line_strength():

	# T = np.arange(1000,20001,100)
	T = np.linspace(1000,20000,1000)
	CaH = np.zeros(T.shape)
	Caabund = 2e-6
	for i in range(0,len(T)):
		NCa = sahabolt_E(T[i],1e2,2,1) # is equal to sahabolt_Ca
		NH = sahabolt_H(T[i],1e2,2)
		CaH[i] = (NCa*Caabund)/NH	

	plt.plot(T,CaH, label=r'Ca$^+$K / H$\alpha$')
	plt.title(r"Strength ratio between Ca$^+$K and H$\alpha$")
	plt.yscale('log')
	plt.grid(linestyle = "--")
	plt.xlabel(r'Temperature [K]')
	plt.ylabel(r'Ca II K / H$\alpha$',)
	plt.legend()
	# plt.savefig("strengthratio.pdf")
	plt.show()

# plot_line_strength()


def plot_T_sens():

	T = np.linspace(2000,12000,101)
	dNCadT = np.zeros(T.shape)
	dNHdT = np.zeros(T.shape)
	dT = 1.
	for i in range(len(T)):
		NCa = sahabolt_E(T[i],1e2,2,1)
		NCa2 = sahabolt_E(T[i]-dT,1e2,2,1)
		dNCadT[i] = (NCa - NCa2)/(dT*NCa)
		NH = sahabolt_H(T[i],1e2,2)
		NH2 = sahabolt_H(T[i]-dT,1e2,2)
		dNHdT[i] = (NH-NH2)/(dT*NH)
	NCa = np.zeros(T.shape)
	NH = np.zeros(T.shape)
	for i in range(len(T)):
		NCa[i] = sahabolt_E(T[i],1e2,2,1)
		NH[i] = sahabolt_H(T[i],1e2,2)


	plt.subplot(2,1,1)
	plt.title(r"Temperature sensitivity, Ca II vs H $\alpha$")
	plt.plot(T,np.absolute(dNHdT), label=r'H $\alpha$, (s=2)')
	plt.plot(T,np.absolute(dNCadT), label=r'Ca II K, (s=1)')
	plt.ylim(1e-7,1)
	plt.grid(ls ="--")
	plt.ylabel(r"$\left| \left( \Delta n(r,s) / \Delta T \right) /  n(r,s) \right|$")
	plt.legend(loc = "lower left")
	plt.yscale('log')

	plt.subplot(2,1,2)
	plt.plot(T,NH/np.amax(NH),  label = r'H $\alpha$')
	plt.plot(T,NCa/np.amax(NCa), label = 'Ca II K')
	plt.grid(ls ="--")	
	plt.xlabel('Temperature [K]')
	plt.ylabel("Population")
	plt.legend()
	# plt.savefig("tempsens.pdf")
	plt.show()	
	
	
# plot_T_sens()


def plot_hot_cool():
	T = np.linspace(1,2e4,10000)

	nH = np.zeros(T.shape)
	for i in range(len(T)):
		nH[i] = sahabolt_H(T[i],1e2,1)

	it = np.where(nH < 0.5 )[0][0]
	plt.title("Hot vs Cool stars")
	plt.plot(T,nH,label = "H I, fraction")
	plt.axvline(x = T[it], ls = "--", color ="black", alpha = 0.5, label = "Transition")
	plt.plot(T[it],nH[it],'o',color = "black", label="T = %d K" %T[it])
	plt.text(10000,0.55,"Hydrogen 50% ionized", fontsize = 14)
	plt.text(1200,0.1, "Cool Stars",fontsize = 13)
	plt.text(15000,0.1, "Hot Stars", fontsize = 13)
	plt.xlabel('Temperature [K]')
	plt.legend()
	plt.ylabel('Neutral hydrogen fraction')
	plt.grid(ls = "--")
	# plt.savefig("hotcool.pdf")
	plt.show()

	# print (T[it])

# plot_hot_cool()

