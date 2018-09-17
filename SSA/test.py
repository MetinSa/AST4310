import numpy as np 
# nrlevels = 100
# g = np.zeros((2,nrlevels))
# chi = np.zeros((2,nrlevels))
# for s in range(nrlevels):
# 	g[0,s] = 2*(s+1)**2
# 	chi[0,s] = 13.598*(1- (1/(s+1)**2))	
# g[1,0] = 1.
# chi[1,0] = 0. 

ratio = np.zeros(3)
g = [2, 8, 18, 32]
chi = [0, 10.20, 12.09, 12.75]
k_eV = 8.61734e-5
k_erg = 1.380658e-16
T = 5777

for i in range(3):

	ratio[i] = (g[i] * np.exp(-chi[i]/(k_eV*T)))/(g[i+1] * np.exp(-chi[i+1]/(k_eV*T)))

# print(ratio)

print((g[1] * np.exp(-chi[1]/(k_eV*T)))/(g[3] * np.exp(-chi[3]/(k_eV*T))))