# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:36:41 2024

@author: angel
"""


u_CO2 = np.load("u_CO2.npy") 
aB = np.load("aB.npy") 
data = np.load("data.npy")
u_CO2 = np.reshape(u_CO2, np.shape(data))
delta_u_z = u_CO2[:, :, 1:] - u_CO2[:, :, :-1] #Concentration difference
flux_z = delta_u_z * aB[:, :, 1:] #flux on each of elements

total_flux_z = [] #total flux on x direction 
for k in range(1, np.shape(data)[2]):
    total_flux_z.append(np.sum(np.sum(flux_z[:, :, k])))
dx = 0.66928e-6 #Scale in X direction
dy = 0.66928e-6 #Scale in Y direction
dz = 0.66928e-6 #Scale in Z direction
#Leaf area
Area = dx * dy * np.shape(data)[0] * np.shape(data)[1]

#flux mol/m2s
Flux = total_flux_z / Area

Flux_upper = -Flux[0]
Flux_lower = Flux[-1]
photo = Flux_upper + Flux_lower
return photo
