# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:36:57 2024

@author: angel
"""

def h_mem(u, T): 
R = 8.314 
D_CO2_25 = 1.81e-6 * np.exp(-16900 / (R * 298.15)) # Frank et al.1996 
D_CO2 = D_CO2_25 * np.exp(16900 / R * (1 / 298.15 - 1 / T)) 
Hen_CO2_25 = 0.83 
Hen_CO2 = Hen_CO2_25 * np.exp(-20256.28 / 8.3144 * (1 / 298.15 - 1 / T)) #Henri constant molCO2l/molCO2g
hm = 3.5e-3 * (D_CO2 / D_CO2_25) * (Hen_CO2 / Hen_CO2_25)
return hm
