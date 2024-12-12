# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:51:03 2024

@author: angel
"""


import numpy as np
import matplotlib.pyplot as plt

u_CO2 = np.load("u_CO2.npy")
data = np.load("data.npy")

conc = u_CO2
T = 298.15 #K
Hen = 0.83 * np.exp(-20256.28/8.3144 * (1/298.15 - 1/T)) #Henri constant, molCO2/L/molCO2g
index = np.where(data > 0)
conc[index] = conc[index] * Hen
min_val = np.min(conc)
max_val = np.max(conc)

for i in range(conc.shape[2]):
    plt.imshow(conc[:,:,i], vmin=min_val, vmax=max_val)
    title = "slice" + str(i)
    plt.title(title)
    plt.axis("equal")
    plt.colorbar()
    plt.pause(0.2)
    plt.clf()

