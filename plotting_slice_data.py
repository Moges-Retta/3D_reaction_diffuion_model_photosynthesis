# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:50:42 2024

@author: angel
"""


import numpy as np
import matplotlib.pyplot as plt

data = np.load("data.npy")

for i in range(data.shape[0]):
    tam = np.zeros((data.shape[1], data.shape[2]))
    tam[:,:] = data[i, :, :]
    plt.imshow(tam)
    plt.title("slice" + str(i + 1))
    plt.axis("equal")
    plt.clim(0, 9)
    plt.colorbar()
    plt.pause(0.2)
    plt.clf()

