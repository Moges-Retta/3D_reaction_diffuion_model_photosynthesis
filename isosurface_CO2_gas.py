# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:37:45 2024

@author: angel
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_surface(data):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.add_collection(Poly3DCollection(data, alpha=.25, facecolor='#800000'))
    plt.show()

def make_surface(data, conc, thresh):
    x, y, z = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]), np.arange(data.shape[2]))
    x = x * 0.33464e-6 * 2
    y = y * 0.33464e-6 * 2
    z = z * 0.33464e-6 * 2
    x, y, z = [arr.astype(np.float32) for arr in (x, y, z)]
    tam = np.zeros(x.shape)
    tam[1:-1, 1:-1, 1:-1] = conc
    conc = tam
    tam[1:-1, 1:-1, 1:-1] = 1 - data
    data = 1 - tam
    data = (conc * (1 - data) * 8.314 * 298 / 101325) * (data < thresh)
    return np.array(np.nonzero(data)).T

data = np.load('data.npy')

epid = (data == 2).astype(int)
meso = ((data == 1) | (data == 3)).astype(int)
chloro = (data == 8).astype(int)
bundle = (data == 5).astype(int)
chloro_bundle = (data == 9).astype(int)
vascular = (data == 7).astype(int)
air = (data == 6).astype(int)
air_stom = np.load('air_stom.npy')

data = epid + meso + chloro + bundle + chloro_bundle + vascular + air + air_stom
data = np.flip(data, 2)
conc = np.load('u_CO2.npy')
conc = np.flip(conc, 2)

surface = make_surface(data, conc, 0)
plot_surface(surface)
