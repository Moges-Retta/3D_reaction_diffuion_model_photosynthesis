# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:38:24 2024

@author: angel
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# load data
data = ...
a = np.where(data == 2)[0]
epid = np.zeros_like(data)
epid[a] = 1

a = np.where(np.isin(data, [1, 3]))[0]
meso = np.zeros_like(data)
meso[a] = 1

a = np.where(data == 8)[0]
chloro = np.zeros_like(data)
chloro[a] = 1

a = np.where(data == 5)[0]
bundle = np.zeros_like(data)
bundle[a] = 1

a = np.where(data == 3)[0]
Bwall = np.zeros_like(data)
Bwall[a] = 1

a = np.where(data == 9)[0]
chloro_bundle = np.zeros_like(data)
chloro_bundle[a] = 1

a = np.where(data == 7)[0]
vascular = np.zeros_like(data)
vascular[a] = 1

a = np.where(data == 6)[0]
air = np.zeros_like(data)
air[a] = 1

data = epid + meso + chloro + bundle + chloro_bundle + vascular

import numpy as np 
import scipy.io
u_CO2 = scipy.io.loadmat('u_CO2') 
conc = u_CO2['u_CO2'] del u_CO2
tam_data = np.zeros(data.shape) 
tam_conc = np.zeros(conc.shape) 
for i in range(data.shape[2]): 
tam_data[:, :, i] = data[:, :, data.shape[2]-i-1]
tam_conc[:, :, i] = conc[:, :, data.shape[2]-i-1] 
data = tam_data 
conc = tam_conc 
del tam_data, tam_conc
x, y, z = np.meshgrid(np.linspace(-1, size(data, 2)-1, size(data, 2)), np.linspace(-1, size(data, 1)-1, size(data, 1)), np.linspace(-1, size(data, 3)-1, size(data, 3))) 
x = x * 0.33464e-6 * 2 y = y * 0.33464e-6 * 2 z = z * 0.33464e-6 * 2
tam = np.zeros(x.shape)
tam[1:data.shape[2]+1, 1:data.shape[1]+1, 1:data.shape[3]+1] = conc
conc = tam 
tam[1:data.shape[2]+1, 1:data.shape[1]+1, 1:data.shape[3]+1] = data 
data = tam
T = 298 Hen = 0.83 * np.exp(-20256.28 / 8.3144 * (1 / 298.15 - 1 / T)) 
conc = conc * Hen * 8.314 * 298 / 101325 
data = 1 - data
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt 
fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
hpatch = ax.plot_surface(x, y, z, data, rstride=1, cstride=1, facecolors=conc)
 ax.set_xlabel('(m)')
 ax.set_ylabel('(m)') 
ax.set_zlabel('(m)')
 ax.set_aspect([1, 4, 4]) 
ax.view_init(elev=30, azim=-37.5) 
plt.show()










Iso surface HCO3
import numpy as np
from scipy.ndimage import label
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# load data
data = np.load("data.npy")

# epid
epid = np.zeros(data.shape)
epid[np.where(data == 2)] = 1

# meso
meso = np.zeros(data.shape)
meso[np.where(data == 1)] = 1

# chloro
chloro = np.zeros(data.shape)
chloro[np.where(data == 8)] = 1

# bundle
bundle = np.zeros(data.shape)
bundle[np.where(data == 5)] = 1

# chloro_bundle
chloro_bundle = np.zeros(data.shape)
chloro_bundle[np.where(data == 9)] = 1

# vascular
vascular = np.zeros(data.shape)
vascular[np.where(data == 7)] = 1

data = epid + meso + chloro + bundle + chloro_bundle + vascular

# load conc
conc = np.load("u_HCO3.npy")

# generate x, y, z
x, y, z = np.meshgrid(
    np.arange(-1, data.shape[2]+1),
    np.arange(-1, data.shape[1]+1),
    np.arange(-1, data.shape[0]+1),
)
x = x * 0.33464e-6 * 2
y = y * 0.33464e-6 * 2
z = z * 0.33464e-6 * 2

# padding
tam = np.zeros(x.shape)
tam[1:data.shape[0]+1, 1:data.shape[1]+1, 1:data.shape[2]+1] = conc
conc = tam
tam[1:data.shape[0]+1, 1:data.shape[1]+1, 1:data.shape[2]+1] = data
data = tam

# define T
T = 298
Hen = 0.83 * np.exp(-20256.28 / 8.3144 * (1 / 298.15 - 1 / T))
conc = conc * Hen

data = 1 - data

# isosurface
def plot_isosurface(x, y, z, data, iso_value):
    vertices, faces, normals, values = measure.marching_cubes_lewiner(
        data, iso_value, spacing=(x[1, 0, 0] - x[0, 0, 0],
                                  y[0, 1, 0] - y[0, 0, 0],
                                  z[0, 0, 1] - z[0, 0, 0]))
    mesh = Poly3DCollection(vertices[faces], alpha=0.5)
    mesh.set_facecolor("red")
    ax.add_collection3d

