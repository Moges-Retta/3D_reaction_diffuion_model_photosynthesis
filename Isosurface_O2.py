# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:39:07 2024

@author: angel
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d, art3d
import matplotlib.pyplot as plt

# Load data
data = ...

# Find data equal to 2 and assign 1 in epid
epid = np.zeros(data.shape)
a = np.where(data == 2)
epid[a] = 1

# Find data equal to 1 and 3 and assign 1 in meso
meso = np.zeros(data.shape)
a = np.where(np.logical_or(data == 1, data == 3))
meso[a] = 1

# Find data equal to 8 and assign 1 in chloro
chloro = np.zeros(data.shape)
a = np.where(data == 8)
chloro[a] = 1

# Find data equal to 5 and assign 1 in bundle
bundle = np.zeros(data.shape)
a = np.where(data == 5)
bundle[a] = 1

# Find data equal to 9 and assign 1 in chloro_bundle
chloro_bundle = np.zeros(data.shape)
a = np.where(data == 9)
chloro_bundle[a] = 1

# Find data equal to 7 and assign 1 in vascular
vascular = np.zeros(data.shape)
a = np.where(data == 7)
vascular[a] = 1

# Find data equal to 6 and assign 1 in air
air = np.zeros(data.shape)
a = np.where(data == 6)
air[a] = 1

# Load air_stom
air_stom = ...
data = chloro + chloro_bundle

# Load u_O2 and assign to conc
conc = ...

# Rotate data and conc upside down
tam_data = np.zeros(data.shape)
tam_conc = np.zeros(conc.shape)
for i in range(data.shape[2]):
    tam_data[:, :, i] = data[:, :, data.shape[2] - i - 1]
    tam_conc[:, :, i] = conc[:, :, data.shape[2] - i - 1]
data = tam_data
conc = tam_conc

# Calculate x, y, z
x, y, z = np.meshgrid(np.arange(-1, data.shape[2] + 1), np.arange(-1, data.shape[1] + 1), np.arange(-1, data.shape[0] + 1))
x = x * 0.33464e-6 * 2
y = y * 0.33464e-6 * 2
z = z * 0.33464e-6 * 2

tam = np.zeros(x.shape)
tam[1:-1, 1:-1, 1:-1] = conc
conc = tam
tam[1:-1, 1:-1, 1:-1] = 1 - data
data = 1 - tam

conc *= (1 - data) * 8.314 * 298 / 101325 / 1000

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

vertex_data = np.array(mesh.points)
tri_indices = np.array(mesh.simplices)
triangles = Poly3DCollection(vertex_data[tri_indices], facecolor='interp', edgecolor='none')
ax.add_collection3d(triangles)

ax.set_xlim3d(-1, size(data,2) + 1)
ax.set_ylim3d(-1, size(data,1) + 1)
ax.set_zlim3d(-1, size(data,3) + 1)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.gca().set_aspect([1, 4, 4])
ax.view_init(elev=30, azim=45)
plt.tight_layout()

plt.show()
