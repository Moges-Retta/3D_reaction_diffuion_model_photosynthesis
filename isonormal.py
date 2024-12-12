# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:39:38 2024

@author: angel
"""


from skimage import measure
import numpy as np

# Creating the isosurface using marching cubes
verts, faces, normals, values = measure.marching_cubes_lewiner(data, 0, step_size=1, allow_degenerate=False)

# Plotting the isosurface
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
mesh = Poly3DCollection(verts[faces])
mesh.set_edgecolor('k')
ax.add_collection3d(mesh)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

