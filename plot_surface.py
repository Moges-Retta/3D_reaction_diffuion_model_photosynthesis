# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:48:33 2024

@author: angel
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_surface(data):
    verts, faces, normals, values = measure.marching_cubes_lewiner(data, 0)
    ax = plt.gca(projection='3d')
    ax.add_collection(Poly3DCollection(verts[faces], 
                      facecolor='red', edgecolor='k', alpha=0.2))
    ax.set_xlim(0, data.shape[0])
    ax.set_ylim(0, data.shape[1])
    ax.set_zlim(0, data.shape[2])
    plt.show()
