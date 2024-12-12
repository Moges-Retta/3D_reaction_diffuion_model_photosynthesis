# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:40:18 2024

@author: angel
"""


import numpy as np
# load data
data = ... # load data

# finding the indexes where data == 2 and marking those as epid
a = np.where(data == 2)
epid = np.zeros(data.shape)
epid[a] = 1

# finding the indexes where data == 1 or data == 3 and marking those as meso
a = np.where(data == 1)
meso = np.zeros(data.shape)
meso[a] = 1

a = np.where(data == 3)
meso[a] = 1

# finding the indexes where data == 8 and marking those as chloro
a = np.where(data == 8)
chloro = np.zeros(data.shape)
chloro[a] = 1

# finding the indexes where data == 5 and marking those as bundle
a = np.where(data == 5)
bundle = np.zeros(data.shape)
bundle[a] = 1

# finding the indexes where data == 3 and marking those as Bwall
a = np.where(data == 3)
Bwall = np.zeros(data.shape)
Bwall[a] = 1

# finding the indexes where data == 9 and marking those as chloro_bundle
a = np.where(data == 9)
chloro_bundle = np.zeros(data.shape)
chloro_bundle[a] = 1

# finding the indexes where data == 7 and marking those as vascular
a = np.where(data == 7)
vascular = np.zeros(data.shape)
vascular[a] = 1

# finding the indexes where data == 6 and marking those as air
a = np.where(data == 6)
air = np.zeros(data.shape)
air[a] = 1

data = chloro + chloro_bundle

# load u_O2
conc = ... # load u_O2

# rotating upside down
tam_data = np.zeros(data.shape)
tam_conc = np.zeros(conc.shape)
for i in range(data.shape[2]):
    tam_data[:, :, i] = data[:, :, data.shape[2]-i-1]
    tam_conc[:, :, i] = conc[:, :, data.shape[2]-i-1]
data = tam_data
conc = tam_conc

x, y, z = np.meshgrid(np.arange(-1, data.shape[2]+1), np.arange(-1, data.shape[1]+1), np.arange(-1, data.shape[0]+1))
x = x * 0.33464e-6 * 2
y = y * 0.33464e-6 * 2
z = z * 0.33464e-6 * 2

tam = np.zeros(x.shape)
tam[1:-1, 1:-1, 1:-1] = conc
conc = tam
tam[1:-1, 1:-1, 1:-1] = data
data = tam

T = 298 #°K
Hen = 0.83 * np.exp(-20256.28 / 8.3144 * (1/298.15 - 1/T)) # Henri constant molCO2l/molCO2g

conc = conc * 1 * 8.314 * 298 / 101325 / 1e6 * 100 # Expressing in %
data = 1 - data
hpatch = mlab.pipeline.surface(mlab.pipeline.iso_surface(x, y, z, scalars=data, contours=[0], vmax=conc))
# vertices

#set properties to obj hpatch
hpatch.actor.property.interpolation = 'phong'
hpatch.actor.property.opacity = 1

# Define the view
mlab.gcf().scene.parallel_projection = True
mlab.gcf().scene.camera.parallel_scale = 4
mlab.view(azimuth=-37.5, elevation=30)  # specifies the az=-37.5 (negative counter clockwise) and elevation to 30 above x, y plane
mlab.gcf().scene.light_manager.light_mode = 'vtk'

mlab.xlabel('(m)')
mlab.ylabel('(m)')
mlab.zlabel('(m)')







