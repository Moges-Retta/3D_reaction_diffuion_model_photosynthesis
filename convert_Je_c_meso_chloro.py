# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:35:39 2024

@author: angel
"""


import numpy as np

dx = 0.66928e-6
dy = 0.66928e-6
dz = 0.66928e-6

data = np.load("data.npy")
J_vol = np.load("J_vol.npy")
meso_ori = np.load("meso_ori.npy")

meso_new = meso_ori
in1 = np.where(data == 8)
chloro = np.zeros(data.shape)
chloro[in1] = np.amax(meso_new) + 1
data2 = chloro + meso_new
data2[in1] = np.amax(meso_new) + 1
Je_c_2 = np.zeros(data.shape)

for j in range(1, int(np.amax(meso_new)) + 1):
    temp = np.where(meso_new == j)
    M = np.zeros(data.shape)
    M[temp] = j
    
    temp2 = np.where(data2 == j)
    M_C = np.zeros(data.shape)
    M_C[temp2] = j
    
    Chl = M - M_C
    index = np.where(Chl == j)
    if index.size == 0:
        a = np.where(data == 8)
        JChl_c = np.sum(JChl) / len(a[0])
        Je_c_2[temp2] = JChl_c
    if index.size != 0:
        JChl = J_vol[index]
        temp3 = np.where(data[temp] == 1)
        cyto = temp3[0].size
        I = np.sum(JChl)
        Je_c_2[temp2] = np.sum(JChl) / cyto

a = np.where(data == 3)
Je_c_2[a] = 0
np.save("Je_c_2.npy", Je_c_2)

