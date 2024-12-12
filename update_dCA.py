# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:49:11 2024

@author: angel
"""

def update_dCA(CA_r, CA_r0, dCA_HCO3, geom, a):
tam = (CA_r - CA_r0)
index = np.where(geom==1)
tam2 = tam[index] / CA_r[index]
tol = np.linalg.norm(tam2) / len(index[0])
if tol > 0.8:
dCA = np.zeros(len(tam))
return dCA
else:
dCA = a * tam
return dCA
