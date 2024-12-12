# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:49:38 2024

@author: angel
"""

def update_drp(rp_r, rp_r0, drp_O2, geom, a):
    tam = (rp_r - rp_r0)
    index = [i for i, x in enumerate(geom) if x == 1]
    tam2 = tam[index] / rp_r[index]
    tol = np.linalg.norm(tam2) / len(index)
    if tol > 0.2:
        drp = np.zeros_like(tam)
        tol
    else:
        drp = a * tam
    return drp
