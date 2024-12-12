# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:33:23 2024

@author: angel
"""
import numpy as np
import math

def coeff_mod_HCO3(i, j, k, i1, j1, k1, D, geom):
    T = 298 #K
    Hen = 0.83 * math.exp(-20256.28 / 8.3144 * (1 / 298.15 - 1 / T)) # Henri constant molCO2l/molCO2g
    resist = 1 / (4.3e-8 * Hen) + 0.2e-6 / (1.17e-9 * Hen)
    dx = 0.66928e-6
    resist = resist / dx # resistance dividing dx
    
    if i == 0 or j == 0 or k == 0 or i == np.shape(geom)[0]+1 or j == np.shape(geom)[1]+1 or k == np.shape(geom)[2]+1:
        a = 0
    else:
        tam1 = D[i,j,k]
        tam2 = D[i1,j1,k1]   
        
        if tam1 == 0 or tam2 == 0:
            a = 0
        else:
            if tam1 != tam2:
                a = 1 / (0.5 / tam1 + 0.5 / tam2 + resist) # Adding the membrane resistance between 2 different materials
            else:
                a = 2 * tam1 * tam2 / (tam1 + tam2)       
    return a
