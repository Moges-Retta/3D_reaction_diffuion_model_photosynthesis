# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:34:05 2024

@author: angel
"""
import numpy as np
import math

def coeff_mod_HCO3_chloro_vac(i, j, k, i1, j1, k1, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid):
    tcw_m = 0.161e-6 # mesophyll cell wall - N1A1 leaves
    tcw_b = 0.188e-6 # BS cell wall

    R = 8.314 # J/mol/K
    Hen_25 = 0.83
    Hen = Hen_25 * math.exp(-20256.28/R * (1/298.15 - 1/T)) #Henry constant

    D_HCO3_25 = 1.17e-9
    D_HCO3 = D_HCO3_25 * math.exp(16900/R * (1/298.15 - 1/T))

    P_HCO3 = 4.3e-8 * (D_HCO3/D_HCO3_25)

    resist = 1/(P_HCO3 * Hen) # resistance of the chloroplast membrane
    dx = 0.66928e-6
    resist = resist/dx # resistance dividing dx

    phi_pd = 0.03 # Area of PD per unit area of M_BS interface
    r_pd = (tcw_m + tcw_b)/(D_HCO3 * Hen * dx * phi_pd)

    a1 = coeff_mod_HCO3(i, j, k, i1, j1, k1, D, geom)
    if i == 0 or j == 0 or k == 0 or i == geom.shape[0]+1 or j == geom.shape[1]+1 or k == geom.shape[2]+1:
        a = a1
    else:
        tam1 = chloro[i][j][k]
        tam2 = chloro[i1][j1][k1]
        if tam1 != tam2:
            a = 1 / (1/a1 + 2 * resist) # Adding the mebrane resistence of chloroplast
        else:
            a = a1
  
        id1 = meso_ori[i][j][k] + epid[i][j][k] + air[i][j][k]
        id2 = meso_ori[i1][j1][k1] + epid[i1][j1][k1] + air[i][j][k]
    
        if id1 != id2:
            if id1 > 0 and id2 == 1 or id1 == 1 and id2 > 0:
                a = 1 / (1/a1 + 2 * resist) # Mc epid interface
            else:
                a = 0 # epid air
      
        # Adding the resitance to the chloroplast layers between two cells
        tam1 = meso_ori[i][j][k] * chloro[i][j][k]
        tam2 = meso_ori[i1][j1][k1] * chloro[i1][j1][k1]  
        if tam1 > 0 and tam2 > 0:
            if tam1 != tam2:
    a = 1 / (1 / a + 4 * resist) # 4 times resistance between two different chloroplasts

id3 = meso_bundle[i][j][k] # Meso bundle vascular
id4 = meso_bundle[i1][j1][k1]
if id3 != id4:
    if id3 == 0 or id4 == 0: # exposed mc or bs
        a = 0
    else:
        if id3 > 100:
            if id4 < 100:
                a = 1 / (1 / a + r_pd) # mc bs interfac
        if id3 < 100:
            if id4 > 100:
                a = 1 / (1 / a + r_pd) # mc bs interfac
        if id4 > 100:
            if id3 > 100:
                a = 1 / (1 / a1 + 2 * resist) # bs-bs and bs-vascu interface

# Adding the resistance to the cell membrane between the vacuole
tam1 = vacuole[i][j][k]
tam2 = vacuole[i1][j1][k1]
if tam1 != tam2:
    a = 1 / (1 / a1 + resist) # Adding the membrane resistance between 2 different materials

