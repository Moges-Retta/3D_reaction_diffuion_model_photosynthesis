# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:32:22 2024

@author: Moges Retta
Determine the coefficient for the problem in the z-direction for 

Inputs
air : geometrical domain of air
chloro : geometrical domain of M chloroplasts
epid : geometrical domain of epidermis
d : array containing voxel dimensions, dx,dy,dz
D : diffusion coefficients assigned to the appropriate domain
geom : geometrical domain of tissue
meso_bundle : geometrical domain of mesophyll and bundle-sheath cells
meso_ori : geometrical domain of mesophyll cells
U : external CO2 concentration (umol/m3)
para_HCO3 : set of parameters for bicarbonate diffusion
vacuole : geometrical domain of vacuole

Outputs
aB : coefficients of the flux leaving the bottom face, z direction
aE : coefficients of the flux leaving the east face, x direction
aN : coefficients of the flux leaving the north face, y direction
aP : coefficients of the flux at the node
aS : coefficients of the flux leaving the south face, y direction
aT : coefficients of the flux leaving the top face, z direction
aW : coefficients of the flux leaving the west face, x direction
Su : constant terms of the source term
Sp : variable terms of the source term

"""
import numpy as np
def coeff_M_z_mod_isolate_2bd_HCO3(geom, T, D, d, R, chloro, vacuole, meso_bundle, meso_ori, air, epid):
    dx, dy, dz = d[0][0], d[0][1], d[0][2]
    aE, aW, aS, aN, aT, aP, Su, Sp = [np.zeros_like(geom) for _ in range(8)]
    
    for k in range(geom.shape[2]):
        for j in range(geom.shape[1]):
            for i in range(geom.shape[0]):
                Sp[i, j, k] = 0
                if k == 0 or k == geom.shape[2]-1:
                    Su_ijk = 0
                else:
                    Su_ijk = 0
                aE[i, j, k] = coeff_mod_HCO3_chloro_vac(i-1, j, k, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dy * dz / dx
                aW[i, j, k] = coeff_mod_HCO3_chloro_vac(i+1, j, k, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dy * dz / dx
                aS[i, j, k] = coeff_mod_HCO3_chloro_vac(i, j-1, k, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dx * dz / dy
                aN[i, j, k] = coeff_mod_HCO3_chloro_vac(i, j+1, k, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dx * dz / dy
                aB[i, j, k] = coeff_mod_HCO3_chloro_vac(i, j, k-1, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dx * dy / dz
                aT[i, j, k] = coeff_mod_HCO3_chloro_vac(i, j, k+1, i, j, k, T, D, geom, chloro, vacuole, meso_ori, meso_bundle, air, epid) * dx * dy / dz
                aP[i, j, k] = aE[i, j, k] + aW[i, j, k] + aS[i, j, k] + aN[i, j, k] + aB[i, j, k] + aT[i, j, k] - Sp[i, j,k]
Su [i, j, k] = Su[i, j, k] +Su_ijk; 
Su=Su+R*dx*dy*dz*geom;
