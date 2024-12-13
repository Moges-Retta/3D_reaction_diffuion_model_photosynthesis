# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:33:03 2024

@author: Moges Retta
Determine the coefficient for the problem in the z-direction for O2

Input
air : geometrical domain of air
chloro : geometrical domain of M chloroplasts
epid : geometrical domain of epidermis
d : array containing voxel dimensions, dx,dy,dz
D : diffusion coefficients assigned to the appropriate domain
geom : geometrical domain of tissue
meso_bundle : geometrical domain of mesophyll and bundle-sheath cells
meso_ori : geometrical domain of mesophyll cells
U : external CO2 concentration (umol/m3)
para_O2 : set of parameters for O2 diffusion
vacuole : geometrical domain of vacuole

Output
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
import coeff_mod_chloro_vacuole

def coeff_M_z_mod_chloro(geom, U, D, d, chloro, air, vacuole, meso_bundle, epid, meso_ori, para_CO2):
    dx, dy, dz = d[0, 0], d[0, 1], d[0, 2]

    # Coefficient initialization
    shape = geom.shape
    aE = np.zeros(shape)
    aW = np.zeros(shape)
    aS = np.zeros(shape)
    aN = np.zeros(shape)
    aT = np.zeros(shape)
    aP = np.zeros(shape)
    Su = np.zeros(shape)
    Sp = np.zeros(shape)
    aB=aE;
    
    for k in range(shape[2]):
        for j in range(shape[1]):
            for i in range(shape[0]):

                Sp[i, j, k] = 0

                if k == 0:
                    Sp[i, j, k] -= 2 * coeff_mod_chloro_vacuole(i, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dy / dz
                    Su_ijk = 2 * coeff_mod_chloro_vacuole(i, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dy * dz / dx * U[0]
                elif k == shape[2] - 1:
                    Sp[i, j, k] -= 2 * coeff_mod_chloro_vacuole(i, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dy / dz
                    Su_ijk = 2 * coeff_mod_chloro_vacuole(i, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dy / dz * U[1]
                else:
                    Sp[i, j, k] = 0
                    Su_ijk = 0

                aE[i, j, k] = coeff_mod_chloro_vacuole(i - 1, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dy * dz / dx
                aW[i, j, k] = coeff_mod_chloro_vacuole(i + 1, j, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dy * dz / dx
                aS[i, j, k] = coeff_mod_chloro_vacuole(i, j - 1, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dz / dy
                aN[i, j, k] = coeff_mod_chloro_vacuole(i, j + 1, k, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dz / dy
                aB[i, j, k] = coeff_mod_chloro_vacuole(i, j, k - 1, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dy / dz
                aT[i, j, k] = coeff_mod_chloro_vacuole(i, j, k + 1, i, j, k, D, geom, chloro, vacuole, air, meso_bundle, epid, meso_ori, para_CO2, dx) * dx * dy / dz

                aP[i, j, k] = aE[i, j, k] + aW[i, j, k] + aS[i, j, k] + aN[i, j, k] + aB[i, j, k] + aT[i, j, k] - Sp[i, j, k]
                Su[i, j, k] += Su_ijk

    return aP, aE, aW, aS, aN, aB, aT, Su, Sp