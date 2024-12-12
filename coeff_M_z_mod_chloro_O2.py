# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:33:03 2024

@author: angel
"""
import numpy as np

def coeff_M_z_mod_chloro_O2(geom, T, U, D, d, R, chloro, air, vacuole, meso_bundle, u_O2, epid, meso_ori):
    U1 = U[0][0]
    U2 = U[0][1]
    dx = d[0][0]
    dy = d[0][1]
    dz = d[0][2]
    
    aE = np.zeros((geom.shape[0], geom.shape[1], geom.shape[2]))
    aW = np.zeros_like(aE)
    aS = np.zeros_like(aE)
    aN = np.zeros_like(aE)
    aT = np.zeros_like(aE)
    aP = np.zeros_like(aE)
    Su = np.zeros_like(aE)
    Sp = np.zeros_like(aE)
    
    for k in range(geom.shape[2]):
        for j in range(geom.shape[1]):
            for i in range(geom.shape[0]):
                Sp[i][j][k] = 0
                if k == 0:
                    Sp[i][j][k] = Sp[i][j][k] - 2 * coeff_mod_chloro_vacuole_O2(i, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dy / dz
                    Su_ijk = 2 * coeff_mod_chloro_vacuole_O2(i, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dy * dz / dx * U1
                elif k == geom.shape[2] - 1:
                    Sp[i][j][k] = Sp[i][j][k] - 2 * coeff_mod_chloro_vacuole_O2(i, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dy / dz
                    Su_ijk = 2 * coeff_mod_chloro_vacuole_O2(i, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dy / dz * U2
                else:
    Sp[i, j, k] = 0
    Su_ijk = 0
    aE[i, j, k] = coeff_mod_chloro_vacuole_O2(i - 1, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dy * dz / dx
    aW[i, j, k] = coeff_mod_chloro_vacuole_O2(i + 1, j, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dy * dz / dx
    aS[i, j, k] = coeff_mod_chloro_vacuole_O2(i, j - 1, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dz / dy
    aN[i, j, k] = coeff_mod_chloro_vacuole_O2(i, j + 1, k, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dz / dy
    aB[i, j, k] = coeff_mod_chloro_vacuole_O2(i, j, k - 1, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dy / dz
    aT[i, j, k] = coeff_mod_chloro_vacuole_O2(i, j, k + 1, i, j, k, D, T, geom, chloro, vacuole, air, meso_bundle, u_O2, epid, meso_ori) * dx * dy / dz

    aP[i, j, k] = aE[i, j, k] + aW[i, j, k] + aS[i, j, k] + aN[i, j, k] + \
        aB[i, j, k] + aT[i, j, k] - Sp[i, j, k]
    Su[i, j, k] = Su[i, j, k] + Su_ijk
Su = Su + R * dx * dy * dz * geom
