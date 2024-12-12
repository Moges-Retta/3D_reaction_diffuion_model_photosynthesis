# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:34:46 2024

@author: angel
"""

import math
def coeff_mod_chloro_vacuole(i, j, k, i1, j1, k1, T, D, geom, chloro, vacuole, air, meso_bundle, u_CO2, epid, meso_ori):
    tcw_m = 0.161e-6
    tcw_b = 0.188e-6
    R = 8.314
    Hen_CO2_25 = 0.83
    Hen = Hen_CO2_25 * math.exp(-20256.28 / R * (1 / 298.15 - 1 / T))
    dx = 0.66928e-6
    D_CO2_25 = 1.81e-6 * math.exp(-16900 / (R * 298.15))
    D_CO2 = D_CO2_25 * math.exp(16900 / R * (1 / 298.15 - 1 / T))
    r_cw = tcw_m / (D_CO2 * Hen * dx * 1.0)
    r_cw2 = (tcw_b + tcw_b) / (D_CO2 * Hen * dx * 0.1)
    r_cw3 = tcw_b / (D_CO2 * Hen * dx * 1e-12)
    phi_pd = 0.03
    r_pd = (tcw_b + tcw_m) / (D_CO2 * Hen * dx * phi_pd)
    h_epid = 6.99e-6 * (D_CO2 / D_CO2_25) * (Hen / Hen_CO2_25)
if (i == 0) or (j == 0) or (k == 0) or (i == len(geom) - 1) or (j == len(geom[0]) - 1) or (k == len(geom[0][0]) - 1):
    a = 0
else:
    D1 = D[i][j][k]
    D2 = D[i1][j1][k1]
    a = 2 * D1 * D2 / (D1 + D2)

    # Adding the cuticle resistence between the epidermis and air
    index1 = air[i][j][k]
    index2 = air[i1][j1][k1]
    id1 = geom[i][j][k]
    id2 = geom[i1][j1][k1]
    if index1 != index2:
        if id1 != id2:
            a = dx / (dx / a + 1 / h_epid)

    # Adding the resistance to membrane and cellwall to the cell of bs and
    # mc
    id1 = meso_bundle[i][j][k]
    id2 = meso_bundle[i1][j1][k1]

    hm = h_mem(u_CO2[i][j][k], T)
    hm1 = h_mem(u_CO2[i1][j1][k1], T)
    h_cm_tam = min(hm, hm1)
    h_cm_tam2 = 1.6e-2
    r_cm = 1 / (h_cm_tam * Hen * dx)  # bundle sheath plasma membrane
    r_cm2 = 1 / (h_cm_tam2 * Hen * dx)  # mesophyll plasma membrane
    if id1 != id2:
        if id1 == 0 or id2 == 0:
            if id1 > 100 or id2 > 100:  # bs exposed to air
                a = 1 / (1 / a + r_cw3 + r_cm)
            elif id1 < 100 or id2 < 100:  # mc exposed to air
                a = 1 / (1 / a + r_cw + r_cm2)
        else:
            if id1 > 100:
                if id2 > 100:
                    a = 1 / (1 / a + r_cw2 + r_cw2 + r_cm + r_cm)  # bs-bs and bs-vascu interface
                else:
                    a = 1 / (1 / a + r_pd)  # mc bs interface
            if id1 < 100:
                if id2 > 100:
                    a = 1 / (1 / a + r_pd)  # mc bs interface
                else:
                    a = 1 / (1 / a + r_cw + r_cw + r_cm2 + r_cm2)  # mc mc interface
    id1 = meso_ori(i, j, k) + epid(i, j, k) + air(i, j, k)
    id2 = meso_ori(i1, j1, k1) + epid(i1, j1, k1) + air(i, j, k)

    if id1 != id2:
        if (id1 == 0 and id2 == 1) or (id1 == 1 and id2 == 0):
            a = 1 / (1 / a + r_cw + r_cm)  # pore epid interface
        elif (id1 > 0 and id2 == 1) or (id1 == 1 and id2 > 0):
            a = 1 / (1 / a + r_cw + r_cm + r_cw + r_cm2)  # Mc epid interface

    tam1 = chloro(i, j, k)
    tam2 = chloro(i1, j1, k1)
    if tam1 != tam2:
        if tam1 > 0:
            hm = h_mem(u_CO2(i, j, k), T)
        if tam2 > 0:
            hm = h_mem(u_CO2(i1, j1, k1), T)
        r_mem = 1 / (hm * Hen * dx)
        a = 1 / (1 / a + 2 * r_mem)

    tam1 = meso_ori(i, j, k) * chloro(i, j, k)
    tam2 = meso_ori(i1, j1, k1) * chloro(i1, j1, k1)
    if tam1 > 0 and tam2 > 0:
        if tam1 != tam2:
            hm = h_mem(u_CO2(i, j, k), T)
            hm1 = h_mem(u_CO2(i1, j1, k1), T)
            r_mem1 = 1 / (hm * Hen * dx)
            r_mem2 = 1 / (hm1 * Hen * dx)
            a = 1 / (1 / a + 2 * r_mem1 + 2 * r_mem2)

    tam1 = vacuole(i, j, k)
    tam2 = vacuole(i1, j1, k1)
    if tam1 != tam2:
        if tam1 > 0:
            hm = h_mem(u_CO2(i, j, k), T)
        if tam2 > 0:
            hm = h_mem(u_CO2(i1, j1, k1), T)
        r_mem = 1 / (hm * Hen * dx)
        a = 1 / (1 / a + r_mem)

    return a

