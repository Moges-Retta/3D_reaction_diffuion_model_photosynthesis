# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:49:58 2024

@author: angel
"""

import numpy as np

def rp_discret_O2(u0_CO2, u0_O2, p_para, Je):
    Vc = p_para[0]
    Kmc = p_para[1]
    gamma = p_para[2]
    KmO = p_para[3]
    
    index = np.where(u0_CO2 < 0)
    u0_CO2[index] = 0
    KM = Kmc * (1 + u0_O2 / KmO)
    
    rp_j = -Je * gamma * (u0_O2) / (3 * u0_CO2 + 7 * gamma * u0_O2)
    rp_c = -Vc * gamma * (u0_O2) / (u0_CO2 + KM)
    div_Sp_rpj = -(3 * u0_CO2 * gamma * Je) / ((3 * u0_CO2 + 7 * gamma * u0_O2) ** 2)
    div_Sp_rpc = -(gamma * Vc * (Kmc + u0_CO2)) / ((u0_CO2 + KM) ** 2)
    
    Ac = np.maximum(rp_c, rp_j)
    index = np.where((rp_c - rp_j) == 0)
    div_Sp_rpc[index] = div_Sp_rpj[index]
    
    rp = 2 * Ac
    div_rp = 2 * div_Sp_rpc
    
    return rp, div_rp
