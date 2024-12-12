# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:50:15 2024

@author: angel
"""

def rp_discret_CO2(u0_CO2, u0_O2, p_para, Je):
    Vc = p_para[0]
    Kmc = p_para[1]
    gamma = p_para[2]
    KmO = p_para[3]
    # to avoid positive photosynthesis
    index = [i for i, x in enumerate(u0_CO2) if x < 0]
    for i in index:
        u0_CO2[i] = 0
    KM = Kmc * (1 + u0_O2 / KmO)
    
    rp_j = -Je * gamma * (u0_O2) / (3 * u0_CO2 + 7 * gamma * u0_O2)
    rp_c = -Vc * gamma * (u0_O2) / (u0_CO2 + KM)
    div_Sp_rpj = -(3 * u0_CO2 + 7 * gamma * u0_O2 - 7 * gamma * u0_O2 * Je) / ((3 * u0_CO2 + 7 * gamma * u0_O2) ** 2)
    div_Sp_rpc = -((u0_CO2 + KM) - gamma * Vc * u0_O2) / ((u0_CO2 + KM) ** 2)
    
    # Combine Ac and Aj
    Ac = [max(rp_c[i], rp_j[i]) for i in range(len(rp_c))]
    index = [i for i, x in enumerate(rp_c - rp_j) if x == 0]
    for i in index:
        div_Sp_rpc[i] = div_Sp_rpj[i]
    
    rp = 2 * Ac
    div_rp = 2 * div_Sp_rpc
    return rp, div_rp
