# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:41:09 2024

@author: angel
"""

def photo_PEP_J(u0_HCO3, p_para, Je):
    Vp, Kp = p_para

    index = [i for i, val in enumerate(u0_HCO3) if val < 0]
    for i in index:
        u0_HCO3[i] = 0

    Vpj = [-Je/2 + u*0 for u in u0_HCO3]
    Vpc = [-Vp*u/(u+Kp) for u in u0_HCO3]

    div_Sp_Vpj = [0 for _ in u0_HCO3]
    div_Sp_Vpc = [-Vp*Kp/((u+Kp)**2) for u in u0_HCO3]

    for i, (vpc, vpj) in enumerate(zip(Vpc, Vpj)):
        if vpc > vpj:
            Vpc[i] = vpc
        else:
            Vpc[i] = vpj
            div_Sp_Vpc[i] = div_Sp_Vpj[i]

    PEP = Vpc
    div_PEP = div_Sp_Vpc

    return PEP, div_PEP
