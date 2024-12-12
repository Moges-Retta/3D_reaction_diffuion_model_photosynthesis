# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:35:58 2024

@author: angel
"""

def diff_coeff_3m(Dg, Dl, D_HCO3_m, D_HCO3c, Depid, epid, meso, BS_cyto, chloro_bundle): 
D = D_HCO3c * chloro_bundle + D_HCO3_m * meso + Depid * epid + Dl * BS_cyto + Dg * (1 - (meso + epid + chloro_bundle + BS_cyto)) geom = epid + meso + chloro_bundle + BS_cyto 
tam = D 
return tam
