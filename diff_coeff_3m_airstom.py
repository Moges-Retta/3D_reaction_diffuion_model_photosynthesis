# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:36:19 2024

@author: angel
"""

def diff_coeff_3m_airstom(Dg, Dl, Dl_m, Dc, Depid, D_st, epid, meso, BS_cyto, chloro_bundle, air_stom): 
D = Dl_m * meso + Dc * chloro_bundle + Dl * BS_cyto + Depid * epid + D_st * air_stom + Dg * (1 - (meso + epid + air_stom + chloro_bundle + BS_cyto)) 
geom = epid + meso + chloro_bundle + BS_cyto 
tam = D 
return tam
