# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:37:24 2024

@author: angel
"""

def h_mem_O2(u, T):
    R = 8.314
    Hen_O2 = 3.2e-2 * (T / 298.15) ** -1700
    Hen_CO2 = 0.83 * (T / 298.15) ** -20256.28 / R
    DO2 = 1.97e-9 * T / 293
    DCO2 = 1.67e-9 * T / 293
    Lmem = 8e-9
    porosity = 0.1
    hm = DO2 * porosity / Lmem
    return hm
