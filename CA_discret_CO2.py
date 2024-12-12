# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:29:21 2024

@author: angel
"""
def CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    Xa = CA * C_CA
    CA_r = -ka * Xa * (u0_CO2 - u0_HCO3 * pro / K_const) / (kC + kC / kHCO3 * u0_HCO3 + u0_CO2)
    div_CA = -ka * Xa * (kC + kC / kHCO3 * u0_HCO3 + u0_HCO3 * pro / K_const) / ((kC + kC / kHCO3 * u0_HCO3 + u0_CO2) ** 2)
    return CA_r, div_CA
