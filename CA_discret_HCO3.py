# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:30:01 2024

@author: angel
"""
def CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    Xa = CA * C_CA
    CA_r = ka * Xa * (u0_CO2 - u0_HCO3 * pro / K_const) / (kC + kC / kHCO3 * u0_HCO3 + u0_CO2)
    div_CA = -ka * Xa * (pro * kC / K_const + pro * u0_CO2 / K_const + u0_CO2 * kC / kHCO3) / ((kC + kC / kHCO3 * u0_HCO3 + u0_CO2) ** 2)
    return CA_r, div_CA
