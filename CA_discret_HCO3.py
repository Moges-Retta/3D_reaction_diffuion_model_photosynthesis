# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:30:01 2024

@author: Moges Retta

Integrate the equation for the net hydration rate of CO2 with respect to bicarbonate

Inputs
CA : a geometrical domain where CA hydration is assumed
ka : turnover rate of carbonic anhydrase enzyme (1/s)
kC : Michaelis-Menten constants of carbonic anhydrase hydration (mol m-3)
k_eq : equilibrium constant for CA
k_HCO3 : Michaelis-Menten constants of carbonic anhydrase dehydration (mol m-3)
H : cocentration of H+ ions (mol/m3), pro
u0_CO2 : CO2 concentration (umol/m3)
u0_HCO3 : bicarbonate concentration (umol/m3)
Xa : concentratiion of CA (mol/m3)

Outputs
CA_r : net hydration of CO2
div_CA : first derivative of CA_r with respect to HCO3

"""
def CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    Xa = CA * C_CA
    CA_r = ka * Xa * (u0_CO2 - u0_HCO3 * pro / K_const) / (kC + kC / kHCO3 * u0_HCO3 + u0_CO2)
    div_CA = -ka * Xa * (pro * kC / K_const + pro * u0_CO2 / K_const + u0_CO2 * kC / kHCO3) / ((kC + kC / kHCO3 * u0_HCO3 + u0_CO2) ** 2)
    return CA_r, div_CA
