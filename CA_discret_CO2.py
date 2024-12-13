# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:29:21 2024

@author: Moges Retta

Integrate the equation for the net hydration rate of CO2 with respect to CO2
% Calculate hydration rate and gradient of hydration at the current iteration
% Source terms are linearized using Picard's method of source term linearization.
% https://www.cfd-online.com/Wiki/Source_term_linearization

Inputs
CA : a geometrical domain where CA hydration is assumed
ka: turnover rate of carbonic anhydrase enzyme (1/s)
kC : Michaelis-Menten constants of carbonic anhydrase hydration (mol m-3)
keq : equilibrium constant for CA
kHCO3: Michaelis-Menten constants of carbonic anhydrase dehydration (mol m-3)
pro : concentration of H+ ions (mol/m3)
u0_CO2 : CO2 concentration at the current iteration (umol/m3)
u0_HCO3 : bicarbonate concentration at the current iteration(umol/m3)
Xa : concentratiion of CA (mol/m3)

Outputs
CA_r: net hydration of CO2
div_CA : first derivative of CA_r with respect to CO2

"""
def CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    Xa = CA * C_CA
    CA_r = -ka * Xa * (u0_CO2 - u0_HCO3 * pro / K_const) / (kC + kC / kHCO3 * u0_HCO3 + u0_CO2)
    div_CA = -ka * Xa * (kC + kC / kHCO3 * u0_HCO3 + u0_HCO3 * pro / K_const) / ((kC + kC / kHCO3 * u0_HCO3 + u0_CO2) ** 2)
    return CA_r, div_CA
