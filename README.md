# Reaction diffusion model of CO2 during C4 photosynthesis

The Python code provided in this repository allows users to simulate the gas 
transport dynamics in a C4 leaf, and to study the effects of different environmental 
conditions on these processes. The code has been extensively tested and validated 
using experimental data from maize plants and can be easily adapted to other 
C4 and C3 plant species.

This coding manual is associated with the publication
"A reaction-diffusion model for photosynthesis gas transport in C4 plant maize" 
(DOI: https://doi.org/10.1093/jxb/erad138), which describes the development 
and validation of the model in detail. 


## Analyze results.py
Calculate the rate of photosynthesis, CO2 concentrations, and CO2 conductances

## diff_coeff_3m_airstom.py
The diffusion coefficients are assigned to respective cells and cell components.

## CA_discret_CO2.py
Integrate the equation for the net hydration rate of CO2 with respect to CO2

## CA_discret_HCO3.py
Integrate the equation for the net hydration rate of CO2 with respect to bicarbonate

## coeff_M_z_mod_chloro_O2.py
Determine the coefficient for the problem in the z-direction for O2

## coeff_mod_chloro_vacuole_O2.py
Determine the terms, Su and Sp, associated with the source terms expressing the resistance of due to chloroplast envelope, plasma membrane, tonoplast, cell wall, and plasmodesma

## discretise_CO2_fun.py
Determine coefficients K and F for the reaction-diffusion model of CO2 in the z-direction

## coeff_M_z_mod_chloro.py
Determine the coefficient for the problem in the z-direction for CO2

## coeff_M_z_mod_chloro.py
Determine the terms Sp and Su associated with the source terms expressing 
the resistance of CO2 due to chloroplast envelope, plasma membrane, tonoplast, 
cell wall, and the plasmodesma.

## coeff_M_z_mod_isolate_2bd_HCO3
Determine the coefficient for the problem in the z-direction for HCO3

## coeff_mod_HCO3_chloro_vac
Determine the terms Sp and Su associated with the source terms expressing 
the resistance of HCO3 due to chloroplast envelop, plasma membrane, tonoplast, and the plasmodesma.

## discretise_O2_fun.py
Determine coefficients K and F for the reaction-diffusion model in the z-direction

## discretise_HCO3
Determine coefficients K and F for the reaction-diffusion model of HCO3 in the z-direction

## discret
Combine the coefficients of the discretized equation for each control volume to form K and F matrices

## flux_calculation_z
calculate the flux of CO2 across the leaf thickness

## h_mem
Membrane permeability to CO2

## h_mem_O2
Membrane permeability to O2

## isosurface_CO2_gas
Make 3-D plots for CO2 profile in gas phase

## isosurface_CO2_liq
Make 3-D plots for CO2 profile in liq phase

## isosurface_O2_gas
Make 3-D plots for O2 profile in gas phase

## isosurface_O2_liq
Make 3-D plots for O2 profile in liq phase

## Main.py
The main program for solving the 3-D model

## Parameters.py
The model parameters consisting of diffusion coefficients, parameters of the photosynthesis model, 
and ATP production rates are defined here.

## photoresp_discret_J_Modified.py
Calculate photorespiration rate at the current iteration

## photo_PEP_J.py
Source terms of PEP carboxylation equation are linearized using Picard's method

## photo_fixation_J.py
Source terms of electron transport limited rate of photosynthesis are linearized using Picard's method

## rp_discret_CO2.py
Discretize the respiration term with respect to CO2

## rp_discret_O2.py
Discretize the respiration term with respect to O2