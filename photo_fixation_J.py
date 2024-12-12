# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:42:22 2024

@author: angel
"""

import numpy as np

def photo_fixation_J(u0_CO2, u0_O2, p_para, Je):
    """
    Calculate photosynthesis rates and their derivatives.

    Parameters:
    u0_CO2 : numpy array
        Initial CO2 concentration.
    u0_O2 : numpy array
        Initial O2 concentration.
    p_para : list or numpy array
        Parameter array containing [Vc, Kmc, gamma, KmO].
    Je : float
        Electron transport rate.

    Returns:
    photo : numpy array
        Photosynthesis rate.
    div_photo : numpy array
        Derivative of photosynthesis rate.
    """
    # Extract parameters
    Vc, Kmc, gamma, KmO = p_para

    # Ensure non-negative concentrations
    u0_CO2 = np.maximum(u0_CO2, 0)
    u0_O2 = np.maximum(u0_O2, 0)

    # Calculate KM
    KM = Kmc * (1 + u0_O2 / KmO)

    # Calculate Aj and Ac
    Aj = -Je * u0_CO2 / (3 * u0_CO2 + 7 * gamma * u0_O2)
    Ac = -Vc * u0_CO2 / (u0_CO2 + KM)

    # Calculate derivatives
    div_Sp_Aj = -Je * 7 * gamma * u0_O2 / ((3 * u0_CO2 + 7 * gamma * u0_O2) ** 2)
    div_Sp_Ac = -Vc * KM / ((u0_CO2 + KM) ** 2)

    # Combine Ac and Aj
    Ac = np.maximum(Ac, Aj)
    index = np.where((Ac - Aj) == 0)
    div_Sp_Ac[index] = div_Sp_Aj[index]

    # Output
    photo = Ac
    div_photo = div_Sp_Ac

    return photo, div_photo

# Example usage
if __name__ == "__main__":
    u0_CO2 = np.array([10, 20, -5, 30])
    u0_O2 = np.array([15, -10, 5, 25])
    p_para = [100, 200, 0.5, 300]
    Je = 50

    photo, div_photo = photo_fixation_J(u0_CO2, u0_O2, p_para, Je)
    print("Photo:", photo)
    print("Div Photo:", div_photo)
