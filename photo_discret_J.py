# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:41:34 2024

@author: angel
"""


def photo_discret_J(u0_CO2, p_para, Je):
    Rd = p_para[1]
    gamma = p_para[2]
    Hen = p_para[3]
    U1 = p_para[4]
    dV = p_para[5]
    Kmc = p_para[6]
    Vcmax = p_para[7]
    Tp = p_para[8]

    index = [i for i, x in enumerate(u0_CO2) if x < 0]
    for i in index:
        u0_CO2[i] = 0

    Aj = -Je * u0_CO2 / (4 * (u0_CO2 + 2 * gamma))
    Ac = -Vcmax * u0_CO2 / (u0_CO2 + Kmc)

    div_Sp_Aj = -1/2 * Je * gamma / ((u0_CO2 + 2 * gamma) ** 2)
    div_Sp_Ac = -Vcmax * Kmc / ((u0_CO2 + Kmc) ** 2)

    Ac = [max(x, y) for x, y in zip(Ac, Aj)]
    div_Sp_Ac = [x if x != y else y for x, y in zip(div_Sp_Ac, div_Sp_Aj)]

    Ap = Ac
    index = [i for i, x in enumerate(u0_CO2) if x > gamma]
    for i in index:
        Ap[i] = -3 * Tp * u0_CO2[i] / (u0_CO2[i] - gamma)
        div_Sp_Ap[i] = 3 * Tp * gamma / ((u0_CO2[i] - gamma) ** 2)

    photo = [max(x, y) for x, y in zip(Ac, Ap)]
    div_Sp = div_Sp_Ac
    div_Sp = [x if x != y else y for x, y in zip(div_Sp, div_Sp_Ap)]

    return photo, div_Sp
