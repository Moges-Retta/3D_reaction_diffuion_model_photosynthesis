# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:40:48 2024

@author: angel
"""

import numpy as np

dx = 0.66928e-6
dy = 0.66928e-6
dz = 0.66928e-6

data = np.load("data")

meso_vascular = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 1)
meso_vascular[a] = 1
a = np.where(data == 5)
meso_vascular[a] = 1
a = np.where(data == 8)
meso_vascular[a] = 1
a = np.where(data == 3)
meso_vascular[a] = 1
a = np.where(data == 7)
meso_vascular[a] = 1
a = np.where(data == 9)
meso_vascular[a] = 1

meso = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 1)
meso[a] = 1

epid = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 2)
epid[a] = 1

bundle = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 5)
bundle[a] = 1

chloro = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 8)
chloro[a] = 1

vacuole = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 3)
vacuole[a] = 1

chloro_bundle = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 9)
chloro_bundle[a] = 1

vascular = np.zeros(data.shape, dtype=np.uint8)
a = np.where(data == 7)
vascular[a] = 1

geom = epid + meso_vascular

Td = [13.5, 18, 25, 30, 34, 39]

E_Jmax = 77900
S_Jmax = 627
H_Jmax = 191929
E_Vcmax = 53400
E_Vpmax = 37000
E_gamma = 27400
E_KmC = 35600
E_KmO = 15100
E_Kp = 27200
S_Vpmax = 663
H_Vpmax = 214150
E_Rd = 41900

Iinc = 1500.092
beta = 0.85
theta = 0.79127
Xm = 2/3
Xb = 1/3
fe = 0.15

Ym = 0.6
Yb = 0.4

Im = Xm * beta * Iinc * (1 - fe) / 2
Ib = Xb * beta * Iinc * (1 - fe)

Jmax_m = Ym * 162
Jmax_b = Yb * 162

Jm = (Im + Jmax_m - (Im + Jmax_m)2 - 4 * theta * Im * Jmax_m)(1/2) / (2 * theta)
Jb = (Ib + Jmax_b - (Ib + Jmax_b)2 - 4 * theta * Ib * Jmax_b)(1/2) / (2 * theta)

Jt = Jm + Jb

T = 25
T = T + 273.15
P = 101325
R = 8.3144621
Hen = 0.83 * np.exp(-20256.28/R * (1/298.15 - 1/T))
U1 = 380 * P / R / T
J = Jt
xx = 0.4
J_PEP = xx * J
J_fix = (1 - xx) * J

Jmax_25 = 158.9588
a1 = E_Jmax * (T - 298.15) / (298.15 * R * T)
a2 = 1 + np.exp((298.15 * S_Jmax - H_Jmax) / (298.15 * R))
a3 = 1 + np.exp((T * S_Jmax - H_Jmax) / (T * R))
Jmax = Jmax_25 * np.exp(a1) * a2 / a3

Rd_25 = 1.70
Rd = Rd_25 * np.exp(E_Rd / R * (1/298.15 - 1/T))
Kmc_25 = 485
Kmc = Kmc_25 * np.exp(E_KmC / R * (1/298.15 - 1/T))

KmO_25 = 146000
KmO = KmO_25 * np.exp(E_KmO / R * (1/298.15 - 1/T))

import numpy as np

# Calculation of volume fractions
f = np.sum(chloro) / np.sum(1 - air)
thickness = np.sum(1 - air) / (air.shape[0] * air.shape[1]) * dz

f_Resp_epi = np.sum(geom) / np.sum(1 - air)
f2 = np.sum(meso) / (np.sum(meso + chloro + vacuole))
f3 = np.sum(bundle) / (np.sum(bundle + chloro_bundle))
f_vs = np.sum(vascular) / np.sum(1 - air)
f_vac = np.sum(vacuole) / np.sum(meso + chloro + vacuole)

f_BS_chloro = np.sum(chloro_bundle) / np.sum(1 - air)
f_MC_cyto = np.sum(meso) / np.sum(1 - air)
f_BS_cyto = np.sum(bundle) / np.sum(1 - air)
f_chloro_cyto = np.sum(chloro) / np.sum(meso)

f_Resp_meso = f_Resp_epi * f2
f_Resp_bundle = f_Resp_epi * f3
f_Resp_vs = f_Resp_epi

# Calculation of volumetric rates
J_PEP = J_PEP / (thickness * f) * (f_chloro_cyto)
J_fix = J_fix / (thickness * f_BS_chloro)
Rd_epid = Rd / (thickness * f_Resp_epi)
Rd_meso = Rd / (thickness * f_Resp_meso)
Rd_bundle = Rd / (thickness * f_Resp_bundle)
Rd_vascu = Rd / (thickness * f_Resp_vs)

Vcmax = Vcmax / (thickness * f_BS_chloro)
Vpmax = Vpmax / (thickness * f_MC_cyto)

p_para = [J_PEP, J_fix, Rd_meso, Rd_bundle, gamma, Hen, U1, dx * dy * dz, Kmc, Vcmax, Vpmax, Kp, Rd_epid, Rd_vascu, Rd, T, KmO]
np.save("p_para", p_para)
