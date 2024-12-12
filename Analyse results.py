# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:31:49 2024

@author: angel
"""
import numpy as np

def flux_calculation_z(data, u_CO2, u_O2, u_HCO3, p_para):
# Constants
dx = 0.66928e-6
dy = 0.66928e-6
dz = 0.66928e-6

# Make array for liquid phase media in meso_vascular
meso_vascular = np.zeros(data.shape)
a = np.where((data == 1) | (data == 5) | (data == 3) | (data == 8) | (data == 7) | (data == 9))
meso_vascular[a] = 1

# Mesophyll cytosol
meso = np.zeros(data.shape)
a = np.where(data == 1)
meso[a] = 1

# Epidermis
epid = np.zeros(data.shape)
a = np.where(data == 2)
epid[a] = 1

# Vacuole
vacuole = np.zeros(data.shape)
a = np.where(data == 3)
vacuole[a] = 1

# Cytosol of bundle sheath cells
 bundle = np.zeros(data.shape)
a = np.where(data == 5)
bundle[a] = 1

    # Chloroplasts of bundle sheath cells
    chloro = np.zeros(data.shape)
    a = np.where(data == 8)
    chloro[a] = 1

    # Chloroplasts of mesophyll cells
    chloro_bundle = np.zeros(data.shape)
    a = np.where(data == 9)
    chloro_bundle[a] = 1

    # Vascular bundles
    vascular = np.zeros(data.shape)
    a = np.where(data == 7)
    vascular[a] = 1

    geom = epid + meso_vascular  # Liquid phase media

    # Load model parameters
    Rd_meso = p_para[2]
    Rd_bundle = p_para[3]
    gamma = p_para[4]
    Hen = p_para[5]
    U1 = p_para[6]
    dV = p_para[7]
    Kmc = p_para[8]
    Vcmax = p_para[9]
    Vpmax = p_para[10]
    Kp = p_para[11]
    Rd_epid = p_para[12]
    Rd_vascu = p_para[13]
    Rd = p_para[14]
    KmO = p_para[16]

    P = 101325  # Pa/atm
    T = p_para[15]  # K
    R = 8.3144621  # m3 Pa/mol/K

    # Load data of the air outside the leaf
    air = load_air()

    # Calculate volume fractions
f = np.sum(chloro)

# Calculate volume fractions
f = np.sum(chloro) / (np.sum(1-air))  # Chloroplasts
thickness = np.sum(1-air) / (air.shape[0] * air.shape[1]) * dz  # Leaf thickness
f_vs = np.sum(vascular) / (np.sum(1-air))  # Vascular bundle
f2 = np.sum(meso) / (np.sum(meso + chloro + vacuole))  # Mesophyll cytosol
f3 = np.sum(bundle) / (np.sum(bundle + chloro_bundle))  # Bundle sheath cytosol
f_BS_chloro = np.sum(chloro_bundle) / (np.sum(1-air))  # Bundle sheath chloroplast
f_MC_cyto = np.sum(meso) / (np.sum(1-air))  # Cytosol of mesophyll
f_BS_cyto = np.sum(bundle) / (np.sum(1-air))  # Bundle sheath cyto
f_epid = np.sum(epid) / (np.sum(1-air))  # Epidermis

# Calculation of fraction of total respiration rate in epidermis, mesophyll, and bundle sheath cells
f_Resp_epi = np.sum(geom) / (np.sum(1-air))  # fraction of epidermis respiration
f_Resp_meso = f_Resp_epi * f2  # fraction of meso respiration
f_Resp_bundle = f_Resp_epi * f3  # fraction of bundle respiration

J_vol = np.load("J_vol.npy")
J_fix = J_vol
Vcmax = Vcmax * (thickness * f_BS_chloro)
Vpmax = Vpmax * (thickness * f_MC_cyto)
J = J_fix * (thickness * f_BS_chloro)
KM = Kmc * (1 + u_O2 / KmO)
Aj = J * (u_CO2) / (3 * (u_CO2 + 7 * gamma * u_O2 / 3))
Ac = Vcmax * (u_CO2) / (u_CO2 + KM)

photo_G = np.minimum(Ac, Aj)
index = np.where(photo_G < -Rd)

Ac = Ac.flatten()
Aj = Aj.flatten()
photo_G = np.maximum(Ac[index], Aj[index])
photo_G = photo_G * chloro_bundle.flatten()
a = np.where(chloro_bundle == 1)
photo_G = photo_G[a]

# Gross photosynthesis
photo_G_result = [np.mean(photo_G), np.min(photo_G), np.max(photo_G)]  # µmol/m2s

Ac = Ac * chloro_bundle.flatten()
a = np.where(chloro_bundle == 1)
import numpy as np

Ac_r = np.mean(np.mean(Ac(a)))

Aj = np.multiply(Aj, chloro_bundle.flatten())
a = np.where(chloro_bundle == 1)
Aj_r = np.mean(Aj[a])

photoresp_Aj = J * gamma * u_O2 / (3 * (u_CO2 + 7 * gamma * u_O2 / 3))
photoresp_Ac = Vcmax * gamma * u_O2 / (u_CO2 + KM)
photoresp_Ac = photoresp_Ac.flatten()
photoresp_Aj = photoresp_Aj.flatten()
photoresp = np.min(np.array([photoresp_Aj, photoresp_Ac]).T, axis=1)
index_chloro = np.where(chloro_bundle == 1)
photoresp_cell = np.mean(photoresp[index_chloro])

J_PEP = np.load("Je_c_2.npy")

J2 = J_PEP * (thickness * f_MC_cyto)

Vpj = J2 / 2 + u_CO2 * 0
Vpc = Vpmax * u_HCO3 / (u_HCO3 + Kp)
PEPc = np.minimum(Vpj.flatten(), Vpc.flatten())
index = np.where(PEPc < -Rd)
PEPc[index] = np.maximum(Vpc[index], Vpj[index])
PEPc = np.multiply(PEPc, meso.flatten())
a = np.where(meso == 1)
PEPc = PEPc[a]
PEPc_result = [np.mean(PEPc), np.min(PEPc), np.max(PEPc)]
Vpc = np.multiply(Vpc, meso.flatten())
a = np.where(meso == 1)
Vpc_r = np.mean(Vpc[a])

Vpj = np.multiply(Vpj, meso.flatten())
a = np.where(meso == 1)
Vpj_r = np.mean(Vpj[a])

data = np.load("air.npy")
a = np.where(data > 0)
cell = np.zeros_like(data)
cell[a] = 1
air_intercellular = 1 - cell - air

Ci_meso = np.multiply(u_CO2, air_intercellular)
tam = np.where(air_intercellular == 1)
Ci_meso = Ci_meso[tam]
Ci_result = [np.mean(Ci_meso), np.min(Ci_meso), np.max(Ci_meso)]
Ci_result = np.multiply(Ci_result, T * R / P)
a = (chloro_bundle == 1).nonzero()[0]
Cc = chloro_bundle[a] * u_CO2[a]
O2 = chloro_bundle[a] * u_O2[a]
Cc_result = [Cc.mean(), Cc.min(), Cc.max()]
Oc_result = [O2.mean(), O2.min(), O2.max()]
CH_B = chloro_bundle[a] * u_HCO3[a]
CH_B_result = [CH_B.mean(), CH_B.min(), CH_B.max()]

BB = Cc_result[0] / CH_B_result[0]

Cc_result = [x * T * R / P for x in Cc_result]
Oc_result = [x * T * R / P for x in Oc_result]

a = (meso == 1).nonzero()[0]
Cm = meso[a] * u_CO2[a]
Om = meso[a] * u_O2[a]
CH = meso[a] * u_HCO3[a]

Cm_result = [Cm.mean(), Cm.min(), Cm.max()]
Om_result = [Om.mean(), Om.min(), Om.max()]
CH_result = [CH.mean(), CH.min(), CH.max()]
MM = CH_result[0] / Cm_result[0]

Cm_result = [x * T * R / P for x in Cm_result]
Om_result = [x * T * R / P for x in Om_result]

A_G = photo_G_result[0]
A_flux = photo
Vp = PEPc_result[0]

Ca = U1 * T * R / P
Ci = Ci_result[0]
Cm = Cm_result[0]
Cc = Cc_result[0]
O2b = Oc_result[0] / 1000
O2m = Om_result[0] / 1000

Rm = Rd_meso * (thickness * f_MC_cyto)
Rs = Rd_bundle * (thickness * f_BS_cyto)
Repi = Rd_epid * (thickness * f_epid)
Rv = Rd_vascu * (thickness * f_vs)

Rd_calc = Rm + Rs + Repi + Rv

A_net = photo_G_result[0] - photoresp_cell - Rd_calc

gs = A_flux / (Ca - Ci)
gm = A_flux / (Ci - Cm)
L = Vp - A_flux - Rm
gbs = L * 1000 / (Cc - Cm)

