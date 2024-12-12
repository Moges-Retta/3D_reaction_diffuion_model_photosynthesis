# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:47:58 2024

@author: angel
"""


import numpy as np
import scipy.io

p_para = scipy.io.loadmat('p_para.mat')
T = p_para['p_para'][15][0]

data = scipy.io.loadmat('data.mat')['data']
air = scipy.io.loadmat('air.mat')['air']
meso_ori = scipy.io.loadmat('meso_ori.mat')['meso_ori']
bundle_ori = scipy.io.loadmat('bundle_ori.mat')['bundle_ori']
u_CO2 = scipy.io.loadmat('u_CO2.mat')['u_CO2']

def discretise_HCO3(data, air, meso_ori, bundle_ori, T):
    aB_HCO3 = np.zeros(data.shape)
    Su_HCO3 = np.zeros(data.shape)
    K_HCO3 = np.zeros(data.shape)
    F_HCO3 = np.zeros(data.shape)
    # calculation of aB_HCO3, Su_HCO3, K_HCO3, F_HCO3 goes here

    return aB_HCO3, Su_HCO3, K_HCO3, F_HCO3

aB_HCO3, Su_HCO3, K_HCO3, F_HCO3 = discretise_HCO3(data, air, meso_ori, bundle_ori, T)

scipy.io.savemat('aB_HCO3.mat', {'aB_HCO3': aB_HCO3})
scipy.io.savemat('Su_HCO3.mat', {'Su_HCO3': Su_HCO3})
scipy.io.savemat('K_HCO3.mat', {'K_HCO3': K_HCO3})
scipy.io.savemat('F_HCO3.mat', {'F_HCO3': F_HCO3})

UO_1 = 0.21 * 101325 / 8.3144 / T * 10**6
u0_O2 = np.zeros(data.shape) + UO_1

def discretise_O2_fun(data, air, meso_ori, bundle_ori, u0_O2, T):
    aB_O2 = np.zeros(data.shape)
    K_O2 = np.zeros(data.shape)
    F_O2 = np.zeros(data.shape)
    # calculation of aB_O2, K_O2, F_O2 goes here

    return aB_O2, K_O2, F_O2

aB_O2, K_O2, F_O2 = discretise_O2_fun(data, air, meso_ori, bundle_ori, u0_O2, T)
scipy.io.savemat('aB_O2.mat', {'aB_O2': aB_O2})

import numpy as np

K_HCO3 = np.load('K_HCO3.npy')
F_HCO3 = np.load('F_HCO3.npy')

data = np.load('data.npy')

a = np.where(data == 1)
meso = np.zeros_like(data)
meso[a] = 1

a = np.where(data == 3)
vacuole = np.zeros_like(data)
vacuole[a] = 1

a = np.where(data == 5)
bundle = np.zeros_like(data)
bundle[a] = 1

a = np.where(data == 7)
vascular = np.zeros_like(data)
vascular[a] = 1

a = np.where(data == 8)
chloro = np.zeros_like(data)
chloro[a] = 1

chloro_bundle = np.zeros_like(data)
a = np.where(data == 9) 
chloro_bundle[a] = 1


a = np.where(data == 2)
epid = np.zeros_like(data)
epid[a] = 1

geom = epid + meso + vascular + chloro + chloro_bundle + bundle

del data

J_vol = np.load('J_vol.npy')
Je = J_vol.ravel()

Je_c_2 = np.load('Je_c_2.npy')
J_meso = Je_c_2.ravel()

JLET_BS = np.load('JLET_BS.npy')
JLET = JLET_BS.ravel()

p_para = np.load('p_para.npy')
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
KmO = p_para[16]
import numpy as np

# Load data
data = np.load("data.npy")
u_HCO3 = np.load("u_HCO3.npy")

# Initial guess for concentrations
u0_HCO3 = u_HCO3.flatten()
UO_1 = 210000 * 101325 / 8.3144 / T
u_O2 = np.zeros_like(data) + UO_1
u0_O2 = u_O2.flatten()

# Enzymatic CA hydration constants
ka_25 = 3e5
R = 8.314
E_CA = 40900
S_CA = 210
H_CA = 64500
a1 = E_CA * (T - 298.15) / (298.15 * R * T)
a2 = 1 + np.exp((298.15 * S_CA - H_CA) / (298.15 * R))
a3 = 1 + np.exp((T * S_CA - H_CA) / (T * R))
ka = ka_25 * np.exp(a1) * a2 / a3

C_CA = 0.7 * 1e6
K_const = 5.6e-7
kC = 2.8 * 1e6 / Hen
kHCO3 = 34 * 1e6 / Hen

# Non-enzymatic CA hydration constants
kn_25 = 0.039
E_Kn = 95000
Kn = kn_25 * np.exp(E_Kn / R * (1 / 298.15 - 1 / T))
k1 = ka * C_CA / kC
proton1 = 10 ** (-7.5)
proton2 = 10 ** (-8)
proton3 = 10 ** (-5)
pro = np.zeros_like(geom.flatten()) + proton1

a = np.where(chloro == 1)
pro[a] = proton2
a = np.where(chloro_bundle == 1)
pro[a] = proton2
a = np.where(vacuole == 1)
pro[a] = proton3

geom = geom + vacuole
import numpy as np

# Initial guess for concentrations are the results from the last iteration

# enzymatic CA hydration constants
ka_25 = 3e5
R = 8.314
E_CA = 40900
S_CA = 210
H_CA = 64500
a1 = E_CA * (T - 298.15) / (298.15 * R * T)
a2 = 1 + np.exp((298.15 * S_CA - H_CA) / (298.15 * R))
a3 = 1 + np.exp((T * S_CA - H_CA) / (T * R))
ka = ka_25 * np.exp(a1) * a2 / a3

C_CA = 0.7 * 1e6
K_const = 5.6e-7
kC = 2.8 * 1e6 / Hen
kHCO3 = 34 * 1e6 / Hen

# non-enzymatic CA hydration constants
kn_25 = 0.039
E_Kn = 95000 # Boyd et al. 2016
Kn = kn_25 * np.exp(E_Kn / R * (1 / 298.15 - 1 / T)) # Temperature dependency
k1 = ka * C_CA / kC
proton1 = 10 ** (-7.5) #mol/L, cytoplasm
proton2 = 10 ** (-8) #mol/L, chloroplast
proton3 = 10 ** (-5) #mol/L, vacuole
pro = np.zeros(geom.shape) + proton1

a = np.where(chloro == 1)[0]
pro[a] = proton2
a = np.where(chloro_bundle == 1)[0]
pro[a] = proton2
a = np.where(vacuole == 1)[0]
pro[a] = proton3

geom = geom + vacuole

CA = chloro * Kn / k1 + meso * 1 + bundle * Kn / k1 + chloro_bundle * Kn / k1 + vacuole * Kn / k1
a = np.where(epid == 1)[0]
CA[a] = Kn / k1 #non-enyzymatic hydration
a = np.where(vascular == 1)[0]
CA[a] = Kn / k1 #non-enyzymatic hydration

u0_CO2 = u_CO2.ravel()
u0_O2 = np.full(data.shape, UO_1)
u0_HCO3 = u_HCO3.ravel()

index_matrix = np.where(geom.ravel() == 1)[0]
K_HCO3 = K_HCO3[index_matrix][:, index_matrix]

rel = 1
tol = 1e-16 # tolerance
i = 1
# relative residuals
relres = 1 
rel_ma_CO2 = []
rel_ma_HCO3 = []
rel_ma_O2 = []

# Matrices for collecting minimum concentrations,tolerance, flag of solvers 
# and number of iterations
min_val = []
tol_result = []
flag = []
iter = []

tol_u_CO2 = []
tol_u_HCO3 = []
tol_u_O2 = []

dCA_CO2 = numpy.zeros(geom.ravel().shape)
dCA_HCO3 = numpy.zeros(geom.ravel().shape)
drp_CO2 = numpy.zeros(geom.ravel().shape)
drp_O2 = numpy.zeros(geom.ravel().shape)

while i < 90 and relres > tol:
tol_HCO3 = 1
j = 1
CA_r, div_CA = CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
CA_r0 = CA_r # calculate initial CA hydration rate
while j < 4 and tol_HCO3 > 1e-4:
CA_r, div_CA = CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
PEP, div_PEP = photo_PEP_J(u0_HCO3, [Vpmax, Kp], J_meso)
Su_HCO3 = (CA_r - div_CA * u0_HCO3) * geom.ravel() * dV + (PEP - u0_HCO3 * div_PEP) * dV * meso.ravel() - dCA_CO2 * geom.ravel() * dV
Sp_HCO3 = div_CA * dV * geom.ravel() + div_PEP * dV * meso.ravel()
F1_HCO3 = F_HCO3 + Su_HCO3
F1_HCO3 = F1_HCO3[index_matrix]
Sp_HCO3 = Sp_HCO3[index_matrix]
K1_HCO3 = K_HCO3 - scipy.sparse.spdiags(Sp_HCO3, 0, sparse.shape[0], sparse.shape[1])
tam = numpy.linalg.norm(F1_HCO3)
# clear Su_HCO3 Sp_HCO3 div_CA
# tic
R_HCO3 = scipy.sparse.linalg.factorized(K1_HCO3.astype(float).tocsc())
# toc
print("finish imcomplete LU factorization for K_HCO3")
tol_HCO3 = 1e-12
# tic
u_HCO3, flag_HCO3, relres_HCO3, iter_HCO3, resvec_HCO3 = scipy.sparse.linalg.cg(K1_HCO3, F1_HCO3, tol=tol, maxiter=200, M=R_HCO3)
tol_HCO3 = 1e-12
start_time = time.time()
u_HCO3, flag_HCO3, relres_HCO3, iter_HCO3, resvec_HCO3 = pcg(K1_HCO3, F1_HCO3, tol, 200, R_HCO3.T, R_HCO3, u0_HCO3[index_matrix])
end_time = time.time()
print("Elapsed time:", end_time - start_time)
print("finish solving HCO3")

relres_HCO3 = np.linalg.norm(F1_HCO3 - np.dot(K1_HCO3, u_HCO3)) / np.linalg.norm(F1_HCO3)

tam = (u_HCO3 - u0_HCO3[index_matrix]) / u0_HCO3[index_matrix]
tol_u_HCO3 = np.append(tol_u_HCO3, np.amax(np.abs(tam)))
tol_HCO3 = np.amax(np.abs(tam))

u_HCO3[u_HCO3 < 0] = 0
u0_HCO3[index_matrix] = u_HCO3
j = j + 1

import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import ichol, pcg

def CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    # Code to calculate CA_r and div_CA
    return CA_r, div_CA

def update_dCA(CA_r, CA_r0, dCA_HCO3, geom, constant):
    # Code to update dCA_HCO3
    return dCA_HCO3

def CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3):
    # Code to calculate CA_r and div_CA
    return CA_r, div_CA

def rp_discret_O2(u0_CO2, u0_O2, parameters, Je):
    # Code to calculate rp_r and div_rp
    return rp_r, div_rp

# Initialize variables
tol_O2 = 1
j = 1
tol_u_O2 = []
rel_ma_O2 = []

# Calculate initial CA hydration rate from CO2
CA_r, div_CA = CA_discret_HCO3(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
dCA_HCO3 = update_dCA(CA_r, CA_r0, dCA_HCO3, geom, 0.5)

CA_r, div_CA = CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
CA_r0 = CA_r

# Calculate initial photorespiration rate
rp_r, div_rp = rp_discret_O2(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)
rp_r0 = rp_r

tol_O2 = 1e-4
j = 0

while j < 4 and tol_O2 > 1e-4:
    rp_r, div_rp = rp_discret_O2(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)
    Su_O2 = (rp_r - div_rp * u0_O2) * chloro_bundle * dV - drp_CO2 * chloro_bundle * dV - \
        Rd_meso * dV * meso - Rd_bundle * dV * bundle + \
        0.25 * dV * J_vol * chloro + 0.25 * dV * JLET * chloro_bundle
    Sp_O2 = div_rp * dV * chloro_bundle
    F1_O2 = F_O2 + Su_O2
    
    K1_O2 = K_O2 - spdiags(Sp_O2, 0, sparse(K_O2.shape[0], K_O2.shape[1]))
    tam = np.linalg.norm(F1_O2)
    del Su_O2, Sp_O2, div_rp
    
    R_O2 = sp.sparse.linalg.spilu(K1_O2, drop_tol=0.2e-3)
    u_O2, _, relres_O2, _, resvec_O2 = sp.sparse.linalg.cg(K1_O2, F1_O2, tol=1e-12, maxiter=200, M=R_O2.solve)
    
    relres_O2 = np.linalg.norm(F1_O2 - K1_O2.dot(u_O2)) / np.linalg.norm(F1_O2)
    tam = (u_O2 - u0_O2) / u0_O2
    tol_u_O2.append(np.max(np.abs(tam)))
    tol_O2 = np.max(np.abs(tam))
    del tam
    np.save('tol_u_O2', tol_u_O2)
    
    u_O2[u_O2 < 0] = 0
    u0_O2 = u_O2
    rel_ma_O2.append(resvec_O2)
    j += 1

def solve_CO2(u0_CO2, u0_O2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3, chloro_bundle, bundle, bundle_ori, geom, PEP, Rd_meso, Rd_bundle, Rd_vascu, Rd_epid, dV, meso, vascular, epid, F, tol, K, s):
    rp_r, div_rp = rp_discret_O2(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)
    drp_O2 = update_drp(rp_r, rp_r0, drp_O2, geom, 0.5)

    rp_r, div_rp = rp_discret_CO2(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)
    rp_r0 = rp_r  # calculate initial rp hydration rate from CO2

    tol_CO2 = 1
    j = 1
    while j < 4 and tol_CO2 > 1e-4:

        photo, div_photo = photo_fixation_J(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)

        photo_resp_cyto = photoresp_discret_J_Modified(u0_CO2, u0_O2, chloro_bundle, bundle, bundle_ori, [Kmc, Vcmax, gamma, KmO], Je)
        PEP_de = -PEP * meso.flatten()
        PEP_de = sum(PEP_de.flatten()) / sum(chloro_bundle.flatten())

        Su = (photo - u0_CO2 * div_photo) * dV * chloro_bundle.flatten() + \
             Rd_meso * dV * meso.flatten() + Rd_bundle * dV * bundle.flatten() + \
             PEP_de * dV * chloro_bundle.flatten() + Rd_vascu * dV * vascular.flatten() + \
             Rd_epid * dV * epid.flatten() + photo_resp_cyto.flatten() * dV * bundle.flatten()

        Sp = div_photo * dV * chloro_bundle.flatten()

        Su = sparse(Su)
        Sp = sparse(Sp)

        CA_r, div_CA = CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
        Su_CO2 = Su + (CA_r - div_CA * u0_CO2) * geom.flatten() * dV - dCA_HCO3 * geom.flatten() * dV - drp_O2 * geom.flatten() * dV
        Sp = Sp + div_CA * geom.flatten() * dV
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np

K1 = K - sp.diags(Sp, 0, (len(K), len(K))) # coupling Sp to K matrix

R_CO2 = spla.splu(K1, permc_spec="NATURAL", diag_pivot_thresh=0.2e-3, options={"ILU_MILU": "SILU"})
print("finish imcomplete LU factorization for K of CO2")

u_CO2, info = spla.cg(K1, F + Su_CO2, x0=u0_CO2, tol=tol, maxiter=200, M=R_CO2.solve)
print("finish solving CO2")

relres_CO2 = np.linalg.norm(F + Su_CO2 - K1 @ u_CO2) / np.linalg.norm(F + Su_CO2)

tam = (u_CO2 - u0_CO2) / u0_CO2
tol_u_CO2 = np.append(tol_u_CO2, np.max(np.abs(tam)))
tol_CO2 = np.max(np.abs(tam))

index = np.where(u_CO2 < 0)
u_CO2[index] = 0
u0_CO2 = u_CO2
rel_ma_CO2 = np.append(rel_ma_CO2, resvec_CO2)
j += 1

CA_r, div_CA = CA_discret_CO2(u0_CO2, u0_HCO3, pro, CA, ka, C_CA, K_const, kC, kHCO3)
dCA_CO2 = update_dCA(CA_r, CA_r0, dCA_CO2, geom, 0.5)

rp_r, div_rp = rp_discret_CO2(u0_CO2, u0_O2, [Vcmax, Kmc, gamma, KmO], Je)
drp_CO2 = update_drp(rp_r, rp_r0, drp_CO2, geom, 0.5)

i += 1

rel_ma_HCO3.append(resvec_HCO3)
min_val.append([min(u0_CO2), max(u0_HCO3)])
tol_result.append([relres_CO2, relres_HCO3, relres_O2])
relres = max(relres_CO2, relres_HCO3)
flag.append([flag_CO2, flag_HCO3])
iter.append([iter_CO2, iter_HCO3])

np.save("u_CO2.npy", u_CO2)
np.save("u_HCO3.npy", u_HCO3)
np.save("u_O2.npy", u_O2)
np.save("flag.npy", flag)
np.save("relres.npy", relres)
np.save("iter.npy", iter)
np.save("min_val.npy", min_val)
np.save("tol_result.npy", tol_result)
np.save("rel_ma_CO2.npy", rel_ma_CO2)
np.save("rel_ma_HCO3.npy", rel_ma_HCO3)

data = np.load("data.npy")
n = data.shape[0]

u_CO2 = np.reshape(u_CO2, n)
np.save("u_CO2.npy", u_CO2)

u_O2 = np.reshape(u_O2, n)
np.save("u_O2.npy", u_O2)

u_HCO3 = u0_HCO3
u_HCO3 = np.reshape(u_HCO3, n)
np.save("u_HCO3.npy", u_HCO3)
