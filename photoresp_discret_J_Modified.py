
def photoresp_discret_J_Modified(u0_CO2, u0_O2, chloro_bundle, bundle, bundle_ori, p_para, Je):
Kmc = p_para[0]
Vcmax = p_para[1]
gamma = p_para[2]
KmO = p_para[3]
index = [i for i in range(len(u0_CO2)) if u0_CO2[i] < 0]
for i in index:
    u0_CO2[i] = 0

KM = Kmc * (1 + u0_O2 / KmO)

photoresp_Aj = Je * gamma * u0_O2 / (3 * (u0_CO2 + 7 * gamma * u0_O2 / 3))
photoresp_Ac = Vcmax * gamma * u0_O2 / (u0_CO2 + KM)
photoresp = [min(photoresp_Aj[i], photoresp_Ac[i]) for i in range(len(u0_CO2))]

n = max(bundle_ori)
photo_resp_cyto = [0 for i in range(len(bundle_ori))]

for i in range(1, n + 1):
    index = [j for j in range(len(bundle_ori)) if bundle_ori[j] == i]
    if len(index) > 0:
        tam = [0 for j in range(len(bundle_ori))]
        for j in index:
            tam[j] = 1
        
        tam_chloro = [chloro_bundle[j] * tam[j] for j in range(len(bundle_ori))]
        tam_cyto = [bundle[j] * tam[j] for j in range(len(bundle_ori))]
        
        index_chloro = [j for j in range(len(tam_chloro)) if tam_chloro[j] == 1]
        u_chloro = sum(u0_CO2[j] for j in index_chloro) / len(index_chloro)
        photoresp_cell = sum(photoresp[j] for j in index_chloro) / len(index_chloro)
        
        index_cyto = [j for j in range(len(tam_cyto)) if tam_cyto[j] == 1]
        if len(index_cyto) > 0:
            f_chloro_cyto = len(index_chloro) / len(index_cyto)
            f_chloro_cyto = min(f_chloro_cyto, 6)
            for j in index_cyto:
                photo_resp_cyto[j] = photoresp_cell * f_chloro_cyto
                
return photo_resp_cyto
