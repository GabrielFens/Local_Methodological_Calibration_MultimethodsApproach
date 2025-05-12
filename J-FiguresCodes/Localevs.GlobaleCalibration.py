#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 23:37:19 2025

@author: gabrielfenisse
"""


#%%Comparaison between global and local calibration
dataComparaison_Cal=pd.read_csv('ComparisonMMethods_Calibration.csv',sep=";",decimal=",")
TANN_MAT_G=dataComparaison_Cal["TANN_MAT_Globale"]
TANN_WAPLS_G=dataComparaison_Cal["TANN_WAPLS_Globale"]
TANN_CREST_G=dataComparaison_Cal["TANN_CREST_Globale"]

ErrTANN_MAT_G=dataComparaison_Cal["TANN_MATSigma_Globale"]
ErrTANN_WAPLS_G=dataComparaison_Cal["TANN_WAPLSSigma_Globale"]
ErrTANN_CREST_G=dataComparaison_Cal["TANN_CRESTSigma_Globale"]

TANN_MAT_L=dataComparaison_Cal["TANN_MAT_Locale"]
TANN_WAPLS_L=dataComparaison_Cal["TANN_WAPLS_Locale"]
TANN_CREST_L=dataComparaison_Cal["TANN_CREST_Locale"]

ErrTANN_MAT_L=dataComparaison_Cal["TANN_MATSigma_Locale"]
ErrTANN_WAPLS_L=dataComparaison_Cal["TANN_WAPLSSigma_Locale"]
ErrTANN_CREST_L=dataComparaison_Cal["TANN_CRESTSigma_Locale"]

#%%Plot TANN
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(20, 7.5))

ax.set_title("a. Locale-Globale LGM TANN reconstructions from multi-methods approach")
ax.axhline(y=0, color='black', linestyle='--', linewidth=1)

x_positions = np.arange(len(clean_Name))
offset = 0.1

ax.errorbar(x_positions - offset, TANN_CREST_L - TANN_CREST_G, 
            yerr=ErrTANN_CREST_G + ErrTANN_CREST_L, color='red', alpha=0.25, markerfacecolor='red', fmt='o', elinewidth=1)
ax.errorbar(x_positions, TANN_MAT_L - TANN_MAT_G, 
            yerr=ErrTANN_MAT_G + ErrTANN_MAT_L, color='blue', alpha=0.25, markerfacecolor='blue', fmt='o', elinewidth=1)
ax.errorbar(x_positions + offset, TANN_WAPLS_L - TANN_WAPLS_G, 
            yerr=ErrTANN_WAPLS_G + ErrTANN_WAPLS_L, color='green', fmt='o', markerfacecolor='green', elinewidth=1, alpha=0.25)

values_crest = TANN_CREST_L - TANN_CREST_G
values_mat = TANN_MAT_L - TANN_MAT_G
values_wapls = TANN_WAPLS_L - TANN_WAPLS_G

ax.scatter(x_positions - offset, values_crest, color='red', label="CREST")
ax.scatter(x_positions, values_mat, color='blue', label="MAT")
ax.scatter(x_positions + offset, values_wapls, color='green', label="WAPLS")

ax.set_ylabel("Biomization effect on LGM TANN (°C)")
ax.set_xticks(range(len(clean_Name)))
ax.set_xticklabels(clean_Name, rotation=90)
ax.legend()

plt.tight_layout()
plt.savefig("Globale_Locale_TANN.pdf")
plt.show()

#%%Legend
plt.scatter(1,np.mean(MAT_TANN),s=150,color='green',label='MAT-Biome',alpha=0.3)
plt.scatter(1,np.mean(MAT_TANN_MB),s=150,color='green',label='MAT-Megabiome')

plt.scatter(1,np.mean(WAPLS_TANN),s=150,color='blue',label='WAPLS-Biome',alpha=0.3)
plt.scatter(1,np.mean(WAPLS_TANN_MB),s=150,color='blue',label="WAPLS-Megabiome")

plt.scatter(1,np.mean(CREST_TANN),s=150,color='red',label='CREST-Biome',alpha=0.3)
plt.scatter(1,np.mean(CREST_TANN_MB),s=150,color='red',label='CREST-Megabiome')
plt.legend()
plt.savefig("ScatterCorrelated_CREST-Legend.pdf")


#%%
slope, intercept, r_value, p_value, std_err = linregress(nylat,TANN_CREST_L - TANN_CREST_G)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

fig = plt.figure(figsize=(10,5))
plt.title("b. Latitudinal dependance of CREST biomization effect")
plt.scatter(nylat,TANN_CREST_L - TANN_CREST_G,color='red')
plt.plot(nylat,slope * nylat + intercept, color='black', label='Biomes', linewidth=2)
plt.xlabel("Latitude (°N)")
plt.ylabel("CREST biomization effect on LGM TANN (°C)")
plt.savefig("Globale_Locale_TANN_LatitudeRelationship.pdf")
plt.show()
