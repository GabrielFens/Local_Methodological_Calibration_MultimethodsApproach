#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 19:39:32 2025

@author: gabrielfenisse
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

Age_LGP2 = pd.read_csv("ModeleAge_LGP2.txt", sep=',', decimal='.')
Age_LGP1 = pd.read_csv("ModeleAge_LGP.txt", sep=',', decimal='.')
Age_LGP=np.append(Age_LGP2["X2"][0:27],Age_LGP1["X2"])
Age_LDB = pd.read_csv("ModeleAge_LDB.txt", sep=',', decimal='.')

climate_SynthesisLGP_MAT = pd.read_csv("Output_Recon_MAT_BiomeWeighted_SynthesisLGP.csv", sep=";", decimal=",")
climate_SynthesisLDB_MAT = pd.read_csv("Output_Recon_MAT_BiomeWeighted_LDB.csv", sep=";", decimal=",")

climate_SynthesisLGP_WAPLS = pd.read_csv("Output_Recon_WAPLS_BiomeWeighted_SynthesisLGP.csv", sep=";", decimal=",")
climate_SynthesisLDB_WAPLS = pd.read_csv("Output_Recon_WAPLS_BiomeWeighted_LDB.csv", sep=";", decimal=",")

climate_SynthesisLGP_CREST = pd.read_csv("Output_Recon_CREST_BiomeWeighted_SynthesisLGP.csv", sep=";", decimal=",")
climate_SynthesisLDB_CREST = pd.read_csv("Output_Recon_CREST_BiomeWeighted_LDB.csv", sep=";", decimal=",")

ModernTANN_LGP=10.7
ModernTANN_LDB=8.3
ModernPANN_LGP=1041.5
ModernPANN_LDB=1001.7
ModernMTWA_LGP=18.2
ModernMTWA_LDB=17.4
ModernMTCO_LGP=2.0
ModernMTCO_LDB=0.2

#%%
NGRIP=pd.read_csv('NGRIP_Records.txt', sep="\t")
plt.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black')
plt.xlim(0,40)
plt.savefig("NGRIP.pdf")
plt.show()

#%%
fig, ax = plt.subplots(2, 2, figsize=(20, 15))

#TANN
ax[0,0].set_title("a. LGP - Mean Annual Temperature Anomaly", fontsize=16)

ax2 = ax[0,0].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)
ax2.legend(bbox_to_anchor=(0.18, 0.18),fontsize=10)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1."]-ModernTANN_LGP, color='blue', label='MAT')
ax[0,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1."] - climate_SynthesisLGP_MAT["predMAT.SEP.boot...1."]-ModernTANN_LGP, 
                  climate_SynthesisLGP_MAT["predMAT.fit...1."] + climate_SynthesisLGP_MAT["predMAT.SEP.boot...1."]-ModernTANN_LGP, color='blue', alpha=0.3)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_WAPLS["pred.fit...1."]-ModernTANN_LGP, color='green', label='WAPLS')
ax[0,0].fill_between(Age_LGP/1E3,climate_SynthesisLGP_WAPLS["pred.fit...1."] - climate_SynthesisLGP_WAPLS["pred.SEP.boot...1."]-ModernTANN_LGP,
                   climate_SynthesisLGP_WAPLS["pred.fit...1."] + climate_SynthesisLGP_WAPLS["pred.SEP.boot...1."]-ModernTANN_LGP, color='green', alpha=0.2)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_CREST['mean']-ModernTANN_LGP, color='red', label='CREST')
ax[0,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_CREST['X0.95_inf']-ModernTANN_LGP, climate_SynthesisLGP_CREST['X0.95_sup']-ModernTANN_LGP, color='red', alpha=0.1)

ax[0,0].set_xlim(0, 40)
ax[0,0].set_xlabel("", fontsize=15)
ax[0,0].set_ylabel("Mean Annual Temperature Anomaly (°C)", fontsize=15)
ax[0,0].legend(loc='lower left',fontsize=10)
ax[0,0].tick_params(axis='both', which='major', labelsize=15)
ax[0,0].set_xticks(np.arange(0, 41, 5))  
ax[0,0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[0,0].axvspan(19, 23, color='grey', alpha=0.3)

ax[0,1].set_title("b. LDB - Mean Annual Temperature Anomaly", fontsize=16)

ax2 = ax[0,1].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=13, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)


ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1."]-ModernTANN_LDB, color='blue', label='MAT')
ax[0,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1."] - climate_SynthesisLDB_MAT["predMAT.SEP.boot...1."]-ModernTANN_LDB, 
                  climate_SynthesisLDB_MAT["predMAT.fit...1."] + climate_SynthesisLDB_MAT["predMAT.SEP.boot...1."]-ModernTANN_LDB, color='blue', alpha=0.3)

ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_WAPLS["pred.fit...1."]-ModernTANN_LDB, color='green', label='WAPLS')
ax[0,1].fill_between(Age_LDB["X2"]/1E3,climate_SynthesisLDB_WAPLS["pred.fit...1."] - climate_SynthesisLDB_WAPLS["pred.SEP.boot...1."]-ModernTANN_LDB,
                   climate_SynthesisLDB_WAPLS["pred.fit...1."] + climate_SynthesisLDB_WAPLS["pred.SEP.boot...1."]-ModernTANN_LDB, color='green', alpha=0.2)

ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['mean']-ModernTANN_LDB, color='red', label='CREST')
ax[0,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['X0.95_inf']-ModernTANN_LDB, climate_SynthesisLDB_CREST['X0.95_sup']-ModernTANN_LDB, color='red', alpha=0.1)

ax[0,1].set_xlim(0, 40)
ax[0,1].set_ylabel("Mean Annual Temperature Anomaly (°C)", fontsize=15)
ax[0,1].tick_params(axis='both', which='major', labelsize=15)

ax[0,1].set_xticks(np.arange(0, 41, 5))  
ax[0,1].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[0,1].axvspan(19, 23, color='grey', alpha=0.3)

#PANN
ax[1,0].set_title("c. LGP - Mean Annual Precipitation Anomaly", fontsize=16)

ax2 = ax[1,0].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)

ax[1,0].plot(Age_LGP/1E3, ((climate_SynthesisLGP_MAT["predMAT.fit...1..1"]/ModernPANN_LGP)-1)*100, color='blue', label='MAT')
ax[1,0].fill_between(Age_LGP/1E3, ((climate_SynthesisLGP_MAT["predMAT.fit...1..1"]/ModernPANN_LGP)-1)*100-np.sqrt(((100/ModernPANN_LGP))**2*(climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..1"])**2), 
                  ((climate_SynthesisLGP_MAT["predMAT.fit...1..1"]/ModernPANN_LGP)-1)*100+np.sqrt(((100/ModernPANN_LGP))**2*(climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..1"])**2), color='blue', alpha=0.3)

ax[1,0].plot(Age_LGP/1E3, ((climate_SynthesisLGP_WAPLS["pred.fit...1..1"]/ModernPANN_LGP)-1)*100, color='green', label='WAPLS')
ax[1,0].fill_between(Age_LGP/1E3,((climate_SynthesisLGP_WAPLS["pred.fit...1..1"]/ModernPANN_LGP)-1)*100- np.sqrt(((100/ModernPANN_LGP))**2*(climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..1"])**2),
                   ((climate_SynthesisLGP_WAPLS["pred.fit...1..1"]/ModernPANN_LGP)-1)*100 + np.sqrt(((100/ModernPANN_LGP))**2*(climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..1"])**2), color='green', alpha=0.2)

ax[1,0].plot(Age_LGP/1E3, ((climate_SynthesisLGP_CREST['mean.1']/ModernPANN_LGP)-1)*100, color='red', label='CREST')
ax[1,0].fill_between(Age_LGP/1E3, ((climate_SynthesisLGP_CREST['X0.95_inf.1']/ModernPANN_LGP)-1)*100, ((climate_SynthesisLGP_CREST['X0.95_sup.1']/ModernPANN_LGP)-1)*100, color='red', alpha=0.1)

ax[1,0].set_xlim(0, 40)
ax[1,0].set_xlabel("", fontsize=15)
ax[1,0].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[1,0].set_ylabel("Mean Annual Precipitation Anomaly (%)", fontsize=15)
ax[1,0].tick_params(axis='both', which='major', labelsize=15)
ax[1,0].set_xticks(np.arange(0, 41, 5))  
ax[1,0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[1,0].axvspan(19, 23, color='grey', alpha=0.3)

ax[1,1].set_title("d. LDB - Mean Annual Precipitation Anomaly", fontsize=16)
ax2 = ax[1,1].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)


ax[1,1].plot(Age_LDB["X2"]/1E3, ((climate_SynthesisLDB_MAT["predMAT.fit...1..1"]/ModernPANN_LDB)-1)*100, color='blue', label='MAT')
ax[1,1].fill_between(Age_LDB["X2"]/1E3, ((climate_SynthesisLDB_MAT["predMAT.fit...1..1"]/ModernPANN_LDB)-1)*100-np.sqrt(((100/ModernPANN_LDB))**2*(climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..1"])**2), 
                  ((climate_SynthesisLDB_MAT["predMAT.fit...1..1"]/ModernPANN_LDB)-1)*100+np.sqrt(((100/ModernPANN_LDB))**2*(climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..1"])**2), color='blue', alpha=0.3)

ax[1,1].plot(Age_LDB["X2"]/1E3, ((climate_SynthesisLDB_WAPLS["pred.fit...1..1"]/ModernPANN_LDB)-1)*100, color='green', label='WAPLS')
ax[1,1].fill_between(Age_LDB["X2"]/1E3,((climate_SynthesisLDB_WAPLS["pred.fit...1..1"]/ModernPANN_LDB)-1)*100- np.sqrt(((100/ModernPANN_LDB))**2*(climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..1"])**2),
                   ((climate_SynthesisLDB_WAPLS["pred.fit...1..1"]/ModernPANN_LDB)-1)*100 + np.sqrt(((100/ModernPANN_LDB))**2*(climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..1"])**2), color='green', alpha=0.2)

ax[1,1].plot(Age_LDB["X2"]/1E3, ((climate_SynthesisLDB_CREST['mean.1']/ModernPANN_LDB)-1)*100, color='red', label='CREST')
ax[1,1].fill_between(Age_LDB["X2"]/1E3, ((climate_SynthesisLDB_CREST['X0.95_inf.1']/ModernPANN_LDB)-1)*100, ((climate_SynthesisLDB_CREST['X0.95_sup.1']/ModernPANN_LDB)-1)*100, color='red', alpha=0.1)

ax[1,1].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[1,1].set_xlim(0, 40)
ax[1,1].set_ylabel("Mean Annual Precipitation Anomaly (%)", fontsize=15)
ax[1,1].tick_params(axis='both', which='major', labelsize=15)
ax[1,1].set_xticks(np.arange(0, 41, 5))  
ax[1,1].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[1,1].axvspan(19, 23, color='grey', alpha=0.3)

plt.tight_layout()
plt.savefig("LGP&LDBReconstruction_AnnualPlot.pdf")
plt.show()

#%%
fig, ax = plt.subplots(2, 2, figsize=(20, 15))

#MTWA
ax[0,0].set_title("a. LGP - MTWA", fontsize=16)
ax2 = ax[0,0].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)
ax2.legend(bbox_to_anchor=(0.20, 0.23),fontsize=10)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..2"]-ModernMTWA_LGP, color='blue',label='MAT')
ax[0,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..2"] - climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..2"]-ModernMTWA_LGP, 
                  climate_SynthesisLGP_MAT["predMAT.fit...1..2"] + climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..2"]-ModernMTWA_LGP, color='blue', alpha=0.3)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_WAPLS["pred.fit...1..2"]-ModernMTWA_LGP, color='green',label='WAPLS')
ax[0,0].fill_between(Age_LGP/1E3,climate_SynthesisLGP_WAPLS["pred.fit...1..2"] - climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..2"]-ModernMTWA_LGP,
                   climate_SynthesisLGP_WAPLS["pred.fit...1..2"] + climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..2"]-ModernMTWA_LGP, color='green', alpha=0.2)

ax[0,0].plot(Age_LGP/1E3, climate_SynthesisLGP_CREST['mean.2']-ModernMTWA_LGP, color='red',label='CREST')
ax[0,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_CREST['X0.95_inf.2']-ModernMTWA_LGP, climate_SynthesisLGP_CREST['X0.95_sup.2']-ModernMTWA_LGP, color='red', alpha=0.1)


ax[0,0].set_xlim(0, 40)
ax[0,0].set_xlabel("", fontsize=15)
ax[0,0].set_ylabel("Mean Monthly Temperature Anomaly (°C)", fontsize=15)
ax[0,0].tick_params(axis='both', which='major', labelsize=15)
ax[0,0].set_xticks(np.arange(0, 41, 5))  
ax[0,0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[0,0].axvspan(19, 23, color='grey', alpha=0.3)
ax[0,0].legend(loc='lower left',fontsize=15)

ax[0,1].set_title("b. LDB - MTWA", fontsize=16)

ax2 = ax[0,1].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)


ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..2"]-ModernMTWA_LDB, color='blue')
ax[0,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..2"] - climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..2"]-ModernMTWA_LDB, 
                  climate_SynthesisLDB_MAT["predMAT.fit...1..2"] + climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..2"]-ModernMTWA_LDB, color='blue', alpha=0.3)

ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_WAPLS["pred.fit...1..2"]-ModernMTWA_LDB, color='green')
ax[0,1].fill_between(Age_LDB["X2"]/1E3,climate_SynthesisLDB_WAPLS["pred.fit...1..2"] - climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..2"]-ModernMTWA_LDB,
                   climate_SynthesisLDB_WAPLS["pred.fit...1..2"] + climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..2"]-ModernMTWA_LDB, color='green', alpha=0.2)

ax[0,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['mean.2']-ModernMTWA_LDB, color='red')
ax[0,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['X0.95_inf.2']-ModernMTWA_LDB, climate_SynthesisLDB_CREST['X0.95_sup.2']-ModernMTWA_LDB, color='red', alpha=0.1)

ax[0,1].set_xlim(0, 40)
ax[0,1].set_ylabel("Mean Annual Temperature Anomaly (°C)", fontsize=15)
ax[0,1].tick_params(axis='both', which='major', labelsize=15)
ax[0,1].axvspan(19, 23, color='grey', alpha=0.3)

ax[0,1].set_xticks(np.arange(0, 41, 5))  
ax[0,1].set_xticklabels(np.arange(0, 41, 5), fontsize=15)

#MTCO
ax[1,0].set_title("c. LGP - MTCO", fontsize=16)
ax2 = ax[1,0].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)

ax[1,0].plot(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..3"]-ModernMTCO_LGP, color='blue')
ax[1,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..3"] - climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..3"]-ModernMTCO_LGP, 
                  climate_SynthesisLGP_MAT["predMAT.fit...1..3"] + climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..3"]-ModernMTCO_LGP, color='blue', alpha=0.3)

ax[1,0].plot(Age_LGP/1E3, climate_SynthesisLGP_WAPLS["pred.fit...1..3"]-ModernMTCO_LGP, color='green')
ax[1,0].fill_between(Age_LGP/1E3,climate_SynthesisLGP_WAPLS["pred.fit...1..3"] - climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..3"]-ModernMTCO_LGP,
                   climate_SynthesisLGP_WAPLS["pred.fit...1..3"] + climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..3"]-ModernMTCO_LGP, color='green', alpha=0.3)

ax[1,0].plot(Age_LGP/1E3, climate_SynthesisLGP_CREST['mean.3']-ModernMTCO_LGP, color='red')
ax[1,0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_CREST['X0.95_inf.3']-ModernMTCO_LGP, climate_SynthesisLGP_CREST['X0.95_sup.3']-ModernMTCO_LGP, color='red', alpha=0.1)

ax[1,0].set_xlim(0, 40)
ax[1,0].set_xlabel("", fontsize=15)
ax[1,0].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[1,0].set_ylabel("Mean Monthly Temperature Anomaly (°C)", fontsize=15)
ax[1,0].tick_params(axis='both', which='major', labelsize=15)
ax[1,0].set_xticks(np.arange(0, 41, 5))  
ax[1,0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[1,0].axvspan(19, 23, color='grey', alpha=0.3)

ax[1,1].set_title("d. LDB - MTCO", fontsize=16)

ax2 = ax[1,1].twinx()
ax2.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black',alpha=0.5,label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y',labelcolor='black', labelsize=15)


ax[1,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..3"]-ModernMTCO_LDB, color='blue')
ax[1,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..3"] - climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..3"]-ModernMTCO_LDB, 
                  climate_SynthesisLDB_MAT["predMAT.fit...1..3"] + climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..3"]-ModernMTCO_LDB, color='blue', alpha=0.3)

ax[1,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_WAPLS["pred.fit...1..3"]-ModernMTCO_LDB, color='green')
ax[1,1].fill_between(Age_LDB["X2"]/1E3,climate_SynthesisLDB_WAPLS["pred.fit...1..3"] - climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..3"]-ModernMTCO_LDB,
                   climate_SynthesisLDB_WAPLS["pred.fit...1..3"] + climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..3"]-ModernMTCO_LDB, color='green', alpha=0.3)

ax[1,1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['mean.3']-ModernMTCO_LDB, color='red')
ax[1,1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['X0.95_inf.3']-ModernMTCO_LDB, climate_SynthesisLDB_CREST['X0.95_sup.3']-ModernMTCO_LDB, color='red', alpha=0.1)

ax[1,1].set_xlim(0, 40)
ax[1,1].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[1,1].set_ylabel("Mean Annual Temperature Anomaly (°C)", fontsize=15)
ax[1,1].tick_params(axis='both', which='major', labelsize=15)
ax[1,1].axvspan(19, 23, color='grey', alpha=0.3)

ax[1,1].set_xticks(np.arange(0, 41, 5))  
ax[1,1].set_xticklabels(np.arange(0, 41, 5), fontsize=15)

plt.tight_layout()
plt.savefig("LGP&LDBReconstruction_MensualPlot.pdf")
plt.show()

#%%Without PANN
fig, ax = plt.subplots(1, 2, figsize=(20, 7.5))  # 1 ligne, 2 colonnes

# --- Panneau a. LGP - MTWA ---
ax[0].set_title("a. LGP - MTWA", fontsize=16)
ax2 = ax[0].twinx()
ax2.plot(NGRIP["Age"]/1E3, NGRIP["180"], color='black', alpha=0.5, label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y', labelcolor='black', labelsize=15)
ax2.legend(bbox_to_anchor=(0.20, 0.23), fontsize=10)

ax[0].plot(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..2"] - ModernMTWA_LGP, color='blue', label='MAT')
ax[0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_MAT["predMAT.fit...1..2"] - climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..2"] - ModernMTWA_LGP, 
                   climate_SynthesisLGP_MAT["predMAT.fit...1..2"] + climate_SynthesisLGP_MAT["predMAT.SEP.boot...1..2"] - ModernMTWA_LGP, color='blue', alpha=0.3)

ax[0].plot(Age_LGP/1E3, climate_SynthesisLGP_WAPLS["pred.fit...1..2"] - ModernMTWA_LGP, color='green', label='WAPLS')
ax[0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_WAPLS["pred.fit...1..2"] - climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..2"] - ModernMTWA_LGP,
                   climate_SynthesisLGP_WAPLS["pred.fit...1..2"] + climate_SynthesisLGP_WAPLS["pred.SEP.boot...1..2"] - ModernMTWA_LGP, color='green', alpha=0.2)

ax[0].plot(Age_LGP/1E3, climate_SynthesisLGP_CREST['mean.2'] - ModernMTWA_LGP, color='red', label='CREST')
ax[0].fill_between(Age_LGP/1E3, climate_SynthesisLGP_CREST['X0.95_inf.2'] - ModernMTWA_LGP, climate_SynthesisLGP_CREST['X0.95_sup.2'] - ModernMTWA_LGP, color='red', alpha=0.1)

ax[0].set_xlim(0, 40)
ax[0].set_xlabel("", fontsize=15)
ax[0].set_ylabel("Mean Monthly Temperature Anomaly (°C)", fontsize=15)
ax[0].tick_params(axis='both', which='major', labelsize=15)
ax[0].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[0].set_xticks(np.arange(0, 41, 5))  
ax[0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[0].axvspan(19, 23, color='grey', alpha=0.3)
ax[0].legend(loc='lower left', fontsize=15)

# --- Panneau b. LDB - MTWA ---
ax[1].set_title("b. LDB - MTWA", fontsize=16)
ax2 = ax[1].twinx()
ax2.plot(NGRIP["Age"]/1E3, NGRIP["180"], color='black', alpha=0.5, label='NGRIP δ18O')
ax2.set_ylabel("NGRIP δ18O (‰)", fontsize=15, color='black')
ax2.tick_params(axis='y', labelcolor='black', labelsize=15)

ax[1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..2"] - ModernMTWA_LDB, color='blue')
ax[1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_MAT["predMAT.fit...1..2"] - climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..2"] - ModernMTWA_LDB, 
                   climate_SynthesisLDB_MAT["predMAT.fit...1..2"] + climate_SynthesisLDB_MAT["predMAT.SEP.boot...1..2"] - ModernMTWA_LDB, color='blue', alpha=0.3)

ax[1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_WAPLS["pred.fit...1..2"] - ModernMTWA_LDB, color='green')
ax[1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_WAPLS["pred.fit...1..2"] - climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..2"] - ModernMTWA_LDB,
                   climate_SynthesisLDB_WAPLS["pred.fit...1..2"] + climate_SynthesisLDB_WAPLS["pred.SEP.boot...1..2"] - ModernMTWA_LDB, color='green', alpha=0.2)

ax[1].plot(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['mean.2'] - ModernMTWA_LDB, color='red')
ax[1].fill_between(Age_LDB["X2"]/1E3, climate_SynthesisLDB_CREST['X0.95_inf.2'] - ModernMTWA_LDB, climate_SynthesisLDB_CREST['X0.95_sup.2'] - ModernMTWA_LDB, color='red', alpha=0.1)

ax[1].set_xlim(0, 40)
ax[1].set_ylabel("Mean Annual Temperature Anomaly (°C)", fontsize=15)
ax[1].tick_params(axis='both', which='major', labelsize=15)
ax[1].set_xlabel("Age cal BP (kyrs)", fontsize=15)
ax[1].axvspan(19, 23, color='grey', alpha=0.3)
ax[1].set_xticks(np.arange(0, 41, 5))  

plt.tight_layout()
plt.savefig("LGP_LDB_MensualPlot_2panels.pdf")
plt.show()

