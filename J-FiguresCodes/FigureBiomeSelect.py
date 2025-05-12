#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 12:07:46 2025

@author: gabrielfenisse
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

dataBiome_Best_LGP = pd.read_csv("Output_Recon_CREST_BiomeBest_SynthesisLGP.csv", sep=";", decimal=",")
dataBiome_Best_LDB = pd.read_csv("Output_Recon_CREST_BiomeBest_LDB.csv", sep=";", decimal=",")
dataBiome_W_LGP = pd.read_csv("Output_Recon_CREST_BiomeWeighted_SynthesisLGP.csv", sep=";", decimal=",")
dataBiome_W_LDB = pd.read_csv("Output_Recon_CREST_BiomeWeighted_LDB.csv", sep=";", decimal=",")

Age_LGP2 = pd.read_csv("ModeleAge_LGP2.txt", sep=",", decimal=".")  
Age_LGP = pd.read_csv("ModeleAge_LGP.txt", sep=",", decimal=".")  
dataBiome_W_LGP["Age"] = np.append(Age_LGP2["X2"][0:27], Age_LGP["X2"]) / 1E3
dataBiome_W_LGP["Age"] = np.append(Age_LGP2["X2"][0:27], Age_LGP["X2"]) / 1E3

Age_LDB = pd.read_csv("ModeleAge_LDB.txt", sep=",", decimal=".") / 1E3
dataBiome_Best_LDB["Age"] = Age_LDB["X2"]
dataBiome_Best_LDB["Age"] = Age_LDB["X2"]

biomeLGP_files = ["CLDE_LGP.csv", "TAIG_LGP.csv", "PION_LGP.csv", "TEDE_LGP.csv", "WAMX_LGP.csv","XERO_LGP.csv","TUND_LGP.csv","COST_LGP.csv","WAST_LGP.csv"]
biomeLGP2_files = ["CLDE_LGP2.csv", "TAIG_LGP2.csv", "PION_LGP2.csv", "TEDE_LGP2.csv", "WAMX_LGP2.csv","XERO_LGP2.csv","TUND_LGP2.csv","COST_LGP2.csv","WAST_LGP2.csv"]
biomeLGPSynthesis_data = {}
n=0

for file in range(len(biomeLGP2_files)):
    df = pd.read_csv("./LGP2_CREST/"+biomeLGP2_files[file], sep=";", decimal=",")
    df.rename(columns={"Unnamed: 1": "Ages", "Unnamed: 3": "mean","Unnamed: 3": "mean.2"}, inplace=True)
    df["Age"] = df["Ages"] / 1E3
    df=df[1:27]
    df2 = pd.read_csv("./LGP_CREST/"+biomeLGP_files[file], sep=";", decimal=",")
    df2.rename(columns={"Unnamed: 1": "Ages", "Unnamed: 3": "mean","Unnamed: 3": "mean.2"}, inplace=True)
    df2["Age"] = df2["Ages"] / 1E3
    df3=pd.concat([df, df2], axis=0, ignore_index=True) 
    biomeLGPSynthesis_data[biomeLGP_files[n].split("_")[0]] = df3 
    n+=1
   
biomeLDB_files = ["CLDE_LDB.csv", "TAIG_LDB.csv", "PION_LDB.csv", "TEDE_LDB.csv", "WAMX_LDB.csv","XERO_LDB.csv","TUND_LDB.csv","COST_LDB.csv","WAST_LDB.csv"]
biomeLDB_data = {}

for file in biomeLDB_files:
    df = pd.read_csv("./LDB_CREST/"+file, sep=";", decimal=",")
    df.rename(columns={"Unnamed: 1": "Ages", "Unnamed: 3": "mean","Unnamed: 3": "mean.2"}, inplace=True)
    df["Age"] = df["Ages"] / 1E3
    biomeLDB_data[file.split("_")[0]] = df 
    
dataBiome_Best_LGP["Ages"]=biomeLGPSynthesis_data["CLDE"]['Ages']
dataBiome_Best_LDB["Ages"]=df["Ages"]
dataBiome_W_LGP["Ages"]=biomeLGPSynthesis_data["CLDE"]['Ages']
dataBiome_W_LDB["Ages"]=df["Ages"]

#%%Plots
fig, ax = plt.subplots(2,1,figsize=(14, 12)) 
colors = {'CLDE':"red",'TAIG':"darkorchid",'PION':"aquamarine",'TEDE':'#117733','WAMX':'#F5F5DC','XERO':'orange','TUND':"#EE82E8",'COST':'powderblue','WAST':'grey'}

ModernTANN_LGP=10.7
ModernTANN_LDB=8.3
ModernPANN_LGP=1141.5
ModernPANN_LDB=1001.7

ax[0].set_title("a. La Grande Pile - CREST", fontsize=15)
for key, df in biomeLGPSynthesis_data.items():
    ax[0].plot(df["Ages"], df["mean"] - ModernTANN_LGP, label=key, color=colors[key], linewidth=2)
ax[0].plot(dataBiome_Best_LGP["Ages"], dataBiome_Best_LGP["mean"] - ModernTANN_LGP, label="Best", linestyle="-", linewidth=4, color="orange")
ax[0].plot(dataBiome_W_LGP["Ages"], dataBiome_W_LGP["mean"] - ModernTANN_LGP, label="Weighted", linestyle="-", linewidth=4, color="red")

ax[0].set_xlim(0, 40)
ax[0].set_ylabel("Temperature anomaly (°C)", fontsize=15)
ax[0].set_xlabel("")
ax[0].set_xticklabels(np.arange(0, 41, 5), fontsize=15)  
ax[0].set_yticklabels(np.arange(-25, 20, 5), fontsize=15)  
ax[0].axvspan(19, 23, color='grey', alpha=0.3)

ax[0].legend()

ax[1].set_title("b. Lac du Bouchet - CREST", fontsize=15)
for key, df in biomeLDB_data.items():
    ax[1].plot(df["Ages"], df["mean"] - ModernTANN_LDB, label=key, color=colors[key], linewidth=2)
ax[1].plot(dataBiome_Best_LDB["Ages"], dataBiome_Best_LDB["mean"] - ModernTANN_LDB , label="Best", linestyle="-", linewidth=4, color="orange")
ax[1].plot(dataBiome_W_LDB["Ages"], dataBiome_W_LDB["mean"] - ModernTANN_LDB, label="Weighted", linestyle="-", linewidth=4, color="red")

ax[1].set_xlim(0, 40)
ax[1].set_ylabel("Temperature anomaly (°C)", fontsize=15)
ax[1].set_xlabel("Age (kyr cal BP)" ,fontsize=15)
ax[1].set_xticklabels(np.arange(0, 41, 5), fontsize=15)
ax[1].set_yticklabels(np.arange(-20, 10, 5), fontsize=15)  
ax[1].axvspan(19, 23, color='grey', alpha=0.3)

plt.tight_layout()  
plt.savefig("WeightedEffect.pdf")
plt.show()

