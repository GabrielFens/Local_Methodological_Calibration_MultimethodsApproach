#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 20:25:03 2024

@author: gabrielfenisse
"""

import pandas as pd
import numpy as np
    
from mpl_toolkits.basemap import Basemap
from scipy.stats import linregress
import matplotlib.ticker as ticker
import netCDF4 as nc4
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
import xarray as xr
from matplotlib.colors import ListedColormap

#%%Map location
#Synthèse_Cooling
data_Pollen=pd.read_csv('PositionInt_CoolingPollen.csv', sep=";")
nylat=data_Pollen["Latitude"]
nylon=data_Pollen["Longitude"]
Num=data_Pollen["SiteNumber"]
RT=data_Pollen["Taxonomic_resolution"]
Elevation=data_Pollen["Elevation"]
Type=data_Pollen["Type"]
DT=data_Pollen["Diversity"]

IntType=np.unique(Type)
Type[Type==IntType[0]]=0
Type[Type==IntType[1]]=1
Type[Type==IntType[2]]=2
Type[Type==IntType[3]]=3
Type[Type==IntType[4]]=4

col_dict = {0: "blue",
            1: "red",
            2: "orange",
            3: "green",
            4: "black"}

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

labels = np.array(['Alluvial Fan', 'Cave', 'Colluvium', 'Lake', 'Peat Bog'])
len_lab = len(labels)

norm_bins = np.sort([*col_dict.keys()]) + 0.5
norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
print(norm_bins)

norm = mcolors.BoundaryNorm(norm_bins, len_lab, clip=True)
fmt = ticker.FuncFormatter(lambda x, pos: labels[norm(x)])

diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff / 2

#%%Places
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l')
plt.title("Fossil site locations")
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)

CS=plt.scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y = m(nylon, nylat)
Pollen_data = plt.scatter(x, y, s=1, color="black")
I=[8,9,35,36,40]
for i, txt in enumerate(Num):
    if txt in I :
        color = 'green'
    else:
        color = 'red'
    
    plt.text(x[i], y[i], str(txt),
             fontsize=6, ha='right', va='top', color=color)

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  

#%%Initilization using iLOVECLIM inputs
from matplotlib.patches import Polygon

Ice[Ice<100]=float('nan')

fig, axs = plt.subplots(1,3,figsize=(18,4))
#fig.suptitle('Pollen data description')

axs[0].set_title("a. Site locations")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0])
plt.title("Fossil site locations")
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)

CS=axs[0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y = m(nylon, nylat)
Pollen_data = axs[0].scatter(x, y, s=1, color="black")
I=[8,9,35,36,40]
for i, txt in enumerate(Num):
    if txt in I :
        color = 'green'
    else:
        color = 'red'
    
    axs[0].text(x[i], y[i], str(txt),
             fontsize=8, ha='right', va='top', color=color)

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  

axs[1].set_title("b. Type of matrix")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')
CS=axs[1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
    
x, y= m(nylon,nylat)
Pollen_data = axs[1].scatter(x,y,c=Type,s=30,linewidths=1,edgecolors='k',cmap=mpl.cm.Accent,norm=norm)    
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
fig.colorbar(Pollen_data,format=fmt,ticks=tickz)

axs[2].set_title("c. Pollen diversity")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2])
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
CS=axs[2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data =axs[2].scatter(x,y,c=DT,s=30,linewidths=1,edgecolors='k',cmap=plt.get_cmap('tab20'),vmin=0,vmax=80)    

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01) 
fig.colorbar(Pollen_data,ticks=np.arange(0,80,8),label="Taxonomic diversity") 

#m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),ax=axs[2])
#lon2, lat2 = np.meshgrid(lon,lat)
#x, y = m(lon2, lat2)
#m.drawcoastlines()
#m.drawmapboundary(fill_color='lightblue')
#m.fillcontinents(color='lightgrey',lake_color='lightblue')
#CS=axs[2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

#x, y= m(nylon,nylat)
#Pollen_data = axs[2].scatter(x,y,c=RT,linewidths=1,edgecolors='k',cmap=plt.get_cmap('Paired',np.max(RT)-np.min(RT)+1),vmin=np.min(RT)-0.5,vmax=np.max(RT)+0.5)    
#m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
#m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
##fig.colorbar(Pollen_data,ticks=np.arange(np.min(RT),np.max(RT)+1),label="Taxonomic resolution")
#plt.tight_layout()

plt.savefig("PollenDescription.pdf")
plt.show()

#%%AP vs. NAP diagram
from collections import OrderedDict
data_Pollen=pd.read_csv('CoolingPollenAP.csv', sep=";", decimal=",")
AP=["Betula","Betula alba","Acacia","Pinus","Asteraceae","Fraxinus ornus","Fraxinus excelsior","Asteraceae (Asteroideae)","Cichorioideae","Compositae","Ilex","Iridaceae","Juglans","Pinus subgen. Diploxylon","Pinus subgen. Haploxylon","Plumbaginaceae","Polemonium","Platanus","Pistacia","Quercus (deciduous)","Quercus (evergreen)","Convolvulaceae","Cornus","Cistus","Cedrus","Buxus","Asteroideae","Asteraceae (Cichorioideae)","Corylus","Ulmus","Alnus","Alnus fruticosa","Alnus glutinosa type","Arbutus","Tilia","Fraxinus","Salix",'Hedera','Acer','Abies',"Carpinus","Fagus","Larix",'Picea','Rosaceae',"Apiaceae"]
Shrub=["Juniperus","Betula nana","Caprifoliaceae","Ephedra","Hippophae","Myrica",'Viscum',"Viburnum","Ephedra fragilis","Euonymus","Myrtus","Nitraria","Euphorbiaceae"]
NAP=["Poaceae","Allium","Cyperaceae","Artemisia","Cannabis","Geraniaceae","Calluna","Dryas","Empetrum","Dipsacaceae",'Frangula',"Crassulaceae","Lonicera","Phillyrea","Pedicularis","Olea","Papaveraceae","Ostrya/Carpinus orientalis","Populus","Polygonaceae","Pterocarya","Onagraceae","Linaria","Liliaceae",'Brassicaceae',"Ranunculus","Rhododendron","Rubiaceae","Rhus","Rhamnus","Rubus","Rubus chamaemorus","Saxifraga","Saxifragaceae","Scrophulariaceae",'Zizyphus',"Urtica",'Vitis','Valerianaceae',"Vaccinium","Tamaricaceae","Taxus","Sanguisorba","Sambucus","Rutaceae","Rumex","Ranunculaceae","Castanea","Ceratonia",'Caryophyllaceae',"Chenopodiaceae","Helianthemum","Thalictrum",'Armeria',"Boraginaceae","Campanulaceae","Ericaceae","Fabaceae","Dipsaceae","Lamiaceae","Gentianaceae",'Centaurea','Lamiaceae',"Plantago","Polygonum","Sanguisorba minor","Urticaceae","Filipendula","Sanguisorba officinalis"]
Name=list(OrderedDict.fromkeys(data_Pollen['SiteName']))
Mean_per_site=data_Pollen.groupby('SiteName').mean()

def rename_columns(df,AP,NAP,Shrub):
    df.rename(columns={col: '1. AP' for col in df.columns if col in AP}, inplace=True)
    df.rename(columns={col: '2. Shrub' for col in df.columns if col in Shrub}, inplace=True)
    df.rename(columns={col: '3. NAP' for col in df.columns if col in NAP}, inplace=True)
    return df

df_renamed = rename_columns(Mean_per_site,AP,NAP,Shrub)
df_summed = df_renamed.groupby(df_renamed.columns, axis=1).sum()

#df_summed.drop(columns=["Unnamed: 120"], inplace=True)
df_normalized = df_summed.div(df_summed.sum(axis=1), axis=0) * 100
df_sorted = df_normalized.loc[Name]
ax = df_sorted.plot(kind='barh', stacked=True, figsize=(14, 10), colormap="tab20")
ax.set_title("d. Percentage of AP/NAP")

#ax.set_yticklabels(Name)
ax.set_xlabel("Fraction (%)")
plt.legend(loc="upper right")
plt.savefig("PolliniqueDiagram.pdf")

#%%NGRIP_Ramusen
NGRIP=pd.read_csv('NGRIP_Records.txt', sep="\t")
plt.plot(NGRIP["Age"]/1E3,NGRIP["180"],color='black')
plt.xlim(0,40)
plt.savefig("NGRIP.pdf")
plt.show()

#%%SynthesePollenSequences
#Output
data_Megabiome_fos=pd.read_csv('BiomeBest_COolingPollen.csv', sep=";")
Biome=data_Megabiome_fos["Biome_Best"]
MegaBiome=data_Megabiome_fos["Megabiome_Best"]
lon_data=data_Megabiome_fos["Longitude"]
lat_data=data_Megabiome_fos["Latitude"]

colors=['#65B2FF',"#117733",'coral',"bisque","red",'#EE82E8']
#%%Biome
Biome[Biome=="TAIG"]=1.
Biome[Biome=="PION"]=2.
Biome[Biome=="CLMX"]=3.
Biome[Biome=="COCO"]=4.
Biome[Biome=="TEDE"]=5.
Biome[Biome=="COMX"]=6.
Biome[Biome=="WAMX"]=7.
Biome[Biome=="COST"]=8.
Biome[Biome=="WAST"]=9.
Biome[Biome=="TUND"]=10.

col_dict = {1: "darkorchid",
            2: 'aquamarine',
            3: 'coral',
            4: "#65B2FF",
            5: "#117733",
            6: "#CAFFCA",
            7: "#F5F5DC",
            8: 'powderblue',
            9: 'grey',
            10: '#EE82E8'}


labels = np.array(['TAIG', 'PION', 'CLMX', 'COCO', 'TEDE','COMX',"WAMX","COST","WAST",'TUND'])
len_lab = len(labels)

norm_bins = np.sort([*col_dict.keys()]) + 0.5
norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
print(norm_bins)

norm = mcolors.BoundaryNorm(norm_bins, len_lab, clip=True)
fmt = ticker.FuncFormatter(lambda x, pos: labels[norm(x)])

diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff / 2

colors=["darkorchid",'aquamarine','coral',"#65B2FF", "#117733","#CAFFCA","#F5F5DC",'powderblue','grey','#EE82E8']
cmap=mcolors.LinearSegmentedColormap.from_list("mycmap1", colors,N=10)

#%%Megabiome
MegaBiome[MegaBiome=="WTFO"]=1.
MegaBiome[MegaBiome=="TEFO"]=2.
MegaBiome[MegaBiome=="BOFO"]=3.
MegaBiome[MegaBiome=="SAVA"]=4.
MegaBiome[MegaBiome=="STEP"]=5.
MegaBiome[MegaBiome=="TUND"]=6.

col_dict = {1: '#65B2FF',
            2: "#117733",
            3: "coral",
            4: 'bisque',
            5: 'red',
            6: '#EE82E8'}

labels_a = np.array(['WTFO', 'TEFO', 'BOFO', 'SAVA', 'STEP','TUND'])
len_lab_a = len(labels_a)

norm_bins_a = np.sort([*col_dict.keys()]) + 0.5
norm_bins_a = np.insert(norm_bins_a, 0, np.min(norm_bins_a) - 1.0)
print(norm_bins)

norm_a = mcolors.BoundaryNorm(norm_bins_a, len_lab_a, clip=True)
fmt1 = ticker.FuncFormatter(lambda x, pos: labels_a[norm_a(x)])

diff_a = norm_bins_a[1:] - norm_bins_a[:-1]
tickz_a = norm_bins_a[:-1] + diff_a / 2

colors_a=['#65B2FF',"#117733","coral", 'bisque','red','#EE82E8']
cmap_a = mcolors.LinearSegmentedColormap.from_list("mycmap1", colors_a,N=6)

#%%Change into c=Biome --> MegaBiome
Ice[Ice<100]=float('nan')

fig, axs = plt.subplots(1,2,figsize=(10,10))

axs[0].set_title("a. Biome")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=axs[0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
               #linewidths=1,edgecolors='k',color="white")
               #levels=[10,50,90],colors=["blue","blue","blue"])
x, y= m(lon_data,lat_data)
Pollen_data = axs[0].scatter(x,y,c=Biome,linewidths=1,edgecolors='k',cmap=cmap,norm=norm)
#data_biome["Best"]
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  

fig.colorbar(Pollen_data,format=fmt,ticks=tickz,label="LGM dominant biome",shrink=0.35)

axs[1].set_title("b. Megabiome")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=axs[1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
               #linewidths=1,edgecolors='k',color="white")
               #levels=[10,50,90],colors=["blue","blue","blue"])
x, y= m(lon_data,lat_data)
Pollen_data = axs[1].scatter(x,y,c=MegaBiome,linewidths=1,edgecolors='k',cmap=cmap_a,norm=norm_a)
#data_biome["Best"]
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
fig.colorbar(Pollen_data, format=fmt1, ticks=tickz_a, label="LGM dominant biome", shrink=0.35)

plt.tight_layout()
plt.savefig("(Mega)biome_Map.pdf")
plt.show()

#%%
#edgecolors intialilization
MegaBiome[MegaBiome!=2]=float('nan')

#%%
Megabiome_fos=pd.read_csv('BiomeBest_COolingPollen_PercentageMegabiomizationProcess.csv', sep=";")

import matplotlib.colors as mcolors

Ice[Ice<100]=float('nan')

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=plt.scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
               #linewidths=1,edgecolors='k',color="white")
               #levels=[10,50,90],colors=["blue","blue","blue"])
x, y = m(lon2, lat2)
cmap_edge = plt.get_cmap('RdBu_r') 

TotalBiome_Extreme=Megabiome_fos["STEP"]+Megabiome_fos["SAVA"]+Megabiome_fos["TUND"]+Megabiome_fos["DESE"]
x, y= m(lon_data,lat_data)
Pollen_data = plt.scatter(x,y,c=TotalBiome_Extreme,s=30,vmin=0,vmax=60,linewidths=1,cmap=plt.get_cmap('Spectral', 6))

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01) 

plt.colorbar(label="Arid megabiomes at the LGM (%)")

is_TEDE = MegaBiome == 2  
plt.scatter(x[is_TEDE], y[is_TEDE], facecolors='none', edgecolors='black', linewidths=1, s=80, label="Dominant TEFO")
plt.legend(fontsize=8)
plt.savefig("AridMegabiomes.pdf")
plt.show()

#%%LGP and LDB position
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=plt.scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
x, y = m(lon2, lat2)
cmap_edge = plt.get_cmap('RdBu_r') 

x, y= m(lon_data[[22,34]],lat_data[[22,34]])
Pollen_data = plt.scatter(x,y,s=30)  

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01) 

plt.savefig("LGP&LDB_Position.pdf")
plt.show()


#%%I-BIOMES COMPARAISONS
#%%Reconstitutions
data_Pollen=pd.read_csv('Test_CoolingPollen.csv', sep=";")
nylat=data_Pollen["Latitude"]
nylon=data_Pollen["Longitude"]

delta_TANN_MAT=data_Pollen["TANN_MAT_Wbiome"]-data_Pollen["TANN_Modern"]
delta_TANN_WAPLS=data_Pollen["TANN_WAPLS_Wbiome"]-data_Pollen["TANN_Modern"]
delta_TANN_CREST=data_Pollen["TANN_CREST_Wbiome"]-data_Pollen["TANN_Modern"]

delta_TANN_Mega_MAT=data_Pollen["TANN_MAT_WMegabiome"]-data_Pollen["TANN_Modern"]
delta_TANN_Mega_WAPLS=data_Pollen["TANN_WAPLS_WMegabiome"]-data_Pollen["TANN_Modern"]
delta_TANN_Mega_CREST=data_Pollen["TANN_CREST_WMegabiome"]-data_Pollen["TANN_Modern"]

delta_PANN_MAT=(data_Pollen["PANN_MAT_Wbiome"]/data_Pollen["PANN_Modern"])-1
delta_PANN_WAPLS=(data_Pollen["PANN_WAPLS_Wbiome"]/data_Pollen["PANN_Modern"])-1
delta_PANN_CREST=(data_Pollen["PANN_CREST_Wbiome"]/data_Pollen["PANN_Modern"])-1

delta_PANN_Mega_MAT=(data_Pollen["PANN_MAT_WMegabiome"]/data_Pollen["PANN_Modern"])-1
delta_PANN_Mega_WAPLS=(data_Pollen["PANN_WAPLS_WMegabiome"]/data_Pollen["PANN_Modern"])-1
delta_PANN_Mega_CREST=(data_Pollen["PANN_CREST_WMegabiome"]/data_Pollen["PANN_Modern"])-1

delta_MTWA_MAT=data_Pollen["MTWA_MAT_Wbiome"]-data_Pollen["MTWA_Modern"]
delta_MTWA_WAPLS=data_Pollen["MTWA_WAPLS_Wbiome"]-data_Pollen["MTWA_Modern"]
delta_MTWA_CREST=data_Pollen["MTWA_CREST_Wbiome"]-data_Pollen["MTWA_Modern"]

delta_MTCO_MAT=data_Pollen["MTCO_MAT_Wbiome"]-data_Pollen["MTCO_Modern"]
delta_MTCO_WAPLS=data_Pollen["MTCO_WAPLS_Wbiome"]-data_Pollen["MTCO_Modern"]
delta_MTCO_CREST=data_Pollen["MTCO_CREST_Wbiome"]-data_Pollen["MTCO_Modern"]

delta_MTWA_Mega_MAT=data_Pollen["MTWA_MAT_WMegabiome"]-data_Pollen["MTWA_Modern"]
delta_MTWA_Mega_WAPLS=data_Pollen["MTWA_WAPLS_WMegabiome"]-data_Pollen["MTWA_Modern"]
delta_MTWA_Mega_CREST=data_Pollen["MTWA_CREST_WMegabiome"]-data_Pollen["MTWA_Modern"]

delta_MTCO_Mega_MAT=data_Pollen["MTCO_MAT_WMegabiome"]-data_Pollen["MTCO_Modern"]
delta_MTCO_Mega_WAPLS=data_Pollen["MTCO_WAPLS_WMegabiome"]-data_Pollen["MTCO_Modern"]
delta_MTCO_Mega_CREST=data_Pollen["MTCO_CREST_WMegabiome"]-data_Pollen["MTCO_Modern"]

#%%Climate reconstruction - Weighted biome
#lon=lon_lgm
#lat=lat_lgm
fig, axs = plt.subplots(2,3,figsize=(12,7))
#fig.suptitle('Annual climate anomalies - Weighted biome')

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. MAT from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=delta_TANN_MAT,linewidths=1,edgecolors='k',vmax=0,vmin=-10,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,1].set_title("b. WAPLS from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=delta_TANN_WAPLS,s=30,linewidths=1,vmax=0,vmin=-10,edgecolors='k',cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,2].set_title("c. CREST from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=delta_TANN_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-10.0,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.54, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data, cax=cbar_ax,label="TANN anomaly (°C)") 
cm_a.set_ticks(np.arange(-10,2,2))


axs[1,0].set_title("d. MAT from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=delta_TANN_Mega_MAT,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-10.0,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,1].set_title("e. WAPLS from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=delta_TANN_Mega_WAPLS,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-10.0,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,2].set_title("f. CREST from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=delta_TANN_Mega_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-10.0,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.13, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data, cax=cbar_ax,label="TANN anomaly (°C)")
cm_a.set_ticks(np.arange(-10,2,2))

plt.savefig("ClimateReconstructionAnnual_W(Mega)Biome.pdf")
plt.show()

#%%SIAnomaly-Mensual

fig, axs = plt.subplots(3,3,figsize=(12,12))
#fig.suptitle('Monthly climate anomalies  - Weighted biome')

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. MAT - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=delta_MTWA_MAT,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,1].set_title("b. WAPLS - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=delta_MTWA_WAPLS,s=30,linewidths=1,vmax=0,vmin=-16,edgecolors='k',cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,2].set_title("c. CREST - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=delta_MTWA_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.66, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="MTWA anomaly (°C)") 
cm_a.set_ticks(np.arange(-16,2,2))


axs[1,0].set_title("d. MAT - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=delta_MTCO_MAT,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,1].set_title("e. WAPLS - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=delta_MTCO_WAPLS,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,2].set_title("f. CREST - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=delta_MTCO_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.40, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="MTCO anomaly (°C)")
cm_a.set_ticks(np.arange(-16,2,2))


axs[2,0].set_title("g. MAT - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,0].scatter(x,y,c=delta_MTWA_MAT-delta_MTCO_MAT,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[2,1].set_title("h. WAPLS - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,1].scatter(x,y,c=delta_MTWA_WAPLS-delta_MTCO_WAPLS,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[2,2].set_title("i. CREST - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,2].scatter(x,y,c=delta_MTWA_CREST-delta_MTCO_CREST,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

cbar_ax = fig.add_axes([0.92, 0.13, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="Seasonal index anomaly (°C)")
cm_a.set_ticks(np.arange(-8,10,2))

plt.savefig("ClimateReconstructionMensual_WBiome.pdf")
plt.show()

#%%SIAnomaly-Mensual

fig, axs = plt.subplots(3,3,figsize=(12,12))
#fig.suptitle('Monthly climate anomalies  - Weighted biome')

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. MAT - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=delta_MTWA_Mega_MAT,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,1].set_title("b. WAPLS - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=delta_MTWA_Mega_WAPLS,s=30,linewidths=1,vmax=0,vmin=-16,edgecolors='k',cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,2].set_title("c. CREST - MTWA")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=delta_MTWA_Mega_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.66, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="MTWA anomaly (°C)") 
cm_a.set_ticks(np.arange(-16,2,2))


axs[1,0].set_title("d. MAT - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=delta_MTCO_Mega_MAT,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,1].set_title("e. WAPLS - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=delta_MTCO_Mega_WAPLS,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,2].set_title("f. CREST - MTCO")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=delta_MTCO_Mega_CREST,s=30,linewidths=1,edgecolors='k',vmax=0,vmin=-16,cmap=plt.get_cmap('Blues_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.40, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="MTCO anomaly (°C)")
cm_a.set_ticks(np.arange(-16,2,2))


axs[2,0].set_title("g. MAT - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,0].scatter(x,y,c=delta_MTWA_Mega_MAT-delta_MTCO_Mega_MAT,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[2,1].set_title("h. WAPLS - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,1].scatter(x,y,c=delta_MTWA_Mega_WAPLS-delta_MTCO_Mega_WAPLS,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[2,2].set_title("i. CREST - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[2,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[2,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[2,2].scatter(x,y,c=delta_MTWA_Mega_CREST-delta_MTCO_Mega_CREST,s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

cbar_ax = fig.add_axes([0.92, 0.13, 0.015, 0.22])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="Seasonal index anomaly (°C)")
cm_a.set_ticks(np.arange(-8,10,2))

plt.savefig("ClimateReconstructionMensual_WMegabiome.pdf")
plt.show()

#%%Mean climate variables from 3 methods
fig, axs = plt.subplots(2,2,figsize=(12,8))
#fig.suptitle('Mean climate anomaly (LGM-Modern Climate) - Weighted biome')

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. TANN anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=np.mean([delta_TANN_MAT,delta_TANN_WAPLS,delta_TANN_CREST],axis=0),linewidths=1,edgecolors='k',s=30,vmax=0,vmin=-10,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="TANN anomaly (°C)")

axs[0,1].set_title("b. Sesonal index anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=np.mean([delta_MTWA_MAT-delta_MTCO_MAT,delta_MTWA_WAPLS-delta_MTCO_WAPLS,delta_MTWA_CREST-delta_MTCO_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="Seasonal Index Anomaly (°C)") 

axs[1,0].set_title("c. TANN anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=np.mean([delta_TANN_Mega_MAT,delta_TANN_Mega_WAPLS,delta_TANN_Mega_CREST],axis=0),linewidths=1,edgecolors='k',s=30,vmax=0,vmin=-10,cmap=plt.get_cmap('Blues_r',10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="TANN anomaly (°C)")

axs[1,1].set_title("d. Sesonal index anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=np.mean([delta_MTWA_Mega_MAT-delta_MTCO_Mega_MAT,delta_MTWA_Mega_WAPLS-delta_MTCO_Mega_WAPLS,delta_MTWA_Mega_CREST-delta_MTCO_Mega_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmax=8,vmin=-8,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="Seasonal Index Anomaly (°C)") 

plt.tight_layout()
plt.savefig("ClimateReconstructionMean_W(Mega)Biome.pdf")
plt.show()

#%%STD
fig, axs = plt.subplots(2,3,figsize=(12,8))
fig.suptitle('STD climate anomaly (LGM - Modern Climate) - Weighted biome')

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. TANN anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=np.std([delta_TANN_MAT,delta_TANN_WAPLS,delta_TANN_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmin=0,vmax=1.4,cmap=plt.get_cmap('Blues_r',7))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="TANN anomaly (°C)")

axs[1,0].set_title("d. TANN anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=np.std([delta_TANN_Mega_MAT,delta_TANN_Mega_WAPLS,delta_TANN_Mega_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmin=0,vmax=1.4,cmap=plt.get_cmap('Blues_r',7))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="TANN anomaly (°C)")

axs[0,1].set_title("b. PANN anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=np.std([delta_PANN_MAT*100,delta_PANN_WAPLS*100,delta_PANN_CREST*100],axis=0),s=30,linewidths=1,vmax=16,vmin=0,edgecolors='k',cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="PANN anomaly (%)")

axs[1,1].set_title("e. PANN anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=np.std([delta_PANN_Mega_MAT*100,delta_PANN_Mega_WAPLS*100,delta_PANN_Mega_CREST*100],axis=0),s=30,linewidths=1,vmax=16,vmin=0,edgecolors='k',cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="PANN anomaly (%)")

axs[0,2].set_title("c. Sesonal Index anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=np.std([delta_MTWA_MAT-delta_MTCO_MAT,delta_MTWA_WAPLS-delta_MTCO_WAPLS,delta_MTWA_CREST-delta_MTCO_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmax=4,vmin=0,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="Seasonal Index Anomaly (°C)") 

axs[1,2].set_title("f. Sesonal Index anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=np.std([delta_MTWA_Mega_MAT-delta_MTCO_Mega_MAT,delta_MTWA_Mega_WAPLS-delta_MTCO_Mega_WAPLS,delta_MTWA_Mega_CREST-delta_MTCO_Mega_CREST],axis=0),s=30,linewidths=1,edgecolors='k',vmax=4,vmin=0,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="Seasonal Index Anomaly (°C)") 

plt.tight_layout()
plt.savefig("ClimateReconstructionSTD_W(Mega)Biome.pdf")
plt.show()

#%%PANN

fig, axs = plt.subplots(2,3,figsize=(12,7))

Ice[Ice<100]=float('nan')

axs[0,0].set_title("a. MAT from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=delta_PANN_MAT*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,1].set_title("b. WAPLS from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=delta_PANN_WAPLS*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,2].set_title("c. CREST from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=delta_PANN_CREST*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.54, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data, cax=cbar_ax,label="PANN anomaly (%)") 
cm_a.set_ticks(np.arange(-100,100,20))


axs[1,0].set_title("d. MAT from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=delta_PANN_Mega_MAT*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,1].set_title("e. WAPLS from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=delta_PANN_Mega_WAPLS*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,2].set_title("f. CREST from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=delta_PANN_Mega_CREST*100,linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.13, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data, cax=cbar_ax,label="PANN anomaly (%)")
cm_a.set_ticks(np.arange(-100,100,20))

plt.savefig("ClimateReconstructionPANNMethods_W(Mega)Biome.pdf")
plt.show()


fig, axs = plt.subplots(1,2,figsize=(15,5))

axs[0].set_title("a. PANN anomaly from biomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0].scatter(x,y,c=np.mean([delta_PANN_MAT*100,delta_PANN_WAPLS*100,delta_PANN_CREST*100],axis=0),linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="PANN anomaly (%)")

axs[1].set_title("b. PANN anomaly from megabiomes")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1].scatter(x,y,c=np.mean([delta_PANN_Mega_MAT*100,delta_PANN_Mega_WAPLS*100,delta_PANN_Mega_CREST*100],axis=0),linewidths=1,s=30,vmax=100,vmin=-100,edgecolors='k',cmap=plt.get_cmap('bwr_r',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
m.colorbar(Pollen_data,label="PANN anomaly (%)")

plt.tight_layout()
plt.savefig("ClimateReconstructionPANN_W(Mega)Biome.pdf")
plt.show()

#%%Difference between biomes and megabiomes - Scatter plot

fig, ax = plt.subplots(2,1,figsize=(20, 12))

ax[0].set_title("a. TANN differences between biomes - megabiomes")
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=1)

x_positions = np.arange(len(clean_Name))
offset = 0.1

values_crest = delta_TANN_CREST-delta_TANN_Mega_CREST
values_mat = delta_TANN_MAT-delta_TANN_Mega_MAT
values_wapls = delta_TANN_WAPLS-delta_TANN_Mega_WAPLS

ax[0].errorbar(x_positions - offset, values_crest, 
            yerr=np.sqrt((data_Pollen["TANN_CREST_Wbiome_Sigma"]**2)+(data_Pollen["TANN_CREST_WMegabiome_Sigma"]**2)), color='red', alpha=0.25, markerfacecolor='red', fmt='o', elinewidth=1)
ax[0].errorbar(x_positions, values_mat, 
            yerr=np.sqrt(data_Pollen["TANN_MAT_Wbiome_Sigma"]**2+data_Pollen["TANN_MAT_WMegabiome_Sigma"]**2), color='blue', alpha=0.25, markerfacecolor='blue', fmt='o', elinewidth=1)
ax[0].errorbar(x_positions + offset, values_wapls, 
            yerr=np.sqrt(data_Pollen["TANN_WAPLS_Wbiome_Sigma"]**2+data_Pollen["TANN_WAPLS_WMegabiome_Sigma"]**2), color='green', fmt='o', markerfacecolor='green', elinewidth=1, alpha=0.25)

ax[0].scatter(x_positions - offset, values_crest, color='red', label="CREST")
ax[0].scatter(x_positions, values_mat, color='blue', label="MAT")
ax[0].scatter(x_positions + offset, values_wapls, color='green', label="WAPLS")

ax[0].set_ylabel("TANN diff (biomes-megabiomes) (°C)")
ax[0].legend()
ax[0].set_xticks([])


ax[1].set_title("b. SI differences between biomes - megabiomes")
ax[1].axhline(y=0, color='black', linestyle='--', linewidth=1)

x_positions = np.arange(len(clean_Name))
offset = 0.1

values_crest = delta_MTWA_CREST-delta_MTCO_CREST-(delta_MTWA_Mega_CREST-delta_MTCO_Mega_CREST)
values_mat = delta_MTWA_MAT-delta_MTCO_MAT-(delta_MTWA_Mega_MAT-delta_MTCO_Mega_MAT)
values_wapls = delta_MTWA_WAPLS-delta_MTCO_WAPLS-(delta_MTWA_Mega_WAPLS-delta_MTCO_Mega_WAPLS)

ax[1].errorbar(x_positions - offset, values_crest, 
            yerr=np.sqrt(data_Pollen["MTWA_CREST_Wbiome_Sigma"]**2+data_Pollen["MTWA_CREST_WMegabiome_Sigma"]**2+data_Pollen["MTCO_CREST_Wbiome_Sigma"]**2+data_Pollen["MTCO_CREST_WMegabiome_Sigma"]**2), color='red', alpha=0.25, markerfacecolor='red', fmt='o', elinewidth=1)
ax[1].errorbar(x_positions, values_mat, 
            yerr=np.sqrt(data_Pollen["MTWA_MAT_Wbiome_Sigma"]**2+data_Pollen["MTWA_MAT_WMegabiome_Sigma"]**2+data_Pollen["MTCO_MAT_Wbiome_Sigma"]**2+data_Pollen["MTCO_MAT_WMegabiome_Sigma"]**2), color='blue', alpha=0.25, markerfacecolor='blue', fmt='o', elinewidth=1)
ax[1].errorbar(x_positions + offset, values_wapls, 
            yerr=np.sqrt(data_Pollen["MTWA_WAPLS_Wbiome_Sigma"]**2+data_Pollen["MTWA_WAPLS_WMegabiome_Sigma"]**2+data_Pollen["MTCO_WAPLS_Wbiome_Sigma"]**2+data_Pollen["MTCO_WAPLS_WMegabiome_Sigma"]**2), color='green', fmt='o', markerfacecolor='green', elinewidth=1, alpha=0.25)

ax[1].scatter(x_positions - offset, values_crest, color='red', label="CREST")
ax[1].scatter(x_positions, values_mat, color='blue', label="MAT")
ax[1].scatter(x_positions + offset, values_wapls, color='green', label="WAPLS")

ax[1].set_ylabel("SI diff (biomes-megabiomes) (°C)")
ax[1].set_xticks(range(len(clean_Name)))
ax[1].set_xticklabels(clean_Name, rotation=90)

plt.tight_layout()
plt.savefig("ClimateReconstruction_Diff_W(Mega)Biome.pdf")
plt.show()

#%%Map diff
fig, axs = plt.subplots(2,3,figsize=(12,8))

axs[0,0].set_title("a. MAT - TANN")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,0].scatter(x,y,c=delta_TANN_MAT-delta_TANN_Mega_MAT,linewidths=1,edgecolors='k',vmax=4,vmin=-4,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,1].set_title("b. WAPLS - TANN")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,1].scatter(x,y,c=delta_TANN_WAPLS-delta_TANN_Mega_WAPLS,s=30,linewidths=1,vmax=4,vmin=-4,edgecolors='k',cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[0,2].set_title("c. CREST - TANN")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[0,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[0,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[0,2].scatter(x,y,c=delta_TANN_CREST-delta_TANN_Mega_CREST,s=30,linewidths=1,edgecolors='k',vmax=4,vmin=-4,cmap=plt.get_cmap('bwr',8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)
cbar_ax = fig.add_axes([0.92, 0.54, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data, cax=cbar_ax,label="TANN diff (biomes-megabiomes) (°C)") 
cm_a.set_ticks(np.arange(-4,6,2))

axs[1,0].set_title("d. MAT - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,0])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,0].scatter(x,y,c=delta_MTWA_MAT-delta_MTCO_MAT-(delta_MTWA_Mega_MAT-delta_MTCO_Mega_MAT),s=30,linewidths=1,edgecolors='k',vmax=4,vmin=-4,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,1].set_title("e. WAPLS - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,1])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,1].scatter(x,y,c=delta_MTWA_WAPLS-delta_MTCO_WAPLS-(delta_MTWA_Mega_WAPLS-delta_MTCO_Mega_WAPLS),s=30,linewidths=1,edgecolors='k',vmax=4,vmin=-4,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

axs[1,2].set_title("f. CREST - SI")

m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l',ax=axs[1,2])
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
axs[1,2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1,2].scatter(x,y,c=delta_MTWA_CREST-delta_MTCO_CREST-(delta_MTWA_Mega_CREST-delta_MTCO_Mega_CREST),s=30,linewidths=1,edgecolors='k',vmax=4,vmin=-4,cmap=plt.get_cmap('bwr', 8))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)

cbar_ax = fig.add_axes([0.92, 0.13, 0.02, 0.33])
cm_a=fig.colorbar(Pollen_data,cax=cbar_ax,label="Seasonal diff (biomes-megabiomes) (°C)")
cm_a.set_ticks(np.arange(-4,6,2))

#plt.tight_layout()
plt.savefig("ClimateReconstruction_Diff_W(Mega)Biome.pdf")
plt.show()


#%%Formatage Name
Name=np.array(Name)
Name[Name=="Navarres-2"]=float('nan')
Name[Name=="Navarres-1"]="Navarres"
clean_Name = Name[Name != 'nan']

#%%Diff between our climate reconstructions (not W and weighted) and those Davis et al., 2024 (MAT)
dataComparaison_Davis=pd.read_csv('DavisComp.csv',sep=";",decimal=",")
lat_data=dataComparaison_Davis["Latitude"]
lon_data=dataComparaison_Davis["Longitude"]
TANN_LGM_Davis=dataComparaison_Davis["TANN_LGM_DAVIS"]
PANN_LGM_Davis=dataComparaison_Davis["PANN_LGM_DAVIS"]
MTWA_LGM_Davis=dataComparaison_Davis["MTWA_LGM_DAVIS"]
MTCO_LGM_Davis=dataComparaison_Davis["MTCO_LGM_DAVIS"]
delta_TANN_MAT_WBiome=dataComparaison_Davis["TANN_Biomization_Weighted"]
delta_PANN_MAT_WBiome=dataComparaison_Davis["PANN_Biomization_Weighted"]
delta_MTWA_MAT_WBiome=dataComparaison_Davis["MTWA_Biomization_Weighted"]
delta_MTCO_MAT_WBiome=dataComparaison_Davis["MTCO_Biomization_Weighted"]

#Cooling_TANN_MAT_Davis=dataComparaison_Davis["LGMCooling_Davis"]
#Cooling_TANN_MAT_WMegaBiome=dataComparaison_Davis["LGMCooling_Megabiomization_Weighted"]
#Cooling_TANN_MAT_WBiome=dataComparaison_Davis["LGMCooling_Biomization_Weighted"]

#%%
fig, axs = plt.subplots(1,3,figsize=(15,10))
#fig.suptitle('LGM climate anomalies using MAT approach between reconstructed values from this study and those of Davis et al., (2024)',fontweight='bold')

bounds = [-10, -5, -2, 0, 2, 5, 10]  
norm = mcolors.BoundaryNorm(bounds, plt.get_cmap('coolwarm').N)

axs[0].set_title("a. TANN")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),ax=axs[0],resolution='l')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=axs[0].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
Diff1=delta_TANN_MAT_WBiome-TANN_LGM_Davis
x, y= m(lon_data,lat_data)
Pollen_data = axs[0].scatter(x,y,c=Diff1,s=40,linewidths=1,edgecolors='k',cmap=plt.get_cmap('coolwarm'),norm=norm)
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
m.colorbar(Pollen_data,ax=axs[0],label="Temperature anomaly (°C)")

axs[1].set_title("b. MTWA")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),ax=axs[1],resolution='l')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=axs[1].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
Diff3=delta_MTWA_MAT_WBiome-MTWA_LGM_Davis
x, y= m(lon_data,lat_data)
Pollen_data = axs[1].scatter(x,y,c=Diff3,s=40,linewidths=1,edgecolors='k',cmap=plt.get_cmap('coolwarm'),norm=norm)
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
m.colorbar(Pollen_data,label="Temperature anomaly (°C)")

axs[2].set_title("c. MTCO")
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),ax=axs[2],resolution='l')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

CS=axs[2].scatter(x,y,c=Ice[:,:],s=1,cmap="binary")
Diff4=delta_MTCO_MAT_WBiome-MTCO_LGM_Davis
x, y= m(lon_data,lat_data)
Pollen_data = axs[2].scatter(x,y,c=Diff4,s=40,linewidths=1,edgecolors='k',cmap=plt.get_cmap('coolwarm'),norm=norm)
#data_biome["Best"]
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0],linewidth=0.01)
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1],linewidth=0.01)  
m.colorbar(Pollen_data,label="Temperature anomaly (°C)")

plt.tight_layout()
plt.savefig("DavisComparison_Map.pdf")
plt.show()

#%%III-Scatter comparison
#Biome
from matplotlib.ticker import  MultipleLocator, FuncFormatter

MAT_TANN=data_Pollen["TANN_MAT_Wbiome"]-data_Pollen["TANN_Modern"]
WAPLS_TANN=data_Pollen["TANN_WAPLS_Wbiome"]-data_Pollen["TANN_Modern"]
CREST_TANN=data_Pollen["TANN_CREST_Wbiome"]-data_Pollen["TANN_Modern"]

MAT_PANN=(data_Pollen["PANN_MAT_Wbiome"]/data_Pollen["PANN_Modern"])-1
WAPLS_PANN=(data_Pollen["PANN_WAPLS_Wbiome"]/data_Pollen["PANN_Modern"])-1
CREST_PANN=(data_Pollen["PANN_CREST_Wbiome"]/data_Pollen["PANN_Modern"])-1

MAT_MTCO=data_Pollen["MTCO_MAT_Wbiome"]-data_Pollen["MTCO_Modern"]
WAPLS_MTCO=data_Pollen["MTCO_WAPLS_Wbiome"]-data_Pollen["MTCO_Modern"]
CREST_MTCO=data_Pollen["MTCO_CREST_Wbiome"]-data_Pollen["MTCO_Modern"]

MAT_MTWA=(data_Pollen["MTWA_MAT_Wbiome"]-data_Pollen["MTWA_Modern"])
WAPLS_MTWA=(data_Pollen["MTWA_WAPLS_Wbiome"]-data_Pollen["MTWA_Modern"])
CREST_MTWA=(data_Pollen["MTWA_CREST_Wbiome"]-data_Pollen["MTWA_Modern"])

fig, axs = plt.subplots(4, 3, figsize=(12, 14),gridspec_kw={'width_ratios': [3,0.5,0.5]},sharey=False)

#TANN
slope, intercept, r_value, p_value, std_err = linregress(MAT_TANN,CREST_TANN)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[0,0].scatter(MAT_TANN, CREST_TANN, color='blue',label="MAT")
axs[0,0].set_title('a. TANN CREST deviation',fontweight='bold')
axs[0,0].plot(MAT_TANN, slope * MAT_TANN + intercept, color='blue', linewidth=2)
axs[0,0].plot(np.linspace(np.min(MAT_TANN),np.max(MAT_TANN),20),np.linspace(np.min(CREST_TANN),np.max(CREST_TANN),20),color='black',label="y=x",linestyle="--",linewidth=2)

slope, intercept, r_value, p_value, std_err = linregress(WAPLS_TANN,CREST_TANN)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[0,0].scatter(WAPLS_TANN, CREST_TANN, color='green',label='WAPLS')
axs[0,0].plot(WAPLS_TANN,slope * WAPLS_TANN + intercept, color='green', linewidth=2)

axs[0,0].set_xlabel('TANN LGM anomaly (MAT or WAPLS) (°C)')
axs[0,0].set_ylabel('CREST TANN LGM anomaly (°C)')

#PANN
slope, intercept, r_value, p_value, std_err = linregress(MAT_PANN*100,CREST_PANN*100)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[1,0].scatter(MAT_PANN*100, CREST_PANN*100, color='blue')
axs[1,0].set_title('b. PANN CREST deviation',fontweight='bold')
axs[1,0].plot(MAT_PANN*100, slope * MAT_PANN*100 + intercept, color='blue', label='Biomes', linewidth=2)
axs[1,0].plot(np.linspace(np.min(MAT_PANN*100),np.max(MAT_PANN*100),20),np.linspace(np.min(CREST_PANN*100),np.max(CREST_PANN*100),20),color='black',linestyle="--",linewidth=2)

slope, intercept, r_value, p_value, std_err = linregress(WAPLS_PANN*100,CREST_PANN*100)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[1,0].scatter(WAPLS_PANN*100, CREST_PANN*100, color='green')
axs[1,0].plot(WAPLS_PANN*100,slope * WAPLS_PANN*100 + intercept, color='green', label='Biomes', linewidth=2)

axs[1,0].set_xlabel('PANN LGM anomaly (MAT or WAPLS) (%)')
axs[1,0].set_ylabel('CREST PANN LGM anomaly (%)')

axs[0,1].set_title('Mean',fontweight='bold')
axs[0,1].scatter(1,np.mean(MAT_TANN),s=150,color='blue',label='MAT-Biome')
axs[0,1].scatter(1,np.mean(WAPLS_TANN),s=150,color='green',label='WAPLS-Biome')
axs[0,1].scatter(1,np.mean(CREST_TANN),s=150,color='red',label='CREST-Biome')
axs[0,1].set_xticklabels([])
axs[0,1].set_xticks([])

axs[0,2].set_title('STD',fontweight='bold')
axs[0,2].scatter(1,np.std(MAT_TANN),s=150,color='blue',label='MAT')
axs[0,2].scatter(1,np.std(WAPLS_TANN),s=150,color='green',label='WAPLS')
axs[0,2].scatter(1,np.std(CREST_TANN),s=150,color='red',label='CREST')
def round_ticks(x, pos):
    return f'{x:.2f}' 
axs[0,2].yaxis.set_major_formatter(FuncFormatter(round_ticks))
axs[0,2].yaxis.set_major_locator(MultipleLocator(0.02))

axs[0,2].set_xticklabels([])
axs[0,2].set_xticks([])

axs[1,1].scatter(1,np.mean(MAT_PANN)*100,s=150,color='blue')
axs[1,1].scatter(1,np.mean(WAPLS_PANN)*100,s=150,color='green')
axs[1,1].scatter(1,np.mean(CREST_PANN)*100,s=150,color='red')
axs[1,1].set_xticklabels([])
axs[1,1].set_xticks([])

axs[1,2].scatter(1,np.std(MAT_PANN)*100,s=150,color='blue')
axs[1,2].scatter(1,np.std(WAPLS_PANN)*100,s=150,color='green')
axs[1,2].scatter(1,np.std(CREST_PANN)*100,s=150,color='red')
def round_ticks(x, pos):
    return f'{x:.1f}'  
axs[1,2].yaxis.set_major_formatter(FuncFormatter(round_ticks))
axs[1,2].yaxis.set_major_locator(MultipleLocator(0.5))

axs[1,2].invert_yaxis()
axs[1,2].set_xticklabels([])
axs[1,2].set_xticks([])

#MTWA
slope, intercept, r_value, p_value, std_err = linregress(MAT_MTWA,CREST_MTWA)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[2,0].scatter(MAT_MTWA, CREST_MTWA, color='blue')
axs[2,0].set_title('c. MTWA CREST deviation',fontweight='bold')
axs[2,0].plot(MAT_MTWA, slope * MAT_MTWA + intercept, color='blue', label='Biomes', linewidth=2)
axs[2,0].plot(np.linspace(np.min(MAT_MTWA),np.max(MAT_MTWA),20),np.linspace(np.min(CREST_MTWA),np.max(CREST_MTWA),20),color='black',linestyle="--",linewidth=2)

slope, intercept, r_value, p_value, std_err = linregress(WAPLS_MTWA,CREST_MTWA)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[2,0].scatter(WAPLS_MTWA, CREST_MTWA, color='green')
axs[2,0].plot(WAPLS_MTWA,slope * WAPLS_MTWA + intercept, color='green', label='Biome', linewidth=2)
axs[2,0].set_xlabel('MTWA LGM anomaly (MAT or WAPLS) (°C)')
axs[2,0].set_ylabel('CREST MTWA LGM anomaly (°C)')

#MTCO
slope, intercept, r_value, p_value, std_err = linregress(MAT_MTCO,CREST_MTCO)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[3,0].scatter(MAT_MTCO, CREST_MTCO, color='blue')
axs[3,0].set_title('d. MTCO CREST deviation',fontweight='bold')
axs[3,0].plot(MAT_MTCO, slope * MAT_MTCO + intercept, color='blue', label='Biomes', linewidth=2)
axs[3,0].plot(np.linspace(np.min(MAT_MTCO),np.max(MAT_MTCO),20),np.linspace(np.min(CREST_MTCO),np.max(CREST_MTCO),20),color='black',linestyle="--",linewidth=2)

slope, intercept, r_value, p_value, std_err = linregress(WAPLS_MTCO,CREST_MTCO)
print(f'Pente (a) : {slope}')
print(f'Ordonnée à l\'origine (b) : {intercept}')
print(f'Coefficient de corrélation (r) : {r_value}')
print(f'P-value : {p_value}')
print(f'Erreur standard de la pente : {std_err}')

axs[3,0].scatter(WAPLS_MTCO, CREST_MTCO, color='green')
axs[3,0].plot(WAPLS_MTCO,slope * WAPLS_MTCO + intercept, color='green', label='Biomes', linewidth=2)

axs[3,0].set_xlabel('MTCO LGM anomaly (MAT or WAPLS) (°C)')
axs[3,0].set_ylabel('CREST MTCO LGM anomaly (°C)')

axs[2,1].scatter(1,np.mean(MAT_MTWA),s=150,color='blue',label='MAT-Biome')
axs[2,1].scatter(1,np.mean(WAPLS_MTWA),s=150,color='green',label='WAPLS-Biome')
axs[2,1].scatter(1,np.mean(CREST_MTWA),s=150,color='red',label='CREST-Biome')
axs[2,1].set_xticklabels([])
axs[2,1].set_xticks([])

axs[2,2].scatter(1,np.std(MAT_MTWA),s=150,color='blue',label='MAT')
axs[2,2].scatter(1,np.std(WAPLS_MTWA),s=150,color='green',label='WAPLS')
axs[2,2].scatter(1,np.std(CREST_MTWA),s=150,color='red',label='CREST')
def round_ticks(x, pos):
    return f'{x:.1f}'  
axs[2,2].yaxis.set_major_formatter(FuncFormatter(round_ticks))
axs[2,2].yaxis.set_major_locator(MultipleLocator(0.2))

axs[2,2].invert_yaxis()
axs[2,2].set_xticklabels([])
axs[2,2].set_xticks([])

axs[3,1].scatter(1,np.mean(MAT_MTCO),s=150,color='blue')
axs[3,1].scatter(1,np.mean(WAPLS_MTCO),s=150,color='green')
axs[3,1].scatter(1,np.mean(CREST_MTCO),s=150,color='red')
axs[3,1].set_xticklabels([])
axs[3,1].set_xticks([])

axs[3,2].scatter(1,np.std(MAT_MTCO),s=150,color='blue')
axs[3,2].scatter(1,np.std(WAPLS_MTCO),s=150,color='green')
axs[3,2].scatter(1,np.std(CREST_MTCO),s=150,color='red')
def round_ticks(x, pos):
    return f'{x:.1f}'  
axs[3,2].yaxis.set_major_formatter(FuncFormatter(round_ticks))
axs[3,2].yaxis.set_major_locator(MultipleLocator(0.1))

axs[3,2].invert_yaxis()
axs[3,2].set_xticklabels([])
axs[3,2].set_xticks([])

plt.savefig("ScatterCorrelated_Synthesis.pdf")
plt.tight_layout()
plt.show()

#%%
plt.figure(figsize=(12, 8))
m=Basemap(projection="mill",llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max())
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='lightgrey',lake_color='lightblue')

lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
CS=axs[1].scatter(x,y,c=Ice[:,:],s=10,cmap="binary")

x, y= m(nylon,nylat)
Pollen_data=axs[1].scatter(x,y,c=delta_PANN*100,s=300,linewidths=1,edgecolors='k',vmax=100,vmin=-100,cmap=plt.get_cmap('bwr_r', 10))
m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1])
axs[1].colorbar()

plt.figure(figsize=(12, 8))
m = Basemap(projection='mill',lat_ts=10.,llcrnrlon=lon.min(),urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='l')
m.drawcoastlines()
m.drawmapboundary(fill_color='white')
lon2, lat2 = np.meshgrid(lon,lat)
x, y = m(lon2, lat2)
CS=axs[3].contour(x,y,Ice[:,:],levels=[10,50,90],colors=["blue","blue","blue"])

x, y= m(nylon,nylat)
Pollen_data=axs[3].scatter(x,y,c=delta_SI,s=300,linewidths=1,edgecolors='k',vmax=10,vmin=-10,cmap=plt.get_cmap('Spectral', 10))

m.drawparallels(np.arange(int(lat.min()),int(lat.max()),10),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(lon.min()),int(lon.max()),20),labels=[0,0,0,1])
axs[3].colorbar()

