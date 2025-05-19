library(stringr)
#Dominent megabiomes-without weighted
biome=c("CLDE","TAIG","PION","CLMX","COCO","TEDE","COMX","WAMX","XERO","TUND","COST","WAST","CODE","HODE")
#biome=c("TEFO","WTFO","BOFO","TUND","STEP","DESE")

OutBiomization=read.csv2("OutLDB_Biomisation_Classe2.csv",sep=';',dec='.')
OutBiomization=OutBiomization
Ex_Biome=read.csv2("CLDE_LDB.csv",sep=';',dec=',')[,-c(1)]
Output_Recon=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=nrow(Ex_Biome)))
names(Output_Recon)=names(Ex_Biome)
row.names(Output_Recon)=rownames(Ex_Biome)
for (n in 1:nrow(OutBiomization))
{Best=OutBiomization$Best[n]
Nom_File=str_c(Best,"_LDB.csv")
Output_Recon[n,]=read.csv2(Nom_File,sep=';',dec=',')[n,-c(1)]}
write.csv2(Output_Recon,"Output_Recon_CREST_Biome_LDB.csv")

#Weighted mean-FossilSeq
biome=c("CLDE","TAIG","PION","CLMX","COCO","TEDE","COMX","WAMX","XERO","TUND","COST","WAST","CODE","HODE")
#biome=c("TEFO","WTFO","BOFO","TUND","STEP")

OutBiomization=read.csv2("OutLDB_Biomisation_Classe2.csv",sep=';',dec=',')
OutBiomization=OutBiomization[,-c(1,16,17,18)]
Ex_Biome=read.csv2("CLDE_LDB.csv",sep=';',dec=',')[,-c(1)]
Output_Recon=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=nrow(Ex_Biome)))
names(Output_Recon)=names(Ex_Biome)
row.names(Output_Recon)=rownames(Ex_Biome)

for (n in 1:nrow(OutBiomization))
{data=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=length(biome)))
for (m in 1:length(biome))
{Nom_File=str_c(biome[m],"_LDB.csv")
data[m,]=read.csv2(Nom_File,sep=';',dec=',')[n,-c(1)]
data[m,]=data[m,]*OutBiomization[n,biome[m]]}
Output_Recon[n,]=colSums(data)/sum(OutBiomization[n,])}
write.csv2(Output_Recon,"Output_Recon_CREST_Biome_Weighted_LDB.csv")


#######Cooling
#Dominent megabiomes-without weighted
biome=c("CLDE","TAIG","PION","CLMX","COCO","TEDE","COMX","WAMX","XERO","TUND","COST","WAST","HODE")
#biome=c("TEFO","WTFO","BOFO","TUND","STEP")
biome=c("TEFO","WTFO","BOFO","TUND","STEP","DESE")


OutBiomization=read.csv2("OutFossil_Biomization_Cooling2.csv",sep=';',dec=',')
#OutBiomization=read.csv2("OutFossil_MegaBiomization_Cooling.csv",sep=',',dec='.')
#Ex_Biome=read.csv2("TEFO_Cooling.csv",sep=';',dec=',')[,-c(1)]#To change!
Ex_Biome=read.csv2("CLDE_Cooling.csv",sep=';',dec=',')[,-c(1)]#To change!
Output_Recon=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=nrow(Ex_Biome)))
names(Output_Recon)=names(Ex_Biome)
row.names(Output_Recon)=rownames(Ex_Biome)
for (n in 1:nrow(OutBiomization))
{Best=OutBiomization$Best[n]
Nom_File=str_c(Best,"_Cooling.csv")
Output_Recon[n,]=read.csv2(Nom_File,sep=';',dec=',')[n,-c(1)]}
write.csv2(Output_Recon,"Output_Recon_MAT_BestBiome_Cooling.csv")

#Weighted mean-Cooling
OutBiomization=read.csv2("OutFossil_Biomization_Cooling2.csv",sep=';',dec=',')
#OutBiomization=read.csv2("OutFossil_MegaBiomization_Cooling.csv",sep=',',dec='.')
OutBiomization=OutBiomization[,-c(1,18)]
#OutBiomization=OutBiomization[,-c(1,5,9,10)]#To change!
#Ex_Biome=read.csv2("TEFO_Cooling.csv",sep=';',dec=',')[,-c(1)]
Ex_Biome=read.csv2("CLDE_Cooling.csv",sep=';',dec=',')[,-c(1)]
Output_Recon=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=nrow(Ex_Biome)))
names(Output_Recon)=names(Ex_Biome)
row.names(Output_Recon)=rownames(Ex_Biome)

for (n in 1:nrow(OutBiomization))
{data=data.frame(matrix(NA,ncol=ncol(Ex_Biome),nrow=length(biome)))
for (m in 1:length(biome))
{Nom_File=str_c(biome[m],"_Cooling.csv")
data[m,]=read.csv2(Nom_File,sep=';',dec=',')[n,-c(1)]
data[m,]=data[m,]*OutBiomization[n,biome[m]]}
Output_Recon[n,]=colSums(data)/sum(OutBiomization[n,])}
write.csv2(Output_Recon,"Output_Recon_MAT_WeightedBiome_Cooling.csv")

