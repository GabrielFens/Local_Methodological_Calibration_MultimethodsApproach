#Set up environment/Traitement
#Biomisation-MAT
#Please, install these packages if you want to process it - rioja and analogue packages.
library(rioja)#1.0-6 version
library(dplyr)

#MAT transfer function reconstruction - SQ chord
#Biomisation process - Load input data
biome=c("CLDE","TAIG","PION","CLMX","COCO","TEDE","COMX","WAMX","XERO","TUND","COST","WAST","CODE","HODE")
biomeSelect=c("./LDB/CLDE_LDB.csv","./LDB/TAIG_LDB.csv","./LDB/PION_LDB.csv","./LDB/CLMX_LDB.csv","./LDB/COCO_LDB.csv","./LDB/TEDE_LDB.csv","./LDB/COMX_LDB.csv","./LDB/WAMX_LDB.csv","./LDB/XERO_LDB.csv","./LDB/TUND_LDB.csv","./LDB/COST_LDB.csv","./LDB/WAST_LDB.csv","./LDB/CODE_LDB.csv","./LDB/HODE_LDB.csv")
biomeSelect=c("./LGP/CLDE_LGP.csv","./LGP/TAIG_LGP.csv","./LGP/PION_LGP.csv","./LGP/CLMX_LGP.csv","./LGP/COCO_LGP.csv","./LGP/TEDE_LGP.csv","./LGP/COMX_LGP.csv","./LGP/WAMX_LGP.csv","./LGP/XERO_LGP.csv","./LGP/TUND_LGP.csv","./LGP/COST_LGP.csv","./LGP/WAST_LGP.csv","./LGP/CODE_LGP.csv","./LGP/HODE_LGP.csv")
biomeSelect=c("./LGP2/CLDE_LGP2.csv","./LGP2/TAIG_LGP2.csv","./LGP2/PION_LGP2.csv","./LGP2/CLMX_LGP2.csv","./LGP2/COCO_LGP2.csv","./LGP2/TEDE_LGP2.csv","./LGP2/COMX_LGP2.csv","./LGP2/WAMX_LGP2.csv","./LGP2/XERO_LGP2.csv","./LGP2/TUND_LGP2.csv","./LGP2/COST_LGP2.csv","./LGP2/WAST_LGP2.csv","./LGP2/CODE_LGP2.csv","./LGP2/HODE_LGP2.csv")
biomeSelect=c("./Furamoos/CLDE_Furamoos.csv","./Furamoos/TAIG_Furamoos.csv","./Furamoos/PION_Furamoos.csv","./Furamoos/CLMX_Furamoos.csv","./Furamoos/COCO_Furamoos.csv","./Furamoos/TEDE_Furamoos.csv","./Furamoos/COMX_Furamoos.csv","./Furamoos/WAMX_Furamoos.csv","./Furamoos/XERO_Furamoos.csv","./Furamoos/TUND_Furamoos.csv","./Furamoos/COST_Furamoos.csv","./Furamoos/WAST_Furamoos.csv","./Furamoos/CODE_Furamoos.csv","./Furamoos/HODE_Furamoos.csv")
biomeSelect=c("./Eifel/CLDE_Eifel.csv","./Eifel/TAIG_Eifel.csv","./Eifel/PION_Eifel.csv","./Eifel/CLMX_Eifel.csv","./Eifel/COCO_Eifel.csv","./Eifel/TEDE_Eifel.csv","./Eifel/COMX_Eifel.csv","./Eifel/WAMX_Eifel.csv","./Eifel/XERO_Eifel.csv","./Eifel/TUND_Eifel.csv","./Eifel/COST_Eifel.csv","./Eifel/WAST_Eifel.csv","./Eifel/CODE_Eifel.csv","./Eifel/HODE_Eifel.csv")
biomeSelect=c("./Cooling/CLDE_Cooling.csv","./Cooling/TAIG_Cooling.csv","./Cooling/PION_Cooling.csv","./Cooling/CLMX_Cooling.csv","./Cooling/COCO_Cooling.csv","./Cooling/TEDE_Cooling.csv","./Cooling/COMX_Cooling.csv","./Cooling/WAMX_Cooling.csv","./Cooling/XERO_Cooling.csv","./Cooling/TUND_Cooling.csv","./Cooling/COST_Cooling.csv","./Cooling/WAST_Cooling.csv","./Cooling/CODE_Cooling.csv","./Cooling/HODE_Cooling.csv")

#Megabiomisation process - Load input data
biome=c("TEFO","WTFO","BOFO","TUND","STEP","DESE")
biomeSelect=c("./LDB_M/TEFO_LDB.csv","./LDB_M/WTFO_LDB.csv","./LDB_M/BOFO_LDB.csv","./LDB_M/TUND_LDB.csv","./LDB_M/STEP_LDB.csv","./LDB_M/DESE_LDB.csv")
biomeSelect=c("./LGP_M/TEFO_LGP.csv","./LGP_M/WTFO_LGP.csv","./LGP_M/BOFO_LGP.csv","./LGP_M/TUND_LGP.csv","./LGP_M/STEP_LGP.csv","./LGP_M/DESE_LGP.csv")
biomeSelect=c("./LGP_2M/TEFO_LGP_2.csv","./LGP_2M/WTFO_LGP_2.csv","./LGP_2M/BOFO_LGP_2.csv","./LGP_2M/TUND_LGP_2.csv","./LGP_2M/STEP_LGP_2.csv","./LGP_2M/DESE_LGP_2.csv")
biomeSelect=c("./Eifel_2M/TEFO_Eifel.csv","./Eifel_2M/WTFO_Eifel.csv","./Eifel_2M/BOFO_Eifel.csv","./Eifel_2M/TUND_Eifel.csv","./Eifel_2M/STEP_Eifel.csv","./Eifel_2M/DESE_Eifel.csv")
biomeSelect=c("./Furamoos_M/TEFO_Furamoos.csv","./Furamoos_M/WTFO_Furamoos.csv","./Furamoos_M/BOFO_Furamoos.csv","./Furamoos_M/TUND_Furamoos.csv","./Furamoos_M/STEP_Furamoos.csv","./Furamoos_M/DESE_Furamoos.csv")
biomeSelect=c("./Cooling_M/TEFO_Cooling.csv","./Cooling_M/WTFO_Cooling.csv","./Cooling_M/BOFO_Cooling.csv","./Cooling_M/TUND_Cooling.csv","./Cooling_M/STEP_Cooling.csv","./Cooling_M/DESE_Cooling.csv")

#Initialisation of climate variables
VariableClim=rep(c(1:4),times=length(unique(biome)))
RMSEP=data.frame(matrix(NA,ncol=4,nrow=length(VariableClim)))
Variable=c("T_ann","P_ann","MTCO","MTWA")

#File names
nameFile_AssF="./CoolingPollen.csv"
nameFile_Pollen="./Name_CoolingPollen_Harmonized.txt"
AgeSelect=seq(1,nrow(AssF))#List
#AssF$Ages=Age for Cooling data
  
#Ages Cal BP
Age_LGP2=read.csv2("ModeleAge_LGP2.txt",sep=',',dec='.')
Age_LGP=read.csv2("ModeleAge_LGP.txt",sep=',',dec='.')
Age_LDB=read.csv2("ModeleAge_LDB.txt",sep=',',dec='.')
Age_Furamoos=read.csv2("ModeleAge_Furamoos.csv",sep=',',dec='.')
Age_Eifel=read.csv2("ModeleAge_Eifel_2.txt",sep=',',dec='.')

m=1
for (b in 1:length(biome))#For each biome
#Preparation before training
{for (h in 1:length(Variable))#For each variable
{#Modern spectrum selection for each biome
  OutBiomization=read.csv2("ModernClimate_Synthesis_BiomeCooling.csv",sep=';',dec=',')#To change into megabiomes
  #OutBiomization=read.csv2("ModernClimate_Synthesis_Megabiome.csv",sep=';',dec=',')#To change into megabiomes
  NotSampled=OutBiomization$SampleName[is.na(OutBiomization$TANN)]
  OutBiomization=subset(OutBiomization,!(OutBiomization$Place %in% NotSampled))
  clim=subset(OutBiomization,OutBiomization$Biome==biome[b])
  Var=clim[,h+5] #Adapt climatic variable
  VarClim=as.numeric(unlist(Var))
  
  #The fossil dataset and EMPD2 assemblage importation
  #Select all taxa from the fossil dataset
  if ((b==1)&(h==1))
  {
    AssCal=read.csv2('EMPD2-Set.txt',sep='\t',dec='.')
    NameAssCal=read.csv2('NamesAssCal.txt',sep='\t',dec='.')
    
    AssF=read.csv2(nameFile_AssF,sep=";",dec=',')
    NameF=read.csv2(nameFile_Pollen,sep=";",dec='.')#To change!
    rownames(AssF)=make.names(AssF$Age,unique=TRUE)#Age
    AssF=AssF[,-c(1)]
    AssF=AssF*100
    colnames(AssF)=NameF$Names
    
    AssF=AssF[,which(colnames(AssF) %in% c(NameF$Names))]
    for (i in unique(NameF$Names))
    {AssF[,i]=rowSums(AssF[grep(i,colnames(AssF),value=TRUE)])}
    AssF=AssF[,which(colnames(AssF) %in% c(unique(NameF$Names)))]
    
    AssCal=subset(AssCal,AssCal$Place %in% OutBiomization$Place)
    rownames(AssCal)=AssCal$Place
    colnames(AssCal)=NameAssCal$Names
    
    AssCal=AssCal[,!colnames(AssCal) %in% grep("Delete",colnames(AssCal),value=TRUE)]
    for (i in names(AssF))
    {AssCal[,i]=rowSums(AssCal[grep(i,colnames(AssCal),value=TRUE)])}
    AssCal=AssCal[,which(colnames(AssCal) %in% c(names(AssF)))]
    
    AssCalF=AssCal
    for (e in 1:nrow(AssCal))
    {AssCal[e,]=100*(AssCal[e,])/sum(AssCalF[e,])}
    
    Age=AgeSelect
    AssCalInt=AssCal
    AssFInt=AssF}
  
  AssCal=subset(AssCalInt,rownames(AssCalInt) %in% OutBiomization$Place[OutBiomization$Biome==biome[b]])
  #Use only taxa with countings (exclude Zeros)
  TaxaSelect=colnames(AssCal)[colSums(AssCal)==0]
  AssCal=AssCal[,!colnames(AssCal) %in% c(TaxaSelect)]
  AssF=AssFInt[,!colnames(AssFInt) %in% c(TaxaSelect)]
  
  taxa=AssCal
  modern_clim=VarClim
  
  #Climate Reconstructions
  modMAT=MAT(taxa,modern_clim,dist.method="sq.chord",lean=FALSE,k=10)#To change for global map calibration
  analogs=which.min(c(performance(modMAT)$object[1:10,1]))
  RMSEP[m,]=performance(modMAT)$object[analogs,][1]
  m=m+1
  
  print(biome[b])
  print(VariableClim[h])
  print(analogs)
  
  modMAT=MAT(taxa,modern_clim,dist.method="sq.chord",lean=FALSE,k=analogs)#Adapt k-analogue!
  predMAT=predict(modMAT,newdata=AssF,sse=TRUE,lean=FALSE,k=analogs,nboot=200)#Adapt k-analogue!
  
  if (VariableClim[h]==1)
  {ReconTANN=data.frame(Age,predMAT$fit[,1],predMAT$SEP.boot[,1])}
  if (VariableClim[h]==2)
  {ReconPANN=data.frame(Age,predMAT$fit[,1],predMAT$SEP.boot[,1])}
  if (VariableClim[h]==3)
  {ReconMTCO=data.frame(Age,predMAT$fit[,1],predMAT$SEP.boot[,1])}
  if (VariableClim[h]==4)
  {ReconMTWA=data.frame(Age,predMAT$fit[,1],predMAT$SEP.boot[,1])
  climateRecon_MAT=data.frame(ReconTANN,ReconPANN,ReconMTCO,ReconMTWA)
  write.csv2(climateRecon_MAT,biomeSelect[b])}
}}

#END