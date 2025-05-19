#Biomisation-CREST/Traitement
##R_Package
library(crestr)
library(rio)

#Biomisation process - Load input data 
biome=c("CLDE","TAIG","PION","CLMX","COCO","TEDE","COMX","WAMX","XERO","TUND","COST","WAST","HODE","CODE")
biomeSelect=c("./LDB/CLDE_LDB.csv","./LDB/TAIG_LDB.csv","./LDB/PION_LDB.csv","./LDB/CLMX_LDB.csv","./LDB/COCO_LDB.csv","./LDB/TEDE_LDB.csv","./LDB/COMX_LDB.csv","./LDB/WAMX_LDB.csv","./LDB/XERO_LDB.csv","./LDB/TUND_LDB.csv","./LDB/COST_LDB.csv","./LDB/WAST_LDB.csv","./LDB/HODE_LDB.csv","./LDB/CODE_LDB.csv")
biomeSelect=c("./LGP/CLDE_LGP.csv","./LGP/TAIG_LGP.csv","./LGP/PION_LGP.csv","./LGP/CLMX_LGP.csv","./LGP/COCO_LGP.csv","./LGP/TEDE_LGP.csv","./LGP/COMX_LGP.csv","./LGP/WAMX_LGP.csv","./LGP/XERO_LGP.csv","./LGP/TUND_LGP.csv","./LGP/COST_LGP.csv","./LGP/WAST_LGP.csv","./LGP/HODE_LGP.csv","./LGP/CODE_LGP.csv")
biomeSelect=c("./LGP2/CLDE_LGP2.csv","./LGP2/TAIG_LGP2.csv","./LGP2/PION_LGP2.csv","./LGP2/CLMX_LGP2.csv","./LGP2/COCO_LGP2.csv","./LGP2/TEDE_LGP2.csv","./LGP2/COMX_LGP2.csv","./LGP2/WAMX_LGP2.csv","./LGP2/XERO_LGP2.csv","./LGP2/TUND_LGP2.csv","./LGP2/COST_LGP2.csv","./LGP2/WAST_LGP2.csv","./LGP2/HODE_LGP2.csv","./LGP2/CODE_LGP.csv")
biomeSelect=c("./Furamoos/CLDE_Furamoos.csv","./Furamoos/TAIG_Furamoos.csv","./Furamoos/PION_Furamoos.csv","./Furamoos/CLMX_Furamoos.csv","./Furamoos/COCO_Furamoos.csv","./Furamoos/TEDE_Furamoos.csv","./Furamoos/COMX_Furamoos.csv","./Furamoos/WAMX_Furamoos.csv","./Furamoos/XERO_Furamoos.csv","./Furamoos/TUND_Furamoos.csv","./Furamoos/COST_Furamoos.csv","./Furamoos/WAST_Furamoos.csv","./Furamoos/HODE_Furamoos.csv","./Furamoos/CODE_Furamoos.csv")
biomeSelect=c("./Eifel/CLDE_Eifel.csv","./Eifel/TAIG_Eifel.csv","./Eifel/PION_Eifel.csv","./Eifel/CLMX_Eifel.csv","./Eifel/COCO_Eifel.csv","./Eifel/TEDE_Eifel.csv","./Eifel/COMX_Eifel.csv","./Eifel/WAMX_Eifel.csv","./Eifel/XERO_Eifel.csv","./Eifel/TUND_Eifel.csv","./Eifel/COST_Eifel.csv","./Eifel/WAST_Eifel.csv","./Eifel/HODE_Eifel.csv","./Eifel/CODE_Eifel.csv")
biomeSelect=c("./Cooling/CLDE_Cooling.csv","./Cooling/TAIG_Cooling.csv","./Cooling/PION_Cooling.csv","./Cooling/CLMX_Cooling.csv","./Cooling/COCO_Cooling.csv","./Cooling/TEDE_Cooling.csv","./Cooling/COMX_Cooling.csv","./Cooling/WAMX_Cooling.csv","./Cooling/XERO_Cooling.csv","./Cooling/TUND_Cooling.csv","./Cooling/COST_Cooling.csv","./Cooling/WAST_Cooling.csv","./Cooling/HODE_Cooling.csv","./Cooling/CODE_Cooling.csv")

#Megabiomisation process - Load input data 
biome=c("TEFO","WTFO","BOFO","TUND","STEP","DESE")
biomeSelect=c("./LDB_M/TEFO_LDB.csv","./LDB_M/WTFO_LDB.csv","./LDB_M/BOFO_LDB.csv","./LDB_M/TUND_LDB.csv","./LDB_M/STEP_LDB.csv","./LDB_M/DESE_LDB.csv")
biomeSelect=c("./LGP_M/TEFO_LGP.csv","./LGP_M/WTFO_LGP.csv","./LGP_M/BOFO_LGP.csv","./LGP_M/TUND_LGP.csv","./LGP_M/STEP_LGP.csv","./LGP_M/DESE_LGP.csv")
biomeSelect=c("./LGP_2M/TEFO_LGP_2.csv","./LGP_2M/WTFO_LGP_2.csv","./LGP_2M/BOFO_LGP_2.csv","./LGP_2M/TUND_LGP_2.csv","./LGP_2M/STEP_LGP_2.csv","./LGP_2M/DESE_LGP_2.csv")
biomeSelect=c("./Eifel_2M/TEFO_Eifel.csv","./Eifel_2M/WTFO_Eifel.csv","./Eifel_2M/BOFO_Eifel.csv","./Eifel_2M/TUND_Eifel.csv","./Eifel_2M/STEP_Eifel.csv","./Eifel_2M/DESE_Eifel.csv")
biomeSelect=c("./Furamoos_M/TEFO_Furamoos.csv","./Furamoos_M/WTFO_Furamoos.csv","./Furamoos_M/BOFO_Furamoos.csv","./Furamoos_M/TUND_Furamoos.csv","./Furamoos_M/STEP_Furamoos.csv","./Furamoos_M/DESE_Furamoos.csv")
biomeSelect=c("./Cooling_M/TEFO_Cooling.csv","./Cooling_M/WTFO_Cooling.csv","./Cooling_M/BOFO_Cooling.csv","./Cooling_M/TUND_Cooling.csv","./Cooling_M/STEP_Cooling.csv","./Cooling_M/DESE_Cooling.csv")

#Ages Cal BP
Age_LGP2=read.csv2("ModeleAge_LGP2.txt",sep=',',dec='.')
Age_LGP=read.csv2("ModeleAge_LGP.txt",sep=',',dec='.')
Age_LDB=read.csv2("ModeleAge_LDB.txt",sep=',',dec='.')
Age_Furamoos=read.csv2("ModeleAge_Furamoos.csv",sep=',',dec='.')
Age_Eifel=read.csv2("ModeleAge_Eifel_2.txt",sep=',',dec='.')

#File names
nameFile_AssF="./CoolingPollen.csv"
nameFile_Pollen="./Name_CoolingPollen_Harmonized.txt"
AgeSelect=seq(1,nrow(AssF))#List

for (b in 8:length(biome))#For each biomeP
{ 
  metadata=read.csv2("ModernClimate_Synthesis_BiomeCooling.csv",sep=';',dec=',')#To change into megabiomes
  #metadata=read.csv2("ModernClimate_Synthesis_Megabiome.csv",sep=';',dec='.')#To change into megabiomes
  df=read.csv2('EMPD2-Set-Pollen_Normalized.txt',sep='\t',dec='.')
  NameAssCal=read.csv2('NamesAssCal.txt',sep='\t',dec='.')
  colnames(df)=NameAssCal$Names
  
  NotSampled=metadata$Place[is.na(metadata$TANN)]
  metadata=subset(metadata,!(metadata$Place %in% NotSampled))
  df=subset(df,(df$Place %in% metadata$Place))
  Place=df$Place

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

  TaxaSelect=colnames(AssF)[colSums(AssF)==0]
  AssF=AssF[,!colnames(AssF) %in% c(TaxaSelect)]
  NameF=TaxaSelect
  
  for (i in names(AssF))
  {df[,i]=rowSums(df[grep(i,colnames(df),value=TRUE)])}
  df=df[,which(colnames(df) %in% c(names(AssF)))]
  df=df[,!colnames(df) %in% grep("Delete",colnames(df),value=TRUE)]

  AssF$Ages=AgeSelect
  df$Place=Place

  full_df <- merge(metadata, df, on=c('place'), all=TRUE)
  metadata_id    <- c(1:2,3:4)
  coordinates_id <- 2:3
  climate_id     <- c(6:9)
  climate_names  <- colnames(full_df)[climate_id]
  taxa_names     <- colnames(full_df)[10:length(full_df)]
  
  #Geographic constraints
  #full_df <- full_df[full_df[, 'Deg_Lat'] < 40, ]
  #Biome selection !!!
  
  print("Attention BIOME")
  full_df <- subset(full_df,(full_df$Biome == biome[b]))
  print(biome[b])
  
  distributions <- full_df[FALSE, c(1,2,coordinates_id,climate_id,2)]
  for(tax in taxa_names) {
    w <- which(full_df[, tax] > 0)
    if (tax  %in% names(AssF)&((length(w)!=0)))
    {{distributions <- rbind(distributions,
                             cbind( 'taxonid' = rep(tax, length(w)),
                                    'ProxyName' = rep(tax, length(w)),
                                    full_df[w, coordinates_id],
                                    full_df[w, climate_id],
                                    'weight'=full_df[w, tax]
                             )
    )
    }}}    
  colnames(distributions)[1:4]=c('taxonid', 'ProxyName','longitude', 'latitude')
  head(distributions)
  
  climate_space <- full_df[, c(coordinates_id, climate_id)]
  colnames(climate_space)[1:2]=c('longitude', 'latitude')
  head(climate_space)
  
  AssF=AssF[,c('Ages',unique(distributions$ProxyName))]
  AssF=na.omit(AssF)

#Reconstitution - climateWithObs
rcnstrctn <- crest.set_modern_data(distributions,weight=TRUE,minGridCells=20,df=AssF,climate_space=climate_space,climate=climate_names)

#For rcnstrct
rcnstrctn <- crest.calibrate(rcnstrctn,geoWeighting=TRUE,climateSpaceWeighting=TRUE,shape=c('normal','lognormal','normal','normal'),bin_width =c(2,10,2,2),npoints=2000,verbose=TRUE)

#Taxa select before rcnstrct
TaxaSelect=data.frame(matrix(NA,ncol=1,nrow=length(rcnstrctn$inputs$taxa.name)))
for (pol in 1:length(rcnstrctn$inputs$taxa.name))
{if((anyNA(list(rcnstrctn$modelling$pdfs[pol]),recursive=TRUE)==TRUE))
{TaxaSelect[pol,1]=rcnstrctn$inputs$taxa.name[pol]}
  if (nrow(unique(data.frame(rcnstrctn$modelling$distributions[pol])))<5)
  {print(rcnstrctn$inputs$taxa.name[pol])
    TaxaSelect[pol,1]=rcnstrctn$inputs$taxa.name[pol]}}
TaxaSelect=na.omit(TaxaSelect)
rcnstrctn <- excludeTaxa(rcnstrctn,taxa=TaxaSelect[,1], climate=rcnstrctn$parameters$climate)

rcnstrctn <- crest.reconstruct(rcnstrctn,presenceThreshold = 0,taxWeight = "normalisation",uncertainties = c(0.5,0.95),verbose=TRUE)
data=data.frame(rcnstrctn$reconstructions$TANN$optima,rcnstrctn$reconstructions$TANN$uncertainties,rcnstrctn$reconstructions$PANN$optima,rcnstrctn$reconstructions$PANN$uncertainties,rcnstrctn$reconstructions$MTCO$optima,rcnstrctn$reconstructions$MTCO$uncertainties,rcnstrctn$reconstructions$MTWA$optima,rcnstrctn$reconstructions$MTWA$uncertainties)
write.csv2(data,biomeSelect[b])}



#Supplements
#Cross-validation - Pollen cooling
rcn=loo(rcnstrctn)
lapply(rcn$reconstructions$TANN$loo, head)
plot_loo(rcn,taxanames=c("Thalictrum","Helianthemum","Brassicaceae","Brassicaceae",'Larix',"Salix","Sanguisorba","Picea","Betula",'Quercus',"Fraxinus","Cyperaceae",'Pinus',"Asteraceae","Corylus"))
##########

