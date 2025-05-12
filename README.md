# Local_Methodological_Calibration_MultimethodsApproach
❗ This is a read-only mirror of the CRAN R package repository (DOI/URL) in R.
(i) rioja package (Analysis of Quaternary Science data, version 1.0-6, R Core Team https://cran.r-project.org/web/packages/rioja/index.html, Juggins, 2020) to process MAT and WAPLS reconstructions. 
(ii) crest package: https://github.com/mchevalier2/crestr (https://doi.org/10.1016/j.earscirev.2020.103384, Chevalier et al., 2020)
(iii) The harmonising plant funtional type distributions for evaluating Earth System Models: https://pure.mpg.de/pubman/faces/ViewItemFullPage.jsp?itemId=item_2574775_13 (May 16, 2024; Dallmeyer
et al., 2019)
(iv) The biomisation algorithm simulating PFT cover fractions: https://doi.org/10.5281/zenodo.7523423 (last access: May 10, 2023; Cao and Tian, 2021) 
(v) LegacyPollen2.0 datasets: https://doi.pangaea.de/10.1594/PANGAEA.965907 (PANGEA, https://doi.org/10.1594/PANGAEA.965907; Li et al.,
2025) 
(vi) The Eurasia Modern Pollen database EMPD: https://empd2.github.io and www.europeanpollendatabase.net/ 


Here,list of material
"SyntheseFiles" – (i) Compilation of pollen data sequences used in this study (Dataset_PollenCooling.xlsx) and (ii) synthesis of multi-method climate reconstructions (Test_CoolingPollen.xlsx)
"ModernCalibrationDatasets" – Modern calibration dataset including climate variables and pollen abundances (CalibrationDataset_EMPD2.xlsx)
"DatasetFass_Int" – Initial pollen spectra from La Grande Pile (de Beaulieu & Reille, 1989; de Beaulieu, 1992; de Beaulieu & Reille, 1992) and Lac du Bouchet (de Beaulieu et al., 1990) (GPILXX_age_model_pollen_count.tab and Lac_du_Bouchet_age_model_pollen_count.tab, respectively)
"Modele_Age" – (i) Predicted ages for the sequences from La Grande Pile and Lac du Bouchet (ModeleAge_LGP.txt; ModeleAge_LGP2.txt and ModeleAge_LDB.txt)
"HamonizationTable" – (i) Harmonization table (harmonization_table.txt), (ii) comparison of taxon names before and after harmonization (HarmonizationTable_Results.xlsx)
"Tools_BiomeScheme" – (i) Taxon-to-PFT biomisation matrix (Taxpftpi_Class1.csv), and (ii) PFT-to-Biome matrix (biopftpi_Class1.csv)
"Tools_MegabiomeScheme" – (i) Pollen abundance weighting for PFT score calculation (pollen_variable.csv), and (ii) matrix linking Taxa-PFTs and PFTs-Megabiomes as described by Li et al., under revision (Taxon_PFT_Europe.csv and PFT_Megabiomes_Europe.csv)
"ReconstructionMethods" – R scripts for pollen-based climate reconstructions using MAT, with biomisation process (MAT-Biomisation.R), WAPLS (WAPLS-Biomisation.R), and CREST (CREST-Biomisation.R)
"TANNGlobaleCalibration_Results" – TANN results of global calibration reconstructions following the three pollen-based reconstruction methods (GlobaleCalibration_MAT.csv, GlobaleCalibration_WAPLS.csv, and GlobaleCalibration_CREST.csv)
