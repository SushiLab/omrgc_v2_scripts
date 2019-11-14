library(readxl)
library(data.table)

# Open the excel file (sheets 6 to 8) 
auxiliary_metaG<-read_excel(path = "Salazar_et_al_2019_Suppl_Info.xlsx",sheet = 6)
auxiliary_metaT<-read_excel(path = "Salazar_et_al_2019_Suppl_Info.xlsx",sheet = 7)
auxiliary_matched<-read_excel(path = "Salazar_et_al_2019_Suppl_Info.xlsx",sheet = 8)

# Save them as tab-delimited files
fwrite(auxiliary_metaG,file = "OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_metaG.tsv",sep = "\t")
fwrite(auxiliary_metaT,file = "OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_metaT.tsv",sep = "\t")
fwrite(auxiliary_matched,file = "OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_matched.tsv",sep = "\t")