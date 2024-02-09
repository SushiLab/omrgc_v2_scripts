# Create a data folder and move there
mkdir data
cd data

# Download released data
wget https://www.ebi.ac.uk/biostudies/files/S-BSST297/_README
wget https://www.ebi.ac.uk/biostudies/files/S-BSST297/OM-RGC_v2_functional_profile_eggNOG.tar.gz
wget https://www.ebi.ac.uk/biostudies/files/S-BSST297/OM-RGC_v2_functional_profile_KEGG.tar.gz
wget https://www.ebi.ac.uk/biostudies/files/S-BSST297/OM-RGC_v2_taxonomic_profiles.tar.gz

# Uncompress data and remove compressed files
tar xvzf OM-RGC_v2_functional_profile_eggNOG.tar.gz
tar xvzf OM-RGC_v2_functional_profile_KEGG.tar.gz
tar xvzf OM-RGC_v2_taxonomic_profiles.tar.gz

rm OM-RGC_v2_functional_profile_eggNOG.tar.gz
rm OM-RGC_v2_functional_profile_KEGG.tar.gz
rm OM-RGC_v2_taxonomic_profiles.tar.gz

# Download Companion Tables
cd ..
wget https://zenodo.org/record/3539258/files/Salazar_et_al_2019_Suppl_Info.xlsx

# Create the folder for the auxiliary data
mkdir OM-RGC_v2_auxiliary_data

# Open the excel file and save the three needed tables as tab-delimited
Rscript lib/Download_auxiliary_data.R
