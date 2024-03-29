# OM-RGC v2

> [Salazar et al., **Gene expression changes and community turnover differentially shape the global ocean metatranscriptome**, *Cell*, 2019](https://doi.org/10.1016/j.cell.2019.10.014)


## Data

- Raw data can be found at the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) with identifiers in [here](https://doi.org/10.5281/zenodo.3473199)
- The Ocean Microbial Reference Gene Catalog v2 (OM-RGC.v2) can be found at [BioStudies](https://www.ebi.ac.uk/biostudies/studies/S-BSST297) and at the [Companion Website](http://ocean-microbiome.org/)

## Scripts

The data used for the production of figures is not synchronized here but can be downloaded with:

```bash
git clone https://github.com/SushiLab/omrgc_v2_scripts
cd omrgc_v2_scripts/analysis/
sh Download_data.sh
```

* **analysis**: Folder containing the scripts and external resources to produce data, tables and figures. 
	- **lib**: Contains some R functions and external resources used by the scripts.
	- **scripts\_for_figures**: Scripts to produce the tables and figures of the paper.
