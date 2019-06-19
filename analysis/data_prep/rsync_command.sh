# From Mac to CDS
rsync -vr /Users/guillemsalazar/polybox/ETH/PROJECTS/Thermogenomics_release_new guillems@earle:/nfs/cds/scratch/guillems/

# From Mac to Euler
rsync -vr /Users/guillemsalazar/polybox/ETH/PROJECTS/Thermogenomics_release_new guillems@euler:/cluster/scratch/guillems/

# From Euler to CDS
rsync -vr /cluster/scratch/guillems/Thermogenomics_release_new /nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/scratch/guillems/

# From Euler to Mac
rsync -vr guillems@euler:/cluster/scratch/guillems/Thermogenomics_release_new /Users/guillemsalazar/polybox/ETH/PROJECTS/

scp guillems@earle:/nfs/cds/scratch/guillems/Thermogenomics_release_new/data/processed/NOG*.gz /Users/guillemsalazar/polybox/ETH/PROJECTS/Thermogenomics_release_new/data/processed
