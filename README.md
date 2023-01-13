# Autoimmune Antigen Features and Classification

**Start here:** [get_data.py](get_data.py)


## Background

Autoimmunity is the process whereby a host's immune system initiates 
an immune response against it's own cells and tissues. The adaptive
immune system reacts with self-proteins via antibodies or T cells 
leading to host damage. This project attempts to explore and find
features/traits of these self-proteins (or autoimmune antigens) which
can then be used to build a model to predict novel antigens.

### Pipeline
1. Run [get_data.py](get_data.py).
    a. Pulls autoimmune data, T cell and B cell assays, from the IEDB API looping through each autoimmune DOID in [autoimmune_disease.json](autoimmune_disease.json).
    b. Concatenates both T cell and B cell assay data into antigens and gets assay, epitope, and reference counts.
    c. Pulls all autoimmune antigen UniProt entries using their API.
    d. Pulls human proteome from UniProt (ID: UP000005640).
    e. Separates autoimmune antigens from the rest of the human proteome into different files (autoimmune_antigens.fasta and non_autoimmune_proteins.fasta).
    f. Combines all data into one file ([combined_data.csv](combined_data.csv)) which has autoimmune indicator and sequence.
2. Run [clustering.sh](clustering.sh).
    a. This requires "cd-hit" folder which can be downloaded [here](https://github.com/weizhongli/cdhit/releases).
    b. Uses CD-HIT to cluster the entire human proteome at 40% sequence identity. This will cluster the human proteome into ~13,000 clusters with range of 1-160 proteins.
    c. The bash script runs [distribute_clusters.py](distribute_clusters.py) which will combine clusters into 5 folds for training and testing. Rewrites ([combined_data.csv](combined_data.csv)).


### GOALS
1. Extract and validate autoimmune antigen data from the [IEDB](https://www.iedb.org/). **DONE**
2. Explore and identify features of autoimmune antigens. **IN PROGRESS**
3. Build a model for autoimmune antigen prediction.
