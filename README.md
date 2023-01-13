# Autoimmune Antigen Features and Classification

## Background

Autoimmunity is the process whereby a host's immune system initiates an immune response against it's own cells and tissues. The adaptive immune system reacts with self-proteins via antibodies or T cells leading to host damage. This project attempts to explore and find features/traits of these autoantigens which can then be used to build a model to predict novel antigens.

The IEDB contains curated immune epitope data from the scientific literature. These include autoimmune epitopes, which are host sequences recognized by the host immune cells, leading to a disease state. We want to first collect these autoimmune epitopes and retrieve their corresponding antigens. This will give us a dataset of known autoimmune antigens that we can then explore for unique features.

## Pipeline
1. Run [get_data.py](get_data.py).
    1. Pulls autoimmune data, T cell and B cell assays, from the IEDB API looping through each autoimmune DOID in [autoimmune_diseases.json](autoimmune_diseases.json).
    2. Concatenates both T cell and B cell assay data into antigens and gets assay, epitope, and reference counts.
    3. Pulls all autoimmune antigen UniProt entries using their API.
    4. Pulls human proteome from UniProt (ID: UP000005640).
    5. Separates autoimmune antigens from the rest of the human proteome into different files (autoimmune_antigens.fasta and non_autoimmune_proteins.fasta).
    6. Combines all data into one file ([combined_data.csv](combined_data.csv)) which has autoimmune indicator and sequence.
2. Run [clustering.sh](clustering.sh).
    1. This requires "cd-hit" folder which can be downloaded [here](https://github.com/weizhongli/cdhit/releases).
    2. Uses CD-HIT to cluster the entire human proteome at 40% sequence identity. This will cluster the human proteome into ~13,000 clusters with range of 1-160 proteins.
    3. The bash script runs [distribute_clusters.py](distribute_clusters.py) which will combine clusters into 5 folds for training and testing. Rewrites ([combined_data.csv](combined_data.csv)).
3. Run [blast_model.py](blast_model.py).
    1. This will generate a model that can be compared to other trained models as a baseline. This script will take each fold and BLAST them against the other folds and take the top hit by sequence identity and use that as a predictor for autoimmune or not autoimmune.
    2. Five files will be generated for each fold (foldX_output.csv) as the BLAST output.
    3. ROC curves are then created for each fold and this can be used to compare to other ROC curves created from other models. 

## Goals
1. Extract and validate autoimmune antigen data from the [IEDB](https://www.iedb.org/). **DONE**
2. Explore and identify features of autoimmune antigens. **IN PROGRESS**
3. Build a model for autoimmune antigen prediction.
