#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import io
import requests
import json
import gzip
import pandas as pd

from Bio import SeqIO

with open('autoimmune_diseases.json' , 'r') as f:
    diseases = json.load(f)

# TODO:
# - Add cancer API pull


#######################################################
################### AUTOIMMUNE DATA ###################
#######################################################


def pull_autoimmune_data(table):
    '''
    Extracts T cell and B cell positive assay data for autoimmunity from the IEDB.
    Parameter table = 'tcell' or 'bcell' to specify. 
    '''
    
    df = pd.DataFrame()
    for doid in diseases.keys(): 

        # first get the total number of assays as first request to loop through API
        url = 'https://query-api.iedb.org/%s_search' % table
        params = {'order': 'structure_id',
                  'qualitative_measure': 'neq.Negative', # select positive assays only
                  'disease_iris': f'cs.{{{"DOID:"+doid}}}'} 
        r = requests.get(url, params=params, headers={'Prefer': 'count=exact'})
        rows = int(r.headers['Content-Range'].split('/')[-1])
      
      # loop through IEDB API pages using requests - read into pandas DataFrame and concat
        for i in range(rows // 10000 + 1): # API limit is 10,000 entries
            params['offset'] = i*10000

            # request API call returning csv formatting using parameters in params
            s = requests.get(url, params=params, headers={'accept': 'text/csv', 'Prefer': 'count=exact'})
            try:
                df = pd.concat([df, pd.read_csv(io.StringIO(s.content.decode('utf-8')))])
            except pd.errors.EmptyDataError:
                continue

    return df

# read in autoimmune IEDB T cell and B cell assay tables
print('Reading in autoimmune T cell assay data...')
tcell = pull_autoimmune_data('tcell')
print('Done.')

print('Reading in autoimmune B cell assay data...')
bcell = pull_autoimmune_data('bcell')
print('Done.')

# select epitopes where host and source organism are identical
tcell = tcell[(tcell['host_organism_iri'] == tcell['parent_source_antigen_source_org_iri']) | # OR
              (tcell['host_organism_iri'] == tcell['r_object_source_organism_iri'])]
bcell = bcell[(bcell['host_organism_iri'] == bcell['parent_source_antigen_source_org_iri']) | # OR  
              (bcell['host_organism_iri'] == bcell['r_object_source_organism_iri'])]

# create a column combining epitope and related object source antigen IDs
tcell['source_antigen_iri'] = tcell['parent_source_antigen_iri'].fillna(tcell['r_object_source_molecule_iri'])
bcell['source_antigen_iri'] = bcell['parent_source_antigen_iri'].fillna(bcell['r_object_source_molecule_iri'])

# create a column combining epitope and related object source antigen names
tcell['source_antigen_name'] = tcell['parent_source_antigen_name'].fillna(tcell['r_object_source_molecule_name'])
bcell['source_antigen_name'] = bcell['parent_source_antigen_name'].fillna(bcell['r_object_source_molecule_name'])

# create a column combining epitope and related object source organism names
tcell['source_organism_name'] = tcell['source_organism_name'].fillna(tcell['r_object_source_organism_name'])
bcell['source_organism_name'] = bcell['source_organism_name'].fillna(bcell['r_object_source_organism_name'])

# get reference, unique epitope, and counts for reference, epitope and assay by unique antigen for T cell
tcell_counts = []
for i, row in tcell.groupby('source_antigen_iri'):
    diseases = row['disease_names']
    counts = []
    counts.append(i)
    counts.append(list(row['source_antigen_name'].dropna())[0])
    counts.append(', '.join(set(diseases)))
    counts.append(len(row['reference_id'].unique()))
    counts.append(len(row['structure_id'].unique()))
    counts.append(len(row))
    counts.append(list(row['source_organism_name'].dropna())[0])

    tcell_counts.append(counts)

# get reference, unique epitope, and counts for reference, epitope and assay by unique antigen for B cell
bcell_counts = []
for i, row in bcell.groupby('source_antigen_iri'):
    diseases = row['disease_names']
    counts = []
    counts.append(i)
    counts.append(list(row['source_antigen_name'].dropna())[0])
    counts.append(', '.join(set(diseases)))
    counts.append(len(row['reference_id'].unique()))
    counts.append(len(row['structure_id'].unique()))
    counts.append(len(row))
    counts.append(list(row['source_organism_name'].dropna())[0])

    bcell_counts.append(counts)

# put counts into pandas DataFrame
tcell_counts = pd.DataFrame(tcell_counts, columns=['Protein ID', 'Protein Name', 'Diseases',
                                     'Reference Count', 'Epitope Count', 'Assay Count', 'Parent Species'])
bcell_counts = pd.DataFrame(bcell_counts, columns=['Protein ID','Protein Name', 'Diseases', 
                                     'Reference Count', 'Epitope Count', 'Assay Count', 'Parent Species'])

# limit autoimmune antigens requiring 2 references
tcell_counts = tcell_counts[tcell_counts['Reference Count'] > 1]
bcell_counts = bcell_counts[bcell_counts['Reference Count'] > 1]

print('Done.')

print('Writing autoimmune data...')

# output dataframes to one file
writer = pd.ExcelWriter('autoimmune_data.xlsx', engine='xlsxwriter')

tcell.to_excel(writer, sheet_name='T Cell Assays')
tcell_counts.to_excel(writer, sheet_name='T Cell Counts')
bcell.to_excel(writer, sheet_name='B Cell Assays')
bcell_counts.to_excel(writer, sheet_name='B Cell Counts')

writer.save()

print('Done.')

# use requests to get all autoimmune antigen sequences
print('Getting autoimmune antigen sequences from UniProt...')

with open('autoimmune_antigens.fasta', 'w') as f:
    for idx in tcell_counts['Protein ID'].append(bcell_counts['Protein ID']):
        i = idx.split(':')[1]
        r = requests.get('https://www.uniprot.org/uniprot/%s.fasta' % i)
        if not r.text:
            continue
        else:
            f.write(r.text)

print('Done.')

#######################################################
##################### CANCER DATA #####################
#######################################################

def pull_cancer_data(table):
    '''
    Extracts T cell and B cell positive assay data for cancer from the IEDB.
    Parameter table = 'tcell' or 'bcell' to specify. 
    '''
    df = pd.DataFrame()

    # first get the total number of assays as first request to loop through API
    url = 'https://query-api.iedb.org/%s_search' % table
    params = {'order': 'structure_id',
              'qualitative_measure': 'neq.Negative',
              'or': '(e_related_object_type.eq.neo-epitope, immunization_description.plfts."Occurrence of cancer")'}
    r = requests.get(url, params=params, headers={'Prefer': 'count=exact'})
    rows = int(r.headers['Content-Range'].split('/')[-1])
      
    # loop through IEDB API pages using requests - read into pandas DataFrame and concat
    for i in range(rows // 10000 + 1): # API limit is 10,000 entries
        params['offset'] = i*10000

        # request API call returning csv formatting using parameters in params
        s = requests.get(url, params=params, headers={'accept': 'text/csv', 'Prefer': 'count=exact'})
        try:
            df = pd.concat([df, pd.read_csv(io.StringIO(s.content.decode('utf-8')))])
        except pd.errors.EmptyDataError:
            continue

    return df

# read in cancer IEDB T cell and B cell assay tables
print('Reading in cancer T cell assay data...')
tcell = pull_cancer_data('tcell')
print('Done.')

print('Reading in cancer B cell assay data...')
bcell = pull_cancer_data('bcell')
print('Done.')

print('Extracting cancer data...')

# select cancer epitopes by neo-epitope marking OR by occurrence of cancer label
c_t_cell_epitopes = t_cell_df[(t_cell_df['Related Object Epitope Relationship'] == 'neo-epitope') |
                              (t_cell_df['1st in vivo Process Type'] == 'Occurrence of cancer') |
                              (t_cell_df['2nd in vivo Process Type'] == 'Occurrence of cancer')]

c_b_cell_epitopes = b_cell_df[(b_cell_df['Related Object Epitope Relationship'] == 'neo-epitope') |
                              (b_cell_df['1st in vivo Process Type'] == 'Occurrence of cancer') |
                              (b_cell_df['2nd in vivo Process Type'] == 'Occurrence of cancer')]

# fill related object antigen name with parent protein THEN by antigen name from the epitope table
# this is because the parent protein has more consistent naming
c_t_cell_epitopes['Related Object Parent Protein'].fillna(c_t_cell_epitopes['Parent Protein'], inplace=True)
c_b_cell_epitopes['Related Object Parent Protein'].fillna(c_b_cell_epitopes['Parent Protein'], inplace=True)
c_t_cell_epitopes['Related Object Parent Organism'].fillna(c_t_cell_epitopes['Parent Species'], inplace=True)
c_b_cell_epitopes['Related Object Parent Organism'].fillna(c_b_cell_epitopes['Parent Species'], inplace=True)
c_t_cell_epitopes['Related Object Parent Protein ID'].fillna(c_t_cell_epitopes['Parent Protein ID'], inplace=True)
c_b_cell_epitopes['Related Object Parent Protein ID'].fillna(c_b_cell_epitopes['Parent Protein ID'], inplace=True)
c_t_cell_epitopes['Related Object Parent Organism ID'].fillna(c_t_cell_epitopes['Parent Species ID'], inplace=True)
c_b_cell_epitopes['Related Object Parent Organism ID'].fillna(c_b_cell_epitopes['Parent Species ID'], inplace=True)

# drop columns used above to fill in nulls
c_t_cell_epitopes.drop(columns=['Parent Protein', 'Parent Species', 'Parent Protein ID', 'Parent Species ID'], inplace=True)
c_b_cell_epitopes.drop(columns=['Parent Protein', 'Parent Species', 'Parent Protein ID', 'Parent Species ID'], inplace=True)

# get reference, unique epitope, and assay counts by unique antigen for T cell
c_t_cell_counts = []
for i, j in c_t_cell_epitopes.groupby('Related Object Parent Protein ID'):
    diseases = list(set(filter(None, j[['1st in vivo Disease State', '2nd in vivo Disease State']].to_numpy(na_value=None).flatten())))
    counts = []
    counts.append(i)
    counts.append(', '.join(diseases))
    counts.append(list(j['Related Object Parent Protein'].dropna())[0])
    counts.append(len(j['Reference ID'].unique()))
    counts.append(len(j['Epitope'].unique()))
    counts.append(len(j))
    counts.append(list(j['Related Object Parent Organism'].dropna())[0])

    if 'virus' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Virus')
    elif 'sapiens' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Animal')
    elif 'musculus' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Animal')
    elif 'boydii' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Bacteria')
    else:
      counts.append('Unknown')

    c_t_cell_counts.append(counts)

# get reference, unique epitope, and assay counts by unique antigen for B cell
c_b_cell_counts = []
for i, j in c_b_cell_epitopes.groupby('Related Object Parent Protein ID'):
    diseases = list(set(filter(None, j[['1st in vivo Disease State', '2nd in vivo Disease State']].to_numpy(na_value=None).flatten())))
    counts = []
    counts.append(i)
    counts.append(', '.join(diseases))
    counts.append(list(j['Related Object Parent Protein'].dropna())[0])
    counts.append(len(j['Reference ID'].unique()))
    counts.append(len(j['Epitope'].unique()))
    counts.append(len(j))
    counts.append(list(j['Related Object Parent Organism'].dropna())[0])

    if 'virus' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Virus')
    elif 'sapiens' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Animal')
    elif 'musculus' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Animal')
    elif 'communis' in list(j['Related Object Parent Organism'].dropna())[0]:
      counts.append('Plant')
    else:
      counts.append('Bacteria')

    c_b_cell_counts.append(counts)

# put counts into Pandas dataframe
c_t_cell_counts = pd.DataFrame(c_t_cell_counts, columns=['Parent Protein ID', 'Diseases', 'Antigen', 'References', 'Epitopes', 'Assays', 'Source Organism', 'Parent'])
c_b_cell_counts = pd.DataFrame(c_b_cell_counts, columns=['Parent Protein ID', 'Diseases', 'Antigen', 'References', 'Epitopes', 'Assays', 'Source Organism', 'Parent'])

print('Writing cancer data...')

# output dataframes to one file
writer = pd.ExcelWriter('cancer_data.xlsx', engine='xlsxwriter')

c_t_cell_epitopes.to_excel(writer, sheet_name='T Cell Epitopes')
c_t_cell_counts.to_excel(writer, sheet_name='T Cell Antigen Counts')
c_b_cell_epitopes.to_excel(writer, sheet_name='B Cell Epitopes')
c_b_cell_counts.to_excel(writer, sheet_name='B Cell Antigen Counts')

writer.save()

print('Done.')

# use requests to get all cancer antigen sequences
print('Getting cancer antigen sequences from UniProt...')

with open('cancer_antigens.fasta', 'w') as f:
    for i in c_t_cell_counts['Parent Protein ID'].append(c_b_cell_counts['Parent Protein ID']):
        r = requests.get('https://www.uniprot.org/uniprot/%s.fasta' % i)
        if not r.text:
            continue
        else:
            f.write(r.text[:-1])

print('Done.')

###############################################################
##################### NORMAL PROTEIN DATA #####################
###############################################################

# get human proteome
r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz', stream=True)
with open('9606.fasta', 'wb') as f:
    f.write(gzip.open(r.raw, 'rb').read())

# read in human proteome with Biopython and remove autoimmune and cancer antigens
human_protein_ids = []
for record in list(SeqIO.parse('9606.fasta', 'fasta')):
    human_protein_ids.append(record.id.split('|')[1])
