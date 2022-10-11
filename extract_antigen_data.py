#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import io
import requests
import json
import gzip
import numpy as np
import pandas as pd

from pepmatch import Preprocessor, Matcher
from collections import Counter
from Bio import SeqIO

from columns import t_cell_data_columns, t_cell_columns, b_cell_data_columns, b_cell_columns

with open('autoimmune_diseases.json' , 'r') as f:
    diseases = json.load(f)

# TODO:
# - Cobmine related object antigens for autoimmune data

def pull_iedb_assay_data(table):
    '''
    Extracts T cell and B cell positive assay data from the IEDB.
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
        pages = int(r.headers['Content-Range'].split('/')[-1])
      
      # loop through IEDB API pages using requests - read into pandas DataFrame and concat
        for i in range(pages // 10000 + 1): # API limit is 10,000 entries
            params['offset'] = i*10000

            # request API call returning csv formatting using parameters in params
            s = requests.get(url, params=params, headers={'accept': 'text/csv', 'Prefer': 'count=exact'})
            try:
                df = pd.concat([df, pd.read_csv(io.StringIO(s.content.decode('utf-8')))])
            except pd.errors.EmptyDataError:
                continue

    return df

# read in positive IEDB T cell and B cell assay tables
print('Reading in positive T cell assay data...')
pos_t_cell = pull_iedb_assay_data('tcell')
print('Done.')

print('Reading in positive B cell assay data...')
pos_b_cell = pull_iedb_assay_data('bcell')
print('Done.')

# select only colums that are needed
t_cell_df = pos_t_cell[t_cell_data_columns]
b_cell_df = pos_b_cell[b_cell_data_columns]

# collapse headers by renaming column names
t_cell_df.columns = t_cell_columns
b_cell_df.columns = b_cell_columns

# -------------------- Begin Data Cleaning and Organizing --------------------
print('Cleaning and organizing data...')
# get assay ids
t_cell_df['Assay ID'] = t_cell_df['Assay ID'].str.split('/').str[-1]
b_cell_df['Assay ID'] = b_cell_df['Assay ID'].str.split('/').str[-1]

# get reference ids
t_cell_df['Reference ID'] = t_cell_df['Reference ID'].str.split('/').str[-1]
b_cell_df['Reference ID'] = b_cell_df['Reference ID'].str.split('/').str[-1]

# get epitope ids
t_cell_df['Epitope ID'] = t_cell_df['Epitope ID'].str.split('/').str[-1]
b_cell_df['Epitope ID'] = b_cell_df['Epitope ID'].str.split('/').str[-1]

# get parent protein ids
t_cell_df['Parent Protein ID'] = t_cell_df['Parent Protein ID'].str.split('/').str[-1]
b_cell_df['Parent Protein ID'] = b_cell_df['Parent Protein ID'].str.split('/').str[-1]

# get parent protein ids
t_cell_df['Parent Species ID'] = t_cell_df['Parent Species ID'].str.split('/').str[-1].str.split('_').str[-1]
b_cell_df['Parent Species ID'] = b_cell_df['Parent Species ID'].str.split('/').str[-1].str.split('_').str[-1]

# get related object parent protein ids
t_cell_df['Related Object Parent Protein ID'] = t_cell_df['Related Object Parent Protein ID'].str.split('/').str[-1]
b_cell_df['Related Object Parent Protein ID'] = b_cell_df['Related Object Parent Protein ID'].str.split('/').str[-1]

# get related object parent organism ids
t_cell_df['Related Object Parent Organism ID'] = t_cell_df['Related Object Parent Organism ID'].str.split('/').str[-1]
b_cell_df['Related Object Parent Organism ID'] = b_cell_df['Related Object Parent Organism ID'].str.split('/').str[-1]

# get host ids
t_cell_df['Host ID'] = t_cell_df['Host ID'].str.split('_').str[-1]
b_cell_df['Host ID'] = b_cell_df['Host ID'].str.split('_').str[-1]

# get 1st in vivo disease ids
t_cell_df['1st in vivo Disease ID'] = t_cell_df['1st in vivo Disease ID'].str.split('/').str[-1].str.split('_').str[-1]
b_cell_df['1st in vivo Disease ID'] = b_cell_df['1st in vivo Disease ID'].str.split('/').str[-1].str.split('_').str[-1]

# get 2nd in vivo disease ids
t_cell_df['2nd in vivo Disease ID'] = t_cell_df['2nd in vivo Disease ID'].str.split('/').str[-1].str.split('_').str[-1]
b_cell_df['2nd in vivo Disease ID'] = b_cell_df['2nd in vivo Disease ID'].str.split('/').str[-1].str.split('_').str[-1]
print('Done.')
# -------------------- End Data Cleaning and Organizing --------------------

#######################################################
################### AUTOIMMUNE DATA ###################
#######################################################

print('Extracting autoimmune data...')

# read in autoimmune diseases
with open('autoimmune_diseases.json' , 'r') as f:
    diseases = json.load(f)

# select autoimmune epitopes by disease labeling
a_t_cell_epitopes = t_cell_df[(t_cell_df['1st in vivo Disease ID'].isin(list(diseases.keys()))) |
                              (t_cell_df['2nd in vivo Disease ID'].isin(list(diseases.keys())))]

a_b_cell_epitopes = b_cell_df[(b_cell_df['1st in vivo Disease ID'].isin(list(diseases.keys()))) |
                              (b_cell_df['2nd in vivo Disease ID'].isin(list(diseases.keys())))]

# parent species is equivalent to host
a_t_cell_epitopes = a_t_cell_epitopes[a_t_cell_epitopes['Parent Species ID'] == a_t_cell_epitopes['Host ID']]
a_b_cell_epitopes = a_b_cell_epitopes[a_b_cell_epitopes['Parent Species ID'] == a_b_cell_epitopes['Host ID']]

# get reference, unique epitope, and counts for reference, epitope and assay by unique antigen for T cell
a_t_cell_counts = []
for i, j in a_t_cell_epitopes.groupby('Parent Protein ID'):
    diseases = list(set(filter(None, j[['1st in vivo Disease State', '2nd in vivo Disease State']].to_numpy(na_value=None).flatten())))
    counts = []
    counts.append(i)
    counts.append(', '.join(diseases))
    counts.append(list(j['Parent Protein'].dropna())[0])
    counts.append(len(j['Reference ID'].unique()))
    counts.append(len(j['Epitope'].unique()))
    counts.append(len(j))
    counts.append(list(j['Parent Species'].dropna())[0])

    a_t_cell_counts.append(counts)

# get reference, unique epitope, and counts for reference, epitope and assay by unique antigen for B cell
a_b_cell_counts = []
for i, j in a_b_cell_epitopes.groupby('Parent Protein ID'):
    diseases = list(set(filter(None, j[['1st in vivo Disease State', '2nd in vivo Disease State']].to_numpy(na_value=None).flatten())))
    counts = []
    counts.append(i)
    counts.append(', '.join(diseases))
    counts.append(list(j['Parent Protein'].dropna())[0])
    counts.append(len(j['Reference ID'].unique()))
    counts.append(len(j['Epitope'].unique()))
    counts.append(len(j))
    counts.append(list(j['Parent Species'].dropna())[0])

    a_b_cell_counts.append(counts)

# put counts into pandas DataFrame
a_t_cell_counts = pd.DataFrame(a_t_cell_counts, columns=['Parent Protein ID', 'Diseases', 'Parent Protein', 'Reference Count',
                                     'Epitope Count', 'Assay Count', 'Parent Species'])
a_b_cell_counts = pd.DataFrame(a_b_cell_counts, columns=['Parent Protein ID', 'Diseases', 'Parent Protein', 'Reference Count',
                                     'Epitope Count', 'Assay Count', 'Parent Species'])

# limit autoimmune antigens requiring 2 references
a_t_cell_counts = a_t_cell_counts[a_t_cell_counts['Reference Count'] > 1]
a_b_cell_counts = a_b_cell_counts[a_b_cell_counts['Reference Count'] > 1]

print('Done.')

print('Writing autoimmune data...')

# output dataframes to one file
writer = pd.ExcelWriter('autoimmune_data.xlsx', engine='xlsxwriter')

a_t_cell_epitopes.to_excel(writer, sheet_name='T Cell Assays')
a_b_cell_epitopes.to_excel(writer, sheet_name='B Cell Assays')
a_t_cell_counts.to_excel(writer, sheet_name='T Cell Counts')
a_b_cell_counts.to_excel(writer, sheet_name='B Cell Counts')

writer.save()

print('Done.')

# TODO: use requests to get all autoimmune antigen sequences
print('Getting autoimmune antigen sequences from UniProt...')

with open('autoimmune_antigens.fasta', 'w') as f:
    for i in a_t_cell_counts['Parent Protein ID'].append(a_b_cell_counts['Parent Protein ID']):
        r = requests.get('https://www.uniprot.org/uniprot/%s.fasta' % i)
        if not r.text:
            continue
        else:
            f.write(r.text[:-1])

print('Done.')

#######################################################
##################### CANCER DATA #####################
#######################################################

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

# TODO: use requests to get all cancer antigen sequences
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
