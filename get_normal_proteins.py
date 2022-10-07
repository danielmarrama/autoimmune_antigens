#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import requests


# read in datasets
ai_t = pd.read_excel('data/autoimmunity/human_autoimmune_data.xlsx', sheet_name = 1)
ai_b = pd.read_excel('data/autoimmunity/human_autoimmune_data.xlsx', sheet_name = 3)

cancer_t = pd.read_excel('cancer/cancer_data.xlsx', index_col = 0, sheet_name = 1)
cancer_b = pd.read_excel('cancer/cancer_data.xlsx', index_col = 0, sheet_name = 3)

# combine T cell and B cell data
ai =     pd.concat([ai_t, ai_b]).dropna(subset='Antigen ID')
cancer = pd.concat([cancer_t, cancer_b]).dropna(subset='Antigen ID')

antigen_ids = {}
antigen_ids['autoimmune'] = list(set(list(ai_t['Antigen ID'].dropna())).union(set(list(ai_b['Antigen ID'].dropna()))))
antigen_ids['cancer'] = list(set(list(cancer_t['Antigen ID'].dropna())).union(set(list(cancer_b['Antigen ID'].dropna()))))

all_proteins = {}
for i in list(SeqIO.parse('9606.fasta', 'fasta')):
  if str(i.id).split('|')[1] in antigen_ids['autoimmune'] or str(i.id).split('|')[1] in antigen_ids['cancer']:
    continue
  else:
    all_proteins[str(i.id).split('|')[1]] = str(i.seq)


count = 0
with open('normal_proteins_4.fasta', 'wb') as f:
  for i in list(all_proteins.keys())[12000:]:
    print(count + 1)
    r = requests.get('https://www.uniprot.org/uniprot/' + i + '.fasta')
    f.write(r.content) 
    count += 1
