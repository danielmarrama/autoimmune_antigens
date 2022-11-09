#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import io
import requests
import json
import gzip
import pandas as pd

from Bio import SeqIO

# TODO:
# - Get all human protein data and remove autoimmune/cancer antigens


def pull_uniprot_antigens(counts, antigen_type):
  # use requests to get all autoimmune antigen sequences
  with open('%s_antigens.fasta' % antigen_type, 'w') as f:
      for idx in counts['Protein ID']:
          print(idx)
          i = idx.split(':')[1]
          r = requests.get('https://www.uniprot.org/uniprot/%s.fasta' % i)
          if not r.text:
              continue
          else:
              f.write(r.text)
  return 0

def write_data_to_file(tcell, bcell, counts, antigen_type):
  # output dataframes to one file
  writer = pd.ExcelWriter('%s_data.xlsx' % antigen_type, engine='xlsxwriter')

  tcell.to_excel(writer, sheet_name='T Cell Assays', index=False)
  bcell.to_excel(writer, sheet_name='B Cell Assays', index=False)
  counts.to_excel(writer, sheet_name='Antigen Counts', index=False)

  writer.save()

  return 0 

def get_counts(tcell, bcell):

  # aggregate data by protein ID and get antigen name, organism, diseases, cell type, and 
  # epitope, assay and reference count (T cell)
  count_map = {}
  for i, row in tcell.groupby('source_antigen_iri'):
    count_map[i] = []
    count_map[i].append(list(row['source_antigen_name'].dropna())[0])
    count_map[i].append(list(row['source_organism_name'].dropna())[0])
    count_map[i].append(len(row['structure_id'].unique()))
    count_map[i].append(len(row))
    count_map[i].append(len(row['reference_id'].unique()))
    count_map[i].append(list(row['disease_names']))
    count_map[i].append('T cell')

  # repeat the same for B cell and add counts if ID was already seen in T cell
  for i, row in bcell.groupby('source_antigen_iri'):
    if i in count_map.keys():
      count_map[i][2] += len(row['structure_id'].unique())
      count_map[i][3] += len(row)
      count_map[i][4] += len(row['reference_id'].unique())
      count_map[i][5] += list(row['disease_names'])
      count_map[i][6] += ' and B cell'
    else:
      count_map[i] = []
      count_map[i].append(list(row['source_antigen_name'].dropna())[0])
      count_map[i].append(list(row['source_organism_name'].dropna())[0])
      count_map[i].append(len(row['structure_id'].unique()))
      count_map[i].append(len(row))
      count_map[i].append(len(row['reference_id'].unique()))
      count_map[i].append(list(row['disease_names']))
      count_map[i].append('B Cell')
          
  # clean up diseases into unique names using set function
  for k, v in count_map.items():
    count_map[k][5] = ', '.join(set(count_map[k][5]))

  # put counts into pandas DataFrame, reset and rename index
  columns = ['Protein Name', 'Source Organism', 'Epitope Count', 'Assay Count',
            'Reference Count', 'Diseases', 'Targeted By']
  counts = pd.DataFrame.from_dict(count_map, 
                                  orient = 'index', 
                                  columns = columns
                                  ).reset_index().rename(
                                  columns={'index': 'Protein ID'})

  return counts

def iterate_api(url, params):
  """
  IEDB API only allows 10,000 entries per request. We use this function to loop through
  all requests using our URL and parameters until we receive no more data.
  """
  df = pd.DataFrame()
  while(True):
    s = requests.get(url, params=params, headers={'accept': 'text/csv', 'Prefer': 'count=exact'})
    try:
      df = pd.concat([df, pd.read_csv(io.StringIO(s.content.decode('utf-8')))])
      params['offset'] += 10000
    except pd.errors.EmptyDataError:
      break
  
  return df

def pull_iedb_data(table, antigen_type):
  '''
  Extracts T cell and B cell positive assay data from the IEDB.
  antigen_type = 'autoimmune' or 'cancer'
  Parameter table = 'tcell' or 'bcell' to specify. 
  '''
  if antigen_type == 'autoimmune':
    df = pd.DataFrame()
    for doid in diseases.keys(): 

      # first get the total number of assays as first request to loop through API
      url = 'https://query-api.iedb.org/%s_search' % table
      params = {'order': 'structure_id',
                'qualitative_measure': 'neq.Negative', # select positive assays only
                'disease_iris': f'cs.{{{"DOID:"+doid}}}',
                'offset': 0}
      
      df = pd.concat([df, iterate_api(url, params)])
    
    # select epitopes where host and source organism are identical
    df = df[(df['host_organism_iri'] == df['parent_source_antigen_source_org_iri']) | # OR
            (df['host_organism_iri'] == df['r_object_source_organism_iri'])]
  
  elif antigen_type == 'cancer':

    # first get the total number of assays as first request to loop through API
    url = 'https://query-api.iedb.org/%s_search' % table
    params = {'order': 'structure_id',
              'qualitative_measure': 'neq.Negative',
              'or': '(e_related_object_type.eq.neo-epitope, immunization_description.plfts."Occurrence of cancer")',
              'offset': 10000}
      
    df = iterate_api(url, params)
  
  else:
    raise(ValueError('Antigen type not specified. Antigen type requested : %s' % antigen_type))
  
  # create a column combining epitope and related object source antigen IDs
  df['source_antigen_iri'] = df['parent_source_antigen_iri'].fillna(df['r_object_source_molecule_iri'])
  # create a column combining epitope and related object source antigen names
  df['source_antigen_name'] = df['parent_source_antigen_name'].fillna(df['r_object_source_molecule_name'])
  # create a column combining epitope and related object source organism names
  df['source_organism_name'] = df['source_organism_name'].fillna(df['r_object_source_organism_name'])
  
  return df

if __name__ == '__main__':
  print('Extracting autoimmune data...')

  with open('autoimmune_diseases.json' , 'r') as f:
      diseases = json.load(f)

  # read in autoimmune IEDB T cell and B cell assay tables
  print('Reading in autoimmune T cell assay data...')
  tcell = pull_iedb_data('tcell', 'autoimmune')
  print('Done.')

  print('Reading in autoimmune B cell assay data...')
  bcell = pull_iedb_data('bcell', 'autoimmune')
  print('Done.')
  
  # count antigens and reduce antigens to those with more than 1 reference
  counts = get_counts(tcell, bcell)
  counts = counts[counts['Reference Count'] > 1]

  print('Writing autoimmune data...')
  write_data_to_file(tcell, bcell, counts, 'autoimmune')
  print('Done.')

  print('Getting autoimmune antigen sequences from UniProt...')
  pull_uniprot_antigens(counts, 'autoimmune')
  print('Done')

  print('Extracting cancer data...')
  print('Reading in cancer T cell assay data...')
  tcell = pull_iedb_data('tcell', 'cancer')
  print('Done.')

  print('Reading in cancer B cell assay data...')
  bcell = pull_iedb_data('bcell', 'cancer')
  print('Done.')

  counts = get_counts(tcell, bcell)

  print('Writing cancer data...')
  write_data_to_file(tcell, bcell, counts, 'cancer')
  print('Done.')

  print('Getting cancer antigen sequences from UniProt...')
  pull_uniprot_antigens(counts, 'cancer')
  print('Done.')


  # get human proteome
  r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz', stream=True)
  with open('9606.fasta', 'wb') as f:
    f.write(gzip.open(r.raw, 'rb').read())

  # read in human proteome with Biopython and remove autoimmune and cancer antigens
  human_protein_ids = []
  for record in list(SeqIO.parse('9606.fasta', 'fasta')):
    human_protein_ids.append(record.id.split('|')[1])
