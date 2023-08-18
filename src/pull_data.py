#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import io
import re
import requests
import json
import gzip
import pickle
import pandas as pd

from Bio import SeqIO

# TODO:
# * use only source antigen column and see how many are missing without filling in with parent proteins



def pull_assay_data(table: str) -> pd.DataFrame:
  """Extracts T cell and B cell positive assay data from the IEDB based on disease ID.

  Args:
    table (str): table to pull data from. Either 'tcell' or 'bcell'.
  """
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
  
  if table == 'tcell':
    df = df[df['linear_sequence'].str.len() >= 5]

  # create a column combining epitope and related object source antigen IDs
  df['source_antigen_iri'] = df['parent_source_antigen_iri'].fillna(df['r_object_source_molecule_iri'])
  # create a column combining epitope and related object source antigen names
  df['source_antigen_name'] = df['parent_source_antigen_name'].fillna(df['r_object_source_molecule_name'])
  # create a column combining epitope and related object source organism names
  df['source_organism_name'] = df['source_organism_name'].fillna(df['r_object_source_organism_name'])
  
  return df


def iterate_api(url: str, params: dict) -> pd.DataFrame:
  """IEDB API only allows 10,000 entries per request. We use this function to loop through
  all requests using our URL and parameters until we receive no more data.

  Args:
    url (str): URL to IEDB API.
    params (dict): parameters to pass to the API. 
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


def pull_uniprot_antigens(antigens: pd.DataFrame) -> None:
  """Pulls antigen sequences from UniProt based on protein ID.
  
  Args:
    antigens (pd.DataFrame): DataFrame of antigens with protein IDs.
  """
  with open('autoimmune_antigens.fasta', 'w') as f:
    for uniprot_id in antigens['Protein ID']:
      r = requests.get(f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta')
      if not r.text: # skip if no sequence
        continue
      else:
        f.write(r.text)


def write_data_to_file(
  tcell: pd.DataFrame, bcell: pd.DataFrame, antigens: pd.DataFrame
) -> None:
  """Writes T cell, B cell, and antigen data to Excel file.
  
  Args:
    tcell (pd.DataFrame): T cell assay data.
    bcell (pd.DataFrame): B cell assay data.
    antigens (pd.DataFrame): antigen data - all antigens sourced from assay data.
  """
  with pd.ExcelWriter('autoimmune_data.xlsx', engine='xlsxwriter') as writer:
    tcell.to_excel(writer, sheet_name='T Cell Assays', index=False)
    bcell.to_excel(writer, sheet_name='B Cell Assays', index=False)
    antigens.to_excel(writer, sheet_name='Antigens', index=False)


def get_antigens(tcell: pd.DataFrame, bcell: pd.DataFrame) -> pd.DataFrame:
  """Concatenate all T cell and B cell assay data to get antigens and
  their counts of epitopes, assays, and references.
  
  Args:
    tcell (pd.DataFrame): T cell assay data.
    bcell (pd.DataFrame): B cell assay data.
  """
  def parse_diseases(x: str) -> str:
    """Parse disease names into comma separated list. Args: x: str of disease names."""
    return ', '.join(re.findall('{"(.*?)"}', x))

  # replace the UniProt IDs that are not the canonical protein ID for each gene
  # this is from the protein tree work
  with open('canonical_protein_mapping.pickle', 'rb') as f:
    canonical_protein_mapping = pickle.load(f)

  # separate actual ID from "UNIPROT:ID"
  tcell['source_antigen_iri'] = tcell['source_antigen_iri'].str.split(':').str[1]
  bcell['source_antigen_iri'] = bcell['source_antigen_iri'].str.split(':').str[1]

  # replace IDs with canonical ID
  tcell['source_antigen_iri'] = tcell['source_antigen_iri'].map(canonical_protein_mapping)
  bcell['source_antigen_iri'] = bcell['source_antigen_iri'].map(canonical_protein_mapping)

  # aggregate data by protein ID and get antigen name, organism, diseases, cell type, and 
  # epitope, assay and reference count (T cell)
  antigen_map = {}
  for i, row in tcell.groupby('source_antigen_iri'):
    antigen_map[i] = []
    antigen_map[i].append(list(row['source_antigen_name'].dropna())[0])
    antigen_map[i].append(list(row['source_organism_name'].dropna())[0])
    antigen_map[i].append(len(row['structure_id'].unique()))
    antigen_map[i].append(len(row))
    antigen_map[i].append(len(row['reference_id'].unique()))
    antigen_map[i].append(list(row['disease_names']))
    antigen_map[i].append('T cell')

  # repeat the same for B cell and add counts if ID was already seen in T cell
  for i, row in bcell.groupby('source_antigen_iri'):
    if i in antigen_map.keys():
      antigen_map[i][2] += len(row['structure_id'].unique())
      antigen_map[i][3] += len(row)
      antigen_map[i][4] += len(row['reference_id'].unique())
      antigen_map[i][5] += list(row['disease_names'])
      antigen_map[i][6] += ' and B cell'
    else:
      antigen_map[i] = []
      antigen_map[i].append(list(row['source_antigen_name'].dropna())[0])
      antigen_map[i].append(list(row['source_organism_name'].dropna())[0])
      antigen_map[i].append(len(row['structure_id'].unique()))
      antigen_map[i].append(len(row))
      antigen_map[i].append(len(row['reference_id'].unique()))
      antigen_map[i].append(list(row['disease_names']))
      antigen_map[i].append('B Cell')
          
  # clean up diseases into unique names using set function
  for k, v in antigen_map.items():
    antigen_map[k][5] = ', '.join(set(antigen_map[k][5]))

  # put antigens into pandas DataFrame, reset and rename index
  columns = ['Protein Name', 'Source Organism', 'Epitope Count', 'Assay Count',
              'Reference Count', 'Diseases', 'Targeted By']
  antigens = pd.DataFrame.from_dict(antigen_map, 
                                    orient = 'index', 
                                    columns = columns
                                    ).reset_index().rename(
                                    columns={'index': 'Protein ID'})

  antigens['Diseases'] = antigens['Diseases'].apply(parse_diseases)

  return antigens


def get_human_proteome() -> None:
  """Pulls human proteome from UniProt."""
  r = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz', stream=True)
  with open('human_proteome.fasta', 'wb') as f:
    f.write(gzip.open(r.raw, 'rb').read())


def remove_autoimmune_antigens_from_proteome() -> None:
  """Removes autoimmune antigens from human proteome file."""
  autoimmune_ids = []
  for record in list(SeqIO.parse('autoimmune_antigens.fasta', 'fasta')):
    autoimmune_ids.append(record.id.split('|')[1])

  records = []
  with open('non_autoimmune_proteins.fasta', 'w') as f:
    for record in SeqIO.parse('human_proteome.fasta', 'fasta'):
      if record.id.split('|')[1] in autoimmune_ids:
        continue
      else:
        records.append(record)
    SeqIO.write(records, f, 'fasta')


def pull_data_from_fasta(proteins: list, category: str) -> list:
  """Pulls data from FASTA file and returns list of lists of data.
  
  Args:
    proteins (list): list of SeqIO records.
    category (str): category of proteins. Either 'autoimmune' or 'non_autoimmune'.
  """
  data = []
  autoimmune = 1 if category == 'autoimmune' else 0
  for protein in proteins:
    protein_data = []
    protein_data.append(autoimmune)
    protein_data.append(protein.id.split('|')[1])
    try:
      protein_data.append(re.search('GN=(.*?) ', protein.description).group(1))
    except AttributeError:
      protein_data.append('')
    protein_data.append(int(re.search('PE=(.*?) ', protein.description).group(1)))
    protein_data.append(str(protein.seq))
    data.append(protein_data)
  
  return data


def combine_data() -> None:
  """Combines autoimmune antigen and non-autoimmune protein data into one dataset."""
  ai_antigens  = list(SeqIO.parse('autoimmune_antigens.fasta', 'fasta'))
  non_ai_proteins = list(SeqIO.parse('non_autoimmune_proteins.fasta', 'fasta'))

  data = pull_data_from_fasta(ai_antigens, 'autoimmune') + pull_data_from_fasta(non_ai_proteins, 'non_autoimmune')
  
  df = pd.DataFrame(data, columns=['autoimmune', 'id', 'gene', 'pe_level', 'sequence'])

  # make sure proteins are at least 10 residues
  df = df[df['sequence'].str.len() >= 11]
  df.to_csv('combined_data.csv', index=False)


if __name__ == '__main__':
  print('Extracting autoimmune data...')

  with open('autoimmune_diseases.json' , 'r') as f:
    diseases = json.load(f)

  # read in autoimmune IEDB T cell and B cell assay tables
  print('Reading in autoimmune T cell assay data...')
  tcell = pull_assay_data('tcell', 'autoimmune')
  print('Done.')

  print('Reading in autoimmune B cell assay data...')
  bcell = pull_assay_data('bcell', 'autoimmune')
  print('Done.')
  
  # count antigens and reduce antigens to those with more than 1 reference
  antigens = get_antigens(tcell, bcell)
  antigens = antigens[antigens['Reference Count'] > 1]

  print('Writing autoimmune data...')
  write_data_to_file(tcell, bcell, antigens, 'autoimmune')
  print('Done.')

  print('Getting autoimmune antigen sequences from UniProt...')
  pull_uniprot_antigens(antigens)
  print('Done')

  print('Getting human proteome from UniProt...')
  get_human_proteome()
  print('Done.')

  print('Removing autoimmune antigens from human proteome...')
  remove_autoimmune_antigens_from_proteome()
  print('Done.')

  print('Combining data into one dataset...')
  combine_data()
  print('Done.')