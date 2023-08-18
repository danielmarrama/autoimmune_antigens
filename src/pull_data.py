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
from pathlib import Path

from Bio import SeqIO

# TODO:
# * use only source antigen column and see how many are missing without filling in with parent proteins


data_path = Path(__file__).parent.parent / 'data'


def pull_assay_data(assay: str) -> pd.DataFrame:
  """Extracts T cell and B cell positive assay data from the IEDB based on disease ID.

  Args:
    assay (str): assay to pull data from. Either 'tcell' or 'bcell'.
  """
  df = pd.DataFrame()
  for doid in diseases.keys():

    # first get the total number of assays as first request to loop through API
    url = f'https://query-api.iedb.org/{assay}_search'
    params = {'order': 'structure_id',
              'qualitative_measure': 'neq.Negative', # select positive assays only
              'disease_iris': f'cs.{{{"DOID:"+doid}}}',
              'offset': 0}
    
    df = pd.concat([df, iterate_api(url, params)])
  
  # select epitopes where host and source organism are identical
  df = df[(df['host_organism_iri'] == df['parent_source_antigen_source_org_iri']) | # OR
          (df['host_organism_iri'] == df['r_object_source_organism_iri'])]  
  
  if assay == 'tcell': # discard epitopes with less than 5 residues
    df = df[df['linear_sequence'].str.len() >= 5]
  
  return clean_dataframe(df, assay)


def iterate_api(url: str, params: dict) -> pd.DataFrame:
  """IEDB API only allows 10,000 entries per request. We use this function to loop through
  all requests using our URL and parameters until we receive no more data.

  Args:
    url (str): URL to IEDB API.
    params (dict): parameters to pass to the API. 
  """
  df = pd.DataFrame()
  while(True):
    
    r = requests.get(url, params=params)
    r.raise_for_status()

    df = pd.concat([df, pd.read_json(r.text)])
    params['offset'] += 10000
    if r.text == '[]':
      break
  
  return df


def clean_dataframe(df: pd.DataFrame, assay: str) -> pd.DataFrame:
  """Clean the dataframe returned from the IEDB API. This includes extracting
  the relevant information from the curated_source_antigen column, and 
  exploding the disease_iris and disease_names columns into separate rows.

  Args:
    df (pd.DataFrame): DataFrame returned after pulling the IEDB API.
    assay (str): assay type. Either 'tcell' or 'bcell'.
  """
  df['curated_source_antigen_accession'] = df['curated_source_antigen'].map(
    lambda x: x['accession'] if isinstance(x, dict) else None)
  
  df['curated_source_antigen_name'] = df['curated_source_antigen'].map(
    lambda x: x['name'] if isinstance(x, dict) else None)
  
  df['curated_source_organism_id'] = df['curated_source_antigen'].map(
    lambda x: x['source_organism_iri'] if isinstance(x, dict) else None)
  
  df['curated_source_organism_name'] = df['curated_source_antigen'].map(
    lambda x: x['source_organism_name'] if isinstance(x, dict) else None)

  # explode the lists in disease_iris and disease_names to make separate rows for each disease
  df = df.explode(column='disease_iris').reset_index(drop=True)
  df['disease_names'] = df.explode(
    column='disease_names')['disease_names'].reset_index(drop=True)

  # drop curated_source_antigen column
  df.drop('curated_source_antigen', axis=1, inplace=True)

  # add column for assay type
  if assay == 'tcell':
    df['assay_type'] = 'T Cell'
  elif assay == 'bcell':
    df['assay_type'] = 'B Cell'
  else:
    raise ValueError('Assay must be either "tcell" or "bcell".')

  return df


def write_assay_data_to_file(tcell: pd.DataFrame, bcell: pd.DataFrame) -> None:
  """Writes T cell and B cell assay data to Excel file.
  
  Args:
    tcell (pd.DataFrame): autoimmune T cell assay data.
    bcell (pd.DataFrame): autoimmune B cell assay data.
  """
  def filter_by_unique_references(group):
    return group['reference_id'].nunique() >= 2
  
  combined_df = pd.concat([tcell, bcell])

  # keep only rows whereby source antigen has at least two references
  grouped = combined_df.groupby('curated_source_antigen_accession')
  assay_data = grouped.filter(filter_by_unique_references)

  assay_data.to_csv(data_path / 'raw_assay_data.tsv', sep='\t', index=False)


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
  antigens = pd.DataFrame.from_dict(
    antigen_map, 
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

  with open(data_path / 'autoimmune_diseases.json' , 'r') as f:
    diseases = json.load(f)

  print('Pulling autoimmune T cell assay data...')
  tcell = pull_assay_data('tcell')
  print('Done.')

  print('Pull autoimmune B cell assay data...')
  bcell = pull_assay_data('bcell')
  print('Done.')

  print('Writing assay data to file...')
  write_assay_data_to_file(tcell, bcell)
  print('Done.')
  

  # print('Getting autoimmune antigen sequences from UniProt...')
  # pull_uniprot_antigens(antigens)
  # print('Done')

  # print('Getting human proteome from UniProt...')
  # get_human_proteome()
  # print('Done.')

  # print('Removing autoimmune antigens from human proteome...')
  # remove_autoimmune_antigens_from_proteome()
  # print('Done.')

  # print('Combining data into one dataset...')
  # combine_data()
  # print('Done.')