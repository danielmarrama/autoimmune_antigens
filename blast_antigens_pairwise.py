#!/usr/bin/env python3

import os, glob

import pandas as pd

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline


def write_antigen_file(antigen_id):
  '''Creates a .fasta file with just the one antigen to query.'''
  with open(antigen_id + '.fasta', 'w') as f:
    SeqIO.write(antigen, f, 'fasta')

def write_all_other_antigens(antigen_id):
  '''Creates a .fasta file with all other antigens minus the one to query.'''
  records = []
  with open('antigens_minus_%s.fasta' % antigen_id, 'w') as f:
    for record in antigens:
      if antigen_id in str(record.id):
        continue
      else:
        records.append(record)
    SeqIO.write(records, f, 'fasta')

def blast_antigen(antigen_id):
  '''
  Command line function from BioPython to run BlastP. 
  '''
  os.system('./blast/makeblastdb' + ' -in ' + 'antigens_minus_%s.fasta' % antigen_id + ' -dbtype prot')
  blastx_cline = NcbiblastpCommandline(
                    cmd='./blast/blastp',
                    query = antigen_id + '.fasta',
                    db = 'antigens_minus_%s.fasta' % antigen_id, 
                    evalue=100, outfmt=7, out='blast_results/%s_BLAST_results.csv' % antigen_id)

  stdout, stderr = blastx_cline()

def remove_files(antigen_id):
  '''
  Removes .fasta files and BLAST DB files created for each antigen.
  '''
  os.remove(antigen_id + '.fasta')
  os.remove('antigens_minus_%s.fasta' % antigen_id)

  for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
    os.remove(glob.glob('*.' + extension)[0])

def create_evalue_matrix(antigens):
  '''
  After running BLAST with one antigen against the rest,
  get the e-values into a matrix of n x n, where n = # of antigens.
  '''
  columns = ['query', 'subject', '% identity', 'align length', 'mismatches', 'gap opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit score']
  index = []
  matrix_map = {} 
  for antigen1 in antigens:
    antigen1_id = antigen1.id.split('|')[1]
    index.append(antigen1_id)
    
    # keep track of evalues for all other antigens in a dict of lists
    matrix_map[antigen1_id] = []

    # read in each antigen BLAST results 
    df = pd.read_csv('../blast_results/%s_BLAST_results.csv' % antigen1_id, skiprows=5, sep='\t', names=columns).iloc[:-1, :]
    
    # this creates a dict from two columns, the matched antigen and it's e-value
    id_to_evalue = dict(zip(df['subject'].str.split('|')[1], df['evalue']))
    for antigen2 in antigens:
      antigen2_id = antigen2.id.split('|')[1]
      if antigen1_id == antigen2_id:
        matrix_map[antigen1_id].append(1)
      else:
        try:
          matrix_map[antigen1_id].append(id_to_evalue[antigen2_id])
        except KeyError:
          matrix_map[antigen1_id].append(0) # if there is no BLAST match
  
  df = pd.DataFrame.from_dict(matrix_map)
  df.index = index
  df.to_csv('evalue_matrix.csv')


if __name__ == '__main__':

  os.mkdir('blast_results')
  antigens = list(SeqIO.parse('autoimmune_antigens.fasta', 'fasta'))

  autoimmune_ids = []
  for record in antigens:
    autoimmune_ids.append(record.id.split('|')[1])

  for antigen in antigens:
    antigen_id = antigen.id.split('|')[1]
    write_antigen_file(antigen_id)
    write_all_other_antigens(antigen_id)
    blast_antigen(antigen_id)
    remove_files(antigen_id)
  
  create_evalue_matrix(antigens)