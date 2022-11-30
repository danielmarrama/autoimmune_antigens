#!/usr/bin/env python3

import os, glob
import pandas as pd

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

def hobohm1(cutoff):
  '''
  This algorithm reduces a list of sequences. It starts by keeping the
  longest sequence in a unique list and BLASTing the next longest 
  sequence against that. If the percentage identity is above a certain
  cutoff, it is removed, otherwise it is add the the list and the next 
  sequence is BLASTed against the unique list.
  '''
  columns = ['query', 'subject', '% identity', 'align length', 'mismatches', 'gap opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit score']
  antigens = list(SeqIO.parse('autoimmune_antigens.fasta', 'fasta'))
  antigens.sort(key=lambda x: -len(x)) # sort antigens by length
  
  unique_list = [antigens[0]] # start with longest antigen in unique list
  
  for antigen in antigens[1:]:
    write_unique_list(unique_list)
    write_single_antigen(antigen)
    blast_antigen(antigen)

    results = pd.read_csv('hobohm_results.csv', skiprows=5, sep='\t', names=columns).iloc[:-1, :]
    
    # check if any matches are >90% identity
    if (results['% identity'] > cutoff).any(): 
      continue
    else:
      unique_list.append(antigen)

    remove_files(antigen)

def write_unique_list(unique_list):
  '''
  Creates a .fasta file with the unique antigen list to BLAST
  the next antigen against.
  '''
  with open('unique_list.fasta', 'w') as f:
    SeqIO.write(unique_list, f, 'fasta')

def write_single_antigen(antigen):
  '''Creates a .fasta file with just the one antigen to query.'''
  with open(str(antigen.id).split('|')[1] + '.fasta', 'w') as f:
    SeqIO.write(antigen, f, 'fasta')

def blast_antigen(antigen):
  '''
  Command line function from BioPython to run BlastP. 
  '''
  os.system('./blast/makeblastdb' + ' -in ' + 'unique_list.fasta' + ' -dbtype prot')
  blastx_cline = NcbiblastpCommandline(
                    cmd='./blast/blastp',
                    query = str(antigen.id).split('|')[1] + '.fasta',
                    db = 'unique_list.fasta', 
                    evalue=100, outfmt=7, out='hobohm_results.csv')

  stdout, stderr = blastx_cline()

def remove_files(antigen):
  '''
  Removes .fasta files and BLAST DB files created for each antigen.
  '''
  os.remove('unique_list.fasta')
  os.remove('hobohm_results.csv')
  os.remove(str(antigen.id).split('|')[1] + '.fasta')

  for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
    os.remove(glob.glob('*.' + extension)[0])


if __name__ == '__main__':
  hobohm1(90.0)