#!/usr/bin/env python3

import os

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline


def run_blast(data):
  for i in range(0, 5):
    seqs = []
    for i, row in data[data['fold'] != 1].iterrows():
      seqs.append(SeqRecord(row['sequence'], id=row['id'], description=''))

    with open('db.fasta', 'w') as f:
      SeqIO.write(seqs, f, 'fasta')
    
    os.system('../blast/makeblastdb' + ' -in ' + 'db.fasta' + ' -dbtype prot')
    blastx_cline = NcbiblastpCommandline(
                cmd='../blast/blastp',
                query = 'fold%i.fasta' % i+1,
                db = 'db.fasta', 
                evalue=100, outfmt=7, out='fold%i_output.csv' % i+1)

    stdout, stderr = blastx_cline()

    remove_files()

def remove_files():
  '''
  Removes .fasta files and BLAST DB files created for each antigen.
  '''
  os.remove('db.fasta')
  for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
    os.remove(glob.glob('*.' + extension)[0])

if __name__ == '__main__':
    data = pd.read_csv('../combined_data.csv')

    run_blast(data)