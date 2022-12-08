#!/usr/bin/env python3

import os, glob

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

def remove_files(fold):
  '''
  Removes db.fasta file and BLAST DB files created for each fold.
  '''
  os.remove('db.fasta')
  os.remove('fold%i.fasta' % fold)
  for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
    os.remove(glob.glob('*.' + extension)[0])

def run_blast(data):
  for fold in range(1, 6):
    
    seqs = []
    for i, row in data[data['fold'] == fold].iterrows():
      seqs.append(SeqRecord(row['sequence'], id=row['id'], description=''))

    with open('fold%i.fasta' % fold, 'w') as f1:
      SeqIO.write(seqs, f1, 'fasta')

    seqs = []
    for i, row in data[data['fold'] != fold].iterrows():
      seqs.append(SeqRecord(row['sequence'], id=row['id'], description=''))

    with open('db.fasta', 'w') as f2:
      SeqIO.write(seqs, f2, 'fasta')
    
    os.system('../blast/makeblastdb' + ' -in ' + 'db.fasta' + ' -dbtype prot')
    blastx_cline = NcbiblastpCommandline(
                cmd='../blast/blastp',
                query = 'fold%i.fasta' % fold,
                db = 'db.fasta', 
                evalue=1, outfmt=10, out='fold%i_output.csv' % fold)

    stdout, stderr = blastx_cline()

    remove_files(fold)

if __name__ == '__main__':
    data = pd.read_csv('../combined_data.csv')

    run_blast(data)