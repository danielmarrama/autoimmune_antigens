#!/usr/bin/env python3

import os, glob

import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score


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

def plot_roc(data):
  category_map = dict(zip(data['id'], data['category']))
  columns = ['query', 'subject', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']

  plt.figure(figsize=(8,8))

  for j in range(1, 6):
    df = pd.read_csv('fold%i_output.csv' % j, names=columns)
    df['query_category'] = df['query'].map(category_map)
    df['subject_category'] = df['subject'].map(category_map)
    
    idx = df.groupby('query')['%_identity'].transform(max) == df['%_identity']
    df = df[idx]

    df['result'] = ''

    for i, row in df.iterrows():
      if row['query_category'] == 'autoimmune' and row['subject_category'] == 'autoimmune':
        df.loc[i, 'result'] = 'TP'
      elif row['query_category'] == 'autoimmune' and row['subject_category'] == 'non_autoimmune':
        df.loc[i, 'result'] = 'FN'
      elif row['query_category'] == 'non_autoimmune' and row['subject_category'] == 'autoimmune':
        df.loc[i, 'result'] = 'FP'
      else:
        df.loc[i, 'result'] = 'TN'
            
    df['trues'] = df['query_category'] == 'autoimmune'
    
    X = df['%_identity']
    y = df['trues']
    fpr, tpr, thresh = roc_curve(y, X)
    roc_auc = metrics.auc(fpr, tpr)
    
    plt.plot(fpr, tpr, label = 'AUC = %0.2f BLAST Model #%i' % (roc_auc, j))
      
  plt.title('ROC Curve for BLAST Models (Each Fold)')
  plt.legend(loc = 'lower right')
  plt.plot([0, 1], [0, 1],'r--')
  plt.xlim([0, 1])
  plt.ylim([0, 1])
  plt.ylabel('True Positive Rate')
  plt.xlabel('False Positive Rate')

  plt.show()

if __name__ == '__main__':
  data = pd.read_csv('combined_data.csv')

  # run_blast(data)
  plot_roc(data)