#!/usr/bin/env python3


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy import stats

aa_map = {'A': 'alanines', 'R': 'arginines', 'N': 'asparagines', 'D': 'aspartates', 'C': 'cysteines',
          'Q': 'glutamines', 'E': 'glutamates', 'G': 'glycines', 'H': 'histidines', 'I': 'isoleucines',
          'L': 'leucines', 'K': 'lysines', 'M': 'methionines', 'F': 'phenylalanines', 'P': 'prolines',
          'S': 'serines', 'T': 'threonines', 'W': 'tryptophans', 'Y': 'tyrosines', 'V': 'valines'}

df = pd.read_csv('master_data.csv', index_col=0)

def select_aa():
  return input('Enter the 1-letter symbol of the amino acid to compare.\n')

def select_feature():
  return input('Select a feature to compare: length (1), amino acid composition (2)\n')

def x_label(feature):
  if feature == 'length':
    return 'Length'
  elif feature in aa_map.values():
    return 'Amino Acid Composition (%)'

def show_distr_plot(feature):
  a = list(df[df['category'] == 'autoimmune'][feature])
  c = list(df[df['category'] == 'cancer'][feature])
  n = list(df[df['category'] == 'normal'][feature])

  bins = 100

  plt.figure(figsize=(10, 8))
  plt.hist(a, bins=bins, weights=np.ones(len(a)) / len(a), alpha=0.5, label='autoimmune');
  plt.hist(c, bins=bins, weights=np.ones(len(c)) / len(c), alpha=0.5, label='cancer');
  plt.hist(n, bins=bins, weights=np.ones(len(n)) / len(n), alpha=0.5, label='normal');
 
  plt.title('Distributions of antigen type by ' + feature + '.')
  plt.xlabel(x_label(feature))
  plt.ylabel('Percentage of Occurrences')
  plt.legend(loc='upper right')

  plt.annotate('Mann-Whitney\nAutoimmune-Cancer:\np=' + 
      str(format(stats.mannwhitneyu(a,c,)[1], '.2e')), (750, 550), xycoords='figure pixels')
  plt.annotate('Mann-Whitney\nAutoimmune-Normal:\np=' + 
      str(format(stats.mannwhitneyu(a,n,)[1], '.2e')), (750, 500), xycoords='figure pixels')
  plt.annotate('Mann-Whitney\nCancer-Normal:\np=' + 
      str(format(stats.mannwhitneyu(c,n,)[1], '.2e')), (750, 450), xycoords='figure pixels')
  plt.annotate('Kruskal-Wallis:\np=' + 
      str(format(stats.kruskal(a,c,n)[1], '.2e')), (750, 400), xycoords='figure pixels')

  plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
  plt.show()

def run():
  f = select_feature()
  
  match f:
    case '1':
      show_distr_plot('length')
    case '2':
      aa = select_aa()
      if aa not in 'ARNDCQEGHILKMFPSTWYV':
        print('Amino acid not valid.')
        return

      show_distr_plot(aa_map[aa])

    case _:
      print('Invalid selection.') 

run()
