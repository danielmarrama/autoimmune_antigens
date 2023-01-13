#!/usr/bin/env python3

from cdhit_reader import read_cdhit

import pandas as pd

df = pd.read_csv('combined_data.csv', index_col=0)

clusters = read_cdhit('human_db40.clstr')

i = 1
k1, k2, k3, k4, k5 = [], [], [], [], []
for cluster in clusters:
  proteins = [member.name.split('|')[1] for member in cluster.sequences]
  if i == 1:
    k1 += proteins
  elif i == 2:
    k2 += proteins
  elif i == 3:
    k3 += proteins
  elif i == 4:
    k4 += proteins
  elif i == 5:
    k5 += proteins
  i+=1
  if i == 6:
    i = 1

df['fold'] = [None for i in range(len(df))]

df.loc[df['id'].isin(k1), 'fold'] = 1
df.loc[df['id'].isin(k2), 'fold'] = 2
df.loc[df['id'].isin(k3), 'fold'] = 3
df.loc[df['id'].isin(k4), 'fold'] = 4
df.loc[df['id'].isin(k5), 'fold'] = 5

print(df['fold'].value_counts())

df.to_csv('combined_data.csv')