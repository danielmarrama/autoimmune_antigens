#!/usr/bin/env python3

import pandas as pd

def aa_composition(data):
  for aa in 'ARNDCQEGHILKMFPSTWYV':
    data['%s_composition' % aa] = (data['sequence'].str.count(aa) / data['sequence'].str.len())
  
  return data


if __name__ == '__main__':
  data = pd.read_csv('combined_data.csv')
  data['autoimmune'] = data['category'] == 'autoimmune'
  data.drop('category', axis=1, inplace=True)

  data = aa_composition(data)

