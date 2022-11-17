#!/usr/bin/env python3

import os, glob

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline


def write_all_other_antigens(antigen_id):
  records = []
  with open('antigens_minus_%s.fasta' % antigen_id, 'w') as f:
    for record in antigens:
      if antigen_id in str(record.id):
        continue
      else:
        records.append(record)
    SeqIO.write(records, f, 'fasta')

if __name__ == '__main__':

  os.mkdir('results')

  antigens = list(SeqIO.parse('autoimmune_antigens.fasta', 'fasta'))

  autoimmune_ids = []
  for record in antigens:
    autoimmune_ids.append(record.id.split('|')[1])

  for antigen in antigens:
    print(antigen)
    antigen_id = antigen.id.split('|')[1]
    with open(antigen_id + '.fasta', 'w') as f:
      SeqIO.write(antigen, f, 'fasta')

    write_all_other_antigens(antigen_id)

    os.system('./blast/makeblastdb' + ' -in ' + 'antigens_minus_%s.fasta' % antigen_id + ' -dbtype prot')
    blastx_cline = NcbiblastpCommandline(cmd='./blast/blastp',
                                        query = antigen_id + '.fasta',
                                        db = 'antigens_minus_%s.fasta' % antigen_id, 
                                        evalue=100, outfmt=10, out='results/%s_BLAST_results.csv' % antigen_id)

    stdout, stderr = blastx_cline()
    
    os.remove(antigen_id + '.fasta')
    os.remove('antigens_minus_%s.fasta' % antigen_id)

    for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
      os.remove(glob.glob('*.' + extension)[0])