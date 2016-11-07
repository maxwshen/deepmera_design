# Stores project-specific functions

import _config

test = 1

def get_grna_fn(gene):
  return _config.DATA_DIR + 'expdesign2/' + gene.capitalize() +
   '_gRNA.txt'