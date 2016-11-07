# Uses twobittofa to extract dna sequences from hg19.2bit

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util


# Default params
DEFAULT_INP_DIR = _config.DATA_DIR
NAME = util.get_fn(__file__)

TOOL = '/cluster/mshen/tools/2bit/twoBitToFa'
TOOL_DB = '/cluster/mshen/tools/2bit/hg19.2bit'
BUFF = _config.d.DEL_LEN + 10
# Functions

def get_genomic_seqs(loc):
  seqs = []
  for i in range(len(loc)):
    chro = loc[i].split(':')[0]
    start = int(loc[i].split(':')[1].split('-')[0]) - BUFF
    end = int(loc[i].split(':')[1].split('-')[1]) + BUFF
  
    ans = subprocess.check_output(TOOL + ' ' + TOOL_DB + ' -noMask stdout -seq=' + chro + ' -start=' + str(start) + ' -end=' + str(end), shell = True)
    ans = ''.join(ans.split('\n')[1:])
    seqs.append(ans.strip())
  return seqs

def get_cut_sites(spacer_num, spacer):
  # cut occurs immediately after the base of the cut site number
  cuts = []
  for i in range(len(spacer)):
    if spacer[i][:2] == 'CC' and spacer[i][-2:] == 'GG':
      spnum = spacer_num[i]
      if spnum in _config.d.AMBI:
        if _config.d.AMBI[spnum] == '+':
          cuts.append( len(spacer[i]) - 7 + BUFF - 1)
        if _config.d.AMBI[spnum] == '-':
          cuts.append(5 + BUFF)
      else:
        print 'Ambiguous', spnum, spacer[i]
        cuts.append(0)
      continue
    if spacer[i][:2] == 'CC':
      cuts.append(5 + BUFF)
    if spacer[i][-2:] == 'GG':
      cuts.append( len(spacer[i]) - 7 + BUFF - 1)
  return [str(s) for s in cuts]

def form_input(inp_fn, out_dir):
  spacer_num, spacer, loc = [], [], []
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i in _config.d.IGNORE:
        continue
      words = line.split()    
      spacer_num.append(words[0])
      spacer.append(words[1])
      loc.append(words[2])

  seqs = get_genomic_seqs(loc)
  cuts = get_cut_sites(spacer_num, spacer)

  out_fn = out_dir + 'data.txt'
  with open(out_fn, 'w') as f:
    for i in range(len(spacer_num)):
      line = '\t'.join([spacer_num[i], cuts[i], seqs[i]])
      f.write(line + '\n')
  return

@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  form_input(_config.LOC_FN, out_dir)
  util.shell_cp(inp_dir + 'Spacer2Sequence.txt', out_dir)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')