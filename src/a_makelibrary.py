# Uses twobittofa to extract dna sequences from hg19.2bit

import _config
import sys, os, fnmatch, datetime, subprocess, random
import numpy as np

from mylib import util

# Default params
DEFAULT_INP_DIR = _config.DATA_DIR
NAME = util.get_fn(__file__)

TOOL = '/cluster/mshen/tools/2bit/twoBitToFa'
TOOL_DB = '/cluster/mshen/tools/2bit/mm10.2bit'

random.seed(7)

# Functions
def get_genomic_seq(chro, start, end):
  ans = subprocess.check_output(TOOL + ' ' + TOOL_DB + ' -noMask stdout -seq=' + chro + ' -start=' + str(start) + ' -end=' + str(end), shell = True)
  ans = ''.join(ans.split('\n')[1:])
  return ans.strip()

def one_muts(seq):
  nms, seqs = [], []
  ctr = 0
  for i in range(len(seq)):
    for j in ['A', 'C', 'G', 'T']:
      if seq[i] != j:
        seqs.append(seq[:i] + j + seq[i+1:])
        nms.append('single_' + str(ctr))
        ctr += 1
  for s in seqs:
    if len(s) != _config.d.RANDOM_LEN:
      print 'ERROR:', s, len(s)
      break
  return nms, seqs

def two_muts(seq):
  nms, seqs = [], []
  ctr = 0
  for i in range(len(seq)):
    for j in range(i + 1, len(seq)):
      for k in ['A', 'C', 'G', 'T']:
        if seq[i] != k:
          for l in ['A', 'C', 'G', 'T']:
            if seq[j] != l:
              seqs.append(seq[:i] + k + seq[i+1:j] + l + seq[j+1:])
              nms.append('pair_' + str(ctr))
              ctr += 1
  for s in seqs:
    if len(s) != _config.d.RANDOM_LEN:
      print 'ERROR:', s, len(s)
      break
  return nms, seqs

def destroy_dna(seq):
  # Mutate a sequence so that all bases are different
  new = ''
  nts = ['A', 'C', 'G', 'T']
  for i in range(len(seq)):
    while True:
      nt = random.choice(nts)
      if nt != seq[i]:
        new += nt
        break
  return new

def four_window_muts(seq):
  num_replicates = 5
  nms, seqs = [], []
  ctr = 0
  for i in range(0, len(seq)-4, 2):
    for rep in range(num_replicates):
      mut = destroy_dna(seq[i:i+4])
      seqs.append(seq[:i] + mut + seq[i+4:])
      nms.append('four_window_' + str(ctr))
      ctr += 1
  for s in seqs:
    if len(s) != _config.d.RANDOM_LEN:
      print 'ERROR:', s, len(s)
      break
  return nms, seqs

def pairs_four_window_muts(seq):
  num_replicates = 4
  nms, seqs = [], []
  ctr = 0
  for i in range(0, len(seq)-4, 2):
    for j in range(i + 4, len(seq)-4, 2):
      for rep in range(num_replicates):
        mut = destroy_dna(seq[i:i+4])
        mut2 = destroy_dna(seq[j:j+4])
        seqs.append(seq[:i] + mut + seq[i+4:j] + mut2 + seq[j+4:])
        nms.append('pairs_four_window_' + str(ctr))
        ctr += 1
  for s in seqs:
    if len(s) != _config.d.RANDOM_LEN:
      print 'ERROR:', s, len(s)
      break
  return nms, seqs

def comb_dna(seq):
  new = ''
  nts = ['A', 'C', 'G', 'T']
  for i in range(len(seq)):
    if i % 2 == 0:
      while True:
        candidate = random.choice(nts)
        if candidate != seq[i]:
          new += candidate
          break
    else:
      new += seq[i]
  return new

def three_comb_muts(seq):
  num_replicates = 5
  nms, seqs = [], []
  ctr = 0
  for i in range(0, len(seq)-5):
    for rep in range(num_replicates):
      mut = comb_dna(seq[i:i+5])
      seqs.append(seq[:i] + mut + seq[i+5:])
      nms.append('three_comb_' + str(ctr))
      ctr += 1
  for s in seqs:
    if len(s) != _config.d.RANDOM_LEN:
      print 'ERROR:', s, len(s)
      break
  return nms, seqs

def make_library(out_dir):
  for i in range(len(_config.d.NAMES)):
    nm, chro, pos = _config.d.NAMES[i], _config.d.CHRMS[i], _config.d.POSS[i]

    start = pos - _config.d.RANDOM_LEN / 2
    end = pos + _config.d.RANDOM_LEN / 2
    if _config.d.RANDOM_LEN % 2 == 1:
      start -= 1

    seq = get_genomic_seq(chro, start, end)

    all_names, all_seqs = [], []
    nms, seqs = one_muts(seq)
    all_names += nms
    all_seqs += seqs

    nms, seqs = two_muts(seq)
    all_names += nms
    all_seqs += seqs

    nms, seqs = four_window_muts(seq)
    all_names += nms
    all_seqs += seqs

    nms, seqs = pairs_four_window_muts(seq)
    all_names += nms
    all_seqs += seqs

    nms, seqs = three_comb_muts(seq)
    all_names += nms
    all_seqs += seqs

    print nm, start, end, end - start, seq
    print len(all_names), len(nms)

    print 'Total oligos :', len(all_seqs)
    if len(all_seqs) == _config.d.MAX_OLIGOS:
      print 'SUCCESS: Matches', _config.d.MAX_OLIGOS
    else:
      print 'WARNING: Fails to match', _config.d.MAX_OLIGOS
      print 'WARNING: Difference is ', _config.d.MAX_OLIGOS - len(all_seqs)

    out_fn = out_dir + nm + '_library.csv'
    with open(out_fn, 'w') as f:
      pass

  return



@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  make_library(out_dir)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')