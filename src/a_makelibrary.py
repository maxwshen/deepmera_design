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

def ensure_crispr_untargetable(all_seqs, nm, start, end, seq):
  # Sufficient criteria for untargetability:
  #   - At least 1 mutation in PAM
  #   - At least 1 mutation in 10bp next to PAM
  #   - At least 2 mutations in 10bp in gRNA distal from PAM
  PAM = _config.d.PAMS[_config.d.NAMES.index(nm)]
  pos = PAM[0] - start
  if PAM[0] + 2 > end or PAM[0] < start:
    print '\tPAM location', PAM, 'is not fully within designed region', start, end
    sys.exit(0)
  edited = 0
  newallseqs = []
  for s in all_seqs:
    ms = [1 if s[i] != seq[i] else 0 for i in range(len(seq))]
    
    # positive orientation
    if PAM[1] == '+':
      if sum(ms[pos - 10 : pos + 3]) > 0:
        pass
      elif sum(ms[pos - 20 : pos - 10]) > 1:
        pass
      else:
        edit = random.choice([1, 2])
        temp_s = list(s)
        temp_s[pos + edit] = random.choice(['A', 'C', 'T'])
        s = ''.join(temp_s)
        edited += 1

    # reverse orientation
    if PAM[1] == '-':
      if sum(ms[pos + 0 : pos + 13]) > 0:
        pass
      elif sum(ms[pos + 13 : pos + 23]) > 1:
        pass
      else:
        edit = random.choice([0, 1])
        temp_s = list(s)
        temp_s[pos + edit] = random.choice(['A', 'G', 'T'])
        s = ''.join(temp_s)
        edited += 1

    newallseqs.append(s)

  print '\tEdited', edited, 'sequences'
  return newallseqs

def make_library(out_dir):
  for _i in range(len(_config.d.NAMES)):
    nm, chro, pos = _config.d.NAMES[_i], _config.d.CHRMS[_i], _config.d.POSS[_i]

    start = pos - _config.d.RANDOM_LEN / 2
    end = pos + _config.d.RANDOM_LEN / 2
    if _config.d.RANDOM_LEN % 2 == 1:
      start -= 1

    seq = get_genomic_seq(chro, start, end)
    print nm, start, end, end - start, seq

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

    all_seqs = ensure_crispr_untargetable(all_seqs, nm, start, end, seq)
    all_seqs = ensure_crispr_untargetable(all_seqs, nm, start, end, seq)

    print '\t', len(all_names), len(nms)

    print '\tTotal oligos :', len(all_seqs)
    if len(all_seqs) == _config.d.MAX_OLIGOS:
      print '\tSUCCESS: Matches', _config.d.MAX_OLIGOS
    else:
      print '\tWARNING: Fails to match', _config.d.MAX_OLIGOS
      print '\tWARNING: Difference is ', _config.d.MAX_OLIGOS - len(all_seqs)

    out_fn = out_dir + nm + '_library.csv'
    with open(out_fn, 'w') as f:
      for i in range(len(all_names)):
        f.write(all_names[i] + ',' + all_seqs[i] + '\n')

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