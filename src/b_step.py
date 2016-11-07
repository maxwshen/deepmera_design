

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util


# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a_step/'
NAME = util.get_fn(__file__)

# Functions



@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
