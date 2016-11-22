# Config parameters
# imported by src/_config

HOMOLOGY_ARM_LEN = 25
RANDOM_LEN = 28
MAX_OLIGOS = 3918

NAMES = ['MshDS1', 'MshDS2', 'Tdgf1DS1', 'Tdgf1DS2']
CHRMS = ['chr17', 'chr17', 'chr9', 'chr9']
# Mark the center of the motif
POSS = [87672305, 87697471, 110946162, 110954373]
PAMS = [(87672295, '-'), (87697482, '+'), (110946149, '-'), (110954379, '+')]

# Used for controls
MOTIF = {'MshDS1': [(87672298, 87672307, 'Klf/Sp')], 
  'MshDS2': [(87697455, 87697464, 'Ets'), (87697472, 87697482, 'E2F')], 
  'Tdgf1DS1': [(110946154, 110946165, 'Klf/Sp'), (110946149, 110946156, 'Tcfap2')], 
  'Tdgf1DS2': [(110954366, 110954378, 'Ctcf')] }

NUM_GFPPOS_CTRL_MUTATIONS = 3