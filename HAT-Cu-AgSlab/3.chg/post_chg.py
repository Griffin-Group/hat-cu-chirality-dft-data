#!/usr/bin/env python3
print("Starting Python")
import numpy as np
import sys
# Indicies
hat = np.r_[103:151]
cu = np.r_[0:3]
slab = np.r_[3:103]
# Load file
chg = np.genfromtxt(sys.argv[1], skip_header=2, usecols=(4,), max_rows=151)
# Charges
print("Total charge: "+str(sum(chg)))
print("HAT charge: "+str(sum(chg[hat])))
print("Cu charge: "+str(sum(chg[cu])))
print("HAT-Cu charge: "+str(sum(chg[cu])+sum(chg[hat])))
print("Ag slab charge: "+str(sum(chg[slab])))
print("Finished Python")
