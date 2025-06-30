
print("Starting Python")
import numpy as np
import sys
# Indicies
hat = np.r_[0:24,27:51]
cu = np.r_[24:27]
# Load file
chg = np.genfromtxt(sys.argv[1], skip_header=2, usecols=(4,), max_rows=51)
# Charges
print("Total charge: "+str(sum(chg)))
print("HAT charge: "+str(sum(chg[hat])))
print("Cu charge: "+str(sum(chg[cu])))
print("Finished Python")
