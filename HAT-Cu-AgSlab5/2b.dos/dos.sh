# Summing dos
sum_dos_np 0 1 3
mv DOS.SUM.[1-9]* DOS.SUM.Cu
# TODO sum_dos_np_lm not on NERSC
sum_dos_np_lm 0 1 3
mv DOS.SUM.[1-9]* DOS.SUM.Cu.lm

sum_dos_np 0 154 201
mv DOS.SUM.[1-9]* DOS.SUM.HAT

sum_dos_np_lm 0 154 201
mv DOS.SUM.[1-9]* DOS.SUM.HAT.lm

sum_dos_np 0 4 153
mv DOS.SUM.[1-9]* DOS.SUM.AgSlab

sum_dos_np_lm 0 4 153
mv DOS.SUM.[1-9]* DOS.SUM.AgSlab.lm
