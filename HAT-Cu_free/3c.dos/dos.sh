# Summing dos
sum_dos_np 0 25 27
mv DOS.SUM.[1-9]* DOS.SUM.Cu

sum_dos_np_lm 0 25 27
mv DOS.SUM.[1-9]* DOS.SUM.Cu.lm

sum_dos_np 1 {1..24} {28..51}
if [ -e DOS.tmp ]
then mv DOS.tmp DOS.SUM.HAT
else mv DOS.SUM.[1-9]* DOS.SUM.HAT
fi

sum_dos_np_lm 1 {1..24} {28..51}
if [ -e DOS.tmp ]
then mv DOS.tmp DOS.SUM.HAT.lm
else mv DOS.SUM.[1-9]* DOS.SUM.HAT.lm
fi
