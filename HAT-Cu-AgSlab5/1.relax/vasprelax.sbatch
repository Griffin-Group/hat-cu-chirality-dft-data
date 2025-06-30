#!/bin/bash
#SBATCH --job-name=VASPRelax
# Partition information
#SBATCH --account=mp149
#SBATCH --qos=regular
#SBATCH --licenses=scratch
#SBATCH --constraint=cpu
# Node information
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
# Job timing
#SBATCH --time=18:00:00
# Output
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=BField@lbl.gov
#SBATCH --output=JobRelax-%j.out

module load cpu
module load vasp/6.4.3-cpu
module list
# If you don't have access to the `vasp` group on NERSC, you can't use their 
# binaries from the module. Uncomment this line to use the Griffin group version
# Please read $CFS/m3349/codes/vasp/README.md if you use this version
#export PATH=$CFS/m3349/codes/vasp/vasp.6.4.3-cpu/bin/:$PATH

# Environment variables for Slingshot Update Sept. 2024
# Large jobs may break, especially during heavy IO, without this
export FI_CXI_DEFAULT_TX_SIZE=16384
export FI_CXI_RDZV_THRESHOLD=65536

export OMP_NUM_THREADS=1

shopt -s extglob

echo "The current job ID is $SLURM_JOB_ID"
echo "Memory allocated (per node): $SLURM_MEM_PER_NODE"
echo "Running on $SLURM_JOB_NUM_NODES nodes"
echo "Node list:"
sacct --format=JobID,Nodelist%100 -j $SLURM_JOB_ID

#This script does a structural relaxation.

echo "Starting relaxation at `date`" | tee -a track
echo 
#Initialise files
cp incar.relax INCAR || exit 1

# Counting number of atoms assumes elements are named.
natom=`awk 'NR==7 {for (i=1; i<=NF; i++) sum+= $i; print sum}' POSCAR`

# To track the history, we'll keep files
# Get the number of the current highest extension (needs extglob)
exti=`ls -1v POSCAR.+([0-9]) 2>/dev/null | sed 's/POSCAR.//' | tail -n 1`
if [ -z $exti ]; then
	exti=0 # If empty, set to 0
	cp POSCAR POSCAR.0
fi

#We repeatedly relax until OSZICAR indicates that we aren't relaxing any more.
rel=10
first_time=1
while [ ${rel} -ne 1 ]; do
	# Increment extension by 1
	exti=$((exti+1))
    if [ ${first_time} -ne 1 ]; then
        # We have not fully relaxed. Do a continuation job
        cp CONTCAR POSCAR
        sed -i "s/ISTART = .*/ISTART = 1/" INCAR
        sed -i "s/ICHARG = .*/ICHARG = 0/" INCAR
    fi
    srun vasp_std
    if [ $? -ne 0 ]; then
	echo "Error while running VASP at `date`. Aborting." | tee -a track
	exit 1
    fi
    rel=`tail -1 OSZICAR | cut -c 1-4`
    first_time=0
    echo "There has been $rel relaxation steps at `date`" | tee -a track
    # Determine the max force
    line=`grep -n TOTAL-FORCE OUTCAR | cut -d : -f1 | tail -n 1`
    # Retrieve the line of the header of the last force table
    max=`awk 'BEGIN {max = -1} NR>='"$line"'+2 && NR<'"${line}+${natom}"+2' {f = sqrt($4*$4+$5*$5+$6*$6); if (f > max) max = f} END {print max}' OUTCAR`
    # Go through each entry of that table and find the greatest magnitude of the force.
    echo "Maximum force is $max eV/Angstrom" | tee -a track
    tail -1 OSZICAR | tee -a track
	cp OUTCAR OUTCAR.$exti
	cp OSZICAR OSZICAR.$exti
	cp CONTCAR POSCAR.$exti
	cp XDATCAR XDATCAR.$exti
done

echo "Finished relaxation at `date`" | tee -a track
tail -1 OSZICAR | tee -a track

echo "Zipping files at `date`"

sed -n "s/^ *TITEL *=//p" POTCAR > POTCAR.spec
zip -mT Outputs.zip POSCAR.* OUTCAR* OSZICAR* XDATCAR.* DOSCAR EIGENVAL vasprun.xml CONTCAR
zip Outputs.zip track
zip Inputs.zip POSCAR KPOINTS INCAR
zip -mT Inputs.zip POTCAR.spec
zip --junk-paths Inputs.zip $0
rm CHG PCDAT IBZKPT REPORT XDATCAR
zip Outputs.zip JobRelax-${SLURM_JOB_ID}.out

date
