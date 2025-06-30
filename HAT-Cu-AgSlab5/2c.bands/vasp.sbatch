#!/bin/bash
#SBATCH --job-name=VASPBand
# Partition information
#SBATCH --account=m3349
#SBATCH --qos=regular
#SBATCH --licenses=scratch
#SBATCH --constraint=cpu
# Node information
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
# Job timing
#SBATCH --time=3:00:00
# Output
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=BField@lbl.gov
#SBATCH --output=Job-%j.out

module load cpu
module load vasp/6.4.3-cpu
module list

# Environment variables for Slingshot Update Sept. 2024
# Large jobs may break, especially during heavy IO, without this
export FI_CXI_DEFAULT_TX_SIZE=16384
export FI_CXI_RDZV_THRESHOLD=65536

export OMP_NUM_THREADS=1

echo "Doing band structure calculation at `date`." | tee -a track

srun vasp_std
if (( $? )); then
    echo "Error while running VASP. Aborting." | tee -a track
    exit 1
fi

echo "Zipping files at `date`." | tee -a track

sed -n 's/^ *TITEL *=//p' POTCAR > POTCAR.spec
rm CHG REPORT PCDAT XDATCAR CONTCAR
zip Main.zip OUTCAR
zip -mT Output.zip DOSCAR EIGENVAL OSZICAR OUTCAR PROCAR
zip Main.zip POSCAR INCAR KPOINTS
zip -mT Main.zip POTCAR.spec
zip --junk-paths Main.zip $0
gzip vasprun.xml
zip Main.zip Job-${SLURM_JOB_ID}.out

echo "Finished execution at `date`." | tee -a track

exit 0

