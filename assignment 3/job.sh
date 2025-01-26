#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=genoa
. /etc/bashrc
. /etc/profile.d/lmod.sh
module load 2022
module load OpenMPI/4.1.4-GCC-11.3.0
APP=$1
if [ -z "$1" ]; then
  echo "Error: No application provided. Usage: $0 <application_path>"
  exit 1
fi

#ARGS=""
# Optional OpenMPI options (commented out here); the example disables `usnic` transport
#OMPI_OPTS="--mca btl ^usnic"
MPI_RUN=mpirun
echo $MPI_RUN
echo "Running application: $APP"
$MPI_RUN $OMPI_OPTS $APP $ARGS
echo "DONE"