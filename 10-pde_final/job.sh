#!/bin/bash
#SBATCH --job-name=test_values_1_prozess
#SBATCH --partition=west
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=job_1_prozess.txt
#SBATCH --error=job_1_prozess.txt

# Slurm gibt die Anzahl Tasks
NTASKS=$SLURM_NTASKS

# initialisation
mpirun -n $NTASKS ./partdiff.x 1 1 0 1 2 0
mpirun -n $NTASKS ./partdiff.x 1 1 0 2 2 0

# Function 1
## iterations
mpirun -n $NTASKS ./partdiff.x 1 1 0 1 2 1
mpirun -n $NTASKS ./partdiff.x 1 1 0 1 2 100
mpirun -n $NTASKS ./partdiff.x 1 1 100 1 2 1
mpirun -n $NTASKS ./partdiff.x 1 1 100 1 2 100

## precision
mpirun -n $NTASKS ./partdiff.x 1 1 0 1 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 1 0 1 1 1e-10
mpirun -n $NTASKS ./partdiff.x 1 1 10 1 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 1 10 1 1 1e-10

# Function 2
## iterations
mpirun -n $NTASKS ./partdiff.x 1 1 0 2 2 1
mpirun -n $NTASKS ./partdiff.x 1 1 0 2 2 100
mpirun -n $NTASKS ./partdiff.x 1 1 100 2 2 1
mpirun -n $NTASKS ./partdiff.x 1 1 100 2 2 100

## precision
mpirun -n $NTASKS ./partdiff.x 1 1 0 2 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 1 0 2 1 1e-10
mpirun -n $NTASKS ./partdiff.x 1 1 10 2 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 1 10 2 1 1e-10

# initialisation
mpirun -n $NTASKS ./partdiff.x 1 2 0 1 2 0
mpirun -n $NTASKS ./partdiff.x 1 2 0 2 2 0

# Function 1
## iterations
mpirun -n $NTASKS ./partdiff.x 1 2 0 1 2 1
mpirun -n $NTASKS ./partdiff.x 1 2 0 1 2 100
mpirun -n $NTASKS ./partdiff.x 1 2 100 1 2 1
mpirun -n $NTASKS ./partdiff.x 1 2 100 1 2 100

## precision
mpirun -n $NTASKS ./partdiff.x 1 2 0 1 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 2 0 1 1 1e-10
mpirun -n $NTASKS ./partdiff.x 1 2 10 1 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 2 10 1 1 1e-10

# Function 2
## iterations
mpirun -n $NTASKS ./partdiff.x 1 2 0 2 2 1
mpirun -n $NTASKS ./partdiff.x 1 2 0 2 2 100
mpirun -n $NTASKS ./partdiff.x 1 2 100 2 2 1
mpirun -n $NTASKS ./partdiff.x 1 2 100 2 2 100

## precision
mpirun -n $NTASKS ./partdiff.x 1 2 0 2 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 2 0 2 1 1e-10
mpirun -n $NTASKS ./partdiff.x 1 2 10 2 1 1e-4
mpirun -n $NTASKS ./partdiff.x 1 2 10 2 1 1e-10