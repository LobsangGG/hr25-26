#!/bin/bash
#SBATCH --job-name=partdiff_jacobi
#SBATCH --output=3jacobi_test_results.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=west

# Modul laden (falls auf deinem Cluster nötig, z.B. OpenMPI)
# module load openmpi

echo "--- Jacobi MPI Testlauf ---"
echo "Datum: $(date)"

# Die Parameter für Jacobi (Methode = 2)
# Aufruf: ./partdiff <num_threads> <method: 2=Jacobi> <interlines> <func: 1 oder 2> <term: 1=prec, 2=iter> <val>

# 1. Test: Initialisierung (Funktion 1 & 2)
echo "Testing Initialization..."
mpiexec -n 3 ./partdiff.x 1 2 0 1 2 0
mpiexec -n 3 ./partdiff.x 1 2 0 2 2 0

# 2. Test: Funktion 1 mit Iterationen (1 und 1000)
echo "Testing Function 1 (Iterations)..."
mpiexec -n 3 ./partdiff.x 1 2 0 1 2 1
mpiexec -n 3 ./partdiff.x 1 2 100 1 2 1000

# 3. Test: Funktion 1 mit Präzision
echo "Testing Function 1 (Precision)..."
mpiexec -n 3 ./partdiff.x 1 2 0 1 1 1e-4
mpiexec -n 3 ./partdiff.x 1 2 10 1 1 1e-10

# 4. Test: Funktion 2 mit Iterationen
echo "Testing Function 2 (Iterations)..."
mpiexec -n 3 ./partdiff.x 1 2 0 2 2 1
mpiexec -n 3 ./partdiff.x 1 2 100 2 2 100

# 5. Test: Funktion 2 mit Präzision
echo "Testing Function 2 (Precision)..."
mpiexec -n 3 ./partdiff.x 1 2 0 2 1 1e-4
mpiexec -n 3 ./partdiff.x 1 2 10 2 1 1e-10

echo "--- Tests abgeschlossen ---"
