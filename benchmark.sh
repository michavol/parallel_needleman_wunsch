#!/bin/bash
# SH parameters
# If no input is given, N = 100
N=${1:-10}
#second input defines stepsize
STEP=${2:-10}
#third input defines number of iterations
RUNS=${3:-1}
# Set number of maximal processors used for parallel implementations
PROCS=${4:-8}
# Maximal block length
BLOCK_LENGTH=${5:-40}
# Initial size
INIT_SIZE=${6:-10}

cd bashScripts

# Run all algorithms and measure time vs sequnece length
echo "\n-----------------------------------------------\n"
echo "Benchmarking time vs. sequence length"
echo "\n-----------------------------------------------\n"
sh seqlength_vs_cycles.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE

# Run all parallel algorithms and measure time vs number of processor
#echo "\n-----------------------------------------------\n"
#echo "Benchmarking time vs. sequence length"
#echo "\n-----------------------------------------------\n"
#sh seqlength_vs_cycles.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE

# Run all parallel algorithms and measure time vs number of processor
#echo "\n-----------------------------------------------\n"
#echo "Benchmarking time vs. number processors"
#echo "\n-----------------------------------------------\n"
#sh procs_vs_cycles.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE

# Run OpenMP algorithm and measure time vs number of threads
#echo "\n-----------------------------------------------\n"
#echo "Benchmarking time vs. number threads"
#echo "\n-----------------------------------------------\n"
#sh threads_vs_cycles.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE

# Run all parallel algorithms and measure process occupation
#echo "\n-----------------------------------------------\n"
#echo "Benchmarking process occupation"
#echo "\n-----------------------------------------------\n"
#sh procs_occupation.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE

# Run all parallel algorithms and measure block_length vs cycles
#echo "\n-----------------------------------------------\n"
#echo "Benchmarking block_length vs. cycles"
#echo "\n-----------------------------------------------\n"
#sh blockLength_vs_cycles.sh $N $STEP $RUNS $PROCS $BLOCK_LENGTH $INIT_SIZE
