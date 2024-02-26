#!/bin/bash
# This file runs all algorithms and times them

# Algorithms used
# alg_00:= sequential baseline template
# alg_0 := optimized sequential baseline
# alg_1 := mpi version 1
# ...

cd ..

# SH parameters
# If no input is given, N = 100
N=${1:-100}
#second input defines stepsize
STEP=${2:-10}
#third input defines number of iterations
RUNS=${3:-10}
# Set number of maximal processors used for parallel implementations
PROCS=${4:-8}
# Block length
BLOCK_LENGTH=${5:-63}
# Initial size
INIT_SIZE=${6:-100}


# Set Input/Output directories
BENCHMARKING_TYPE=threads_vs_cycles
INPUT=data/sequences
OUTPUT=benchmarking/runtimes_${BENCHMARKING_TYPE}
OUTPUT_ALIGN=data/alignments_${BENCHMARKING_TYPE}

# Define algorithm directory names
ALG_13=parallelOpenMp_13
ALG_14=parallelOpenMpOpt_14

# Clean up and create folders if not existent 
mkdir -p $OUTPUT_ALIGN #Create if it does not exist
rm -rf $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN/$ALG_13
mkdir $OUTPUT_ALIGN/$ALG_14

# Clean up
mkdir -p $OUTPUT #Create if it does not exist
rm -rf $OUTPUT
mkdir $OUTPUT
mkdir $OUTPUT/$ALG_13
mkdir $OUTPUT/$ALG_14

# Make sequenceAlignment
cd build
make 
cd ..

# Generate sequence files
# cd bashScripts
# bash generateSeq.sh $N $STEP $INIT_SIZE
# cd ..

cd $OUTPUT

# Parallel openmp baseline algorithm 
rm -rf $ALG_13
mkdir $ALG_13

# Parallel openmp data type algorithm 
rm -rf $ALG_14
mkdir $ALG_14

cd ../..

for num in $(seq 1 1 $PROCS) # If we start at 1 it will work indefinitely 
do
    echo "Number Threads: $num"

    # ALG_13
    PATH=$PATH:./build/src/parallelNW/$ALG_13
    alg_13  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_13/alignments_${BENCHMARKING_TYPE}_${ALG_13}.txt     $OUTPUT/$ALG_13/runtimes_${BENCHMARKING_TYPE}_${ALG_13}_${N}.txt   $RUNS   $num 
    
    # ALG_14
    PATH=$PATH:./build/src/parallelNW/$ALG_14
    alg_14  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_14/alignments_${BENCHMARKING_TYPE}_${ALG_14}.txt     $OUTPUT/$ALG_14/runtimes_${BENCHMARKING_TYPE}_${ALG_14}_${N}.txt   $RUNS   $num 
    
    echo "Done"

done
