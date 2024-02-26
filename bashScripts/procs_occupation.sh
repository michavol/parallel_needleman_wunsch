#!/bin/bash
# This file runs all algorithms and times them

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
BENCHMARKING_TYPE=procs_occupation
INPUT=data/sequences
OUTPUT=benchmarking/runtimes_${BENCHMARKING_TYPE}
OUTPUT_ALIGN=data/alignments_${BENCHMARKING_TYPE}

# Define algorithm directory names
ALG_10=parallelBaseline_10
ALG_12=parallelShared_12

# Clean up and create folders if not existent 
mkdir -p $OUTPUT_ALIGN #Create if it does not exist
rm -rf $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN/$ALG_10
mkdir $OUTPUT_ALIGN/$ALG_12

# Clean up
mkdir -p $OUTPUT #Create if it does not exist
rm -rf $OUTPUT
mkdir $OUTPUT
mkdir $OUTPUT/$ALG_10
mkdir $OUTPUT/$ALG_12

# Make sequenceAlignment
cd build
make 
cd ..

# Generate sequence files
#cd bashScripts
#bash generateSeq.sh $N $STEP $INIT_SIZE
#cd ..

cd $OUTPUT

# Parallel baseline algorithm 
rm -rf $ALG_10
mkdir $ALG_10

# Parallel shared memory algorithm
rm -rf $ALG_12
mkdir $ALG_12

cd ../..

# Trash text for overall runtime - not needed for this type of measurment (proc_occupations)
touch "trash.txt"

# ALG_1
PATH=$PATH:./build/src/parallelNW/$ALG_10
mpirun -n $PROCS alg_10  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_10/alignments_${BENCHMARKING_TYPE}_${ALG_10}.txt     "./trash.txt"   $RUNS   $PROCS  $BLOCK_LENGTH     $OUTPUT/$ALG_10/runtimes_${BENCHMARKING_TYPE}_${ALG_10}_${N}.txt

# ALG_1
PATH=$PATH:./build/src/parallelNW/$ALG_12
mpirun -n $PROCS alg_12  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_12/alignments_${BENCHMARKING_TYPE}_${ALG_12}.txt     "./trash.txt"   $RUNS   $PROCS  $BLOCK_LENGTH     $OUTPUT/$ALG_12/runtimes_${BENCHMARKING_TYPE}_${ALG_12}_${N}.txt
rm "./trash.txt"
