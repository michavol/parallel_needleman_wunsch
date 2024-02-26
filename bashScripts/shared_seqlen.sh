#!/bin/bash
# This file runs the specified algorithms and times them

cd ..

# SH Arguments
#first input defines maximal size of sequences
N=${1:-120000}
#second input defines stepsize
STEP=${2:-10000}
#third input defines number of iterations
RUNS=${3:-20}
# Set number of processors used for parallel implementations
PROCS=${4:-48}
# Block length for parallel implementation
BLOCK_LENGTH=${5:-40}
# Initial size
INIT_SIZE=${6:-10000}

# Set Input/Output directories
INPUT=data/sequences
OUTPUT=benchmarking/runtimes
OUTPUT_ALIGN=data/scores

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
make alg_10
make alg_12
cd ..

# Generate sequence files
# cd bashScripts
# bash generateSeq.sh $N $STEP $INIT_SIZE
# cd ..

for n in $(seq $INIT_SIZE $STEP $N) # If we start at 1 it will work indefinitely 
do
    echo "Sequence Length: $n"

    # ALG_10
    LOC=./build/src/parallelNW/$ALG_10
    mpirun -n $PROCS $LOC/alg_10  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_10/scores_${ALG_10}.txt     $OUTPUT/$ALG_10/runtime_${ALG_10}_${N}.txt      $RUNS   $PROCS     $BLOCK_LENGTH
    
    # ALG_12
    LOC=./build/src/parallelNW/$ALG_12
    mpirun -n $PROCS $LOC/alg_12  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_12/scores_${ALG_12}.txt     $OUTPUT/$ALG_12/runtime_${ALG_12}_${N}.txt      $RUNS   $PROCS     $BLOCK_LENGTH

    echo "Done"
done




