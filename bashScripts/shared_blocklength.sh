#!/bin/bash
# This file runs all algorithms and times them

# Algorithms used
# alg_00:= sequential baseline template
# alg_0 := optimized sequential baseline
# alg_1 := mpi version 1
# ...

cd ..

# SH parameters
# If no input is given, N = 10000
N=${1:-10000}
#second input defines stepsize
STEP=${2:-10000}
#third input defines number of iterations
RUNS=${3:-20}
# Set number of maximal processors used for parallel implementations
PROCS=${4:-36}
# Maximal block length
BLOCK_LENGTH=${5:-100}
# Initial size of block length
INIT_SIZE=${6:-10}


# Set Input/Output directories
BENCHMARKING_TYPE=shared_blocklength
INPUT=data/sequences
OUTPUT=benchmarking/runtimes_${BENCHMARKING_TYPE}
OUTPUT_ALIGN=data/scores_${BENCHMARKING_TYPE}

# Define algorithm directory names
#ALG_10=parallelBaseline_10
ALG_12=parallelShared_12

# Clean up and create folders if not existent 
mkdir -p $OUTPUT_ALIGN #Create if it does not exist
rm -rf $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN
#mkdir $OUTPUT_ALIGN/$ALG_10
mkdir $OUTPUT_ALIGN/$ALG_12

# Clean up
mkdir -p $OUTPUT #Create if it does not exist
rm -rf $OUTPUT
mkdir $OUTPUT
#mkdir $OUTPUT/$ALG_10
mkdir $OUTPUT/$ALG_12

# Make sequenceAlignment
cd build
make alg_12
cd ..

# # Generate sequence files
# cd bashScripts
# bash generateSeq.sh $N $STEP $INIT_SIZE
# cd ..

cd $OUTPUT

# Parallel baseline algorithm 
#rm -rf $ALG_10
#mkdir $ALG_10

# Parallel shared memory algorithm
#rm -rf $ALG_12
#mkdir $ALG_12

cd ../..

#PROC_NUM="36"
#STEP_SIZE=1
#for num in $PROC_NUM
#do

echo "start computing"

for block_length in $(seq $INIT_SIZE $STEP $BLOCK_LENGTH)
do

        echo "Block Length: $block_length"
        echo "Processor Number: $PROCS"

        # ALG_10
        #PATH=$PATH:./build/src/parallelNW/$ALG_10
        #mpirun -n $num alg_10  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_10/alignments_${BENCHMARKING_TYPE}_${ALG_10}.txt     $OUTPUT/$ALG_10/runtimes_${BENCHMARKING_TYPE}_${ALG_10}_${N}.txt   $RUNS   $num     $block_length 
                
        # ALG_12
        LOC=./build/src/parallelNW/$ALG_12
        mpirun -n $PROCS $LOC/alg_12  $INPUT/seq${N}.txt   $OUTPUT_ALIGN/$ALG_12/scores_${BENCHMARKING_TYPE}_${ALG_12}.txt     $OUTPUT/$ALG_12/runtimes_${BENCHMARKING_TYPE}_${ALG_12}_${N}.txt   $RUNS   $PROCS     $block_length 
                
        echo "Done"

done
#done
  
