#!/bin/bash
# This file runs all algorithms and times them

# Algorithms used
# alg_00:= sequential baseline template
# alg_01 := optimized sequential baseline
# alg_02 := diagonal sequential baseline
# alg_10 := mpi version 1
# ...

cd ..

# SH Arguments
# generate data 
# If no input is given, N = 100
N=${1:-100}
#second input defines stepsize
STEP=${2:-10}
#third input defines number of iterations
RUNS=${3:-10}
# Set number of processors used for parallel implementations
PROCS=${4:-2}
# Block length for parallel implementation
BLOCK_LENGTH=${5:-63}
# Initial size
INIT_SIZE=${6:-100}

# Set Input/Output directories
INPUT=data/sequences
OUTPUT=benchmarking/runtimes
OUTPUT_ALIGN=data/alignments

# Define algorithm directory names
ALG_00=sequentialBaseline_00
ALG_01=sequentialOpt_01
ALG_02=sequentialDiag_02
ALG_03=sequentialDiagOpt_03
ALG_04=sequentialDiagSimd_04
ALG_05=sequentialDiagOptNested_05
ALG_06=sequentialDiagSimdPragma_06
ALG_10=parallelBaseline_10
ALG_11=parallelBaselineOpt_11
ALG_12=parallelShared_12
ALG_13=parallelOpenMp_13
ALG_14=parallelOpenMpOpt_14

# Clean up and create folders if not existent 
mkdir -p $OUTPUT_ALIGN #Create if it does not exist
rm -rf $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN
mkdir $OUTPUT_ALIGN/$ALG_00
mkdir $OUTPUT_ALIGN/$ALG_01
mkdir $OUTPUT_ALIGN/$ALG_02
mkdir $OUTPUT_ALIGN/$ALG_03
mkdir $OUTPUT_ALIGN/$ALG_04
mkdir $OUTPUT_ALIGN/$ALG_05
mkdir $OUTPUT_ALIGN/$ALG_06
mkdir $OUTPUT_ALIGN/$ALG_10
mkdir $OUTPUT_ALIGN/$ALG_11
mkdir $OUTPUT_ALIGN/$ALG_12
mkdir $OUTPUT_ALIGN/$ALG_13
mkdir $OUTPUT_ALIGN/$ALG_14

# Clean up
mkdir -p $OUTPUT #Create if it does not exist
rm -rf $OUTPUT
mkdir $OUTPUT
mkdir $OUTPUT/$ALG_00
mkdir $OUTPUT/$ALG_01
mkdir $OUTPUT/$ALG_02
mkdir $OUTPUT/$ALG_03
mkdir $OUTPUT/$ALG_04
mkdir $OUTPUT/$ALG_05
mkdir $OUTPUT/$ALG_06
mkdir $OUTPUT/$ALG_10
mkdir $OUTPUT/$ALG_11
mkdir $OUTPUT/$ALG_12
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

#sequentialBaseline_00 algorithm 
rm -rf $ALG_00
mkdir $ALG_00

#sequentialOpt_01 algorithm 
rm -rf $ALG_01
mkdir $ALG_01

#sequentialDiag_02 algorithm 
rm -rf $ALG_02
mkdir $ALG_02

#sequentialDiagOpt_03 algorithm 
rm -rf $ALG_03
mkdir $ALG_03

#sequentialDiagSimd_04 algorithm 
rm -rf $ALG_04
mkdir $ALG_04

#sequentialDiagOptNested_05 algorithm 
rm -rf $ALG_05
mkdir $ALG_05

#sequentialDiagSimdPragma_06 algorithm 
rm -rf $ALG_06
mkdir $ALG_06

#parallelBaseline_10 algorithm
rm -rf $ALG_10
mkdir $ALG_10

#parallelBaseline_11 algorithm
rm -rf $ALG_11
mkdir $ALG_11

#parallelShared_12 algorithm
rm -rf $ALG_12
mkdir $ALG_12

#parallelOpenMp_13 algorithm
rm -rf $ALG_13
mkdir $ALG_13

#parallelOpenMpOpt_14 algorithm
rm -rf $ALG_14
mkdir $ALG_14
cd ../..

for n in $(seq $INIT_SIZE $STEP $N) # If we start at 1 it will work indefinitely 
do
    echo "Sequence Length: $n"

    # ALG_00
    #echo "ALG_00"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_00
    #alg_00                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_00/alignments_${ALG_00}.txt     $OUTPUT/$ALG_00/runtime_${ALG_00}_${N}.txt      $RUNS   

    # ALG_01
    #echo "ALG_01"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_01
    #alg_01                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_01/alignments_${ALG_01}.txt     $OUTPUT/$ALG_01/runtime_${ALG_01}_${N}.txt      $RUNS

    # ALG_02
    #echo "ALG_02"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_02
    #alg_02                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_02/alignments_${ALG_02}.txt   $OUTPUT/$ALG_02/runtime_${ALG_02}_${N}.txt    $RUNS

    # ALG_03
    #echo "ALG_03"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_03
    #alg_03                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_03/alignments_${ALG_03}.txt   $OUTPUT/$ALG_03/runtime_${ALG_03}_${N}.txt    $RUNS

    # ALG_04
    #echo "ALG_04"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_04
    #alg_04                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_04/alignments_${ALG_04}.txt   $OUTPUT/$ALG_04/runtime_${ALG_04}_${N}.txt    $RUNS

    # ALG_05
    #echo "ALG_05"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_05
    #alg_05                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_05/alignments_${ALG_05}.txt   $OUTPUT/$ALG_05/runtime_${ALG_05}_${N}.txt    $RUNS

    # ALG_06
    #echo "ALG_06"
    #PATH=$PATH:./build/src/sequentialNW/$ALG_06
    #alg_06                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_06/alignments_${ALG_06}.txt   $OUTPUT/$ALG_06/runtime_${ALG_06}_${N}.txt    $RUNS

    # ALG_10
    #echo "ALG_10"
    #PATH=$PATH:./build/src/parallelNW/$ALG_10
    #mpirun -n $PROCS alg_10  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_10/alignments_${ALG_10}.txt     $OUTPUT/$ALG_10/runtime_${ALG_10}_${N}.txt      $RUNS   $PROCS     $BLOCK_LENGTH
    
     # ALG_11
    echo "ALG_11"
    PATH=$PATH:./build/src/parallelNW/$ALG_11
    mpirun -n $PROCS alg_11  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_11/alignments_${ALG_11}.txt     $OUTPUT/$ALG_11/runtime_${ALG_11}_${N}.txt      $RUNS   $PROCS     $BLOCK_LENGTH
    
    # ALG_12
    #echo "ALG_12"
    #PATH=$PATH:./build/src/parallelNW/$ALG_12
    #mpirun -n $PROCS alg_12  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_12/alignments_${ALG_12}.txt     $OUTPUT/$ALG_12/runtime_${ALG_12}_${N}.txt      $RUNS   $PROCS     $BLOCK_LENGTH

    # ALG_13
    #echo "ALG_13"
    #PATH=$PATH:./build/src/parallelNW/$ALG_13
    #alg_13                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_13/alignments_${ALG_13}.txt   $OUTPUT/$ALG_13/runtime_${ALG_13}_${N}.txt    $RUNS   $PROCS

    # ALG_14
    #echo "ALG_14"
    #PATH=$PATH:./build/src/parallelNW/$ALG_14
    #alg_14                  $INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_14/alignments_${ALG_14}.txt   $OUTPUT/$ALG_14/runtime_${ALG_14}_${N}.txt    $RUNS   $PROCS

    echo "Done"
done




