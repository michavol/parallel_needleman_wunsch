#!/bin/bash

cd ..

# Recompile source code
cd build
make 
cd ..

# Define in and out streams
INPUT_FILES=testing/data/NWData
OUTPUT_ALIGN=testing/data
HELPER_SCRIPTS=testing/helperScripts
RUNS=

# Test all sequences found in INPUT_FILES
for FILE in $INPUT_FILES/*
do
  echo "Testing" $FILE
  python3 $HELPER_SCRIPTS/generateAlignment.py $FILE

  rm $OUTPUT_ALIGN/OutSeq_0/out.txt

  PATH=$PATH:./build/parallelBaseline_10 # Acess to alg_1
  #$INPUT/seq${n}.txt   $OUTPUT_ALIGN/$ALG_10/alignments_${ALG_10}.txt     $OUTPUT/$ALG_10/runtime_${ALG_10}_${N}.txt      $RUNS
  mpirun -n 2 alg_10 $FILE $OUTPUT_ALIGN/OutSeq_0/out.txt

  python3 $HELPER_SCRIPTS/formatAlignment.py 
  python3 $HELPER_SCRIPTS/compareFiles.py $OUTPUT_ALIGN/OutBP/alignmentOut.txt $OUTPUT_ALIGN/OutSeq_0/out_format.txt
done


