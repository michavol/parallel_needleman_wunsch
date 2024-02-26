#!/bin/bash

cd ..

cd build
make 
cd ..

# Define in and out streams
INPUT_FILES=testing/data/NWData
OUTPUT_ALIGN=testing/data
HELPER_SCRIPTS=testing/helperScripts

count = 0

for FILE in $INPUT_FILES/*
do
  counter=$((counter+1))
  echo "Testing" $FILE
  python3 $HELPER_SCRIPTS/generateAlignment.py  $FILE

  rm $OUTPUT_ALIGN/OutSeq_0/out.txt

  PATH=$PATH:./build/sequentialBaseline_00
  alg_00 $FILE $OUTPUT_ALIGN/OutSeq_0/out.txt >> $OUTPUT_ALIGN/OutSeq_0/timing.txt

  python3 $HELPER_SCRIPTS/formatAlignment.py 

  python3 $HELPER_SCRIPTS/compareFiles.py $OUTPUT_ALIGN/OutBP/alignmentOut.txt $OUTPUT_ALIGN/OutSeq_0/out_format.txt
done

echo "Passed" $counter "tests of" 
ls $OUTPUT_ALIGN/NWData | wc -l

#rm $OUTPUT_ALIGN/OutSeq_0/timing.txt