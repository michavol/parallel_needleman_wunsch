#!/bin/bash

cd bashScripts

# Test sequential implementation
echo "Testing sequential implementation"
sh testSeq.sh

# Test parallel implementation
echo "Testing parallel implementation"
sh testPar.sh