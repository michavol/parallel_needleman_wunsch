#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include <nwSolver.hpp>
#include "../../timer/timer.h"

#ifndef RUNS
#define RUNS 10
#endif

int main(int argc, const char *argv[])
{
   
   //setup input stream
   std::ifstream inStream;
   inStream.open(argv[1],std::ios::in);
   
   //setup output stream for alignments and timing
   std::ofstream out_alignments(argv[2]);
   std::ofstream out_benchmarking(argv[3],std::fstream::app);

   // Read sequences and store in seq variables
   NW_SEQ_OPT::seq_t seq1, seq2;
   NW_SEQ_OPT::readSequences(inStream,seq1,seq2);
   //tstart_ = std::chrono::high_resolution_clock::now();

   //set number of runs for timing
   int runs = atoi(argv[4]);

   // Inilialise score matrix
   const NW_SEQ_OPT::len_t m = seq1.size() + 1;
   const NW_SEQ_OPT::len_t n = seq2.size() + 1;

   // Define matrix size: m rows and n columns
   NW_SEQ_OPT::NWMatrix nw_matrix(m, n);

   for(unsigned run = 0; run < runs+4; ++run){

      myInt64 T_start = start_tsc();

      //find optimal alignment of the two sequences
      NW_SEQ_OPT::findOptimalAlignment(nw_matrix,seq1,seq2);

      myDouble64 T_stop = stop_tsc(T_start); 

      if(run > 3) {
         out_benchmarking << m-1 <<"\t" << T_stop << std::endl;
      }

   }

   //NW_SEQ_OPT::printScoreMatrix(score_matrix)
   // Backtracking to find optimal alignment
   NW_SEQ_OPT::seq_t aligned_seq1, aligned_seq2;
   NW_SEQ_OPT::score_t opt_score = NW_SEQ_OPT::backtrackSolution(nw_matrix, seq1, seq2, aligned_seq1, aligned_seq2);

   // Print alignment
   NW_SEQ_OPT::printAlignment(out_alignments, opt_score, aligned_seq1, aligned_seq2);
   
   return 0;
}