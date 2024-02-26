#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include <nwSolver.hpp>
#include "../../timer/timer.h"

int main(int argc, const char *argv[])
{
   // Open sequences file
   std::ifstream inStream;
   inStream.open(argv[1],std::ios::in);
   
   // Setup output stream for alignments
   std::ofstream out_alignments(argv[2]);
   std::ofstream out_benchmarking(argv[3], std::fstream::app);

   //set number of runs for timing
   int runs = atoi(argv[4]);

   // Read sequences and store in seq variables
   NW_SEQ::seq_t seq1, seq2;
   NW_SEQ::readSequences(inStream,seq1,seq2);

   // Inilialise score matrix
   const NW_SEQ::len_t m = seq1.size() + 1;
   const NW_SEQ::len_t n = seq2.size() + 1;

   // Define matrix size: m rows and n columns
   NW_SEQ::score_mat_t score_matrix(m,NW_SEQ::score_line_t(n));

   for(unsigned run = 0; run < runs+4; ++run){

      //tstart_ = std::chrono::high_resolution_clock::now();
      myInt64 T_start = start_tsc();

      //find optimal alignment of the two sequences
      NW_SEQ::findOptimalAlignment(score_matrix, seq1, seq2);

      myDouble64 T_stop = stop_tsc(T_start); 

      if(run > 3) {
         out_benchmarking << m-1 <<"\t" << T_stop << std::endl;
      }

   }
   //tend_ = std::chrono::high_resolution_clock::now();
   //std::cout << m-1 << "\t" << std::chrono::duration<double>(tend_ - tstart_).count() << "\n";
   //printScoreMatrix(score_matrix);

   // Backtracking to find optimal alignment
   NW_SEQ::seq_t aligned_seq1, aligned_seq2;
   NW_SEQ::score_t opt_score = NW_SEQ::backtrackSolution(score_matrix, seq1, seq2, aligned_seq1, aligned_seq2);

   // Print alignment
   NW_SEQ::printAlignment(out_alignments, opt_score, aligned_seq1, aligned_seq2);

   return 0;
}