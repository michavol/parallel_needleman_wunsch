#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdint>
#include <omp.h>
#include <string>
#include <immintrin.h>
#include <smmintrin.h>

#include <userFunctions.hpp>
#include "../../../timer/timer.h"

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
   std::ofstream out_benchmarking(argv[3], std::fstream::app);
   std::ofstream out_benchmarking_threads(argv[6], std::fstream::app);
   //Check for timing type
   bool time_occupation = (argc == 7);

   // Read sequences and store in seq variables
   NW_SEQ_OPT::seq_t seq1, seq2;
   NW_SEQ_OPT::readSequences(inStream,seq1,seq2);
   //tstart_ = std::chrono::high_resolution_clock::now();

   //set number of runs for timing
   int runs = atoi(argv[4]);

   //set number of n_threads
   const int n_threads = atoi(argv[5]);

   // Inilialise score matrix
   const u_int64_t m = seq1.size() + 1;
   const u_int64_t n = seq2.size() + 1;

   // Define matrix size: m rows and n columns
   NW_SEQ_OPT::NWMatrix nw_matrix(m, n);

   for(unsigned run = 0; run < runs+5; ++run){

      myInt64 T_start = start_tsc();

      //find optimal alignment of the two sequences
      //initalise ScoreMatrix
        const u_int64_t m = seq1.size() + 1;
        const u_int64_t n = seq2.size() + 1;
        const NW_SEQ_OPT::len_t init_gap_length = 0;
        
        NW_SEQ_OPT::initialiseScoreMatrix(nw_matrix, m, n, init_gap_length);

        // Fill score matrix
        u_int64_t nw_matrix_size = m*n;

        // Counting in which diagonal and element in diagonal we are in order to compute indices of possible predecessors
        u_int64_t diag = 1;

        // Indices for predecessors
        u_int64_t idx_up;
        u_int64_t idx_diag;
        u_int64_t idx_left;

        // Corresponding row and column in the score matrix
        u_int32_t row;
        u_int32_t col;

        // Fill upper half of matrix
        u_int64_t indx = 1;
        for(; diag < n; ++diag){

            #pragma omp parallel for shared(nw_matrix) num_threads(n_threads)
            //iterate through diagonal
            for(u_int64_t diag_elem = 0; diag_elem <= diag; ++diag_elem){
                  //get each thread a private timing variable
                  // if(time_occupation) myInt64 T_start = start_tsc();

                  u_int64_t i = indx + diag_elem;

                  // get row and col
                  row = diag - diag_elem;
                  col = diag - row;

                  // get idx of up,diag and left
                  idx_up = i - diag;
                  idx_diag = idx_up - diag;
                  idx_left = idx_up - 1;

                  
                  // skip first row and column (already initialized)
                  if(row > 0 && col > 0){
                     //cell up -> gap penalty
                     const NW_SEQ_OPT::len_t gap_up = nw_matrix.gapMatrix[idx_up] + 1;
                     const NW_SEQ_OPT::score_t score_up = nw_matrix.scoreMatrix[idx_up] + NW_SEQ_OPT::gapPenalty(gap_up);

                     //cell diagonal -> similarity function
                     const NW_SEQ_OPT::score_t score_diag = nw_matrix.scoreMatrix[idx_diag] + NW_SEQ_OPT::similarity(seq1[row-1],seq2[col-1]);

                     //cell left -> gap penalty
                     const NW_SEQ_OPT::len_t gap_left = nw_matrix.gapMatrix[idx_left] + 1;
                     const NW_SEQ_OPT::score_t score_left = nw_matrix.scoreMatrix[idx_left] + NW_SEQ_OPT::gapPenalty(gap_left);

                     //find max score and corresponding predecessor cell
                     //if two cells have the same score we decide on the following order on which one to choose: UP - DIAGONAL - LEFT
                     u_int8_t pred;
                     NW_SEQ_OPT::score_t score;
                     NW_SEQ_OPT::len_t gapLength;

                     NW_SEQ_OPT::score_t max_score = std::max(score_up,std::max(score_diag,score_left));

                     //upper cell is best choice -> gapPenalty
                     if(max_score == score_up){
                        pred = 0;
                        gapLength = gap_up;
                     }
                     //diagnoal cell is best choice -> similarityFunction
                     else if(max_score == score_diag){
                        pred = 1;
                        gapLength = 0;
                     } 
                     //left cell is best choice -> gapPenalty
                     else{
                        pred = 2;
                        gapLength = gap_left;
                     }

                     //current cell
                     nw_matrix.scoreMatrix[i] = max_score;
                     nw_matrix.gapMatrix[i] = gapLength;
                     nw_matrix.predMatrix[i] = pred;
                  }

               // if(time_occupation){
               //    #pragma omp critical 
               //    {
               //       myDouble64 T_stop = stop_tsc(T_start); 
               //       if(run != 0) {
               //          int thread_id = omp_get_thread_num();
               //          out_benchmarking_threads << m-1 << "\t" << T_stop << "\t" << thread_id << std::endl;
               //       } 
               //    } 
               // } 
            }
            
            indx += diag + 1;
        }

        // total number of diagnoals in our matrix
        u_int64_t diag_tot = 2 * n - 1;

        // precompute differnce makes it easier to update in loop
        u_int64_t idx_up_diff = diag_tot - diag;

        uint32_t diag_len = diag_tot - diag;

        // Fill lower half of matrix
        for(; diag < diag_tot; ++diag){

            #pragma omp parallel for shared(nw_matrix) num_threads(n_threads)
            for(uint32_t diag_elem = 0; diag_elem < diag_len; ++diag_elem){
               // if(time_occupation) myInt64 T_start = start_tsc();

               u_int64_t i = indx + diag_elem;

               // get row and col
               col = n - (diag_len - diag_elem);
               row = diag - col;

               // get idx of up,diag and left
               idx_up = i - idx_up_diff;
               idx_diag = idx_up - (idx_up_diff + 2);
               idx_left = idx_up - 1;

               // special case near diagonal
               if(diag == n) idx_diag += 1;

               //cell up -> gap penalty
               const NW_SEQ_OPT::len_t gap_up = nw_matrix.gapMatrix[idx_up] + 1;
               const NW_SEQ_OPT::score_t score_up = nw_matrix.scoreMatrix[idx_up] + NW_SEQ_OPT::gapPenalty(gap_up);

               //cell diagonal -> similarity function
               const NW_SEQ_OPT::score_t score_diag = nw_matrix.scoreMatrix[idx_diag] + NW_SEQ_OPT::similarity(seq1[row-1],seq2[col-1]);

               //cell left -> gap penalty
               const NW_SEQ_OPT::len_t gap_left = nw_matrix.gapMatrix[idx_left] + 1;
               const NW_SEQ_OPT::score_t score_left = nw_matrix.scoreMatrix[idx_left] + NW_SEQ_OPT::gapPenalty(gap_left);

               //find max score and corresponding predecessor cell
               //if two cells have the same score we decide on the following order on which one to choose: UP - DIAGONAL - LEFT
               u_int8_t pred;
               NW_SEQ_OPT::score_t score;
               NW_SEQ_OPT::len_t gapLength;

               NW_SEQ_OPT::score_t max_score = std::max(score_up,std::max(score_diag,score_left));

               //upper cell is best choice -> gapPenalty
               if(max_score == score_up){
                  pred = 0;
                  gapLength = gap_up;
               }
               //diagnoal cell is best choice -> similarityFunction
               else if(max_score == score_diag){
                  pred = 1;
                  gapLength = 0;
               } 
               //left cell is best choice -> gapPenalty
               else{
                  pred = 2;
                  gapLength = gap_left;
               }

               //current cell
               nw_matrix.scoreMatrix[i] = max_score;
               nw_matrix.gapMatrix[i] = gapLength;
               nw_matrix.predMatrix[i] = pred;

               // if(time_occupation){
               //    #pragma omp critical 
               //    {
               //       myDouble64 T_stop = stop_tsc(T_start); 
               //       if(run != 0) {
               //          int thread_id = omp_get_thread_num();
               //          out_benchmarking_threads << m-1 << "\t" << T_stop << "\t" << thread_id << std::endl;
               //       } 
               //    } 
               // } 
            }

            indx += diag_len;
            idx_up_diff--;
            --diag_len;
        }

      myDouble64 T_stop = stop_tsc(T_start); 

      if(run > 4) {
         out_benchmarking << m-1 <<"\t" << T_stop << "\t" << n_threads << std::endl;
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