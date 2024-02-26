#include <immintrin.h>
#include <smmintrin.h>

#include <nwSolver.hpp>
#include <userFunctions.hpp>


namespace NW_SEQ_OPT{
    void findOptimalAlignment(NWMatrix& nw_matrix,const seq_t& seq1, const seq_t& seq2)
    {
        //initalise ScoreMatrix
        const u_int64_t m = seq1.size() + 1;
        const u_int64_t n = seq2.size() + 1;
        const len_t init_gap_length = 0;
        
        initialiseScoreMatrix(nw_matrix, m, n, init_gap_length);

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
        u_int64_t i = 1;
        for(; diag < n; ++diag){

            #pragma omp for simd
            #pragma simd_level(10)
            //iterate through diagonal
            for(u_int64_t diag_elem = 0; diag_elem <= diag; ++diag_elem){

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
                    const len_t gap_up = nw_matrix.gapMatrix[idx_up] + 1;
                    const score_t score_up = nw_matrix.scoreMatrix[idx_up] + gapPenalty(gap_up);

                    //cell diagonal -> similarity function
                    const score_t score_diag = nw_matrix.scoreMatrix[idx_diag] + similarity(seq1[row-1],seq2[col-1]);

                    //cell left -> gap penalty
                    const len_t gap_left = nw_matrix.gapMatrix[idx_left] + 1;
                    const score_t score_left = nw_matrix.scoreMatrix[idx_left] + gapPenalty(gap_left);

                    //find max score and corresponding predecessor cell
                    //if two cells have the same score we decide on the following order on which one to choose: UP - DIAGONAL - LEFT
                    u_int8_t pred;
                    score_t score;
                    len_t gapLength;

                    score_t max_score = std::max(score_up,std::max(score_diag,score_left));

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

                // keep index up-to-date
                ++i;

            }
        }

        // total number of diagnoals in our matrix
        u_int64_t diag_tot = 2 * n - 1;

        // precompute differnce makes it easier to update in loop
        u_int64_t idx_up_diff = diag_tot - diag;

        uint32_t diag_len = diag_tot - diag;

        // Fill lower half of matrix
        for(; diag < diag_tot; ++diag){

            #pragma omp for simd
            #pragma simd_level(10)
            for(uint32_t diag_elem = 0; diag_elem < diag_len; ++diag_elem){

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
                const len_t gap_up = nw_matrix.gapMatrix[idx_up] + 1;
                const score_t score_up = nw_matrix.scoreMatrix[idx_up] + gapPenalty(gap_up);

                //cell diagonal -> similarity function
                const score_t score_diag = nw_matrix.scoreMatrix[idx_diag] + similarity(seq1[row-1],seq2[col-1]);

                //cell left -> gap penalty
                const len_t gap_left = nw_matrix.gapMatrix[idx_left] + 1;
                const score_t score_left = nw_matrix.scoreMatrix[idx_left] + gapPenalty(gap_left);

                //find max score and corresponding predecessor cell
                //if two cells have the same score we decide on the following order on which one to choose: UP - DIAGONAL - LEFT
                u_int8_t pred;
                score_t score;
                len_t gapLength;

                score_t max_score = std::max(score_up,std::max(score_diag,score_left));

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

                ++i;
            }
            idx_up_diff--;
            --diag_len;
        }
    }
} //namespace NW_SEQ_OPT