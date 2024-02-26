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
        u_int64_t firstloop_size= (nw_matrix_size + n)/2;

        // Counting in which diagonal and element in diagonal we are in order to compute indices of possible predecessors
        u_int64_t diag_counter = 1;
        u_int64_t diag = 1;

        // Indices for predecessors
        u_int64_t idx_up;
        u_int64_t idx_diag;
        u_int64_t idx_left;

        // Corresponding row and column in the score matrix
        u_int32_t row;
        u_int32_t col;

        matrix_score_t up_gap(8);
        matrix_score_t left_gap(8);

        //matrix_score_t up_penalty(8);
        //matrix_score_t diag_similarity(8);
        //matrix_score_t left_penalty(8);

        matrix_score_t up_score(8);
        matrix_score_t diag_score(8);
        matrix_score_t left_score(8);

        //matrix_score_t up_score_new(8);
        //matrix_score_t diag_score_new(8);
        //matrix_score_t left_score_new(8);
        matrix_score_t vec_max_score(8);

        u_int64_t global_idx;
        score_t max_score_j;

        // Fill upper half of matrix
        for(u_int64_t i = 1; i < firstloop_size; ++i){

            // get row and col
            row = diag_counter;
            col = diag - row;

            // get idx of up,diag and left
            idx_up = i - diag;
            idx_diag = idx_up - diag;
            idx_left = idx_up - 1;

            // Vectorize filling elements in score vector if possible
            // Diagonal has more than 8 values to be computed
            // Avoid first and last element of diagoanl
            if(diag_counter > 8 && diag_counter < diag){
                for(u_int64_t j = 0; j < 8; ++j){

                    up_gap[j] = nw_matrix.gapMatrix[idx_up + j] + 1;
                    left_gap[j] = nw_matrix.gapMatrix[idx_left + j] + 1;

                    //up_penalty[j] = gapPenalty(up_gap[j]);
                    //diag_similarity[j] = similarity(seq1[row-j-1],seq2[col+j-1]);
                    //left_penalty[j] = gapPenalty(left_gap[j]);

                    up_score[j] = nw_matrix.scoreMatrix[idx_up + j] + gapPenalty(up_gap[j]);
                    diag_score[j] = nw_matrix.scoreMatrix[idx_diag + j] + similarity(seq1[row-j-1],seq2[col+j-1]);
                    left_score[j] = nw_matrix.scoreMatrix[idx_left + j] + gapPenalty(left_gap[j]);
                }

                // mmAdd(&up_penalty.front(),&up_score.front(),up_score_new);
                // mmAdd(&diag_similarity.front(),&diag_score.front(),diag_score_new);
                // mmAdd(&left_penalty.front(),&left_score.front(),left_score_new);

                mmMax(&up_score.front(),&diag_score.front(),&left_score.front(), vec_max_score);

                global_idx = i;
                for(unsigned j = 0; j < 8; j++){
                    max_score_j = vec_max_score[j];
                    //assign evaluated new score to score matrix
                    nw_matrix.scoreMatrix[global_idx] = max_score_j;

                    //upper cell is best choice -> gapPenalty
                    if(max_score_j == up_score[j]){
                        nw_matrix.predMatrix[global_idx] = 0;
                        nw_matrix.gapMatrix[global_idx] = up_gap[j];
                    }

                    //diagnoal cell is best choice -> similarityFunction
                    else if(max_score_j == diag_score[j]){
                        nw_matrix.predMatrix[global_idx] = 1;
                        nw_matrix.gapMatrix[global_idx] = 0;
                    } 

                    //left cell is best choice -> gapPenalty
                    else{
                        nw_matrix.predMatrix[global_idx] = 2;
                        nw_matrix.gapMatrix[global_idx] = left_gap[j];
                    }

                    global_idx++;
                }

                diag_counter -= 8;
                i += 7;
                continue;
            }

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

            // end of our current diagonal -> need to go into next one
            if(diag_counter == 0){
                diag++;
                diag_counter = diag;
            }else{
                diag_counter--; // next element of diagonal is visited
            }
            
        }

        // total number of diagnoals in our matrix
        u_int64_t diag_tot = 2 * n - 1;

        // precompute differnce makes it easier to update in loop
        u_int64_t idx_up_diff = diag_tot - diag;

        // for bottom half of matrix
        diag_counter = diag_tot - diag -1;

        // Fill lower half of matrix
        for(u_int64_t i = firstloop_size; i < nw_matrix_size; ++i){

            // get row and col
            col = n - 1 - diag_counter;
            row = diag - col;

            // get idx of up,diag and left
            idx_up = i - idx_up_diff;
            idx_diag = idx_up - (idx_up_diff + 2);
            idx_left = idx_up - 1;

            // special case near diagonal
            if(diag == n) idx_diag += 1;

            // Vectorize filling elements in score vector if possible
            // Diagonal has more than 8 values to be computed
            // Because we are in the lower triangular part of the matrix the condition is looser
            if(diag_counter > 7){

                for(u_int64_t j = 0; j < 8; ++j){

                    up_gap[j] = nw_matrix.gapMatrix[idx_up + j] + 1;
                    left_gap[j] = nw_matrix.gapMatrix[idx_left + j] + 1;

                    //up_penalty[j] = gapPenalty(up_gap[j]);
                    //diag_similarity[j] = similarity(seq1[row-j-1],seq2[col+j-1]);
                    //left_penalty[j] = gapPenalty(left_gap[j]);

                    up_score[j] = nw_matrix.scoreMatrix[idx_up + j] + gapPenalty(up_gap[j]);
                    diag_score[j] = nw_matrix.scoreMatrix[idx_diag + j] + similarity(seq1[row-j-1],seq2[col+j-1]);
                    left_score[j] = nw_matrix.scoreMatrix[idx_left + j] + gapPenalty(left_gap[j]);

                }

                //mmAdd(&up_penalty.front(),&up_score.front(),up_score_new);
                //mmAdd(&diag_similarity.front(),&diag_score.front(),diag_score_new);
                //mmAdd(&left_penalty.front(),&left_score.front(),left_score_new);

                mmMax(&up_score.front(),&diag_score.front(),&left_score.front(), vec_max_score);

                u_int64_t global_idx = i;
                for(unsigned j = 0; j < 8; j++){
                    max_score_j = vec_max_score[j];
                    //assign evaluated new score to score matrix
                    nw_matrix.scoreMatrix[global_idx] = max_score_j;

                    //upper cell is best choice -> gapPenalty
                    if(max_score_j == up_score[j]){
                        nw_matrix.predMatrix[global_idx] = 0;
                        nw_matrix.gapMatrix[global_idx] = up_gap[j];
                    }

                    //diagnoal cell is best choice -> similarityFunction
                    else if(max_score_j == diag_score[j]){
                        nw_matrix.predMatrix[global_idx] = 1;
                        nw_matrix.gapMatrix[global_idx] = 0;
                    } 

                    //left cell is best choice -> gapPenalty
                    else{
                        nw_matrix.predMatrix[global_idx] = 2;
                        nw_matrix.gapMatrix[global_idx] = left_gap[j];
                    }

                    global_idx++;
                }
            }
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

            // end of our current diagonal -> need to go into next one
            if(diag_counter == 0){
                diag++;
                diag_counter = diag_tot - diag -1;
                idx_up_diff--;
            }else{
                diag_counter--; // next element of diagonal is visited
            }
        }
    }
} //namespace NW_SEQ_OPT