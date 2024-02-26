#include <nwHelper.hpp>
#include <userFunctions.hpp>
#include <chrono>

namespace NW_SEQ{
    void findOptimalAlignment(score_mat_t& score_matrix, const seq_t& seq1, const seq_t& seq2)
    {
        const len_t m = score_matrix.size();
        const len_t n = score_matrix[0].size();
        const len_t init_gap_length = 0; 

        initialiseScoreMatrix(score_matrix, m, n, init_gap_length);

        // Fill score matrix
        for(unsigned i = 1; i < m; ++i){
            for(unsigned j = 1; j < n; ++j){

                //cell up -> gap penalty
                Cell cell_up = score_matrix[i-1][j];
                const len_t gap_up = cell_up.gapLength + 1;
                const score_t score_up = cell_up.score + gapPenalty(gap_up);

                //cell diagonal -> similarity function
                Cell cell_diag = score_matrix[i-1][j-1];
                const score_t score_diag = cell_diag.score + similarity(seq1[i-1],seq2[j-1]);

                //cell left -> gap penalty
                Cell cell_left = score_matrix[i][j-1];
                const len_t gap_left = cell_left.gapLength + 1;
                const score_t score_left = cell_left.score + gapPenalty(gap_left);

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
                score_matrix[i][j] = {max_score,gapLength,pred};
            }
        }
    }
} //namespace NW_SEQ