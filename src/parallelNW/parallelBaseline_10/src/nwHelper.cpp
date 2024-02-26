#include "../include/nwHelper.hpp"

// implementation of functions to perform NW algorithm

namespace NW_PAR {

    void readSequences(std::istream& inStream,
            seq_t& seq1, seq_t& seq2)
    {
        std::vector<std::string> sequences(2);  // Array to store two DNA sequences

        if (inStream.good()){   //checking whether the stream is available
            std::getline(inStream, sequences[0]);
            std::getline(inStream, sequences[1]);

        } else {
            std::cerr << "unable to read sequences\n";
            sequences[0] = "";
            sequences[1] = "";
        }

        const unsigned m = sequences[0].size();
        const unsigned n = sequences[1].size();

        // have to guarantee that first sequence is not longer than second sequence
        unsigned int first_index = 0;
        unsigned int second_index = 1;
        if (m > n) {
            second_index = 0;
            first_index = 1;
        }

        for(char letter : sequences[first_index]){
            switch(letter){
                case 'A': seq1.push_back(base_t::A); break;
                case 'T': seq1.push_back(base_t::T); break;
                case 'C': seq1.push_back(base_t::C); break;
                case 'G': seq1.push_back(base_t::G); break;
                default : std::cerr <<"invalid DNA sequence"<<std::endl; break;
            }
        }

        for(char letter : sequences[second_index]) {
            switch(letter){
                case 'A': seq2.push_back(base_t::A); break;
                case 'T': seq2.push_back(base_t::T); break;
                case 'C': seq2.push_back(base_t::C); break;
                case 'G': seq2.push_back(base_t::G); break;
                default : std::cerr<<"invalid DNA sequence"<<std::endl; break;
            }
        }

    }

    void initialiseScoreMatrix(score_mat_t& score_matrix,
            const len_t num_rows,
            const len_t num_cols,
            const len_t initial_gap_length)
    {
        // fill first column
        len_t gap_length = initial_gap_length;

        score_t score = 0;
        // calculate initial score
        for (unsigned i = 1; i < initial_gap_length; ++i) {
            score += gapPenalty(i);
        }

        for (unsigned i = 0; i < num_rows; ++i) {
            score += gapPenalty(gap_length);
            score_matrix[i][0] = (Cell){score, gap_length, pre_t::UP};
            ++gap_length;

        }
    }

    void fillRowWithGaps(score_line_t& row,
            const len_t start_index,
            const len_t num_entries,
            const len_t initial_gap_length)
    {

        len_t index = start_index;
        len_t gap_length = initial_gap_length;
        score_t score = 0;

        // calculate initial score
        for (unsigned i = 1; i < initial_gap_length; ++i) {
            score += gapPenalty(i);
        }

        //fill first entry
        if (index == 0) {
            row[index] = (Cell){score, gap_length, pre_t::UNDEF};
        } else {
            row[index] = (Cell){score, gap_length, pre_t::LEFT};
        }
        ++index;
        ++gap_length;

        //fill row
        for (unsigned i = 1; i < num_entries; ++i) {
            score += gapPenalty(gap_length);
            row[index] = (Cell){score, gap_length, pre_t::LEFT};
            ++index;
            ++gap_length;
        }


    }

    Cell computeCell(const Cell& up, const Cell& diag, const Cell& left, base_t a, base_t b) {
        const len_t gap_up = up.gapLength + 1;
        const score_t score_up = up.score + gapPenalty(gap_up);

        //cell diagonal -> similarity function
        const score_t score_diag = diag.score + similarity(a, b);

        //cell left -> gap penalty
        const len_t gap_left = left.gapLength + 1;
        const score_t score_left = left.score + gapPenalty(gap_left);

        score_t max_score = std::max(score_up,std::max(score_diag,score_left));

        //upper cell is best choice -> gapPenalty
        if(max_score == score_up){
            return {max_score, gap_up, pre_t::UP};
        }
        //diagnoal cell is best choice -> similarityFunction
        else if(max_score == score_diag){
            return {max_score, 0, pre_t::DIAG};
        }
        //left cell is best choice -> gapPenalty
        else{
            return {max_score, gap_left, pre_t::LEFT};
        }

    }

    void computeScoreBlock(const seq_t& seq1, const seq_t& seq2,
            const score_line_t& row_above,
            score_mat_t& score_matrix,
            const len_t num_rows,
            const len_t start_col,
            const len_t num_cols)
    {
        // first row
        len_t row = 0;
        base_t row_base = seq1[row];
        for (len_t col = start_col; col < start_col + num_cols; ++col) {
            score_matrix[row][col] = computeCell(row_above[col], row_above[col - 1],
                    score_matrix[row][col - 1], row_base, seq2[col - 1]);
        }
        // other rows
        for (row = 1; row < num_rows; ++row) {
            base_t row_base = seq1[row];
            for (len_t col = start_col; col < start_col + num_cols; ++col) {
                score_matrix[row][col] = computeCell(score_matrix[row - 1][col], score_matrix[row - 1][col - 1],
                        score_matrix[row][col - 1], row_base, seq2[col - 1]);
            }
        }
    }

    score_t backtrackSolution(const score_mat_t& score_matrix,
            const seq_t& seq1,
            const seq_t& seq2,
            seq_t& aligned_seq1,
            seq_t& aligned_seq2)
    {
        // Initialize dimensions
        const len_t m = score_matrix.size();
        const len_t n = score_matrix[0].size();

        //indices for backtracking
        int row = m-1;
        int col = n-1;

        //extract score
        score_t total_score = score_matrix[row][col].score;

        // backtrack until top row is passed
        while(row >= 0) {

            // Retrieve predecssor of current cell
            const uint32_t pred = score_matrix[row][col].predecessor;

            // Fill aligned sequences accordingly
            // Because of initialization row and col, our original sequence is misaligned with our score_matrix rows and cols. Hence -1.
            switch(pred) {
                case pre_t::UP: { //up
                    aligned_seq1.push_back(seq1[row]);
                    aligned_seq2.push_back(base_t::GAP);
                    --row;
                    break;
                }
                case pre_t::DIAG: { //diag
                    aligned_seq1.push_back(seq1[row]);
                    aligned_seq2.push_back(seq2[col-1]);
                    --row;
                    --col;
                    break;
                }
                case pre_t::LEFT: { //left
                    aligned_seq1.push_back(base_t::GAP);
                    aligned_seq2.push_back(seq2[col-1]);
                    --col;
                    break;
                }
                default: {
                    std::cerr << "Invalid predecessor entry!";
                    return 0;
                }
            }
        }

        //go through top row
        while (col > 0) {
            aligned_seq1.push_back(base_t::GAP);
            aligned_seq2.push_back(seq2[col-1]);
            --col;
        }

        return total_score;
    }


    void printAlignment(std::ostream& outStream,
            const score_t total_score,
            const seq_t& aligned_seq1,
            const seq_t& aligned_seq2)
    {
        len_t l = aligned_seq1.size();

        //sanity check: both sequences must have same length
        if (aligned_seq2.size() != l) {
            std::cerr << "Aligned sequences must be of same length\n";
            return;
        }

        // Print alignment, start from back since order is reversed in aligned_seq
        for(signed i = l - 1; i >= 0; i--){
            base_t base = aligned_seq1[i];
            switch(base){
                case base_t::A: outStream << 'A'; break;
                case base_t::T: outStream << 'T'; break;
                case base_t::C: outStream << 'C'; break;
                case base_t::G: outStream << 'G'; break;
                case base_t::GAP: outStream << '-'; break;
                default: std::cerr << "Invalid base at position "<< l-i << "\n";
            }
        }
        outStream << "\n";

        for(signed i = l - 1; i >= 0; i--){
            base_t base = aligned_seq2[i];
            switch(base){
                case base_t::A: outStream << 'A'; break;
                case base_t::T: outStream << 'T'; break;
                case base_t::C: outStream << 'C'; break;
                case base_t::G: outStream << 'G'; break;
                case base_t::GAP: outStream << '-'; break;
                default: std::cerr << "Invalid base at position "<< l-i << "\n";
            }
        }

        // print score
        outStream << "\n" << total_score << "\n";
    }


    void printScoreMatrix(const score_mat_t& score_matrix){

       unsigned m = score_matrix.size();
       unsigned n = score_matrix[0].size();

       for (unsigned row = 0; row < m; row++){
          for (unsigned col = 0; col < n; col++){
             Cell cell = score_matrix[row][col];
             std::cout << cell.score << ":" << cell.gapLength << ':' << cell.predecessor << '\t';
          }
          std::cout << std::endl;
       }
    }

    void printSequence(const seq_t& sequence) {
        for (base_t b : sequence) {
            switch (b) {
                case base_t::A:
                    std::cout << 'A';
                    break;
                case base_t::G:
                    std::cout << 'G';
                    break;
                case base_t::T:
                    std::cout << 'T';
                    break;
                case base_t::C:
                    std::cout << 'C';
                    break;
            }
        }
        std::cout << std::endl;
    }

} // namespace NW_PAR

