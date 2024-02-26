#include <nwHelper.hpp>
// implementation of functions to perform NW algorithm

namespace NW_PAR {



    void stringToSequence(const std::string& s, seq_t& seq) {
        len_t size = s.size();
        seq.size = size;
        seq.data = new base_t[size];
        for (len_t i = 0; i < size; ++i) {
            base_t base;
            const char letter = s[i];
            switch(letter) {
                case 'A': base = base_t::A; break;
                case 'T': base = base_t::T; break;
                case 'C': base = base_t::C; break;
                case 'G': base = base_t::G; break;
                default : std::cerr << "invalid DNA sequence" << std::endl; break;
            }
            seq.data[i] = base;
        }

    }

    void readSequences(std::istream& inStream, seq_t& seq1, seq_t& seq2){

        std::vector<std::string> sequences(2);  // Array to store two DNA sequences

        if (inStream.good()){   //checking whether the stream is available
            std::getline(inStream, sequences[0]);
            std::getline(inStream, sequences[1]);

        } else {
            std::cerr << "unable to read sequences\n";
            sequences[0] = "";
            sequences[1] = "";
        }

        // have to guarantee that first sequence is not longer than second sequence
        if (sequences[0].size() <= sequences[1].size()) {
            stringToSequence(sequences[0], seq1);
            stringToSequence(sequences[1], seq2);
        } else {
            stringToSequence(sequences[1], seq1);
            stringToSequence(sequences[0], seq2);
        }

    }


    void setupScoreMatrix(NwMatrix& nw_matrix){
        //fill first row for first process
        if (nw_matrix.global_row_index == 0) {
            //fill gaps
            NwLine*& first_row = nw_matrix.prev_row;
            std::iota(first_row->gaps, first_row->gaps + first_row->size, 0);
            // set pres
            *(first_row->pres) = pre_t::UNDEF;
            std::fill(first_row->pres + 1, first_row->pres + first_row->size, pre_t::LEFT);
            //calculate scores for all gap lengths
            std::transform(first_row->gaps, first_row->gaps + first_row->size, first_row->scores, gapPenalty);
            std::partial_sum(first_row->scores, first_row->scores + first_row->size, first_row->scores);

        }
        //Fill first column
        score_t score = 0;
        len_t gap =  nw_matrix.global_row_index + 1;
        pre_t pre = pre_t::UP;
        len_t pos = 0;

        for (len_t i = 0; i < gap; ++i) {
            score += gapPenalty(i);
        }

        for (len_t i = 0; i < nw_matrix.rows; ++i, ++gap, pos += nw_matrix.cols) {
            score += gapPenalty(gap);
            nw_matrix.scores[pos] = score;
            nw_matrix.gaps[pos] = gap;
            nw_matrix.pres[pos] = pre;
        }
    }

    void computeBlock(NwMatrix& nw_matrix, const seq_t& seq1, const seq_t& seq2,  const len_t start, const len_t end) {

        len_t segment_size = end - start;
        // set pointers for first row
        score_t* above_score = nw_matrix.prev_row->scores;
        len_t* above_gap = nw_matrix.prev_row->gaps;
        pre_t* above_pre = nw_matrix.prev_row->pres;

        score_t* score = nw_matrix.scores;
        len_t* gap = nw_matrix.gaps;
        pre_t* pre = nw_matrix.pres;

        base_t* first_seq = seq1.data;
        base_t* second_seq = seq2.data;



        for (len_t i = 0; i < nw_matrix.rows; ++i) {
            for (len_t j = start;  j < end; ++j) {
                const len_t gap_up = above_gap[j] + 1;
                const score_t score_up = above_score[j] + gapPenalty(gap_up);

                const score_t score_diag = above_score[j-1] + similarity(first_seq[i], second_seq[j-1]);

                const len_t gap_left = gap[j - 1] + 1;
                const score_t score_left = score[j - 1] + gapPenalty(gap_left);

                const score_t max_score = std::max(score_up, std::max(score_diag,score_left));


                //upper cell is best choice -> gapPenalty
                if(max_score == score_up){
                    pre[j] = pre_t::UP;
                    gap[j] = gap_up;
                }
                //diagnoal cell is best choice -> similarityFunction
                else if(max_score == score_diag){
                    pre[j] = pre_t::DIAG;
                    gap[j] = 0;
                }
                //left cell is best choice -> gapPenalty
                else{
                    pre[j] = pre_t::LEFT;
                    gap[j] = gap_left;
                }

                score[j] = max_score;

            }
            //advance one row
            above_score = score;
            above_gap = gap;
            above_pre = pre;
            score += nw_matrix.cols;
            gap += nw_matrix.cols;
            pre += nw_matrix.cols;
        }

    }



    score_t backtrackSolution(const NwMatrix& nw_matrix,
            const seq_t& seq1, const seq_t& seq2,
            seq_t& aligned_seq1,
            seq_t& aligned_seq2){

        // Initialize dimensions

        // While we're not at the end point
        // (matrix[0] is the only cell with direction == undef)
        // Stop loop when score_matrix[0][0] has been reached.
            // Retrieve predecssor of current cell

            // Fill aligned sequences accordingly
            // Because of initialization row and col, our origninal sequence is misaligned with our score_matrix rows and cols. Hence -1.

        // Extract optimal score and return it
        return 0;
    }

    void printAlignment(std::ostream& outStream,
            const score_t total_score,
            const seq_t& aligned_seq1,
            const seq_t& aligned_seq2){

        len_t l = aligned_seq1.size;

        // Print alignment
        for(signed i = l - 1; i >= 0; i--){
            outStream << convertToChar(aligned_seq1.data[i]);

        }
        outStream << "\n";

        for(signed i = l - 1; i >= 0; i--){
            outStream << convertToChar(aligned_seq2.data[i]);
        }

        // Store Score
        outStream << "\n";
        outStream << total_score;
    }

    void printSequence(std::ostream& outStream, const seq_t& seq) {
        for (unsigned i = 0; i < seq.size; ++i) {
            outStream << convertToChar(seq.data[i]);
        }
        outStream << "\n";
    }

    void printScoreMatrix(std::ostream& outStream, const NwMatrix& nw_matrix) {
        for (unsigned i = 0; i < nw_matrix.prev_row->size; ++i) {
            outStream << nw_matrix.prev_row->scores[i] << " ";
        }
        outStream << "\n";
        score_t* it = nw_matrix.scores;
        for (unsigned i = 0; i < nw_matrix.rows; ++i) {
            for(unsigned j = 0; j < nw_matrix.cols; ++j) {
                outStream << *it++ << " ";
            }
            outStream << "\n";
        }

    }

    void printGapMatrix(std::ostream& outStream, const NwMatrix& nw_matrix) {
        for (unsigned i = 0; i < nw_matrix.prev_row->size; ++i) {
            outStream << nw_matrix.prev_row->gaps[i] << " ";
        }
        outStream << "\n";
        len_t* it = nw_matrix.gaps;
        for (unsigned i = 0; i < nw_matrix.rows; ++i) {
            for(unsigned j = 0; j < nw_matrix.cols; ++j) {
                outStream << *it++ << " ";
            }
            outStream << "\n";
        }

    }

    void printPreMatrix(std::ostream& outStream, const NwMatrix& nw_matrix) {
        for (unsigned i = 0; i < nw_matrix.prev_row->size; ++i) {
            outStream << convertToChar(nw_matrix.prev_row->pres[i]) << " ";
        }
        outStream << "\n";
        pre_t* it = nw_matrix.pres;
        for (unsigned i = 0; i < nw_matrix.rows; ++i) {
            for(unsigned j = 0; j < nw_matrix.cols; ++j) {
                outStream << convertToChar(*it++) << " ";
            }
            outStream << "\n";
        }

    }


    char convertToChar(base_t base)
    {
        switch(base){
            case base_t::A: return 'A';
            case base_t::T: return 'T';
            case base_t::C: return 'C';
            case base_t::G: return 'G';
            case base_t::GAP: return '-';
            default: return ' ';
        }
    }

    char convertToChar(pre_t pre) {
        switch(pre) {
            case pre_t::UP: return 'u';
            case pre_t::LEFT: return 'l';
            case pre_t::DIAG: return 'd';
            case pre_t::UNDEF: return '-';
            default: return ' ';
        }
    }

} // namespace NW_PAR


