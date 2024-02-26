#include <nwHelper.hpp>
#include <userFunctions.hpp>
// implementation of functions to perform NW algorithm

namespace NW_SEQ {
   void readSequences(std::ifstream& inStream, seq_t& seq1, seq_t& seq2){

      std::vector<std::string> sequences;  // Array to store two DNA sequences
      
      if (inStream.is_open()){   //checking whether the file is open
         std::string seq;//store one sequence
         while(std::getline(inStream, seq)){  //read data from file object and put it into string.
            sequences.push_back(seq);  
            //std::cout << seq << "\n";   //print the data of the string
         }

      inStream.close();
      }

      unsigned m = sequences[0].size();
      unsigned n = sequences[1].size();

      for(unsigned i = 0; i < m; ++i){
         unsigned char letter = sequences[0][i];
         switch(letter){
            case 'A': seq1.push_back(base_t::A); break;
            case 'T': seq1.push_back(base_t::T); break;
            case 'C': seq1.push_back(base_t::C); break;
            case 'G': seq1.push_back(base_t::G); break;
            default : std::cerr<<"invalid DNA sequence"<<std::endl; break;
         }
      }

      for(unsigned j = 0; j < n; ++j){
         unsigned char letter = sequences[1][j];
         switch(letter){
            case 'A': seq2.push_back(base_t::A); break;
            case 'T': seq2.push_back(base_t::T); break;
            case 'C': seq2.push_back(base_t::C); break;
            case 'G': seq2.push_back(base_t::G); break;
            default : std::cerr<<"invalid DNA sequence"<<std::endl; break;
         }
      }
   }

   void initialiseScoreMatrix(score_mat_t& score_matrix, const len_t num_rows, const len_t num_cols, const len_t init_gap_length){
      //Fill first Cell
      score_matrix[0][0] = {0,init_gap_length,3};
      
      // Fill first row
      for(unsigned i = 1; i < num_cols; ++i){
         len_t gap_length = score_matrix[0][i-1].gapLength + 1;
         score_t score = score_matrix[0][i-1].score + gapPenalty(gap_length);
         score_matrix[0][i] = {score,gap_length,2};
      }

      //Fill first column
      for(unsigned j = 1; j < num_rows; ++j){
         len_t gap_length = score_matrix[j-1][0].gapLength + 1;
         score_t score = score_matrix[j-1][0].score + gapPenalty(gap_length);
         score_matrix[j][0] = {score,gap_length,0};
      }
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

   score_t backtrackSolution(const score_mat_t& score_matrix,
                              const seq_t& seq1, const seq_t& seq2,
                              seq_t& aligned_seq1,
                              seq_t& aligned_seq2){
      
      // Initialize dimensions
      unsigned m = score_matrix.size();
      unsigned n = score_matrix[0].size();
      unsigned row = m-1;
      unsigned col = n-1;

      // While we're not at the end point
      // (matrix[0][0] is the only cell with direction == 3 (undefined))
      u_int8_t pred;
      while(pred != 3){ // Stop loop when score_matrix[0][0] has been reached.
         // Retrieve predecssor of current cell
         pred = score_matrix[row][col].predecessor;
         
         // Fill aligned sequences accordingly
         // Because of initialization row and col, our origninal sequence is misaligned with our score_matrix rows and cols. Hence -1.
         switch(pred){
            case 0: { //up 
               aligned_seq1.push_back(seq1[row-1]); 
               aligned_seq2.push_back(base_t::GAP);
               row--; 
               break;
            }
            case 1: { //diag
               aligned_seq1.push_back(seq1[row-1]);
               aligned_seq2.push_back(seq2[col-1]);
               row--; col--; 
               break;
            }
            case 2: { //left
               aligned_seq1.push_back(base_t::GAP);
               aligned_seq2.push_back(seq2[col-1]);
               col--;
               break;
            }
            case 3: break; //that it does not print error for last iteration of while loop

            default: {
               std::cerr << "Invalid predecessor entry!";
            }
         }
      }
      
      // Extract optimal score and return it
      score_t opt_score = score_matrix[m-1][n-1].score;
      return opt_score;
   }

   void printAlignment(std::ofstream& outStream,
                        const score_t total_score,
                        const seq_t& aligned_seq1,
                        const seq_t& aligned_seq2){

      len_t l = aligned_seq1.size();

      // Print alignment 
      for(signed i = l - 1; i >= 0; i--){
         base_t base = aligned_seq1[i];
         switch(base){
            case base_t::A: outStream << 'A'; break;
            case base_t::T: outStream << 'T'; break;
            case base_t::C: outStream << 'C'; break;
            case base_t::G: outStream << 'G'; break;
            case base_t::GAP: outStream << '-'; break;
            default: std::cerr << "Invalid base!";
         }
      }
      outStream << std::endl;

      for(signed i = l - 1; i >= 0; i--){
         base_t base = aligned_seq2[i];
         switch(base){
            case base_t::A: outStream << 'A'; break;
            case base_t::T: outStream << 'T'; break;
            case base_t::C: outStream << 'C'; break;
            case base_t::G: outStream << 'G'; break;
            case base_t::GAP: outStream << '-'; break;
            default: std::cerr << "Invalid base!";
         }
      }

      // Store Score
      outStream << std::endl;
      outStream << total_score;
   }
    
} // namespace NW_SEQ


