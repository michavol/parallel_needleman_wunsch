#include <nwHelper.hpp>
#include <userFunctions.hpp>
// implementation of functions to perform NW algorithm

namespace NW_SEQ_OPT {
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

   void initialiseScoreMatrix(NWMatrix& nw_matrix, const len_t num_rows, const len_t num_cols, const len_t init_gap_length){
      //Fill first Cell
      len_t indx_row = 0;
      len_t indx_col = 0;
      nw_matrix.scoreMatrix[0] = 0;
      nw_matrix.gapMatrix[0] = init_gap_length;
      nw_matrix.predMatrix[0] = 3;

      // Fill first column
      for(unsigned row = 1; row < num_rows; ++row){
         // gap length is equal to row
         score_t score = nw_matrix.scoreMatrix[indx_row] + gapPenalty(row);

         // quick way of converting row to index
         indx_row += row;
         nw_matrix.scoreMatrix[indx_row] = score;
         nw_matrix.gapMatrix[indx_row] = row;
         nw_matrix.predMatrix[indx_row] = 0;
      }
      
      //Fill first row
      for(unsigned col = 1; col < num_cols; ++col){
         // gap length is equal to col
         score_t score = nw_matrix.scoreMatrix[indx_col] + gapPenalty(col);

         // quick way of converting col to index
         indx_col += col + 1;
         nw_matrix.scoreMatrix[indx_col] = score;
         nw_matrix.gapMatrix[indx_col] = col;
         nw_matrix.predMatrix[indx_col] = 2;
      }
   }

   void printScoreMatrix(const NWMatrix& nw_matrix){

      // unsigned m = nw_matrix.m;
      // unsigned n = nw_matrix.n;

      // for (unsigned row = 0; row < m; row++){
      //    for (unsigned col = 0; col < n; col++){
      //       len_t indx = nw_matrix.indxMatrix[row][col];
      //       std::cout << nw_matrix.scoreMatrix[indx] << ":" << nw_matrix.gapMatrix[indx]  << ':' << nw_matrix.predMatrix[indx]  << '\t';
      //    }
      //    std::cout << std::endl;
      // }
   }

   score_t backtrackSolution(const NWMatrix& nw_matrix,
                              const seq_t& seq1, const seq_t& seq2,
                              seq_t& aligned_seq1,
                              seq_t& aligned_seq2){
      
      // '// Initialize dimensions
      unsigned m = nw_matrix.m;
      unsigned n = nw_matrix.n;
      // unsigned row = m-1;
      // unsigned col = n-1;

      // // While we're not at the end point
      // // (matrix[0][0] is the only cell with direction == 3 (undefined))
      // u_int8_t pred;
      // while(pred != 3){ // Stop loop when score_matrix[0][0] has been reached.
      //    // Retrieve predecssor of current cell
      //    len_t idx = nw_matrix.indxMatrix[row][col];
      //    pred = nw_matrix.predMatrix[idx];
         
      //    // Fill aligned sequences accordingly
      //    // Because of initialization row and col, our origninal sequence is misaligned with our score_matrix rows and cols. Hence -1.
      //    switch(pred){
      //       case 0: { //up 
      //          aligned_seq1.push_back(seq1[row-1]); 
      //          aligned_seq2.push_back(base_t::GAP);
      //          row--; 
      //          break;
      //       }
      //       case 1: { //diag
      //          aligned_seq1.push_back(seq1[row-1]);
      //          aligned_seq2.push_back(seq2[col-1]);
      //          row--; col--; 
      //          break;
      //       }
      //       case 2: { //left
      //          aligned_seq1.push_back(base_t::GAP);
      //          aligned_seq2.push_back(seq2[col-1]);
      //          col--;
      //          break;
      //       }
      //       case 3: break; //that it does not print error for last iteration of while loop

      //       default: {
      //          std::cerr << "Invalid predecessor entry!";
      //       }
      //    }
      // }
      

      // // Extract optimal score and return it
      // len_t idx_last = nw_matrix.indxMatrix[m-1][n-1];'
      score_t opt_score = nw_matrix.scoreMatrix[m*n - 1];
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
    
} // namespace NW_SEQ_OPT


