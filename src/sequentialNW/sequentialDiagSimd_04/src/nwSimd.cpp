/** 
 *  @file   userFunctions.hpp 
 *  @brief  Functions to calculate gap penalties and similarity between bases
 *  @author Sequence Aligners 
 *  @date   2021-11-15 
 ***********************************************/
#include <nwSimd.hpp>


namespace NW_SEQ_OPT {
    void print_vec(const std::vector<u_int32_t>& vec){
        unsigned size = vec.size();
        for (unsigned i = 0; i < size; i++){
            std::cout << vec[i] << '\t';
        }
        std::cout << std::endl;
    }

    void print_mmvec(const __m256i_u vec){
        unsigned size = 8;
        std::vector<int32_t> x(8);
        x[0] = _mm256_extract_epi32(vec,0);
        x[1] = _mm256_extract_epi32(vec,1);
        x[2] = _mm256_extract_epi32(vec,2);
        x[3] = _mm256_extract_epi32(vec,3);
        x[4] = _mm256_extract_epi32(vec,4);
        x[5] = _mm256_extract_epi32(vec,5);
        x[6] = _mm256_extract_epi32(vec,6);
        x[7] = _mm256_extract_epi32(vec,7);
        for (unsigned i = 0; i < size; i++){
            std::cout << x[i] << '\t';
        }
        std::cout << std::endl;
    }

    void mmExtract(const __m256i_u &mmvec, std::vector<u_int32_t> &vec){
        if (vec.size() < 8) std::cerr << "Vector argument too small";
        vec[0] = _mm256_extract_epi32(mmvec,0);
        vec[1] = _mm256_extract_epi32(mmvec,1);
        vec[2] = _mm256_extract_epi32(mmvec,2);
        vec[3] = _mm256_extract_epi32(mmvec,3);
        vec[4] = _mm256_extract_epi32(mmvec,4);
        vec[5] = _mm256_extract_epi32(mmvec,5);
        vec[6] = _mm256_extract_epi32(mmvec,6);
        vec[7] = _mm256_extract_epi32(mmvec,7);
    }

    // Vectorized retrieval of 
    void mmMax( const u_int32_t *vecup_p, 
                const u_int32_t *vecdiag_p,
                const u_int32_t *vecdown_p,
                std::vector<u_int32_t> &max){

        __m256i_u mmvecup_p = _mm256_loadu_si256((__m256i_u *)vecup_p);
        __m256i_u mmvecdiag_p = _mm256_loadu_si256((__m256i_u *)vecdiag_p);
        __m256i_u mmvecdown_p = _mm256_loadu_si256((__m256i_u *)vecdown_p);

        mmvecup_p = _mm256_max_epi32(mmvecup_p, mmvecdiag_p);
        mmvecdown_p = _mm256_max_epi32(mmvecup_p, mmvecdown_p);

        mmExtract(mmvecdown_p, max);
    }

    void mmAdd( const u_int32_t *vec_penalty, 
                const u_int32_t *vec_score_old,
                std::vector<u_int32_t> &vec_score_new
                ){
        
        __m256i_u mmvec_penalty = _mm256_loadu_si256((__m256i_u *)vec_penalty);
        __m256i_u mmvec_score_old = _mm256_loadu_si256((__m256i_u *)vec_score_old);

        mmvec_score_old = _mm256_add_epi32(mmvec_penalty, mmvec_score_old);

        mmExtract(mmvec_score_old, vec_score_new);
    }

}   // namespace NW_SEQ_OPT