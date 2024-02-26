#include <immintrin.h>
#include <vector>
#include <iostream>

#pragma once


namespace NW_SEQ_OPT {
    void print_vec(const std::vector<u_int32_t>& vec);

    void print_mmvec(const __m256i_u vec);

    void mmExtract(const __m256i_u &mmvec, std::vector<u_int32_t> &vec);

    // Vectorized retrieval of 
    void mmMax( const u_int32_t *vecup_p, 
                const u_int32_t *vecdiag_p,
                const u_int32_t *vecdown_p,
                std::vector<u_int32_t> &max);

    void mmAdd( const u_int32_t *vec_penalty, 
                const u_int32_t *vec_score_old,
                std::vector<u_int32_t> &vec_score_new);
}   // namespace NW_SEQ_OPT