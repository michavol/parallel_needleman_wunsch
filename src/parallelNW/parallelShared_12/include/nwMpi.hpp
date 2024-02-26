#pragma once

/**
 *  @file   nwMpi.hpp
 *  @brief  Helper Functions to use MPI communication in Sequence Alignment
 *  @author Sequence Aligners
 *  @date   2021-11-18
 ***********************************************/

#include <mpi.h>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "nwHelper.hpp"

namespace NW_PAR {


    /**
     * @struct MpiEnvironment
     * @brief  a struct containing all relevant information on mpi communicators
     *
     * uccessor_rank is only important for last process in each node
     *
     * arrays and mpi objects must be freed separately
     */
    struct MpiEnvironment {
            static const int RANK_UNDEF = -1;
            //communicators
            MPI_Comm global;
            MPI_Comm node;

            MPI_Win global_matrix; // has displacement 1
            MPI_Win node_matrix; // has displacement 1
            MPI_Win node_sequences; // has displacement sizeof(base_t)

            MPI_Request send_request;


            void* global_mat_base;
            void* node_mat_base; // actual allocated memory is 1 byte lagrger (flag)
            void* node_seq_base;

            size_t global_mat_size;
            size_t node_mat_size;
            size_t node_seq_size;

            len_t* predecessor_end;
            len_t* local_end;

            len_t block_size;

            int global_rank;
            int node_rank;
            int global_size;
            int node_size;


            //int global_partner_rank;

            int global_predecessor;
            int global_successor;
            int node_predecessor;
            int node_successor;
    };

    struct NodeGroup {
            int leader;
            int last;
            int size;
    };
    using NodeList = std::vector<NodeGroup>;

    /**
     * use type aliases for mpi types
     * base_t   =   MPI_UNSIGNED_CHAR;
     * pre_t    =   MPI_UNSIGNED_CHAR;
     * score_t  =   MPI_INT;
     * len_t    =   MPI_UNSIGNED;
     */


    bool hasSuccessorOnNode(MpiEnvironment& mpi_env);

    bool hasPredecessorOnNode(MpiEnvironment& mpi_env);

    bool hasSuccessorOnOtherNode(MpiEnvironment& mpi_env);

    bool hasPredecessorOnOtherNode(MpiEnvironment& mpi_env);

    /**
     * @brief
     *
     * @pre MPI has been initialised, if there are several nodes, then on each node must be at least two processes
     * @post proesses know each other and where they are located
     * @param mpi_env
     * @param nodes only relevant for i/o rank
     *
     * create communicator on each node
     *
     * get ranks in all communicators
     *
     * setup for communication:
     * - create groups for neighbouring processes
     * - know where to get intermediate results
     * - know where to get sequences
     */
    void prepareMpi(MpiEnvironment& mpi_env, NodeList& nodes);



    /**
     * @brief
     *
     * @pre global process 0 has read in sequences
     * @post all nodes have copies of the squeunces available to all, sizo of score matrix is known
     * @param mpi_env
     * @param seq1
     * @param seq2
     *
     * IMPORTANT: second sequence is stored **in front of** first sequence
     */
    void shareSequences(MpiEnvironment& mpi_env, NodeList& nodes, NwMatrix& nw_matrix, seq_t& seq1, seq_t& seq2);

    /**
     * @brief
     *
     * @pre all processes know their segments of the sequences
     * @post
     * @param mpi_env
     * @param nw_matrix
     */
    void shareMatrix(MpiEnvironment& mpi_env, NwMatrix& nw_matrix);


    void calculateBlockSize(MpiEnvironment& mpi_env, NwMatrix& nw_matrix);

    /**
     * @brief
     *
     * @pre
     * @post
     * @param mpi_env
     * @param nw_matrix
     * @param next_position past the end index of calculated columns
     */
    void getPrevRow(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, const len_t start, const len_t end);



    /**
     * @brief
     *
     * @pre previous access epoch to the corresponding window has been completed
     * @post
     * @param mpi_env
     * @param nw_matrix
     * @param next_position index of the first row that is not yet computed, >= block_size
     *
     *
     */
    void sendLastRow(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, const len_t start, const len_t end);

    /**
     * @brief delete communicators, close windows and free memory
     *
     * @pre
     * @post
     */
    void cleanMpi(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, seq_t& seq1, seq_t& seq2);



} // namespace NW_PAR
