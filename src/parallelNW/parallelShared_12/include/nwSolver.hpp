#pragma once

#include "nwHelper.hpp"
#include "nwMpi.hpp"
#include "nwSolver.hpp"

namespace NW_PAR {


    /**
     * @brief
     *
     * @pre one process has read sequences, processes know where they are
     * @post sequences and relevant parts of (empty) score matrix are available for all processes
     * @param nw_matrix
     * @param mpi_env
     * @param nodes
     * @param seq1
     * @param seq2
     */
    void setup(NwMatrix& nw_matrix, MpiEnvironment& mpi_env, NodeList& nodes, seq_t& seq1, seq_t& seq2);
    /**
     * @brief
     *
     * @pre preceeding process finished the corresponding row in its score matrix
     * @post a block of the score matrix is computed an available for the succeeding process
     *
     * when accessing first row of succeeding process, lock the window for that process
     * -> other process is locked only when he also calls corresponding function
     * -> locking guarantees that changes made in the locked region are visible  after the unlock, but not  before
     *
     * notation: lock referes to MPI functions "post, complete" and "start, wait"
     */
    void completeBlock();


    void collectMatrix();
}//namespace NW_PAR
