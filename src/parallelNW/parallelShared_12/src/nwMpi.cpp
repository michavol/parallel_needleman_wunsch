#include "../include/nwMpi.hpp"

namespace NW_PAR {

    bool hasSuccessorOnNode(MpiEnvironment& mpi_env) {
        return mpi_env.node_rank < mpi_env.node_size - 1;
    }

    bool hasPredecessorOnNode(MpiEnvironment& mpi_env) {
        return mpi_env.node_rank > 0;
    }

    bool hasSuccessorOnOtherNode(MpiEnvironment& mpi_env) {
        return !hasSuccessorOnNode(mpi_env) && mpi_env.global_successor != MpiEnvironment::RANK_UNDEF;
    }

    bool hasPredecessorOnOtherNode(MpiEnvironment& mpi_env) {
        return mpi_env.node_rank == 0 && mpi_env.global_rank != 0;

    }

    void prepareMpi(MpiEnvironment& mpi_env, NodeList& nodes) {
        // rank and size in world
        MPI_Comm_dup(MPI_COMM_WORLD, &(mpi_env.global));
        MPI_Comm_rank(mpi_env.global, &(mpi_env.global_rank));
        MPI_Comm_size(mpi_env.global, &(mpi_env.global_size));

        //initialise send request
        mpi_env.send_request = MPI_REQUEST_NULL;

        // split to nodes
        const int split_key = 0;
        // if all processes use same key, they are ordered as in global
        // the ordering might be assumed in other communications
        // (for example to determine where to send data designated for another node)
        MPI_Info split_info = MPI_INFO_NULL; //we don't provide any information
        MPI_Comm_split_type(mpi_env.global, MPI_COMM_TYPE_SHARED, split_key, split_info, &(mpi_env.node));

        // rank and size in node
        MPI_Comm_rank(mpi_env.node, &(mpi_env.node_rank));
        MPI_Comm_size(mpi_env.node, &(mpi_env.node_size));


        // create groups for neighbouring pairs on node
        if (hasSuccessorOnNode(mpi_env)) {
             mpi_env.node_successor = mpi_env.node_rank + 1;
             mpi_env.global_successor = MpiEnvironment::RANK_UNDEF;
        } else {
            mpi_env.node_successor = MpiEnvironment::RANK_UNDEF;

        }
        if (hasPredecessorOnNode(mpi_env)) {
            mpi_env.node_predecessor = mpi_env.node_rank - 1;
            mpi_env.global_predecessor = MpiEnvironment::RANK_UNDEF;
        } else {
            mpi_env.node_predecessor = MpiEnvironment::RANK_UNDEF;

        }


        // prepare communication between nodes
        if (mpi_env.global_rank == 0) {
            //set predecessor to undefined
            mpi_env.global_predecessor = MpiEnvironment::RANK_UNDEF;


            // get information from own nodes:
            int last_rank = 0;
            if (hasSuccessorOnNode(mpi_env)) {
                /*
                 * int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                 *      int source, int tag, MPI_Comm comm, MPI_Status *status)
                 */
                MPI_Recv(&last_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_env.node, MPI_STATUS_IGNORE);
            }

            NodeGroup group;
            group.leader = mpi_env.global_rank;
            group.last = last_rank;
            group.size = mpi_env.node_size;
            nodes.push_back(group);
            // get information from other nodes:
            int known_processes = mpi_env.node_size;

            while (known_processes < mpi_env.global_size) {
                int node_info[3];
                MPI_Recv(node_info, 3, MPI_INT, MPI_ANY_SOURCE, 0, mpi_env.global, MPI_STATUS_IGNORE);
                group.leader = node_info[0];
                group.last = node_info[1];
                group.size= node_info[2];
                nodes.push_back(group);

                known_processes += group.size;
            }

            // all nodes and their respective leaders are known to rank 0

            // send rank of leader of successor node to each node
            // in order to send intermediate results
            unsigned num_other_nodes = nodes.size() - 1;

            if (num_other_nodes > 0) {

                int neighbour_info[nodes.size()][2];

                std::vector<MPI_Request> sends(num_other_nodes);

                for (unsigned i = 1; i < num_other_nodes; ++i) {
                    neighbour_info[i][0] = nodes[i + 1].leader; // global rank of successor node
                    neighbour_info[i][1] = nodes[i - 1].last; // global rank of predecessor node

                    MPI_Isend(neighbour_info[i], 2, MPI_INT, nodes[i].leader, 0,
                            mpi_env.global, &(sends[i - 1]));
                }
                // inform last that it has no successor
                neighbour_info[num_other_nodes][0] = MpiEnvironment::RANK_UNDEF; // last node has no successor
                neighbour_info[num_other_nodes][1] = nodes[num_other_nodes - 1].last; // global rank of predecessor node

                /*
                 * int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
                 *      int tag, MPI_Comm comm, MPI_Request *request)
                 */
                MPI_Isend(neighbour_info[num_other_nodes], 2, MPI_INT, nodes[num_other_nodes].leader, 0,
                        mpi_env.global, &(sends[num_other_nodes - 1]));

                // inform last process on node where to send data
                /*
                 * int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                 *      int tag, MPI_Comm comm)
                 */
                MPI_Send(&(nodes[1].leader), 1, MPI_INT, nodes[0].last, 0, mpi_env.global);

                //wait until al sends are finished
                MPI_Waitall(num_other_nodes, sends.data(), MPI_STATUSES_IGNORE);

            } else {
                // inform last process that there is no other node
                int no_successor = MpiEnvironment::RANK_UNDEF;

                if (hasSuccessorOnNode(mpi_env)) {
                    MPI_Send(&no_successor, 1, MPI_INT, nodes[0].last, 0, mpi_env.global);
                } else {
                    mpi_env.global_successor = no_successor;
                }
            }

        } else if (mpi_env.node_rank == 0) {
            // send information on node group
            int last_rank = mpi_env.global_rank;
            if (hasSuccessorOnNode(mpi_env)) {
                MPI_Recv(&last_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_env.node, MPI_STATUS_IGNORE);
            }

            int node_info[3] = {mpi_env.global_rank, last_rank, mpi_env.node_size};
            /*
             * int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
             *      int tag, MPI_Comm comm)
             */
            MPI_Send(node_info, 3, MPI_INT, 0, 0, mpi_env.global);

            // get leader of successor node and last of predecessor node
            int neighbour_info[2];
            /*
             * int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
             *      int source, int tag, MPI_Comm comm, MPI_Status *status)
             */
            MPI_Recv(neighbour_info, 2, MPI_INT, 0, 0, mpi_env.global, MPI_STATUS_IGNORE);

            // set global partner to predecessor
            mpi_env.global_predecessor = neighbour_info[1];

            if (hasSuccessorOnNode(mpi_env)) {
                MPI_Send(neighbour_info, 1, MPI_INT, last_rank, 0, mpi_env.global);
            } else {
                mpi_env.global_successor = neighbour_info[0];
            }


        } else if (!hasSuccessorOnNode(mpi_env)) {
            // send rank to leader
            MPI_Send(&(mpi_env.global_rank), 1, MPI_INT, 0, 0, mpi_env.node);

            // get information on where to send results
            MPI_Recv(&(mpi_env.global_successor), 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_env.global, MPI_STATUS_IGNORE);

        }

    }

    void shareSequences(MpiEnvironment& mpi_env, NodeList& nodes, NwMatrix& nw_matrix, seq_t& seq1, seq_t& seq2) {
        //create accessible windows to distribute sequences
        MPI_Win global_sequences;
        int global_seq_size;
        void* global_seq_base;

        if (mpi_env.global_rank == 0) {
            global_seq_size = seq1.size + seq2.size;

        } else {
            global_seq_size = 0;
        }
        /*
         * int MPI_Win_allocate (MPI_Aint size, int disp_unit, MPI_Info info,
         *            MPI_Comm comm, void *baseptr, MPI_Win *win)
         */

        MPI_Win_allocate(global_seq_size*sizeof(base_t), sizeof(base_t),
                MPI_INFO_NULL, mpi_env.global, &global_seq_base, &global_sequences);


        // ==================================
        // distribute sequences between nodes
        // ==================================


        if (mpi_env.node_rank == 0) {
            int first_size;     // total size of seq 1 on corresponding node
            int second_size;    // total size of seq 2 on node
            int first_start;    // position in the full seq 1 where local part of seq 1 starts

            if (mpi_env.global_rank == 0) {

                //copy sequences to window
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, global_sequences);

                base_t* current_pos = (base_t*)global_seq_base;
                std::copy(seq2.data, seq2.data + seq2.size, current_pos);
                current_pos += seq2.size;
                std::copy(seq1.data, seq1.data + seq1.size, current_pos);

                MPI_Win_sync(global_sequences);

                MPI_Win_unlock(0, global_sequences);


                // distribute seq1 evenly between processes

                // entries in seq1 to be distributed to nodes
                const int min_rows = seq1.size / mpi_env.global_size;
                int remaining_rows = seq1.size % mpi_env.global_size;

                //data for first node
                second_size = seq2.size;
                first_start = 0;

                int additional =  std::min(mpi_env.node_size, remaining_rows);
                remaining_rows -= additional;

                first_size = min_rows * mpi_env.node_size + additional;

                int global_row_index = first_size; //position where next block starts

                int num_other_nodes = nodes.size() - 1;

                int sequence_info[nodes.size()][3];

                std::vector<MPI_Request> sends(num_other_nodes);

                for (unsigned i = 1; i <= num_other_nodes; ++i) {

                    additional =  std::min(nodes[i].size, remaining_rows);
                    int rows = min_rows * nodes[i].size + additional;
                    remaining_rows -= additional;

                    // send information
                    sequence_info[i][0] = rows;            // size of corresponding part of seq 1
                    sequence_info[i][1] = second_size;     // size of seq 2
                    sequence_info[i][2] = global_row_index;// position of part of seq 1 in global seq 1

                    /*
                     * int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
                     *      int tag, MPI_Comm comm, MPI_Request *request)
                     */
                    MPI_Isend(sequence_info[i], 3, MPI_INT, nodes[i].leader, 0, mpi_env.global, sends.data() + i - 1);

                    global_row_index += rows; // advance position

                }

                //wait until al sends are finished
                MPI_Waitall(num_other_nodes, sends.data(), MPI_STATUSES_IGNORE);

            } else {
                //get information about sequences
                int sequence_info[3];
                MPI_Recv(sequence_info, 3, MPI_INT, 0, 0, mpi_env.global, MPI_STATUS_IGNORE);
                first_size = sequence_info[0];
                second_size = sequence_info[1];
                first_start = sequence_info[2];

            }

            // create window on node
            mpi_env.node_seq_size = first_size + second_size;
            MPI_Win_allocate_shared(mpi_env.node_seq_size * sizeof(base_t), sizeof(base_t), MPI_INFO_NULL,
                    mpi_env.node, &(mpi_env.node_seq_base), &(mpi_env.node_sequences));

            // copy data to shared window
            if (mpi_env.global_rank == 0) {
                //start access epoch - not needed since no RMA accesses are done on shared memory
                /*
                 * int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win)
                 */
//                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, mpi_env.node_sequences);
                base_t* current_pos = (base_t*)mpi_env.node_seq_base;
                std::copy(seq2.data, seq2.data + seq2.size, current_pos);
                current_pos += seq2.size;
                std::copy(seq1.data + first_start, seq1.data + first_start + first_size, current_pos);

                /*
                 * int MPI_Win_sync (MPI_Win win)
                 */
                MPI_Win_sync(mpi_env.node_sequences);
//
//                MPI_Win_unlock(0, mpi_env.node_sequences);

                delete[] seq1.data;
                delete[] seq2.data;
            } else {
                // get sequences from global window and store in node window
                /*
                 * MPI_Get(void *origin_addr, int origin_count, MPI_Datatype
                 *      origin_datatype, int target_rank, MPI_Aint target_disp,
                 *      int target_count, MPI_Datatype target_datatype, MPI_Win win)
                 */
                // start access epoch
                MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, global_sequences);
//                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, mpi_env.node_sequences);


                // second sequence comes first
                base_t* seq_ptr = (base_t*)mpi_env.node_seq_base;
                MPI_Get(seq_ptr, second_size, MPI_UNSIGNED_CHAR, 0, 0, second_size,
                        MPI_UNSIGNED_CHAR, global_sequences);

                // part of first sequence
                seq_ptr += second_size;
                MPI_Get(seq_ptr, first_size, MPI_UNSIGNED_CHAR, 0, second_size + first_start, first_size,
                        MPI_UNSIGNED_CHAR, global_sequences);

                //end access
                MPI_Win_flush_local(0, global_sequences);

                MPI_Win_sync(mpi_env.node_sequences);

                MPI_Win_unlock(0, global_sequences);

            }

            // ============================
            // distribute sequences on node
            // ============================

            // set known values for leader
            nw_matrix.global_row_index = first_start;
            seq2.data = (base_t*)mpi_env.node_seq_base;
            seq2.size = second_size;
            seq1.data = (base_t*)mpi_env.node_seq_base + second_size;

            // distribute seq1 evenly between processes
            const int min_rows = first_size / mpi_env.node_size;
            int remaining_rows = first_size % mpi_env.node_size;


            int rows = min_rows;
            if (remaining_rows > 0) {
                ++rows;
                --remaining_rows;
            }

            seq1.size = rows;

            int local_row_index = rows;

            int sequence_info[mpi_env.node_size - 1][4];

            std::vector<MPI_Request> sends(mpi_env.node_size - 1);

            for (int i = 1; i < mpi_env.node_size; ++i) {
                rows = min_rows;

                if (remaining_rows > 0) {
                    ++rows;
                    --remaining_rows;
                }

                // send information
                sequence_info[i-1][0] = local_row_index; // position of corresponding part of seq 1 on node seq 1
                sequence_info[i-1][1] = rows;            // size of corresponding part of seq 1
                sequence_info[i-1][2] = second_size;     // size of seq 2
                sequence_info[i-1][3] = local_row_index + first_start;   // position of correspoinding part of seq 1 in global seq 1

                MPI_Isend(sequence_info[i-1], 4, MPI_INT, i, 0, mpi_env.node, sends.data() + i - 1);

                local_row_index += rows;

            }

            MPI_Waitall(sends.size(), sends.data(), MPI_STATUSES_IGNORE);

        } else {
            /*
             * MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info,
             *               MPI_Comm comm, void *baseptr, MPI_Win *win)
             */
            mpi_env.node_seq_size = 0;
            MPI_Win_allocate_shared(mpi_env.node_seq_size * sizeof(base_t), sizeof(base_t), MPI_INFO_NULL,
                    mpi_env.node,  &(mpi_env.node_seq_base), &(mpi_env.node_sequences));

            int sequence_info[4];
            MPI_Recv(sequence_info, 4, MPI_INT, 0, 0, mpi_env.node, MPI_STATUS_IGNORE);

            int node_row_index = sequence_info[0];
            seq1.size = sequence_info[1];
            seq2.size = sequence_info[2];
            nw_matrix.global_row_index = sequence_info[3];

            /*
             * int MPI_Win_shared_query (MPI_Win win, int rank, MPI_Aint *size,
             *            int *disp_unit, void *baseptr)
             */
            void* window_ptr;
            int displacement;
            MPI_Aint window_size;
            MPI_Win_shared_query(mpi_env.node_sequences, 0, &window_size, &displacement, &window_ptr);

            seq2.data = (base_t*)window_ptr;
            seq1.data = seq2.data + seq2.size + node_row_index;

        }

        //set information for matrix
        nw_matrix.rows = seq1.size;
        nw_matrix.cols = seq2.size + 1;

        //free global window
        MPI_Win_free(&global_sequences);

    }

    void shareMatrix(MpiEnvironment& mpi_env, NwMatrix& nw_matrix) {
        // determine size of local matrix
        // add one in the size for an additional byte used as flag
        mpi_env.node_mat_size = nw_matrix.rows * nw_matrix.cols + 1;
        // create window for matrix on node
        /*
         * MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info,
         *               MPI_Comm comm, void *baseptr, MPI_Win *win)
         */

        const int displacement = sizeof(score_t) + sizeof(len_t) + sizeof(pre_t);
        MPI_Win_allocate_shared(mpi_env.node_mat_size*displacement, 1, MPI_INFO_NULL,
                mpi_env.node, &(mpi_env.node_mat_base), &(mpi_env.node_matrix));

        // set pointers to matrix
        nw_matrix.scores = (score_t*)mpi_env.node_mat_base;
        nw_matrix.gaps = (len_t*)(nw_matrix.scores + mpi_env.node_mat_size - 1);
        nw_matrix.pres = (pre_t*)(nw_matrix.gaps + mpi_env.node_mat_size - 1);

        mpi_env.local_end = (len_t*)(nw_matrix.pres + mpi_env.node_mat_size - 1);
        *(mpi_env.local_end) = 0;

        // create global window for last rows
        if (hasSuccessorOnOtherNode(mpi_env)) {
            mpi_env.global_mat_size = nw_matrix.cols;
        } else {
            mpi_env.global_mat_size = 0;
        }
        /*
         * int MPI_Win_allocate (MPI_Aint size, int disp_unit, MPI_Info info,
         *            MPI_Comm comm, void *baseptr, MPI_Win *win)
         */
        MPI_Win_allocate(mpi_env.global_mat_size*displacement, 1, MPI_INFO_NULL,
                mpi_env.global, &(mpi_env.global_mat_base), &(mpi_env.global_matrix));



        //create object for last row of predecessor
        NwLine* other_row = new NwLine();

        // query matrix of successor
        /*
         * int MPI_Win_shared_query (MPI_Win win, int rank, MPI_Aint *size,
         *            int *disp_unit, void *baseptr)
         */
        if (hasPredecessorOnNode(mpi_env)) {
            void* other_base;
            int other_disp;
            MPI_Aint other_total_size;
            MPI_Win_shared_query(mpi_env.node_matrix, mpi_env.node_rank - 1,
                    &other_total_size, &other_disp, &other_base);

            // size of one row
            len_t row_size = nw_matrix.cols;

            // number of entries in the matrix without last row, ie. index of first entry in last row
            // subtract one because of the flag entry
            len_t last_row_start = other_total_size / displacement - row_size - 1;

            // set pointer in prev_row
            other_row->size = row_size;
            //last row of score matrix
            other_row->scores = (score_t*)other_base + last_row_start;

            // one row of scores, rest is gaps
            other_row->gaps = (len_t*) (other_row->scores + row_size) + last_row_start;
            other_row->pres = (pre_t*) (other_row->gaps + row_size) + last_row_start;

            // set pointer to predecessor flag
            mpi_env.predecessor_end = (len_t*)(other_row->pres + row_size);

        } else {
            other_row->size = nw_matrix.cols;

            //allocate memory
            void* other_base = std::malloc(displacement*other_row->size);

            other_row->scores = (score_t*)other_base;
            other_row->gaps = (len_t*) (other_row->scores + other_row->size);
            other_row->pres = (pre_t*) (other_row->gaps + other_row->size);

            // set value to o
           mpi_env.predecessor_end = new len_t(0);

        }

        nw_matrix.prev_row = other_row;
    }

    void calculateBlockSize(MpiEnvironment& mpi_env, NwMatrix& nw_matrix) {
        mpi_env.block_size = 63;
    }

    void getPrevRow(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, const len_t start, const len_t end) {
        if (hasPredecessorOnOtherNode(mpi_env)) {
            // predecessor is not on same node

            const len_t segment_size = end - start;
            const len_t row_size = nw_matrix.cols;
            MPI_Aint target_start_position = sizeof(score_t)*start;

            //wait as long as predecessor has not completed the row above
            while(*(mpi_env.predecessor_end) < end) {
                /*
                 * int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                 *      int source, int tag, MPI_Comm comm, MPI_Status *status)
                 */
                MPI_Recv(mpi_env.predecessor_end, 1, MPI_INT, mpi_env.global_predecessor, 0, mpi_env.global, MPI_STATUS_IGNORE);
            }
            /*
             * MPI_Get(void *origin_addr, int origin_count, MPI_Datatype
             *      origin_datatype, int target_rank, MPI_Aint target_disp,
                   int target_count, MPI_Datatype target_datatype, MPI_Win win)
             */

            MPI_Win_lock(MPI_LOCK_SHARED, mpi_env.global_predecessor, 0, mpi_env.global_matrix);
            // get scores from other window
            void* score_insert = (void*)(nw_matrix.prev_row->scores + start);
            MPI_Get(score_insert, segment_size, MPI_INT, mpi_env.global_predecessor, target_start_position,
                    segment_size, MPI_INT, mpi_env.global_matrix);

            // get gaps form other window
            void* gap_insert = (void*)(nw_matrix.prev_row->gaps + start);
            target_start_position = sizeof(score_t)*row_size + sizeof(len_t)*start;
            MPI_Get(gap_insert, segment_size, MPI_UNSIGNED, mpi_env.global_predecessor, target_start_position,
                    segment_size, MPI_UNSIGNED, mpi_env.global_matrix);

            // get pres
            void* pre_insert = (void*)(nw_matrix.prev_row->pres + start);
            target_start_position = (sizeof(score_t) + sizeof(len_t))*(row_size) + sizeof(pre_t)*start;
            MPI_Get(pre_insert, segment_size, MPI_UNSIGNED_CHAR, mpi_env.global_predecessor, target_start_position,
                    segment_size, MPI_UNSIGNED_CHAR, mpi_env.global_matrix);

            //complete transfer to local memory
            MPI_Win_flush_local(mpi_env.global_predecessor, mpi_env.global_matrix);

            MPI_Win_unlock(mpi_env.global_predecessor, mpi_env.global_matrix);


        } else if (hasPredecessorOnNode(mpi_env)) {
            // predecessor is on same node
            //spin until predecessor completed the row above
            while (*mpi_env.predecessor_end < end) {
                __asm__("nop");
            }

        }

    }



    void sendLastRow(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, const len_t start, const len_t end) {

       if (hasSuccessorOnOtherNode(mpi_env)) {
            //if successor is not on same node


            //copy data of last row to window

            /*
             * int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win)
             */
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mpi_env.global_rank, 0, mpi_env.global_matrix);

            const len_t segment_size = end - start;
            const len_t row_size = nw_matrix.cols;
            const len_t segment_start = (nw_matrix.rows - 1) * nw_matrix.cols + start;

            score_t* score_insert = (score_t*)mpi_env.global_mat_base + start;
            score_t* scores = nw_matrix.scores + segment_start;
            // copy scores to win
            std::copy(scores, scores + segment_size, score_insert);

            len_t* gap_insert = (len_t*) (score_insert - start + row_size) + start;
            len_t* gaps = nw_matrix.gaps + segment_start;
            // copy gaps to win
            std::copy(gaps, gaps + segment_size, gap_insert);

            pre_t* pre_insert = (pre_t*) (gap_insert - start + row_size) + start;
            pre_t* pres = nw_matrix.pres + segment_start;
            //copy pres to win
            std::copy(pres, pres + segment_size, pre_insert);


            /*
             * int MPI_Win_sync (MPI_Win win)
             */
            MPI_Win_sync(mpi_env.global_matrix);

            MPI_Win_unlock(mpi_env.global_rank, mpi_env.global_matrix);

//            MPI_Wait(&(mpi_env.send_request), MPI_STATUS_IGNORE);
            *(mpi_env.local_end) = end;


            /*
             * int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
             *      int tag, MPI_Comm comm, MPI_Request *request)
             */
            MPI_Isend(mpi_env.local_end, 1, MPI_INT, mpi_env.global_successor, 0, mpi_env.global, &(mpi_env.send_request));


        } else if (hasSuccessorOnNode(mpi_env)){
            // start access epoch
            /*
             * int MPI_Win_post(MPI_Group group, int assert, MPI_Win win)
             */
            *(mpi_env.local_end) = end;
            MPI_Win_sync(mpi_env.node_matrix);

        }

    }

    void cleanMpi(MpiEnvironment& mpi_env, NwMatrix& nw_matrix, seq_t& seq1, seq_t& seq2) {
        //wait for other to complete access on local matrix
        /*
         * int MPI_Wait(MPI_Request *request, MPI_Status *status)
         */
        MPI_Wait(&(mpi_env.send_request), MPI_STATUS_IGNORE);

        //free windows
        MPI_Win_free(&(mpi_env.node_sequences));
        MPI_Win_free(&(mpi_env.node_matrix));
        MPI_Win_free(&(mpi_env.global_matrix));

        // free communicators
        MPI_Comm_free(&(mpi_env.global));
        MPI_Comm_free(&(mpi_env.node));


        if (mpi_env.node_rank == 0) {
            std::free(nw_matrix.prev_row->scores);
        }

        mpi_env.node_mat_base = nullptr;
        mpi_env.node_seq_base  = nullptr;

    }
} // namespace NW_PAR
