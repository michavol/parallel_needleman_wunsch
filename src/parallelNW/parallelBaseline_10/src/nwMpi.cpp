#include "../include/nwMpi.hpp"

namespace NW_PAR {
    /*
     void createMpiCellType(MPI_Datatype* mpi_cell_t) {
        // assumes that total size of Cell is size of some MPI integer type
        // tested in a simplified setting :
        // MacBook Pro, Dual-Core Intel i5
        // using Open MPI 4.1.1, Apple Clang 12.0.0
        MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(Cell), mpi_cell_t);
    }
    */
    
    /*
    void createMpiBaseType(MPI_Datatype* mpi_base_t) {
        // assumes that total size of base_t is size of some MPI integer type
        // tested in a simplified setting :
        // MacBook Pro, Dual-Core Intel i5
        // using Open MPI 4.1.1, Apple Clang 12.0.0
        MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(base_t), mpi_base_t);
    }
    */
    
    void distributeSequences(const seq_t& seq1, const seq_t& seq2,
                len_t& local_block_height,
                len_t& other_block_height,
                MPI_Comm communicator,
                int local_rank, int comm_size)
    {
        // decide block height
        len_t total_rows = seq1.size();
        len_t total_cols = seq2.size();
        other_block_height = total_rows / comm_size;
        local_block_height = other_block_height + total_rows % comm_size;

        //optimise height distribution
        if (local_block_height > 1.5 * other_block_height && comm_size < local_block_height) {
            ++other_block_height;
            local_block_height -= (comm_size  - 1);
        }

        len_t global_row_index = local_block_height;
        std::vector<MPI_Request> requests;
        requests.reserve(3*(comm_size - 1));
        for (int i = 0; i < comm_size; ++i) {
            if (i == local_rank) continue;
            // send size information to other processes
            unsigned int info[3] = {global_row_index, other_block_height, total_cols};

            MPI_Request current_request;
            MPI_Isend(&info, 3, MPI_INT, i, 0, communicator, &current_request);
            requests.push_back(current_request);

            // send sequences
            // non_blocking send can be used since sequences exist for the entire duration of the program

            // part of first sequence
            MPI_Isend(&(seq1[global_row_index]), other_block_height, MPI_UNSIGNED_CHAR, i, 1, communicator, &current_request);
            requests.push_back(current_request);

            // full second sequence
            MPI_Isend(seq2.data(), total_cols, MPI_UNSIGNED_CHAR, i, 2, communicator, &current_request);
            requests.push_back(current_request);

            global_row_index += other_block_height;
        }

        // wait until all processes got data
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }

    void recieveSequences(seq_t& seq1, seq_t& seq2,
               len_t& global_row_index,
               MPI_Comm communicator,
               int sender_rank, int local_rank)
    {
        // get size information
        unsigned int info[3];
        MPI_Recv(&info, 3, MPI_INT, sender_rank, 0, communicator, MPI_STATUS_IGNORE);

        global_row_index = info[0];
        len_t num_rows = info[1];
        len_t num_cols = info[2];

        // reserve space for sequences
        seq1.resize(num_rows);
        seq2.resize(num_cols);

        // get first sequence
        MPI_Recv(seq1.data(), num_rows, MPI_UNSIGNED_CHAR, sender_rank, 1, communicator, MPI_STATUS_IGNORE);

        // get second sequence
        MPI_Recv(seq2.data(), num_cols, MPI_UNSIGNED_CHAR, sender_rank, 2, communicator, MPI_STATUS_IGNORE);
    }

    void sendSegmentOfRow(const score_mat_t& score_matrix,
            const len_t row_index,
            const len_t start_index,
            const len_t num_entries,
            unsigned int tag,
            MPI_Comm communicator,
            int local_rank, int recieve_rank)
    {
        MPI_Request ignore;
        MPI_Isend(&(score_matrix[row_index][start_index]), num_entries,
                MPI_LONG_LONG_INT, recieve_rank, tag, communicator, &ignore);
    }

    void recieveSegmentOfRow(score_line_t& row_above,
               const len_t start_index,
               const len_t num_entries,
               unsigned int tag,
               MPI_Comm communicator,
               int sender_rank, int local_rank)
    {
        MPI_Recv(&(row_above[start_index]), num_entries, MPI_LONG_LONG_INT, sender_rank,
                tag, communicator, MPI_STATUS_IGNORE);
    }

    void sendScoreMatrix(const score_mat_t& score_matrix,
                const len_t num_rows,
                const len_t num_cols,
                const len_t global_row_index,
                MPI_Comm communicator,
                int local_rank, int recieve_rank)
    {
        // use blocking send since matrix is deleted afterwards
        // send all rows independently
        for (unsigned i = 0; i < num_rows; ++i) {
            MPI_Send(score_matrix[i].data(), num_cols, MPI_LONG_LONG_INT, recieve_rank, i, communicator);
        }


    }

    void recieveScoreMatrices(score_mat_t& score_matrix,
                const len_t total_rows,
                const len_t total_cols,
                const len_t local_block_height,
                const len_t other_block_height,
                MPI_Comm communicator,
                int local_rank, int comm_size)
    {
        //calculate position for insertions
        len_t row_index = local_block_height;

        // store al recieve requests
        std::vector<MPI_Request> requests;
        requests.reserve((total_cols-local_block_height)*(comm_size-1));

        //recieve data from all processes, append in score matrix
        for (unsigned process = 0; process < comm_size; ++process) {
            if (process == local_rank) continue;
            MPI_Request req;
            for (unsigned i = 0; i < other_block_height; ++i) {
                MPI_Irecv(score_matrix[row_index + i].data(), total_cols, MPI_LONG_LONG_INT, process, i, communicator, &req);
                requests.push_back(req);
            }
            row_index += other_block_height;
        }
        // wait until all processes wrote data
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    }

} // namespace NW_PAR
