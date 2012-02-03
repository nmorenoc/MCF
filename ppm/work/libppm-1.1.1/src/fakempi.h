integer, parameter :: MPI_REAL = 1275069468
integer, parameter :: MPI_INTEGER = 1275069467
integer, parameter :: MPI_DOUBLE_PRECISION = 1275070495
integer, parameter :: MPI_LOGICAL = 1275069469
integer, parameter :: MPI_DOUBLE_COMPLEX = 1275072546
integer, parameter :: MPI_STATUS_SIZE = 5
integer, parameter :: MPI_LOR = 1476395015
integer, parameter :: MPI_COMPLEX = 1275070494

integer, parameter :: MPI_MAX = 1476395009
integer, parameter :: MPI_MIN = 1476395010
integer, parameter :: MPI_SUM = 1476395011
integer, parameter :: MPI_REQUEST_NULL = 738197504



external :: MPI_Abort
external :: MPI_Initialized
external :: MPI_Comm_Size
external :: MPI_Comm_Rank
external :: MPI_Allreduce
external :: MPI_Alltoall
external :: MPI_Gather
external :: MPI_Irecv
external :: MPI_Isend
external :: MPI_Recv
external :: MPI_Reduce
external :: MPI_Send
external :: MPI_SendRecv
external :: MPI_Wtime


