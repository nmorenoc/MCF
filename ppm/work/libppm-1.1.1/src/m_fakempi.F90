#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, &
     DATATYPE,  OP,  COMM, ierror) 

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_ALLREDUCE

subroutine MPI_ALLTOALL(SENDBUF,  SENDCOUNT,  SENDTYPE, &
     RECVBUF, RECVCOUNT,  RECVTYPE, &
     COMM, IERROR) 

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_ALLTOALL

subroutine  MPI_Gather(SENDBUF, SENDCOUNT , SENDTYPE, RECVBUF, &
     RECVTYPE, ROOT, COMM, ierror)

  integer, intent(out) :: ierror
  ierror = 0 
end subroutine 

subroutine MPI_ISEND( BUF, COUNT,  DATATYPE,  SOURCE, &
     TAG,  COMM, REQUEST, IERROR)

   integer, intent(out) :: ierror
   ierror = 0
end subroutine MPI_ISEND


subroutine MPI_IRECV( BUF, COUNT,  DATATYPE,  SOURCE, &
     TAG,  COMM, REQUEST, IERROR)

   integer, intent(out) :: ierror
   ierror = 0
end subroutine MPI_IRECV

subroutine MPI_IRSEND( BUF, COUNT,  DATATYPE,  DEST, &
     TAG, COMM, REQUEST, IERROR)

      integer, intent(out) :: ierror
      ierror = 0
end subroutine MPI_IRSEND

subroutine MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, &
      TAG, COMM,  STATUS,  IERROR)

      integer, intent(out) :: ierror
      ierror = 0
end subroutine MPI_RECV

subroutine MPI_REDUCE( SENDBUF,  RECVBUF, COUNT, &
     DATATYPE, OP,  ROOT,  COMM, &
     IERROR)

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_REDUCE

subroutine MPI_SEND( BUF, COUNT, DATATYPE, DEST, &
       TAG,  COMM, IERROR)

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_SEND

subroutine MPI_SENDRECV(SENDBUF,  SENDCOUNT,  SENDTYPE, &
       DEST,  SENDTAG,  RECVBUF,  RECVCOUNT, &
       RECVTYPE,  SOURCE,  RECVTAG,  COMM, &
       STATUS,  IERROR)

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_SENDRECV

subroutine MPI_INITIALIZED(FLAG,  IERROR)

  integer, intent(out) :: ierror
  integer, intent(out)   ::  flag
  flag = 1
  ierror = 0
end subroutine MPI_INITIALIZED


subroutine MPI_COMM_RANK(COMM, RANK, IERROR)

  integer, intent(out) :: rank 
  integer, intent(out) :: ierror
  rank = 0
  ierror = 0
end subroutine MPI_COMM_RANK


subroutine MPI_COMM_SIZE(COMM, SIZE, IERROR)

  integer, intent(out) :: ierror
  integer, intent(out) :: size 
  size = 1
  ierror = 0
end subroutine MPI_COMM_SIZE

subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, &
     COMM,  IERROR)

  integer, intent(out) :: ierror
  ierror = 0
end subroutine MPI_BCAST

function  MPI_WTIME()
  real :: MPI_WTIME
  MPI_WTIME = 0.0

end function MPI_WTIME

subroutine MPI_ABORT(COMM,  ERRORCODE,  IERROR)
  print *, __FILE__, ':', __LINE__, "MPI_ABORT called"
  stop 20
end subroutine MPI_ABORT

