module MPIParameters
#ifndef _DEBUG
  use MPI
#endif
implicit none
  integer, save :: irank_mpi, ierr_mpi, iprocess
contains
  subroutine MPIInitialization
#ifdef _DEBUG
    irank_mpi = 0
    ierr_mpi  = 0
    iprocess  = 1
#else
    call MPI_INIT(ierr_mpi)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, iprocess , ierr_mpi)
    call MPI_COMM_RANK(MPI_COMM_WORLD, irank_mpi, ierr_mpi)
#endif
  end subroutine
end module
