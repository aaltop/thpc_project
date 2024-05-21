module globals
    use mpi
    implicit none
    integer, parameter :: real_kind=8
    integer, parameter :: mpi_r=MPI_REAL8
    integer, parameter :: int_kind=4
    integer, parameter :: mpi_i=MPI_INTEGER4

end module globals