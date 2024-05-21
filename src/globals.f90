module globals
    use mpi
    implicit none
    ! seemingly best kept at 8, 4 seems to cause some issues
    ! which have yet to be looked at more closely
    integer, parameter :: real_kind=8
    integer, parameter :: mpi_r=MPI_REAL8
    ! depending on the MPI installation, using other integer kinds
    ! could also cause issues
    integer, parameter :: int_kind=4
    integer, parameter :: mpi_i=MPI_INTEGER4

end module globals