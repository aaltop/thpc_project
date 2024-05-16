program tsp_ga
    use mpi
    use globals
    use utils
    use tsp
    implicit none

    integer :: io, i, iostat, num_cities, num_considered, generations
    integer, allocatable :: idx(:), routes(:,:)

    real(kind=real_kind), allocatable :: random_val(:), weights(:), &
        distances(:,:)

    ! MPI variables
    ! ------------------------------------
    integer, parameter :: tag = 50

    integer :: id, ntasks, source_id, dest_id, rc, nlen, left_neighbor, right_neighbor

    real(real_kind), allocatable :: real_recv(:), real_send(:)

    integer,dimension(mpi_status_size) :: status

    character(len=MPI_MAX_PROCESSOR_NAME) :: host

    ! ======================================

    ! MPI Initialisation
    ! ----------------------------------

    ! rc is the error message variable
    call mpi_init(rc)
    if (rc/=mpi_success) then
        print *,'MPI initialization failed.'
        stop
    end if

    ! Determines the size of the group associated with a communicator
    call mpi_comm_size(mpi_comm_world,ntasks,rc)

    ! Determines the rank of the calling process in the communicator
    ! (rank in this case is probably just a simple numbering for the process)
    call mpi_comm_rank(mpi_comm_world,id,rc)

    ! The name returned should identify a particular piece of hardware; 
    ! the exact format is implementation defined.
    ! (In this case, this would be the name of the computer/node/whatever?)
    call mpi_get_processor_name(host,nlen,rc)

    ! ======================================

    if ( 0 == id ) print "(a,i3)", "the number of tasks is", ntasks
    ! read in the array of distances form a file
    !
    ! the file should be formatted such that the first line contains
    ! the number of cities, and after that there is the array of
    ! cities' distances from each other
    ! ---------------------------------------------------------------
    io = 1234
    open(newunit=io, file="wg59_dist_array.txt", status="old", action="read", iostat=iostat)
    read(io, *) num_cities
    allocate(distances(num_cities, num_cities))


    iostat = 0
    do i = 1,num_cities
        read(io, *, iostat=iostat) distances(1:num_cities, i)
        if (0 /= iostat) then
            print *, "something went wrong"
            stop
        end if
    end do
    close(io)
    ! ====================================================================

    num_cities = 11
    generations = 1000

    ! allocate(random_val(num_cities))
    allocate(weights(generations))

    ! allocate(idx(generations))
    ! call shuffle(weights, idx)


    allocate(routes(num_cities, generations))

    call find_optimal_route(distances(1:num_cities, 1:num_cities), 10, 15, real(0.95, real_kind), generations, routes)


    ! Calculate the distances gained from the breeding algorithm
    io = 1234
    open(io, file="breed.txt", status="replace", action="write")
    do i = 1, generations
        
        ! print "(11i3)", routes(:,i)
        call calculate_total_distance(routes(:,i), distances, weights(i))
        write(io, *) weights(i)

    end do
    close(io)

    ! Do a random search
    io = 1234
    open(io, file="random_search.txt", status="replace", action="write")
    do i = 1, generations
        call new_route(routes(:,1))
        call calculate_total_distance(routes(:,1), distances, weights(1))
        write(io, *) weights(1)
    end do
    close(io)

    ! Do a "random breed", whereby two random routes are bred
    ! and the child is the new tested route
    io = 1234
    open(io, file="random_breed.txt", status="replace", action="write")
    do i = 1, generations
        call new_route(routes(:,1))
        call new_route(routes(:,2))
        call breed(routes(:,1), routes(:,2), distances(1:num_cities, 1:num_cities), routes(:,3))
        call calculate_total_distance(routes(:,3), distances(1:num_cities, 1:num_cities), weights(1))
        write(io, *) weights(1)
    end do
    close(io)

    ! Do a "better random breed", which is the same as above, but
    ! now the best performing route out of the three is picked each 
    ! generation
    io = 1234
    open(io, file="random_breed_better.txt", status="replace", action="write")
    do i = 1, generations
        call new_route(routes(:,1))
        call calculate_total_distance(routes(:,3), distances(1:num_cities, 1:num_cities), weights(1))
        call new_route(routes(:,2))
        call calculate_total_distance(routes(:,3), distances(1:num_cities, 1:num_cities), weights(2))
        call breed(routes(:,1), routes(:,2), distances(1:num_cities, 1:num_cities), routes(:,3))
        call calculate_total_distance(routes(:,3), distances(1:num_cities, 1:num_cities), weights(3))
        write(io, *) minval(weights(1:3))
    end do
    close(io)

    call mpi_finalize(rc)
    
end program tsp_ga