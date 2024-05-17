program tsp_ga
    use mpi
    use globals
    use utils
    use tsp
    implicit none

    integer :: io, i, iostat, num_cities, num_considered, generations, &
        t0, t1, clock_rate, idx
    integer, allocatable :: routes(:,:)

    real(kind=real_kind) :: shortest_distance
    real(kind=real_kind), allocatable :: random_val(:), weights(:), &
        distances(:,:)

    character(len=80) :: file_name
    character(len=:), allocatable :: folder_command

    ! MPI variables
    ! ------------------------------------
    integer, parameter :: tag = 50

    integer :: id, ntasks, source_id, dest_id, rc, nlen, mpi_neighbor

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

    ! make folder for generated data
    if ( 0 == id ) then
        print "(a)", "Creating folder for data..."
        write(*, "(4x)", advance="no")
        folder_command = "mkdir generated_data"
        call system(folder_command)
        print *
    end if

    if ( 0 == id ) print "(a,i3)", "the number of processes is", ntasks
    call system_clock(t0, clock_rate)
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

    num_cities = 15
    generations = 10000

    ! allocate(random_val(num_cities))
    allocate(weights(generations))

    ! allocate(idx(generations))
    ! call shuffle(weights, idx)


    allocate(routes(num_cities, generations))

    call parallel_find_optimal_route( &
    distances(1:num_cities, 1:num_cities), &
    10, 15, real(0.95, real_kind), generations, 3, 5, routes)

    call system_clock(t1)
    print '(a,x,g0,x,g16.8,a)', 'Wall clock time for process', id, real(t1-t0,real_kind)/clock_rate, ' seconds'

    write(file_name,"(g0)") id
    file_name = "generated_data/parallel_breed" // trim(file_name) // ".txt"
    
    io = 1234 + id
    open(io, file=file_name, status="replace", action="write")
    do i = 1, generations
        
        call calculate_total_distance(routes(:,i), distances, weights(i))

        ! track best route
        if (1 == i) then
            shortest_distance = weights(i)
            idx = 1
        else
            if (weights(i) < shortest_distance) then
                shortest_distance = weights(i)
                idx = i
            end if
        end if

        write(io, *) weights(i)

    end do
    close(io)

    print "(a,g0,a)", "Best route in process ", id, ":"
    print *, routes(:,idx)
    print "(a)", "Route distance:"
    print *, weights(idx)

    call mpi_finalize(rc)

contains

! Attempts to find the optimal route in the Travelling Salesman Problem
! setting using a genetic algorithm.
!
! - distances(i,j) contains the distance between city "i" and city "j".
! It is assumed that the total number of cities is the size of the first
! dimension of this array.
!
! - <num_candidates> is the number of parents that can partake in breeding
! each generation.
!
! - <num_bred> is the number of bred children each generation.
!
! - <mutation_chance> is the probability of a child mutating.
!
! - <generations> is the number of generations that should pass.
! 
! - <num_migrators> is the number of migrating specimen each generation.
!
! - <migration_freq> tells how many generations should pass between
! migrations.
!
! - The optimal routes each generation. The number of rows should equal
! the number of rows of <distances>, while the number of columns should
! equal the number of generations.
subroutine parallel_find_optimal_route( &
    distances, num_candidates, num_bred, mutation_chance, generations, &
    num_migrators, migration_freq, optimal_route &
    )
    implicit none

    real(kind=real_kind), intent(in) :: distances(:,:), mutation_chance
    integer(kind=int_kind), intent(in) :: num_candidates, num_bred, generations, num_migrators, migration_freq
    integer(kind=int_kind), intent(out) :: optimal_route(size(distances, dim=1), generations)

    integer(kind=int_kind) :: & 
        possible_partners(2, num_candidates*(num_candidates-1)), &
        candidates(size(distances, dim=1), num_candidates), &
        children(size(distances, dim=1), num_bred), &
        idx(num_candidates*(num_candidates-1)), &
        migrators(size(distances, dim=1), num_migrators), &
        i, j, gen, n

    real(kind=real_kind) :: &
        weights(2,num_candidates*(num_candidates-1)), &
        random_val


    if ( num_bred > num_candidates*(num_candidates-1) ) then
        print *, "Number of children should not exceed the number of possible combinations of parents"
        stop
    end if

    if (num_migrators > num_candidates) then
        print *, "Number of migrators should not exceed the number of candidate parents"
        stop
    end if

    possible_partners = partner_permutations(num_candidates)

    ! get random routes, calculate distances
    do i = 1, num_candidates
            
        call new_route(candidates(:,i))
        call calculate_total_distance(candidates(:,i), distances, weights(1,i))

    end do
    call calculate_fitness(weights(1,:num_candidates))

    do gen = 1, generations
    
        ! calculate a combined fitness for partners by summing
        ! their individuals fitnesses
        do i = 1, size(possible_partners, dim=2)
            weights(2,i) = sum(weights(1, possible_partners(:,i)))
        end do

        ! Randomly choose indices of partners, weighted by
        ! their total fitness, then breed these partners to make
        ! children
        call shuffle(weights(2,:), idx)
        do i = 1, num_bred
            call breed( &
                candidates(:,possible_partners(1,idx(i))), &
                candidates(:, possible_partners(2, idx(i))), &
                distances, &
                children(:, i) &
            )


            ! mutation
            call random_number(random_val)
            if (random_val < mutation_chance) call mutate(children(:,i))

            ! keep track of best route
            call calculate_total_distance(children(:,i), distances, weights(1,i))
            if (1 == i) then
                optimal_route(:, gen) = children(:,1)
            else
                if (weights(1,i) < weights(1,i-1)) optimal_route(:,gen) = children(:,i)
            end if
        end do


        ! calculate fitness for children
        call calculate_fitness(weights(1,1:num_bred))
        ! cull randomly based on fitness
        call shuffle(weights(1,1:num_bred), idx(1:num_bred))
        candidates(:,:) = children(:,idx(1:num_candidates))

        ! send migrators to neighboring process, accept from other
        ! neighboring process
        if (1 < ntasks .and. 0 == modulo(gen, migration_freq)) then
            migrators(:,:) = candidates(:, 1:num_migrators)

            if ( 0 == id ) then
                call mpi_recv( &
                candidates(:, 1:num_migrators), &
                size(candidates(:, 1:num_migrators)), &
                MPI_INTEGER, modulo(id-1,ntasks), tag, &
                MPI_COMM_WORLD, status, rc&
                )

                call mpi_send( &
                migrators, size(migrators), MPI_INTEGER, &
                modulo(id+1,ntasks), tag, MPI_COMM_WORLD, rc &
                )
            else
                call mpi_send( &
                migrators, size(migrators), MPI_INTEGER, &
                modulo(id+1,ntasks), tag, MPI_COMM_WORLD, rc &
                )

                call mpi_recv( &
                candidates(:, 1:num_migrators), &
                size(candidates(:, 1:num_migrators)), &
                MPI_INTEGER, modulo(id-1,ntasks), tag, &
                MPI_COMM_WORLD, status, rc&
                )
            end if

            ! recalculate the weights (could also send, but seems easier
            ! this way)
            do i = 1, num_candidates
                call calculate_total_distance(candidates(:,i), distances, weights(1,i))
            end do
            call calculate_fitness(weights(1, 1:num_candidates))

            ! shuffle again to get new indices (maybe not the best,
            ! but its also what's easiest right now)
            call shuffle(weights(1,1:num_candidates), idx(1:num_candidates))
            candidates(:,:) = children(:,idx(1:num_candidates))
        end if
        ! set the living children's (new parents') fitnesses
        weights(1,1:num_candidates) = weights(1,idx(1:num_candidates))

    end do


end subroutine
    
end program tsp_ga