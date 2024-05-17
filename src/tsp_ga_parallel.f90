program tsp_ga
    use mpi
    use globals
    use utils
    use tsp
    implicit none

    integer :: io, i, iostat, num_cities, num_considered, generations, &
        t0, t1, clock_rate, idx, ia, num_bred, num_candidates, &
        num_migrators, migration_freq
    integer, allocatable :: routes(:,:)

    real(kind=real_kind) :: shortest_distance, mutation_chance
    real(kind=real_kind), allocatable :: distances(:,:), weights(:)

    character(len=80) :: file_name
    character(len=200) :: distances_file, arg

    ! MPI variables
    ! ------------------------------------
    integer, parameter :: tag = 50

    integer :: id, ntasks, nlen, rc

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


    ! distances(1:num_cities, 1:num_cities), &
    ! 10, 20, real(0.95, real_kind), generations, 2, 5, routes

    ! read in parameters from command line
    ! -------------------------------------
    ia=command_argument_count()
    if (7 /= ia) then
        call get_command_argument(0,arg)
        write(6,'(a,a,a)')   'usage: ',trim(arg),&
            ' distances_file num_candidates num_bred mutation_chance generations num_migrators migration_freq'

        write(6,'(a)')       '    distances_file = file to read city distances matrix from'
        write(6,'(a)')       '    num_candidates  = number of candidate parents in each generation'
        write(6,'(a)')       '    num_bred = number of children bred by the parents each generation'
        write(6,'(a)')       '    mutation_chance = chance of mutation'
        write(6,'(a)')       '    generations = number of generations to run the solver for'
        write(6,'(a)')       '    num_migrators = how many should migrate'
        write(6,'(a)')       '    migration_freq = how many generations should pass between migrations'
        stop
    end if

    call get_command_argument(1,arg); read(arg,*) distances_file
    call get_command_argument(2,arg); read(arg,*) num_candidates
    call get_command_argument(3,arg); read(arg,*) num_bred
    call get_command_argument(4,arg); read(arg,*) mutation_chance
    call get_command_argument(5,arg); read(arg,*) generations
    call get_command_argument(6,arg); read(arg,*) num_migrators
    call get_command_argument(7,arg); read(arg,*) migration_freq
    ! =============================================

    ! make folder for generated data
    if ( 0 == id ) then
        write(*, "(4x)", advance="no")
        call system("mkdir generated_data 2> /dev/null")
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
    do i = 1, num_cities
        read(io, *, iostat=iostat) distances(1:num_cities, i)
        if (0 /= iostat) then
            print *, "something went wrong"
            stop
        end if
    end do
    close(io)
    ! ====================================================================


    ! TSP solve with migration
    ! ---------------------------------
    generations = 10000
    allocate(routes(num_cities, generations))
    allocate(weights(generations))

    call breed_statistics(10)
    call mpi_finalize(rc)
    stop

    call parallel_find_optimal_route( &
    distances, &
    num_candidates, num_bred, mutation_chance, generations, num_migrators, migration_freq, routes)

    call system_clock(t1)
    print '(a,x,g0,x,g16.8,a)', 'Wall clock time for process', id, real(t1-t0,real_kind)/clock_rate, ' seconds'

    ! write distances to file
    write(file_name,"(g0)") id
    file_name = "generated_data/parallel_breed1_" // trim(file_name) // ".txt"
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
        write(io, *) weights(i),  routes(:, i)

    end do
    close(io)

    ! write(*, "(a,g0,a)", advance="no") "Best route in process ", id, ": "
    print *, "Best route of distance", weights(idx), "in process", id, ":", routes(:,idx)
    ! ==========================================


    ! TSP solve without migration
    ! ----------------------------------------
    call find_optimal_route( &
    distances(1:num_cities, 1:num_cities), &
    num_candidates, num_bred, mutation_chance, generations, routes)

    ! write distances to file
    write(file_name,"(g0)") id
    file_name = "generated_data/parallel_breed2_" // trim(file_name) // ".txt"
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

        write(io, *) weights(i), routes(:,i)

    end do
    close(io)

    ! =============================================

    call mpi_finalize(rc)

contains

! For <repeat> runs of the TSP algorithm, collect stats of minimum- ,
! mean- , maximum distance, and time, and print out the (minimum) mean of these stats
subroutine breed_statistics(repeat)
    implicit none

    integer(kind=int_kind), intent(in) :: repeat

    integer(kind=int_kind) :: i, j

    real(kind=real_kind) :: stats(4,repeat), recv_stats(4,repeat)

    ! collect stats
    do i = 1, repeat
        if (0 == id) print "(a,g0,a,g0)", "round ", i, " out of ", repeat
        call system_clock(t0, clock_rate)
        call parallel_find_optimal_route( &
        distances, &
        num_candidates, num_bred, mutation_chance, generations, &
        num_migrators, migration_freq, routes)
        call system_clock(t1)
        ! time
        stats(4,i) = real(t1-t0,real_kind)/clock_rate

        do j = 1, generations
            call calculate_total_distance(routes(:, j), distances, weights(j))
        end do
        ! min, mean, max
        stats(1,i) = minval(weights)
        stats(2,i) = sum(weights)/size(weights)
        stats(3,i) = maxval(weights)
        
    end do


    ! collect the stats from each process
    if ( 1 < ntasks ) then
        
        if ( 0 == id ) then
            
            do i = 1, ntasks-1
                
                call mpi_recv( &
                recv_stats, &
                size(recv_stats), &
                MPI_REAL, i, tag, &
                MPI_COMM_WORLD, status, rc&
                )

                ! the best distance metrics
                where (recv_stats(1:3,:) < stats(1:3,:))
                    stats(1:3,:) = recv_stats(1:3,:)
                end where

                ! the worst times
                where (recv_stats(4,:) > stats(4,:))
                    stats(4,:) = recv_stats(4,:)
                end where

            end do

            print *, "min", sum(stats(1,:))/size(stats(1,:))
            print *, "mean", sum(stats(2,:))/size(stats(2,:))
            print *, "max", sum(stats(3,:))/size(stats(3,:))
            print *, "time", sum(stats(4,:))/size(stats(4,:))

        else

            call mpi_send( &
            stats, size(stats), MPI_REAL, &
            0, tag, MPI_COMM_WORLD, rc &
            )

        end if

    end if

end subroutine breed_statistics


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