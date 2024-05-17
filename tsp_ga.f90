program tsp_ga
    use globals
    use utils
    use tsp
    implicit none


    integer :: io, i, iostat, num_cities, num_considered, generations, &
        t0, t1, clock_rate, idx, ia, num_bred, num_candidates
    integer, allocatable :: routes(:,:)

    real(kind=real_kind) :: shortest_distance, mutation_chance
    real(kind=real_kind), allocatable :: random_val(:), weights(:), &
        distances(:,:)

    character(len=200):: distances_file, arg


    ! read in parameters from command line
    ! -------------------------------------
    ia=command_argument_count()
    if (5 /= ia) then
        call get_command_argument(0,arg)
        write(6,'(a,a,a)')   'usage: ',trim(arg),' distances_file num_candidates num_bred mutation_chance generations'
        write(6,'(a)')       '    distances_file = file to read city distances matrix from'
        write(6,'(a)')       '    num_candidates  = number of candidate parents in each generation'
        write(6,'(a)')       '    num_bred = number of children bred by the parents each generation'
        write(6,'(a)')       '    mutation_chance = chance of mutation'
        write(6,'(a)')       '    generations = number of generations to run the solver for'
        stop
    end if

    call get_command_argument(1,arg); read(arg,*) distances_file
    call get_command_argument(2,arg); read(arg,*) num_candidates
    call get_command_argument(3,arg); read(arg,*) num_bred
    call get_command_argument(4,arg); read(arg,*) mutation_chance
    call get_command_argument(5,arg); read(arg,*) generations
    ! =============================================

    ! make folder for generated data
    print "(a)", "Creating folder for data..."
    write(*, "(4x)", advance="no")
    call system("mkdir generated_data")
    print *

    call system_clock(t0, clock_rate)
    ! read in the array of distances form a file
    !
    ! the file should be formatted such that the first line contains
    ! the number of cities, and after that there is the array of
    ! cities' distances from each other
    ! ---------------------------------------------------------------
    io = 1234
    open(newunit=io, file=trim(distances_file), status="old", action="read", iostat=iostat)
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

    call find_optimal_route(distances(1:num_cities, 1:num_cities), num_candidates, num_bred, mutation_chance, generations, routes)

    call system_clock(t1)
    print '(a,g16.8,a)', 'Wall clock time: ',real(t1-t0,real_kind)/clock_rate,' seconds'

    ! Calculate the distances gained from the breeding algorithm
    io = 1234
    open(io, file="generated_data/breed.txt", status="replace", action="write")
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

    print "(a,g0,a)", "Best route:"
    print *, routes(:,idx)
    print "(a)", "Route distance:"
    print *, weights(idx)

    ! Do a random search
    io = 1234
    open(io, file="generated_data/random_search.txt", status="replace", action="write")
    do i = 1, generations
        call new_route(routes(:,1))
        call calculate_total_distance(routes(:,1), distances, weights(1))
        write(io, *) weights(1)
    end do
    close(io)

    ! Do a "random breed", whereby two random routes are bred
    ! and the child is the new tested route
    io = 1234
    open(io, file="generated_data/random_breed.txt", status="replace", action="write")
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
    open(io, file="generated_data/random_breed_better.txt", status="replace", action="write")
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
    
end program tsp_ga