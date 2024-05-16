program tsp_ga
    implicit none

    integer, parameter :: real_kind=8
    integer, parameter :: int_kind=4

    integer :: io, i, iostat, num_cities, num_considered, generations
    integer, allocatable :: idx(:), routes(:,:)

    real(kind=real_kind), allocatable :: random_val(:), weights(:), &
        distances(:,:)

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


    io = 1234
    open(io, file="breed.txt", status="replace", action="write")
    do i = 1, generations
        
        print "(11i3)", routes(:,i)
        call calculate_total_distance(routes(:,i), distances, weights(i))
        write(io, *) weights(i)

    end do
    close(io)

    io = 1234
    open(io, file="random_search.txt", status="replace", action="write")
    do i = 1, generations
        call new_route(routes(:,1))
        call calculate_total_distance(routes(:,1), distances, weights(1))
        write(io, *) weights(1)
    end do
    close(io)

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


    ! allocate(idx(5))
    ! call merge_argsort(distances(1:5, 1), idx)
    ! print *, distances(1:5, 1)
    ! print *, idx(5:1:-1)
    

    contains

    subroutine find_optimal_route(distances, num_candidates, num_bred, mutation_chance, generations, optimal_route)
        implicit none

        real(kind=real_kind), intent(in) :: distances(:,:), mutation_chance
        integer(kind=int_kind), intent(in) :: num_candidates, num_bred, generations
        integer(kind=int_kind), intent(out) :: optimal_route(size(distances, dim=1), generations)

        integer(kind=int_kind) :: & 
            possible_partners(2, num_candidates*(num_candidates-1)), &
            candidates(size(distances, dim=1), num_candidates), &
            children(size(distances, dim=1), num_bred), &
            idx(num_candidates*(num_candidates-1)), &
            i, j, gen, n

        real(kind=real_kind) :: &
            weights(2,num_candidates*(num_candidates-1)), &
            random_val


        if ( num_bred > num_candidates*(num_candidates-1) ) then
            print *, "Number of children should not exceed the number of possible combinations of parents"
            stop
        end if

        possible_partners = partner_permutations(num_candidates)

        ! get random routes, calculate distances
        do i = 1, num_candidates
                
            call new_route(candidates(:,i))
            call calculate_total_distance(candidates(:,i), distances, weights(1,i))

        end do


        ! As "fitness", use this
        ! weights(1,:num_candidates) = minval(weights(1,:num_candidates))/weights(1,:num_candidates)
        weights(1,:num_candidates) = 1 - weights(1,:num_candidates)/maxval(weights(1,:num_candidates) + 1)


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

            ! print "(11i3)", children
            ! print "(f5.0)", weights(1,1:num_bred)
            ! stop


            ! calculate fitness for children
            ! print "(f10.3)", weights(1,1:num_bred)
            ! weights(1,1:num_bred) = minval(weights(1,1:num_bred))/weights(1,1:num_bred)
            weights(1,1:num_bred) = 1 - weights(1,1:num_bred)/(maxval(weights(1,1:num_bred))+1)
            weights(1,1:num_bred) = weights(1,1:num_bred)/sum(weights(1,1:num_bred))
            ! print *, sum(weights(1,1:num_bred))/size(weights(1,1:num_bred)), maxval(weights(1,1:num_bred))
            ! print "(f6.3)", weights(1,1:num_bred)
            ! cull randomly based on fitness
            call shuffle(weights(1,1:num_bred), idx(1:num_bred))
            ! print *
            ! print "(i3)", idx(1:num_bred)
            ! print *
            ! print "(f6.2)", weights(1,1:num_bred) - weights(1,idx(1:num_bred))
            ! stop
            candidates(:,:) = children(:,idx(1:num_candidates))
            ! set the living children's (new parents') fitnesses
            weights(1,1:num_candidates) = weights(1,idx(1:num_candidates))
            ! print *, sum(weights(1,1:num_candidates))/size(weights(1,1:num_candidates)), maxval(weights(1,1:num_candidates))

        end do


    end subroutine

    ! for <num_candidates>, calculates all possible permutations
    ! of two partners, such that no candidate is paired with itself.
    function partner_permutations(num_candidates)
        implicit none

        integer(kind=int_kind), intent(in) :: num_candidates
        integer(kind=int_kind) :: &
            partner_permutations(2, num_candidates*(num_candidates-1)), &
            idx, i, n
        i = 0
        idx = 1
        do while (idx <= num_candidates*(num_candidates-1))

            n = modulo(i, num_candidates)+1
            if (n == (i/num_candidates+1)) then
                i = i + 1
                cycle
            end if

            partner_permutations(1,idx) = i/num_candidates + 1
            partner_permutations(2,idx) = n

            i = i + 1
            idx = idx + 1
        end do

    end function partner_permutations

    ! switch the places of two cities on <route>.
    subroutine mutate(route)
        implicit none

        integer(kind=int_kind), intent(inout) :: route(:)

        integer(kind=int_kind) :: idx(size(route)-1)
        real(kind=real_kind) :: weights(size(route)-1)

        if ( 1 == size(route) .or. 2 == size(route) .or. 3 == size(route) ) return

        weights = 1.0
        call shuffle(weights, idx)
        ! always keep one (1) as first element
        idx = idx + 1
        idx(3) = route(idx(1))
        route(idx(1)) = route(idx(2))
        route(idx(2)) = idx(3)

    end subroutine

    ! Breed <route1> and <route2> to create a new route <child>.
    ! <distances> contain distances between cities, with distances(i,j)
    ! being the distance being city "i" and city "j".
    subroutine breed(route1, route2, distances, child)
        implicit none

        integer(kind=int_kind), intent(in) :: route1(:), route2(size(route1))
        integer(kind=int_kind), intent(out) :: child(size(route1))
        real(kind=real_kind), intent(in) :: distances(:,:)

        integer(kind=int_kind) :: not_added(size(route1)-1), n, i, idx(1)
        real(kind=real_kind) :: random_val

        n = size(route1)

        child = 0
        child(1) = 1
        do i = 2, n
            not_added(i-1) = i
        end do
        
        do i = 2, n


            ! keep track of cities not yet added
            if ( 2 < i ) then

                idx = findloc(not_added, child(i-1))

                if ( idx(1) == n+2-i) then
                    ! no need to move things around if the previously added value was
                    ! the last in the not_added array
                    continue
                else
                    ! if the added value was not last, move all the ones
                    ! that come after it one space to the "left"
                    not_added(idx(1):n+1-i) = not_added(idx(1)+1:n+2-i)
                end if
                
            end if

            ! add a new city on from the parents based on which parents'
            ! city is closer to the current last child city.
            if ( distances(child(i-1), route1(i)) < distances(child(i-1), route2(i)) ) then
                if ( .not. any(child == route1(i)) ) then
                    child(i) = route1(i)
                    cycle
                end if
            else if ( .not. any(child == route2(i)) ) then
                child(i) = route2(i)
                cycle
            end if

            ! if a suitable candidate is not found in either parent,
            ! randomly sample from the values that have not been added
            ! yet
            call random_number(random_val)
            idx(1) = ceiling(random_val*(n+1-i))
            child(i) = not_added(idx(1))

        end do
    end subroutine breed


    ! calculate the round-trip distance of <route>, where the elements
    ! of <route> match the indices of <distances>, and the element
    ! distances(i,j) contains the distance between city "i" and city "j".
    !
    ! Return distance in <dist>.
    subroutine calculate_total_distance(route, distances, dist)
        implicit none

        integer(kind=int_kind), intent(in) :: route(:)
        real(kind=real_kind), intent(in) :: distances(:,:)
        real(kind=real_kind), intent(out) :: dist

        integer(kind=int_kind) :: n, i

        n = size(route)

        dist = 0.0
        do i = 1, n
            dist = dist + distances(route(i),route(modulo(i,n)+1))
        end do

    end subroutine calculate_total_distance


    ! Returns a new route in <route>. Because the order does not
    ! matter for the cities, the city "1" is always first.
    subroutine new_route(route)
        implicit none
        integer(kind=int_kind), intent(inout) :: route(:)

        integer(kind=int_kind) :: num_cities
        real(kind=real_kind) :: weights(size(route)-1)

        num_cities = size(route)

        route = 0
        weights = 1.0
        ! this basically means that city "1" is always first
        call shuffle(weights, route(2:num_cities))
        route = route + 1

    end subroutine new_route

    ! Get shuffled indices, emphasise picking based on <weights>.
    ! Return shuffled indices <idx> corresponding to <weights>.
    subroutine shuffle(weights, idx)
        implicit none

        real(kind=real_kind), intent(in), dimension(:) :: weights

        integer(kind=int_kind), intent(out) :: idx(size(weights))

        real(kind=real_kind) :: random_val(size(weights))
        integer(kind=int_kind) :: i, idx_nonzero(size(weights)), nonzero_weights, &
            idx_zero(size(weights)), zero_weights, n

        n = size(weights)

        ! find the weights' indices that are non-zero and zero.
        ! (Would probably be more sensible to just not allow non-positive
        ! weights.)
        ! -------------------------------------
        nonzero_weights = 0
        zero_weights = 0
        do i = 1, n
            
            if (0 /= weights(i)) then
                idx_nonzero(nonzero_weights+1) = i
                nonzero_weights = nonzero_weights + 1
            else
                idx_zero(zero_weights+1) = i
                zero_weights = zero_weights + 1
            end if
                
        end do
        ! =============================================

        ! sample based on okay weights
        call random_number(random_val)
        call merge_argsort( &
        1/weights(idx_nonzero(1:nonzero_weights))*log(random_val(1:nonzero_weights)), & ! see Weighted Random Sampling (2005; Efraimidis, Spirakis)
        idx(1:nonzero_weights) &
        )

        ! sample based on zero weights (uniform sample for all zero weights)
        call merge_argsort(random_val(nonzero_weights+1:n), idx(nonzero_weights+1:n))
        ! actually sort zero indices based on the sample order
        do i = 0, zero_weights-1
            idx(n-i) = idx_zero(idx(n-i))
        end do




    end subroutine shuffle

    ! from https://github.com/Astrokiwi/simple_fortran_argsort
    ! sorts in descending order
    subroutine merge_argsort(r,d)
        implicit none

        real(kind=real_kind), intent(in), dimension(:) :: r
        integer(kind=int_kind), intent(out), dimension(size(r)) :: d
      
        integer(kind=int_kind), dimension(size(r)) :: il

        integer(kind=int_kind) :: stepsize
        integer(kind=int_kind) :: i,j,n,left,k,ksize
      
        n = size(r)
      
        do i=1,n
            d(i)=i
        end do
      
        if ( n==1 ) return
      
        stepsize = 1
        do while (stepsize<n)
            do left=1,n-stepsize,stepsize*2
                i = left
                j = left+stepsize
                ksize = min(stepsize*2,n-left+1)
                k=1
          
                do while ( i<left+stepsize .and. j<left+ksize )
                    if ( r(d(i))>r(d(j)) ) then
                        il(k)=d(i)
                        i=i+1
                        k=k+1
                    else
                        il(k)=d(j)
                        j=j+1
                        k=k+1
                    endif
                enddo
          
                if ( i<left+stepsize ) then
                    ! fill up remaining from left
                    il(k:ksize) = d(i:left+stepsize-1)
                else
                    ! fill up remaining from right
                    il(k:ksize) = d(j:left+ksize-1)
                endif
                d(left:left+ksize-1) = il(1:ksize)
            end do
            stepsize=stepsize*2
        end do

        return
    end subroutine


    
end program tsp_ga