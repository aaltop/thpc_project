module tsp
    use globals
    use utils
    implicit none
    
    
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
! - <optimal_route> will contain the optimal routes each generation. The number of rows should equal
! the number of rows of <distances>, while the number of columns should
! equal the number of generations.
!
! - <route_lengths> will contain the lengths of the corresponding routes in
! <optimal_route>, and should therefore be of length <generations>.
subroutine find_optimal_route(distances, num_candidates, num_bred, mutation_chance, generations, optimal_route, route_lengths)
    implicit none

    real(kind=real_kind), intent(in) :: distances(:,:), mutation_chance
    integer(kind=int_kind), intent(in) :: num_candidates, num_bred, generations
    integer(kind=int_kind), intent(out) :: optimal_route(size(distances, dim=1), generations)
    real(kind=real_kind), intent(out) :: route_lengths(generations)

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

    if (num_bred < num_candidates) then
        print *, "Number of children should be greater than or equal to the number of candidate parents"
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
                route_lengths(gen) = weights(1,1)
            else
                if (weights(1,i) < weights(1,i-1)) then
                    optimal_route(:,gen) = children(:,i)
                    route_lengths(gen) = weights(1,i)
                end if
            end if
        end do


        ! calculate fitness for children
        call calculate_fitness(weights(1,1:num_bred))
        ! cull randomly based on fitness
        call shuffle(weights(1,1:num_bred), idx(1:num_bred))
        candidates(:,:) = children(:,idx(1:num_candidates))
        ! set the living children's (new parents') fitnesses
        weights(1,1:num_candidates) = weights(1,idx(1:num_candidates))

    end do


end subroutine

subroutine calculate_fitness(weights)
    implicit none

    real(kind=real_kind), intent(inout) :: weights(:)

    weights = 1 - weights/(maxval(weights)+1)

end subroutine calculate_fitness

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
        ! don't pair with self
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

    ! pointless to switch in these cases
    if ( 1 == size(route) .or. 2 == size(route) .or. 3 == size(route) ) return

    ! NOTE: practically speaking, calling the shuffle is currently quite
    ! expensive compared to what it is actually needed for here,
    ! which is a uniform sample of two indices. Still, it does the
    ! job.
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
    
end module tsp