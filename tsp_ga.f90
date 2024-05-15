program tsp_ga
    implicit none

    integer, parameter :: real_kind=8
    integer, parameter :: int_kind=4

    integer :: io, i, iostat, num_cities, num_considered
    integer, allocatable :: idx(:), routes(:,:)

    real(kind=real_kind), allocatable :: random_val(:), weights(:), &
        distances(:,:)

    ! read in the array of distances form a file
    !
    ! the file should be formatted such that the first line contains
    ! the number of cities, and after that there is the array of
    ! cities' distances from each other
    ! ---------------------------------------------------------------
    open(newunit=io, file="wg59_dist_array.txt", status="old", action="read")
    read(io, *) num_cities
    allocate(distances(num_cities, num_cities))

    iostat = 0
    do i = 1,num_cities
        read(io, *, iostat=iostat) distances(1:num_cities, i)
        if (0 /= iostat) then
            print *, "something went wrong"
            exit
        end if
    end do
    close(io)
    ! ====================================================================

    allocate(random_val(num_cities))
    allocate(weights(num_cities))
    weights = 0.0
    weights(1) = 1.0

    allocate(idx(num_cities))
    call shuffle(weights, idx)

    num_considered = 10

    allocate(routes(11, num_considered))
    do i = 1, 100
        call new_route(routes(:, 1))
        call new_route(routes(:, 2))
        call breed(routes(:,1), routes(:,2), distances, routes(:,3))

        call calculate_total_distance(routes(:,1), distances, random_val(1))
        call calculate_total_distance(routes(:,2), distances, random_val(2))
        call calculate_total_distance(routes(:,3), distances, random_val(3))
        print "(3f5.0)", random_val(1:3)

    end do

    ! allocate(idx(5))
    ! call merge_argsort(distances(1:5, 1), idx)
    ! print *, distances(1:5, 1)
    ! print *, idx(5:1:-1)
    

    contains

    ! Breed <route1> and <route2> to create a new route <child>.
    ! <distances> contain distances between cities, with distances(i,j)
    ! being the distance being city "i" and city "j".
    subroutine breed(route1, route2, distances, child)
        implicit none

        integer(kind=int_kind), intent(in) :: route1(:), route2(size(route1))
        integer(kind=int_kind), intent(out) :: child(size(route1))
        real(kind=real_kind), intent(in) :: distances(:,:)

        integer(kind=int_kind) :: not_added(size(route1)-1), n, i, idx
        real(kind=real_kind) :: random_val

        n = size(route1)

        child = 0
        child(1) = 1
        do i = 2, n
            not_added(i-1) = i
        end do
        
        do i = 2, n

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
            idx = ceiling(random_val*(n+1-i))
            child(i) = not_added(idx)
            ! no need to move things around if the sampled value was
            ! the last in the not_added array
            if ( idx == n+1-i) cycle

            ! if the sampled value was not last, move all the ones
            ! that come after it one space to the "left"
            not_added(idx:n-i) = not_added(idx+1:n+1-i)

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