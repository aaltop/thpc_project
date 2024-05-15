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
    do i = 1, 10000000
        call new_route(routes(:, 1))
        call calculate_total_distance(routes(:,1), distances, random_val(2))
        if (1 == i) then
            random_val(1) = random_val(2)
        else
            if (random_val(1) > random_val(2)) random_val(1) = random_val(2)
        end if

    end do

    print *, random_val(1)

    ! allocate(idx(5))
    ! call merge_argsort(distances(1:5, 1), idx)
    ! print *, distances(1:5, 1)
    ! print *, idx(5:1:-1)
    

    contains

    ! subroutine breed(route1, route2, child)
    !     implicit none

    !     integer(kind=int_kind), intent(in) :: route1(:), route2(size(route1))
    !     integer(kind=int_kind), intent(out) :: child(size(route1))

    !     integer(kind=int_kind) :: not_added(size(route1)), n, i

    !     n = size(route1)

    !     do i = 1, n
    !         not_added(i) = i
    !     end do



    ! end subroutine breed


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