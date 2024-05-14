program tsp_ga
    implicit none

    integer, parameter :: real_kind=8
    integer, parameter :: int_kind=4

    integer :: io, i, iostat, num_cities
    integer, allocatable :: idx(:)

    real(kind=real_kind), allocatable :: random_val(:), weights(:)
    real(kind=real_kind), allocatable :: distances(:,:)

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
    weights = 1.0
    weights(1) = 10.0
    weights = weights/sum(weights)

    allocate(idx(num_cities))
    call sample(distances(1:num_cities,1), weights, idx)
    print *, idx(1:5)

    ! allocate(idx(5))
    ! call merge_argsort(distances(1:5, 1), idx)
    ! print *, distances(1:5, 1)
    ! print *, idx(5:1:-1)
    

    contains

    ! sample <num_to_pick> values for <values> randomly, emphasise
    ! picking based on <weights> which should correspond to <values>.
    ! Return indices corresponding to <values>.
    !
    ! Based on Weighted Random Sampling (2005; Efraimidis, Spirakis)
    subroutine sample(values, weights, idx)
        implicit none

        real(kind=real_kind), intent(in), dimension(:) :: values
        real(kind=real_kind), intent(in), dimension(size(values)) :: weights

        integer(kind=int_kind), intent(out) :: idx(size(values))

        real(kind=real_kind) :: random_val(size(values)), &
            keys(size(values))


        call random_number(random_val)
        ! values in weights could be zero, but at least my compiler
        ! does allow this by using infinities, and frankly, I'm fine
        ! with that: it results ultimately in the correct behaviour,
        ! and having to bypass the division by zero just to end
        ! up in the same (correct) behaviour seems stupid. Just sort your
        ! compilers out or use a more sensible language.
        keys = 1/weights*log(random_val)
        call merge_argsort(keys, idx)




    end subroutine sample

    ! from https://github.com/Astrokiwi/simple_fortran_argsort
    ! sorts in descending order
    subroutine merge_argsort(r,d)
        implicit none

        real(kind=real_kind), intent(in), dimension(:) :: r
        integer, intent(out), dimension(size(r)) :: d
      
        integer, dimension(size(r)) :: il

        integer :: stepsize
        integer :: i,j,n,left,k,ksize
      
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