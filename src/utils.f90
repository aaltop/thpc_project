module utils
    use globals
    implicit none
    
contains

! Get shuffled indices, emphasise picking based on <weights>.
! Return shuffled indices <idx> corresponding to <weights>.
subroutine shuffle(weights, idx)
    implicit none

    real(kind=real_kind), intent(in), dimension(:) :: weights

    integer(kind=int_kind), intent(out) :: idx(size(weights))

    real(kind=real_kind) :: random_val(size(weights))
    integer(kind=int_kind) :: i, n

    if ( any(weights <= 0) ) then
        print *, "Value error in shuffle: Weights should be positive"
        stop
    end if

    n = size(weights)

    ! sample based on okay weights
    call random_number(random_val)
    call merge_argsort( &
    1/weights*log(random_val), & ! see Weighted Random Sampling (2005; Efraimidis, Spirakis)
    idx &
    )

end subroutine shuffle

! ==================================================================
! Licensed under the LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.html
! from https://github.com/Astrokiwi/simple_fortran_argsort (David Williamson)
! ==================================================================
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
end subroutine merge_argsort
    
end module utils