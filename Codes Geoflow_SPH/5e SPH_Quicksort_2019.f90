module quicksort_2019
!
! this module uses quick sorting method to resort any array in one dimension 
! input :array, index of this array, number of elements in this array   
! output: resorted array, resorted index, number of elements in this array   
!    
! Chuan nov.2 2017   
!    
implicit none
public :: sort
private :: Partition

contains

recursive subroutine sort(A,b,npoin)
  real,    intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: b        !
  integer :: npoin    !
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, b, iq)
     call sort(A(:iq-1), b(:iq-1), npoin)
     call sort(A(iq:), b(iq:), npoin)
  endif
end subroutine sort

subroutine Partition(A, b, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: b    !
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp, temp2
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp  = A(i)
        A(i)  = A(j)
        A(j)  = temp
        temp2 = b(i)
        b(i)  = b(j)
        b(j)  = temp2
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module quicksort_2019