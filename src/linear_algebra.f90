module linear_algebra
!===============================================================================================
!  Module for numerical linear algebra routines
!===============================================================================================

!Global constants
use constants, only: &
   eps2, &
   i4, &
   r8, &
   r4

implicit none

contains 

!===============================================================================================
!   The following routines were taken from iModel  https://github.com/pedrospeixoto/iModel
!===============================================================================================

function cross_product(a,b)
    !-----------------------------------------------------------------------
    !  CROSS_PRODUCT
    !
    !  Returns the right-handed vector cross product of two 3-vectors:  
    !				C = A x B.
    !-----------------------------------------------------------------------
    implicit none

    real (r8), intent(in):: a(1:3)
    real (r8), intent(in):: b(1:3)
    real (r8):: cross_product(1:3)

    cross_product(1) = a(2)*b(3) - a(3)*b(2)                                    
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)

    return
end function cross_product

function det(p1, p2, p3)
    !-----------------------------------------------------------------------
    !  DET
    !
    !  Returns the determinant of the matrix made of the 3 points p1, p2, p3
    !   as columns
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: p1(1:3)
    real (r8), intent(in) :: p2(1:3)
    real (r8), intent(in) :: p3(1:3)
    real (r8):: det

    det=dot_product(cross_product(p1,p2),p3)

    return
end function det


function norm(p)
    !-----------------------------------------
    ! NORM
    ! Calculates the euclidian norm of a vector
    !----------------------------------------------
    real(r8), intent(in) :: p(:)
    real(r8):: norm

    norm=dot_product( p, p)
    norm=dsqrt(norm)

    return
end function norm

function error_norm_max_rel(f, g)
    !-------------------------------------------
    !Calculates the maximum absolute value of
    !  f-g csgrid array divided by the max of g
    !-------------------------------------------
    real (r8),  intent(in) :: f(:,:,:)
    real (r8),  intent(in) :: g(:,:,:)
    real (r8):: error_norm_max_rel
    real (r8):: maxfg
    real (r8):: maxg

    maxfg=maxval(abs(f-g))
    maxg=maxval(abs(g))

    if(abs(maxg) <= eps2 )then
       !print*, "error_norm_max_rel error: division by zero"
       maxg = 1._r8
    end if

    error_norm_max_rel=maxfg/maxg
    return
end function error_norm_max_rel

function error_norm_2_rel(f, g)
    !-------------------------------------------
    !Calculates relative L2 error
    !  That is, the square root of the
    !  the sum of the squares of (f-g) vector
    !  divided by the norm 2 of g
    !-------------------------------------------
    real (r8),  intent(in) :: f(:,:,:)
    real (r8),  intent(in) :: g(:,:,:)
    real (r8):: error_norm_2_rel

    real (r8):: sum_sq
    real (r8):: sum_sq_g

    sum_sq = sum((f-g)**2)
    sum_sq_g = sum(g**2)

    if(abs(sum_sq_g) <= eps2 )then
       !print*, "error_norm_2_rel error: division by zero"
       sum_sq_g = size(f)
    end if
    error_norm_2_rel=dsqrt(sum_sq/sum_sq_g)

    return
end function error_norm_2_rel

function error_norm_1_rel(f, g)
    !-------------------------------------------
    !Calculates the
    !  the sum of the absolute values of (f-g) vector
    !-------------------------------------------
    real (r8), intent(in) :: f(:,:,:)
    real (r8), intent(in) :: g(:,:,:)
    real (r8):: error_norm_1_rel

    real (r8):: sumfg
    real (r8):: sumg

    sumfg = sum(abs(f-g))
    sumg = sum(abs(g))
    if(abs(sumg) <= eps2 )then
       sumg = size(f)
    end if
    error_norm_1_rel = sumfg/sumg
    return
end function error_norm_1_rel

function norm_max(f)
    !-------------------------------------------
    ! max norm
    !-------------------------------------------
    real (r8),  intent(in) :: f(:,:,:)
    real (r8):: norm_max

    norm_max = maxval(abs(f))
    return
end function norm_max

function norm_2(f)
    !-------------------------------------------
    ! L2 norm
    !-------------------------------------------
    real (r8),  intent(in) :: f(:,:,:)
    real (r8):: norm_2

    norm_2 = dsqrt(sum((f)**2))
    return
end function norm_2

function norm_1(f)
    !-------------------------------------------
    ! L1 norm
    !-------------------------------------------
    real (r8), intent(in) :: f(:,:,:)
    real (r8):: norm_1

    norm_1 = sum(abs(f))
    return
end function norm_1


end module linear_algebra 
