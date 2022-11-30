module linear_algebra
  !===============================================================================================
  !  Module for numerical linear algebra routines
  !===============================================================================================

  !Global constants
  use constants, only: &
       eps, &
       i4, &
       r8, &
       r4, &
       r16

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

  function robdet(p1, p2, p3)
    !-----------------------------------------------------------------------
    !  ROBDET
    !
    !  Returns the robust determinant of the matrix made of the 3 points p1, p2, p3
    !   as columns - The robust part is to ensure that the sign is correct if
    !   it is near zero.
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: p1(1:3)
    real (r8), intent(in) :: p2(1:3)
    real (r8), intent(in) :: p3(1:3)
    real (r8):: robdet

    real (r16):: a(1:6)
    real (r16):: robdettmp
    real(r16), dimension(1:3) :: q1
    real(r16), dimension(1:3) :: q2
    real(r16), dimension(1:3) :: q3

    q1=p1
    q2=p2
    q3=p3
    a(1)=q1(1)*q2(2)*q3(3)
    a(2)=q2(1)*q3(2)*q1(3)
    a(3)=q3(1)*q1(2)*q2(3)
    a(4)=-q3(1)*q2(2)*q1(3)
    a(5)=-q1(1)*q3(2)*q2(3)
    a(6)=-q2(1)*q1(2)*q3(3)

    robdettmp=a(1)+a(2)+a(3)+a(4)+a(5)+a(6)
    robdet=real(robdettmp, r8)

    return
  end function robdet

  function solve3x3(l1, l2, l3, b)
    !-----------------------------------------------------------------------
    !  3x3 linear system solution
    !  li = line i
    !  b - right hand side vector
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: l1(1:3)
    real (r8), intent(in) :: l2(1:3)
    real (r8), intent(in) :: l3(1:3)
    real (r8), intent(in) :: b(1:3)
    real (r8):: solve3x3(1:3) !Solution
    real (r8):: detdiv !Matrix determinant
    real (r8):: v(1:3) !Auxiliar vector

    detdiv=det(l1,l2,l3)

    !Check for null determinant
    if(abs(detdiv)<eps)then
       print*, "    solve3x3 error: null determinant"
       stop
    end if

    v(1) = (l1(2) * l2(3) * b(3) - l1(2) * b(2) * l3(3) + l1(3) * l3(2) * &
         b(2) - l1(3) * b(3) * l2(2) + b(1) * l2(2) * l3(3) - b(1) * l2(3) * l3(2))
    v(2) = (-l1(1) * l2(3) * b(3) + l1(1) * b(2) * l3(3) + l2(1) * l1(3) * &
         b(3) - l2(1) * b(1) * l3(3) + l2(3) * l3(1) * b(1) - b(2) * l3(1) * l1(3))
    v(3) =  (l3(1) * l1(2) * b(2) + l3(2) * l2(1) * b(1) - l3(2) * l1(1) * b(2) - &
         l3(1) * b(1) * l2(2) - b(3) * l2(1) * l1(2) + b(3) * l1(1) * l2(2))

    solve3x3=v/detdiv

    return
  end function solve3x3

  function solvelintri(a, n)
    !---------------------------------------------
    !Solve a nxn Triangular linear system with weighted columns
    !The A matrix be of format (n+1) x (n+1)
    !  X 0 0 0 0
    !  X X 0 0 0
    !  X X X 0 0
    !  X X X X 0
    !  Y Y Y Y Y   -> Right Hand Side
    !-----------------------------------------------
    !Matrix size
    integer (i4):: n

    !Augmented Matrix for the least square problem
    real (r8), intent(in) :: a(1:n+1,1:n+1)

    !Solution
    real (r8):: solvelintri(1:n)

    !Aux vars
    real (r8), allocatable :: x(:)
    real (r8):: tmp
    integer:: i
    integer:: j

    allocate(x(1:n))

    do i=n, 1, -1
       tmp=0
       do j=n, i+1, -1
          tmp=tmp+a(j,i)*x(j)
       end do
       x(i)=(a(n+1,i)-tmp)/a(i,i)
    end do
    solvelintri(1:n)=x(1:n)

    return
  end function solvelintri

  function rot_point (p, theta)
    !-------------------------------------------------------
    !  ROT_POINT
    !
    !   This subroutine applies a rotation 'gimble like'
    ! around x, y, z-axis respectively using angles theta(x, y, z),
    ! of the point p=(x, y, z)
    !
    ! On input:
    !       p=(x,y,z) = coordinates of a point on the unit sphere.
    !       theta_? = angles of rotation in radians
    !
    ! On output:
    !       pr=(xr,yr,zr) = rotated point
    !---------------------------------------------------------
    real (r8):: p(1:3)
    real (r8):: theta(1:3)
    real (r8):: pr(1:3, 1)
    real (r8):: rot_point(1:3)
    real(r8):: R(1:3, 1:3)

    !Save rotated point
    pr(1:3,1)=p(1:3)

    !Rotate around x axis
    R(1,1)=1; R(1,2)=0;             R(1,3) =0
    R(2,1)=0; R(2,2)=dcos(theta(1)); R(2,3) = -dsin(theta(1))
    R(3,1)=0; R(3,2)=dsin(theta(1)); R(3,3) = dcos(theta(1))
    pr=matmul(R, pr)

    !Rotate around y axis
    R(1,1)=dcos(theta(2));  R(1,2)=0; R(1,3) = dsin(theta(2))
    R(2,1)=0;              R(2,2)=1; R(2,3) =0
    R(3,1)=-dsin(theta(2)); R(3,2)=0; R(3,3) = dcos(theta(2))
    pr=matmul(R, pr)

    !Rotate around z axis
    R(1,1)=dcos(theta(3));  R(1,2)=-dsin(theta(3)); R(1,3) = 0
    R(2,1)=dsin(theta(3)); R(2,2)=dcos(theta(3)); R(2,3) = 0
    R(3,1)=0;              R(3,2)=0; R(3,3) = 1
    pr=matmul(R, pr)

    rot_point(1:3)=pr(1:3, 1)
    return
  end function rot_point

  subroutine choleskydecomp(A, n, L)
    !-----------------------------------------
    ! Cholesky Decomposition (A=L D L')
    !    A is assumed to be symetric
    !    A must be positive definite
    !    A needs only values on lower triangular part (i>=j)
    !    D is diagonal
    !    L is lower triangular with ones in diagonal
    ! As output, D is places in the diagonal os L, only for
    !   storing purposes
    !
    ! Wikipedia version
    !------------------------------------------
    !Dimension of the matrix
    integer(i4), intent(in) :: n

    !Matrix to be decomposed
    real(r8), intent(inout) :: A(1:n,1:n)

    !Decomposed Matrix
    !  Contains lower triangular values and
    !   on its diagonal the D matriz values
    real(r8), intent(out) :: L(1:n,1:n)

    !Diagonal Matrix
    real(r8):: D(1:n)
    real(r8):: D2(1:n,1:n)
    real(r8):: D3(1:n,1:n)
    real(r8)::  T(1:n,1:n)

    !indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k

    !Check ?
    logical:: check
    real(r8):: error
    character(len=16):: fmt

    !Do the LDL' decomp
    L(1:n,1:n)=0._r8
    D(1:n)=0._r8
    do j=1,n
       D(j)=A(j,j)
       do k=1,j-1
          D(j)=D(j)-D(k)*L(j,k)**2
       end do
       if(D(j)<-eps/1000)then
          !print*, "Warningcholeskydecomp:"
          !print*, "   Matrix not positive definite, diagonal up to now="
          !print*, D(1:j)
          !print*
          !elseif(abs(D(j))<eps/1000000)then
       elseif(abs(D(j))==0)then
          !print*, "Error on choleskydecomp:"
          !print*, "   Matrix not positive definite, diagonal with zero"
          !print*, D(1:j)
          !print*
          !stop
       end if

       L(j,j)=1._r8
       do i=j+1,n
          L(i,j)=A(i,j)
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
          end do
          if(abs(D(j))==0)then
             L(i,j)=0
          else
             L(i,j)=L(i,j)/D(j)
          end if
       end do
    end do

    !Debugging checks
    check=.false.
    if(check)then
       !Fill in symetric part of the matrix A (received blank)
       do i=1, n
          do j=i+1, n
             A(i,j)=A(j,i)
          end do
       end do
       fmt="(xxf16.8)"
       write(fmt(2:3),'(i2)') n
       print*, "A"
       print fmt,transpose(A)
       print*
       print*,"L"
       print fmt,transpose(L)
       print*
       print *, "D"
       print fmt,D
       print*
       do i=1,n
          do j=1,n
             if(i==j)then
                D2(i,j)=D(i)
                D3(i,j)=dsqrt(D(i))
             else
                D2(i,j)=0._r8
                D3(i,j)=0._r8
             end if
          end do
       end do
       T=matmul(matmul(L, D2), transpose(L))-A
       error=maxval(T(1:n,1:n))
       print*, "T=LDL'-A"
       print fmt, transpose(T)
       print*, "Max error (max T):", error
       if(abs(error)>eps)then
          print*,"Error in choleskydecomp:"
          print*,"   Decomposition with too much error, error=", error
          stop
       end if
    end if

    !Put diagonal D in the diagonal of L
    do i=1, n
       L(i,i)=D(i)
    end do

    return
  end subroutine choleskydecomp

  function condnumest(L, n)
    !---------------------------------------------------------
    ! MATRI_CONDNUMEST
    !Calculate the matrix condition number lower estimate
    !
    !Uses maximum value diveded by minimun value of the diagonal
    !  matrix generated by a cholesky LDLt decomposition
    !
    !Receives a matrix decomposed by cholesky LDLt, that is,
    !  a lower triangular matrix L and a diagonal D stored
    !  in the diagonal of L, and tha size 'n' of the matrix
    !
    !Returns tha max/min value of the diagonal
    !---------------------------------------------------------
    !Condition number estimative
    real (r8):: condnumest

    !Size the matrix
    integer (i4):: n

    !Matrix
    real (r8):: L(1:n, 1:n)

    !Auxiliar variables
    integer (i4):: i
    real (r8):: dmin
    real (r8):: dmax

    condnumest=0

    !Get max and min of diagonal
    dmin=100000.
    dmax=0.
    do i=1, n
       dmin=min(abs(L(i,i)), dmin)
       dmax=max(abs(L(i,i)), dmax)
    end do

    if(dmin==0)then
       print*, "Warning condnumest: dmin=0"
       condnumest=9999999999999.
    else
       condnumest=dmax/dmin
    end if

    return
  end function condnumest

  subroutine choleskysolve(L, x, b, n)
    !--------------------------------------------------
    !  solve cholesky decomposed system (L D L')
    !
    !  L is lower diagonal and has D on its diagonal
    !  x is the solution of L D L' x = b
    !
    !--------------------------------------------------
    !Dimension of the matrix
    integer(i4), intent(in) :: n

    !Lower triangular matrix
    real(r8), intent(inout) :: L(1:n,1:n)

    !RHS vetor
    real(r8), intent(in) :: b(1:n)

    !Solution
    real(r8), intent(out) :: x(1:n)

    !auxiliar Matrix
    real(r8):: y(1:n)
    real(r8):: D(1:n,1:n)
    real(r8)::  T(1:n,1:n)

    !indexes
    integer (i4):: i
    integer (i4):: j

    !Check ?
    logical:: check
    real(r8):: res(1:n)
    real(r8):: normres
    character(len=16):: fmt

    !Solve Ly=b
    do i=1,n
       y(i)=b(i)
       do j=1,i-1
          y(i)=y(i)-y(j)*L(i,j)
       end do
    end do

    !Solve Dz=y, and put z-> y
    do i=1,n
       if(L(i,i)==0)then
          y(i)=0
       else
          y(i)=y(i)/L(i,i)
       end if
    end do

    !Solve L'x=y
    do i=n,1, -1
       x(i)=y(i)
       do j=n,i+1,-1
          x(i)=x(i)-x(j)*L(j,i)
       end do
    end do

    !Debugging checks
    check= .false. !.true.
    if(check)then
       print*, "CholeskySolve Debug"
       fmt="(xxf16.8)"
       write(fmt(2:3),'(i1)') n
       print*, "x:"
       print fmt, x
       D(1:n,1:n)=0._r8
       do i=1,n
          D(i,i)=L(i,i)
          L(i,i)=1._r8
       end do
       T=matmul(matmul(L,D), transpose(L))
       res=b-matmul(T,x)
       print*,"Residual:"
       print fmt, res
       normres=norm(res)
       print*, "Residual norm : ", normres
       if(normres>eps)then
          print*,"Warning in choleskysolve:"
          print*,"   Residual too large, norm(residual)=", normres
       end if
    end if

    return
  end subroutine choleskysolve


  !===============================================================================================
  !    ERROR NORMS
  !===============================================================================================

  function positive(x)
    !-----------------------------------------
    ! POSITIVE
    ! If x negative return 0, else return x
    !----------------------------------------------
    real(r8), intent(in) :: x
    real(r8):: positive

    if(x<0)then
       positive=0._r8
    else
       positive=x
    end if

    return
  end function positive

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

  subroutine normalize(p)
    !-----------------------------------------
    ! NORMALIZE
    ! Normalizes a vector
    !----------------------------------------------
    real(r8), intent(inout):: p(:)
    real(r8):: normp

    normp=norm(p) !sqrt(dot_product ( p(1:3), p(1:3) ))
    p=p/normp

    return
  end subroutine normalize

  function distance(x,y)
    !------------------------------------------------------------
    ! DISTANCE
    !	Calculates the distance between the points x e y in R3
    !
    !	Pedro Peixoto - Dec 2010
    !---------------------------------------------------------------

    real (r8), dimension(3)::x
    real (r8), dimension(3)::y
    real (r8), dimension(3)::z
    real (r8):: distance

    z=x-y
    distance=norm(z)

    return
  end function distance

  function error_norm_max(f, g, n)
    !-------------------------------------------
    !Calculates the maximum absolute value of
    !  f-g vector both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_max

    error_norm_max=maxval(abs(f(1:n)-g(1:n)))


    return
  end function error_norm_max

  function error_norm_max_rel(f, g, n)
    !-------------------------------------------
    !Calculates the maximum absolute value of
    !  f-g vector divided by the max of g
    !   both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_max_rel
    real (r8):: maxfg
    real (r8):: maxg

    maxfg=maxval(abs(f(1:n)-g(1:n)))
    maxg=maxval(abs(g(1:n)))
    error_norm_max_rel=maxfg/maxg

    return
  end function error_norm_max_rel

  function error_norm_2(f, g, n)
    !-------------------------------------------
    !Calculates the square root of the
    !  the sum of the squares of (f-g) vector
    !  both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_2

    real (r8):: sum_sq

    sum_sq=dot_product((f-g),(f-g))

    error_norm_2=dsqrt(sum_sq/real(n, r8))

    return
  end function error_norm_2

  function error_norm_2_rel(f, g)
    !-------------------------------------------
    !Calculates relative L2 error
    !  That is, the square root of the
    !  the sum of the squares of (f-g) vector
    !  divided by the norm 2 of g
    !  both having size 1:n
    !-------------------------------------------
    real (r8),  intent(in) :: f(:)
    real (r8),  intent(in) :: g(:)
    real (r8):: error_norm_2_rel

    real (r8):: sum_sq
    real (r8):: sum_sq_g

    sum_sq=dot_product((f-g),(f-g))
    sum_sq_g=dot_product(g,g)

    if(abs(sum_sq_g) == 0 )then
       print*, "error_norm_2_rel error: division by zero"
       stop
    end if
    error_norm_2_rel=dsqrt(sum_sq/sum_sq_g)

    return
  end function error_norm_2_rel

  function error_norm_1(f, g, n)
    !-------------------------------------------
    !Calculates the
    !  the sum of the absolute values of (f-g) vector
    !  both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_1

    error_norm_1=sum(abs(f-g))

    return
  end function error_norm_1

end module linear_algebra 
