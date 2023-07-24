module duogrid_interpolation
!========================================================================
!
! Module for Lagrange interpolation at duogrid
!
! Luan Santos 2022
! Routines based on iModel (https://github.com/pedrospeixoto/iModel)
!========================================================================

!Global constants
use constants, only: &
    griddir, &
    i4, &
    pardir, &
    r8, &
    showonscreen, &
    simulcase, &
    nbfaces, &
    erad, &
    pi, pio4, deg2rad,&
    i0, iend, &
    j0, jend, &
    n0, nend, &
    nghost

!Data structures
use datastruct, only: &
    cubedsphere, &
    scalar_field, &
    vector_field, &
    simulation, &
    lagrange_poly_cs

! CS mapping
use sphgeo, only: &
    inverse_equiangular_gnomonic_map2

implicit none

contains 

subroutine gethalodata(Q)
    !---------------------------------------------------
    ! Fill the center ghost cell values of Q that are needed
    ! for duogrid interpolation
    !--------------------------------------------------
    type(scalar_field), intent(inout) :: Q
    integer(i4) :: p, east, north, south, west


    ! --------------------- Panel 1 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 1
    north = 5
    south = 6
    east  = 2
    west  = 4
 
    ! Data of panel 1 from east
    Q%f(iend+1:,j0:jend,p) = Q%f(i0:i0+nghost, j0:jend, east) ! Panel 2

    ! Data of panel 1 from  west
    Q%f(:i0-1,j0:jend,p) = Q%f(iend-nghost:iend, j0:jend, west) ! Panel 4

    ! Data of panel 1 from north
    Q%f(i0:iend,jend+1:,p) = Q%f(i0:iend, j0:j0+nghost, north) ! Panel 5

    ! Data of panel 1 from south
    Q%f(i0:iend,:j0-1,p)  = Q%f(i0:iend, jend-nghost:jend, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 2 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 2
    north = 5
    south = 6
    east  = 3
    west  = 1

    ! Data of panel 2 from east
    Q%f(iend+1:,j0:jend,p) = Q%f(i0:i0+nghost, j0:jend, east) ! Panel 3

    ! Data of panel 2 from west
    Q%f(:i0-1,j0:jend,p) = Q%f(iend-nghost:iend, j0:jend, west) ! Panel 1

    ! Data of panel 2 from north
    Q%f(i0:iend,jend+1:,p) = transpose(Q%f(jend:jend-nghost:-1, i0:iend, north)) ! Panel 5

    ! Data of panel 2 from south
    Q%f(i0:iend,:j0-1,p)  = transpose(Q%f(jend:jend-nghost:-1, iend:i0:-1,south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 3 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 3
    north = 5
    south = 6
    east  = 4
    west  = 2

    ! Data of panel 3 from east
    Q%f(iend+1:,j0:jend,p) = Q%f(i0:i0+nghost, j0:jend,east) ! Panel 4

    ! Data of panel 3 from west
    Q%f(:i0-1,j0:jend,p) = Q%f(iend-nghost:iend, j0:jend, west) ! Panel 2

    ! Data of panel 3 from north
    Q%f(i0:iend,jend+1:,p) = Q%f(iend:i0:-1, jend:jend-nghost:-1, north) ! Panel 5

    ! Data of panel 3 from south
    Q%f(i0:iend,:j0-1,p)  = Q%f(iend:i0:-1, j0+nghost:j0:-1, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 4 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 4
    north = 5
    south = 6
    east  = 1
    west  = 3

    ! Data of panel 4 from east
    Q%f(iend+1:,j0:jend,p) = Q%f(i0:i0+nghost, j0:jend, east) ! Panel 1

    ! Data of panel 4 from west
    Q%f(:i0-1,j0:jend,p) = Q%f(iend-nghost:iend, j0:jend, west) ! Panel 3

    ! Data of panel 4 from north
    Q%f(i0:iend,jend+1:,p) = transpose(Q%f(i0:i0+nghost, jend:jend-nghost:-1, north)) ! Panel 5

    ! Data of panel 4 from south
    Q%f(i0:iend,:j0-1,p)   = transpose(Q%f(i0+nghost:i0:-1, j0:jend, south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 5 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 5
    north = 3
    south = 1
    east  = 2
    west  = 4

    ! Data of panel 5 from east
    Q%f(iend+1:,j0:jend,p) = transpose(Q%f(i0:iend, jend:jend-nghost:-1, east)) ! Panel 2

    ! Data of panel 5 from west
    Q%f(:i0-1,j0:jend,p) = transpose(Q%f(iend:i0:-1,jend-nghost:jend, west)) ! Panel 4

    ! Data of panel 5 from north
    Q%f(i0:iend,jend+1:,p) = Q%f(iend:i0:-1, jend:jend-nghost:-1, north) ! Panel 3

    ! Data of panel 5 from south
    Q%f(i0:iend,:j0-1,p)   = Q%f(i0:iend, jend-nghost:jend, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE

    ! --------------------- Panel 6 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, i0, iend, j0, jend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 6
    north = 1
    south = 3
    east  = 2
    west  = 4

    ! Data of panel 6 from east
    Q%f(iend+1:,j0:jend,p) = transpose(Q%f(iend:i0:-1, j0:j0+nghost, east)) ! Panel 2

    ! Data of panel 6 from west
    Q%f(:i0-1,j0:jend,p) = transpose(Q%f(i0:iend,j0+nghost:j0:-1, west)) ! Panel 4

    ! Data of panel 6 from north
    Q%f(i0:iend,jend+1:,p) = Q%f(i0:iend, j0:j0+nghost, north) ! Panel 3

    ! Data of panel 6 from south
    Q%f(i0:iend,:j0-1,p)   = Q%f(iend:i0:-1, j0+nghost:j0:-1, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE


end subroutine gethalodata



subroutine dg_interp(Q, mesh, simul)
    !---------------------------------------------------
    !   duogrid interpolation of scalar field Q
    ! (ghost cells are defined at cell centers)
    !--------------------------------------------------
    type(cubedsphere), intent(in):: mesh
    type(scalar_field), intent(inout) :: Q
    type(simulation), intent(in) :: simul

    ! Fill ghost cell centers
    call gethalodata(Q)
end subroutine dg_interp



subroutine compute_lagrange_cs(L, mesh)
    !---------------------------------------------------
    ! compute the lagrange polynomials at ghost cell
    ! points of the cubed-sphere
    !--------------------------------------------------
    type(cubedsphere), intent(inout):: mesh
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4) ::  j, p, g, d

    if(mesh%kind .ne. 'equiangular') then
        print*, 'ERROR in compute_lagrange_cs: invalid mesh kind, ', mesh%kind
        stop
    end if

    if(L%pos==1) then
        ! Compute the nodes
        ! Init local coordinates grid
        L%y_support(n0) = -pio4 - mesh%halosize*mesh%dx*0.5
        do j = n0+1, nend
            L%y_support(j) = L%y_support(j-1) + mesh%dx
        end do

        ! Invert node points using the cs mapping from panel 2
        p = 2
        do g = 1, nghost
            do j = n0, nend
                ! invert it
                call inverse_equiangular_gnomonic_map2(L%x_nodes(j,g), L%y_nodes(j,g), p, &
                mesh%pc(iend+g,j,1)%p, mesh)
            end do
        end do

        ! Compute stencils
        do g = 1, nghost
            do j = n0, nend
                L%kend(j,g) = (L%y_nodes(j,g) - L%y_support(n0))/mesh%dx
                L%kend(j,g) = L%kend(j,g) + ceiling(L%order*0.5_r8)
                L%k0(j,g)   = L%kend(j,g) - L%order + 1

                if (j>=j0 .and. j<=jend)then
                    if(L%kend(j,g)>jend)then
                        L%kend(j,g) = jend-1
                        L%k0(j,g)   = L%kend(j,g) - L%order + 1 
                    else if (L%k0(j,g)<j0)then
                        L%k0(j,g)   = j0
                        L%kend(j,g) = L%k0(j,g) + L%order - 1 
                    end if

                else if (j>=jend+1) then
                    if (L%kend(j,g)>=nend+1)then
                        L%kend(j,g) = nend
                        L%k0(j,g)   = L%kend(j,g) - L%order + 1
                    end if
                else !i<i0
                    if (L%k0(j,g)<n0) then
                        L%k0(j,g)   = n0
                        L%kend(j,g) = L%k0(j,g) + L%order - 1
                    end if
                end if
                !print*,j, L%k0(j,g), L%kend(j,g), j0, jend
            end do
            !read(*,*)
        end do



        ! Store in f the nearest support points used in Lagrange interpolation
        do g = 1, nghost
            do j = n0, nend
                L%f_nodes(j,g,:) = L%y_support(L%k0(j,g):L%kend(j,g))
            end do
        end do

        ! Compute the Lagrange nodes at halo region
        do g = 1, nghost
            do j = n0, nend
                do d = 1, L%order
                    call lagrange_basis(L%y_nodes(j,g), L%f_nodes(j,g,:), L%degree, d, L%p_nodes(j,g,d))
                end do
            end do
        end do

        do g = 1, nghost
            do j = n0, nend
                !print*, abs(sum(L%p_nodes(j,g,:))-1._r8)
            end do
            !read(*,*)
        end do
    end if
end subroutine compute_lagrange_cs



subroutine lagrange_basis(x, x_support, N, j, Lj)
    !---------------------------------------------------
    ! Compute the jth Lagrange polynomial of degree N
    !--------------------------------------------------
    integer(i4), intent(in) :: N, j
    real(r8), intent(in) :: x, x_support(1:N+1)
    real(r8), intent(inout) :: Lj
    integer(i4) :: i

    Lj = 1._r8
    do i = 1, N+1
        if (i .ne. j) then
            Lj = Lj*(x-x_support(i))/(x_support(j)-x_support(i))
        end if
    end do
end subroutine lagrange_basis

end module duogrid_interpolation
