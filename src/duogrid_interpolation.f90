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
    pi, deg2rad,&
    i0, iend, &
    j0, jend, &
    n0, nend, &
    nghost

!Data structures
use datastruct, only: &
    cubedsphere, &
    scalar_field, &
    vector_field, &
    simulation

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




end module duogrid_interpolation
