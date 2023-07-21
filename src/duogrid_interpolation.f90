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

    ! --------------------- Panel 0 ----------------------------
    p = 1
    north = 5
    south = 6
    east  = 2
    west  = 4
 
    ! Data of panel 0 from east
    Q%f(iend+1:,j0:jend,p) = Q%f(i0:i0+nghost,j0:jend,east) ! Panel 1

    ! Data of panel 0 from  west
    Q%f(:i0-1,j0:jend,p) = Q%f(iend-nghost:iend,j0:jend, west) ! Panel 3

    ! Data of panel 0 from north
    Q%f(i0:iend,jend+1:,p) = Q%f(i0:iend,j0:j0+nghost, north) ! Panel 4

    ! Data of panel 0 from south
    Q%f(i0:iend,:j0-1,p)  = Q%f(i0:iend,jend-nghost:jend, south) ! Panel 5
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
