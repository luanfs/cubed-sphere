module datastruct
!====================================================================
!
!
! Module for cubed-sphere mesh structuring.
!
! Reference: "C. Ronchi, R. Iacono, P.S. Paolucci, The “Cubed Sphere”:
! A New Method for the Solution of Partial Differential Equations
! in Spherical Geometry."
!
! Luan Santos 2022
! (luan.santos@usp.br)
!
! Based on module 'datastruct' from iModel
! https://github.com/pedrospeixoto/iModel
!=====================================================================

!Use global constants and kinds
use constants, only: i4, r8

!-------------------------------------------------
! Simple 3d cartesian vector
!------------------------------------------------
type vector
    real (r8), dimension(1:3) :: v
end type vector

!-------------------------------------------------
! Simple 2x2 matrix
!------------------------------------------------
type matrix
    real (r8), dimension(1:2,1:2) :: M
end type matrix

!-------------------------------------------------
!General sphere point structure
! This point is not necessarily a grid point 
!-------------------------------------------------
type point_structure
    !Vector form of cartesian coordinates
    real (r8), dimension(1:3) :: p

    !Spherical/geographic coordinates of node in radians
    !  lat in [-pi/2, pi/2] , lon in [-pi, pi[   
    real (r8) :: lat, lon
end type point_structure

!---------------------------------
! Cubed sphere structure
!---------------------------------
type cubedsphere
    !---------------------------------------------------------
    !  The quadrilateral points are labeled as below
    !
    !  po-------pv--------po
    !  |                  |
    !  |                  |
    !  |                  |
    !  pu       pc        pu
    !  |                  |
    !  |                  |
    !  |                  |
    !  po--------pv-------po
    !---------------------------------------------------------

    ! Center of quadrilaterals 
    type(point_structure), allocatable :: pc(:,:,:)

    ! Midpoint of quadrilateral edges (y direction)
    type(point_structure), allocatable :: pu(:,:,:)

    ! Midpoint of quadrilateral edges (x direction)
    type(point_structure), allocatable :: pv(:,:,:)

    ! Vertices of quadrilaterals
    type(point_structure), allocatable :: po(:,:,:)

    ! Lat/lon grid (for plotting)
    type(point_structure), allocatable :: ll(:,:)

    ! Tangent vector at pu in x direction
    type(vector), allocatable :: tgx_pu(:,:,:)

    ! Tangent vector at pu in y direction
    type(vector), allocatable :: tgy_pu(:,:,:)

    ! Tangent vector at pv in x direction
    type(vector), allocatable :: tgx_pv(:,:,:)

    ! Tangent vector at pv in y direction
    type(vector), allocatable :: tgy_pv(:,:,:)

    ! Tangent vector at pc in x direction
    type(vector), allocatable :: tgx_pc(:,:,:)

    ! Tangent vector at pc in y direction
    type(vector), allocatable :: tgy_pc(:,:,:)

    ! latlon to contravariant conversion at pu
    type(matrix), allocatable :: ll2contra_pu(:,:,:)

    ! latlon to contravariant conversion at pv
    type(matrix), allocatable :: ll2contra_pv(:,:,:)

    ! Contravariant to latlon conversion at pu
    type(matrix), allocatable :: contra2ll_pu(:,:,:)

    ! Contravariant to latlon conversion at pv
    type(matrix), allocatable :: contra2ll_pv(:,:,:)

    ! Metric tensor at pc
    real(r8), allocatable:: mt_pc(:,:,:)

    ! Metric tensor at pu
    real(r8), allocatable:: mt_pu(:,:,:)

    ! Metric tensor at pu
    real(r8), allocatable:: mt_pv(:,:,:)

    ! Minimum/Maximum geodesic distance between vertice points in radians
    real(r8):: mindist, maxdist, meandist

    ! Minimum/Maximum geodesic areas
    real(r8):: minarea, maxarea, meanarea

    ! Local coordinates uniform grid size
    real(r8):: dx, dy

    ! Latlon grid size (for plotting)
    real(r8):: dlon, dlat

    ! Latlon grid points nearest neighbours indexes
    integer(i4), allocatable:: ix_ll(:,:), jy_ll(:,:), panels_ll(:,:)

    ! Number of cells along a coordinate axis
    integer(i4):: n

    ! Total number of cells (6*n*n)
    integer(i4):: nbcells

    ! size of halo region 
    integer(i4):: halosize

    ! Number of cells along a coordinate axis (including ghost cells)
    integer(i4):: ntotal

    ! Lat/lon grid size
    integer(i4):: nlon, nlat

    ! Interior of each panel grid indexes
    integer(i4):: i0, iend, j0, jend

    ! Cell start/end indexes
    integer(i4) :: n0, nend

    !Flag for loadable grid (1- yes, 0- no)
    integer (i4) :: loadable

    !Mesh kind
    !equiangular
    character(len=16) :: kind

    !Grid name, used for file names as outputs
    character (len=128) :: name

    !File where grid may be loadable
    character (len=128) :: filename

end type cubedsphere 


!---------------------------------------------------------
! Variable for scalar values on cubed-sphere grid
!---------------------------------------------------------
type scalar_field
    ! Values array, ordered in the same sequence as the
    ! cubed sphere points
    real (r8), allocatable  :: f(:,:,:)

    ! Position of the values relative to a mesh
    !   0 - Centers (pc)
    !   1 - Vertices (po) 
    !   2 - Midpoint at u position (pu)
    !   3 - Midpoint at v position (pv)
    integer (i4) :: pos

    !Variable name - long name - detailed name
    ! This is used to for filenames of this variable
    character (len=256) :: name

end type scalar_field 

!---------------------------------------------------------
! Variable for vector values on cubed-sphere grid
!---------------------------------------------------------
type vector_field
    ! Geographical coordinates 
    type(scalar_field) :: u
    type(scalar_field) :: v

    ! Contravariant components
    type(scalar_field) :: ucontra
    type(scalar_field) :: vcontra

    ! Contravariant time-averaged winds
    type(scalar_field) :: ucontra_av
    type(scalar_field) :: vcontra_av

    ! Covariant components
    type(scalar_field) :: ucovari
    type(scalar_field) :: vcovari

    ! Position of the values relative to a mesh
    !   0 - Centers (pc)
    !   1 - Vertices (po) 
    !   2 - Midpoint at u position (pu)
    !   3 - Midpoint at v position (pv)
    integer (i4) :: pos

    !Variable name - long name - detailed name
    ! This is used to for filenames of this variable
    character (len=256) :: name
end type vector_field 

!---------------------------------------------------------
! Data structure for simulation class
!---------------------------------------------------------
type simulation
    ! variables for errors 
    real(r8):: linf_error, l1_error, l2_error

    ! variables for mass
    real(r8):: mass, mass0, mass_variation

    ! Time step
    real(r8):: dt

    ! Initial condition
    integer(i4) :: ic

    ! Velocity field
    integer(i4) :: vf

    ! Logical for exact solution
    logical :: exactsolution

    ! Initial condition
    character(len=16) :: ic_name

    ! Velocity field
    character(len=16) :: vf_name

    ! One dimensional flux scheme
    character(len=16) :: recon1d

    ! Two dimensional splitting scheme
    character(len=16) :: opsplit

    ! Metric tensor scheme
    character(len=16) :: mt

    !Variable name - long name - detailed name
    ! This is used to for filenames of this variable
    character (len=256) :: name

end type simulation

!---------------------------------------------------------
! Variable for ppm parabola
!---------------------------------------------------------
type ppm_parabola
    ! parabola coefficients ! Notation from Colella and  Woodward 1984
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    real (r8), allocatable  :: q_L(:,:,:)
    real (r8), allocatable  :: q_R(:,:,:)
    real (r8), allocatable  :: dq(:,:,:)
    real (r8), allocatable  :: q6(:,:,:)

    ! parabola fluxes
    real (r8), allocatable  :: f_L(:,:,:) ! flux from left
    real (r8), allocatable  :: f_R(:,:,:) ! flux from right
    real (r8), allocatable  :: f_upw(:,:,:)  ! upwind flux

    ! Divergence of flux
    real (r8), allocatable  :: df(:,:,:)  ! flux divergence

    ! field to be reconstructed
    type(scalar_field) :: Q

    ! Direction of reconstruction
    !   1 - x direction
    !   2 - y direction 
    integer (i4) :: dir

    ! N - number of cells along a coordinate axis
    integer (i4) :: N

    ! One dimensional reconstruction scheme
    character(len=16) :: recon

    ! Metric tensor scheme
    character(len=16) :: mt


end type ppm_parabola
 
end module datastruct
