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
use constants, only: i4

!-------------------------------------------------
! Simple 3d cartesian vector
!------------------------------------------------
type vector
    real (kind=8), dimension(1:3) :: v
end type vector

!-------------------------------------------------
! Simple 2x2 matrix
!------------------------------------------------
type matrix
    real (kind=8), dimension(1:2,1:2) :: M
end type matrix

!-------------------------------------------------
!General sphere point structure
! This point is not necessarily a grid point 
!-------------------------------------------------
type point_structure
    !Vector form of cartesian coordinates
    real (kind=8), dimension(1:3) :: p

    !Spherical/geographic coordinates of node in radians
    !  lat in [-pi/2, pi/2] , lon in [-pi, pi[   
    real (kind=8) :: lat, lon
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

    !--------------------------------------------------------------------
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

    ! Tangent vector at po in x direction
    type(vector), allocatable :: tgx_po(:,:,:)

    ! Tangent vector at po in y direction
    type(vector), allocatable :: tgy_po(:,:,:)

    !--------------------------------------------------------------------
    ! latlon to contravariant conversion at pu
    type(matrix), allocatable :: ll2contra_pu(:,:,:)

    ! latlon to contravariant conversion at pv
    type(matrix), allocatable :: ll2contra_pv(:,:,:)

    ! latlon to contravariant conversion at pc
    type(matrix), allocatable :: ll2contra_pc(:,:,:)

    ! latlon to contravariant conversion at po
    type(matrix), allocatable :: ll2contra_po(:,:,:)

    !--------------------------------------------------------------------
    ! Contravariant to latlon conversion at pu
    type(matrix), allocatable :: contra2ll_pu(:,:,:)

    ! Contravariant to latlon conversion at pv
    type(matrix), allocatable :: contra2ll_pv(:,:,:)

    ! Contravariant to latlon conversion at pc
    type(matrix), allocatable :: contra2ll_pc(:,:,:)

    ! Contravariant to latlon conversion at po
    type(matrix), allocatable :: contra2ll_po(:,:,:)

    !--------------------------------------------------------------------
    ! latlon to covari conversion at pu
    type(matrix), allocatable :: ll2covari_pu(:,:,:)

    ! latlon to covariant conversion at pv
    type(matrix), allocatable :: ll2covari_pv(:,:,:)

    ! latlon to covariant conversion at pc
    type(matrix), allocatable :: ll2covari_pc(:,:,:)

    ! latlon to covariant conversion at po
    type(matrix), allocatable :: ll2covari_po(:,:,:)

    !--------------------------------------------------------------------
    ! covariant to latlon conversion at pu
    type(matrix), allocatable :: covari2ll_pu(:,:,:)

    ! covariant to latlon conversion at pv
    type(matrix), allocatable :: covari2ll_pv(:,:,:)

    ! covariant to latlon conversion at pc
    type(matrix), allocatable :: covari2ll_pc(:,:,:)

    ! covariant to latlon conversion at po
    type(matrix), allocatable :: covari2ll_po(:,:,:)

    !--------------------------------------------------------------------
    ! Contravariant to covariant conversion at pu
    type(matrix), allocatable :: contra2covari_pu(:,:,:)

    ! Contravariant to covariant conversion at pv
    type(matrix), allocatable :: contra2covari_pv(:,:,:)

    ! Contravariant to covariant conversion at pc
    type(matrix), allocatable :: contra2covari_pc(:,:,:)

    ! Contravariant to covariant conversion at po
    type(matrix), allocatable :: contra2covari_po(:,:,:)

    !--------------------------------------------------------------------
    ! Covariant to contravariant conversion at pu
    type(matrix), allocatable :: covari2contra_pu(:,:,:)

    ! Covariant to contravariant conversion at pv
    type(matrix), allocatable :: covari2contra_pv(:,:,:)

    ! Covariant to contravariant conversion at pc
    type(matrix), allocatable :: covari2contra_pc(:,:,:)

    ! Covariant to contravariant conversion at po
    type(matrix), allocatable :: covari2contra_po(:,:,:)


    !--------------------------------------------------------------------

    ! Metric tensor at pc
    real(kind=8), allocatable:: mt_pc(:,:,:)

    ! Metric tensor at pu
    real(kind=8), allocatable:: mt_pu(:,:,:)

    ! Metric tensor at pu
    real(kind=8), allocatable:: mt_pv(:,:,:)

    ! Metric tensor at po
    real(kind=8), allocatable:: mt_po(:,:,:)
    !--------------------------------------------------------------------

    ! Minimum/Maximum geodesic distance between vertice points in radians
    real(kind=8):: mindist, maxdist, meandist

    ! Minimum/Maximum geodesic areas
    real(kind=8):: minarea, maxarea, meanarea

    ! Local coordinates uniform grid size
    real(kind=8):: dx, dy

    ! Latlon grid size (for plotting)
    real(kind=8):: dlon, dlat

    ! Sphere radius
    real(kind=8):: radius
 
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
    real (kind=8), allocatable  :: f(:,:,:)

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
type velocity_field
    ! Geographical coordinates 
    type(scalar_field) :: u
    type(scalar_field) :: v

    ! Contravariant components
    type(scalar_field) :: ucontra
    type(scalar_field) :: vcontra

    ! Contravariant components from previous time step
    type(scalar_field) :: ucontra_old
    type(scalar_field) :: vcontra_old

    ! Contravariant time-averaged winds
    type(scalar_field) :: ucontra_time_av
    type(scalar_field) :: vcontra_time_av

    ! Contravariant time-averaged centered at time
    type(scalar_field) :: ucontra_time_centered
    type(scalar_field) :: vcontra_time_centered

    ! Covariant components
    type(scalar_field) :: ucovari
    type(scalar_field) :: vcovari

    ! Covariant components from previous time step
    type(scalar_field) :: ucovari_old
    type(scalar_field) :: vcovari_old

    ! Position of the values relative to a mesh
    !   0 - Centers (pc)
    !   1 - Vertices (po) 
    !   2 - Midpoint at u position (pu)
    !   3 - Midpoint at v position (pv)
    integer (i4) :: pos

    !Variable name - long name - detailed name
    ! This is used to for filenames of this variable
    character (len=256) :: name
end type velocity_field 


!---------------------------------------------------------
! Data structure for simulation class
!---------------------------------------------------------
type simulation
    ! variables for errors 
    real(kind=8):: linf_error_h    , l1_error_h    , l2_error_h
    real(kind=8):: linf_error_div  , l1_error_div  , l2_error_div
    real(kind=8):: linf_error_rv   , l1_error_rv   , l2_error_rv
    real(kind=8):: linf_error_av   , l1_error_av   , l2_error_av
    real(kind=8):: linf_error_av_pu, l1_error_av_pu, l2_error_av_pu
    real(kind=8):: linf_error_av_pv, l1_error_av_pv, l2_error_av_pv
    real(kind=8):: linf_error_h_po
    real(kind=8):: linf_error_gradh_pu
    real(kind=8):: linf_error_gradh_pv
    real(kind=8):: linf_error_Ku_po
    real(kind=8):: linf_error_Kv_po
    real(kind=8):: linf_error_K_po
    real(kind=8):: linf_error_ucovari_po
    real(kind=8):: linf_error_vcovari_po


    ! variables for mass
    real(kind=8):: mass, mass0, mass_variation

    ! Time step
    real(kind=8):: dt
    real(kind=8):: dto2
    real(kind=8):: t, tf
    integer(i4) :: nsteps, nplot, n, plotcounter

    ! CFL
    real(kind=8):: cfl

    ! var used in mass fixer
    real(kind=8):: a2

    ! Initial condition
    integer(i4) :: ic

    ! Velocity field
    integer(i4) :: vf

    ! Degree of interpolation at ghost cells
    integer(i4) :: id

    ! Degree of interpolation from D grid to A grid (1 or 3)
    integer(i4) :: avd

    ! Logical for exact solution
    logical :: exactsolution

    ! Initial condition
    character(len=16) :: ic_name

    ! Velocity field
    character(len=16) :: vf_name

    ! Degree of interpolation at ghost cells
    character(len=16) :: id_name

    ! Degree of grid averaging
    character(len=16) :: avd_name

    ! One dimensional flux scheme
    character(len=16) :: recon1d

    ! Two dimensional splitting scheme
    character(len=16) :: opsplit

    ! Metric tensor scheme
    character(len=16) :: mt

    ! Mass fixer scheme
    character(len=16) :: mf

    ! Departure point scheme
    character(len=16) :: dp

    ! Edge treatment
    character(len=16) :: et

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
    real (kind=8), allocatable  :: q_L(:,:,:)
    real (kind=8), allocatable  :: q_R(:,:,:)
    real (kind=8), allocatable  :: dq(:,:,:)
    real (kind=8), allocatable  :: q6(:,:,:)

    ! parabola fluxes
    real (kind=8), allocatable  :: f_upw(:,:,:)  ! upwind flux

    ! Divergence of flux
    real (kind=8), allocatable  :: df(:,:,:)  ! flux divergence

    ! field to be reconstructed
    type(scalar_field) :: Q

    ! Direction of reconstruction
    !   1 - x direction
    !   2 - y direction 
    integer (i4) :: dir

    ! reference point (pc or po)
    !   1 - pc
    !   2 - po
    integer (i4) :: point


    ! N - number of cells along a coordinate axis
    integer (i4) :: N

    ! One dimensional reconstruction scheme
    character(len=16) :: recon

    ! Metric tensor scheme
    character(len=16) :: mt

    ! Edge treatment
    character(len=16) :: et


end type ppm_parabola

!---------------------------------------------------------
! Variable Lagrange polynomials at ghost cells on the
! cubed-sphere
!---------------------------------------------------------
type lagrange_poly_cs
    real (kind=8), allocatable  :: y_support(:)  ! Support points
    real (kind=8), allocatable  :: f_support(:,:)  ! Value of the function at the support points
    real (kind=8), allocatable  :: x_nodes(:,:) ! Nodes where we perform interpolation
    real (kind=8), allocatable  :: y_nodes(:,:) ! Nodes where we perform interpolation
    real (kind=8), allocatable  :: p_nodes(:,:,:) ! Lagrange polynomials at nodes
    real (kind=8), allocatable  :: f_nodes(:,:,:) ! Support values used by node stenil
    real (kind=8), allocatable  :: halodata_east(:,:,:) ! var to store the needed halo data
    real (kind=8), allocatable  :: halodata_west(:,:,:) ! var to store the needed halo data
    real (kind=8), allocatable  :: halodata_north(:,:,:) ! var to store the needed halo data
    real (kind=8), allocatable  :: halodata_south(:,:,:) ! var to store the needed halo data

    ! stencil
    integer(i4), allocatable :: k0(:,:), kend(:,:)

    ! degree
    integer (i4) :: degree

    ! order
    integer (i4) :: order

    ! position 1-center, 2-edges
    integer (i4) :: pos

end type lagrange_poly_cs
 
end module datastruct
