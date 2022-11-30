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
  !General sphere point structure
  ! This point is not necessarily a grid point 
  !-------------------------------------------------
  type point_structure

     !Vector form of cartesian coordinates
     real (r8), dimension(1:3) :: p

     !Spherical/geographic coordinates of node in radians
     !  lat in [-pi/2, pi/2] , lon in [-pi, pi[   
     real (r8) :: lat, lon

     !Cubed - sphere local coordinates of node
     real (r8) :: x, y

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

   ! Geodesical areas
   real(r8), allocatable:: area(:,:,:)

   ! Metric tensor at pc
   real(r8), allocatable:: gc(:,:,:)

   ! Minimum/Maximum geodesic distance between vertice points in radians
   real(r8):: mindist, maxdist, meandist

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

   ! Number of ghost cells on the left
   integer(i4):: nbgl

   ! Number of ghost cells on the right
   integer(i4):: nbgr

   ! Number of ghost cells on the right and on the left
   integer(i4):: nbg

   ! Number of cells along a coordinate axis (including ghost cells)
   integer(i4):: ntotal

   ! Lat/lon grid size
   integer(i4):: nlon, nlat

   ! Interior of each panel grid indexes
   integer(i4):: i0, iend, j0, jend

   !Flag for loadable grid (1- yes, 0- no)
   integer (i4) :: loadable

   !Mesh kind
   !equiangular
   character(len=16) :: kind

   !Mesh midpoints
   !geo (midpoints on geodesics)
   !local (midpoints on local coordinate system)
   character(len=16) :: midpos

   !Mesh resolution
   !unif (uniform)
   !var  (variable)
   character(len=16) :: resolution


   !Grid name, used for file names as outputs
   character (len=128) :: name

  end type cubedsphere 

 
end module datastruct
