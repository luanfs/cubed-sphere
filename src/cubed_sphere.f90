module cubed_sphere  
  !=================================================================================
  !    MAIN ROUTINES FOR THE  CUBED-SPHERE MESH CONSTRUCTION
  !
  ! Routines based on iModel https://github.com/pedrospeixoto/iModel
  !
  !=================================================================================

  !Global constants
  use constants, only: &
       datadir, &
       griddir, &
       i4, &
       nbfaces, &
       pardir, &
       pi, &
       pi2, &
       pio2, &
       pio4, &
       unitspharea, &
       acube, &
       rad2deg, &
       r8, &
       showonscreen

  !Data structures
  use datastruct, only: &
      cubedsphere

  !Input routines
  use input, only: &
      meshload, &
      meshread

  !Output routines
  use output, only: &
      meshstore
 
  !Data allocation
  use allocation, only: &
      meshallocation, &
      r8_1darray_allocation

  !Spherical geometry
  use sphgeo, only: &
      equidistant_gnomonic_map, &
      norm, &
      midpoint, &
      sphquadarea, &
      inverse_equidistant_gnomonic_map, &
      sph2cart, &
      cart2sph, &
      binary_search

  implicit none

  contains

  subroutine meshbuild(mesh)
    !---------------------------------------------------------------------
    !
    !	MESHBUILD
    !
    !	Creates or loads a mesh structure with parameters given implicitly by
    !
    !	mesh%kind="equiangular", "equiedge"
    !
    !	mesh%n  -> number of cells along a coordinate axis
    !
    !	mesh%nbcells  -> total number of cells
    !
    !	mesh%loadable -> 1 or 0, where 1 indicates to load and 0 not to load
    !          the grid structure from files
    !
    !	mesh%mid -> geo or local (midpoint based on geodesic or local coordinate system) 
    !
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    !Auxiliar variables
    character (len=256):: header
    logical:: ifile

    !Give a name to the grid
    call namegrid(mesh)

    !Check if file to be read exists
    if(trim(mesh%kind)=="read")then
      header=trim(griddir)//trim(mesh%name)
      inquire(file=header, exist=ifile)
      if(.not.ifile)then
        print*, "meshbuild error: cannot find mesh data file to read"
        print*, trim(mesh%name)
        stop
      end if
    end if

    !Set header filename and check existence
    header=trim(griddir)//trim(mesh%name)//".dat"
    inquire(file=header, exist=ifile)

    ! Define number of ghost cells
    call initghostcells(mesh)

    ! Allocate the mesh
    call meshallocation(mesh)

    if(mesh%loadable==1 .and. ifile) then
      !------------------------------------------------
      !Load the whole mesh
      !------------------------------------------------
      call meshload(mesh, header)

    else
      !------------------------------------------------
      !Generate mesh
      !------------------------------------------------
      select case(trim(mesh%kind))
      case("equiangular") !Equiangular grid
        call equiangular_cubedsphere_generation(mesh)

      case("read") !Read nodes from file
         call meshread(mesh)

      case default
         print*, "MESH BUILD ERROR: Invalid mesh kind : ", mesh%kind
         stop
         return
      end select

      ! Generate the grid properties
      call cubedsphere_properties(mesh)

      ! Create lat/lon grid
      call latlon_grid(mesh)

      ! Store the grid
      call meshstore(mesh, header)
  end if

    print*
    print*, "Mesh created or loaded: ", mesh%name
    print*,"-------------------------------------------------"
    print*

    return
  end subroutine meshbuild

  subroutine namegrid(mesh)
    !---------------------------------------------------------------------
    !
    ! NAMEGRID
    !
    ! Sets a name for the grid bases on type, number of cells, ...
    ! Name is returned in mesh%name
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    character (len=60):: an

    !Set named based on number of cells along a coordinate axis
    write(an,'(i8)') mesh%n

    if(trim(mesh%kind)=="read")then
       mesh%name = trim(mesh%name)
    else
       mesh%name = trim(mesh%kind)//"_"//trim(adjustl(an))
    end if

    mesh%name = trim(mesh%name)//"_"//trim(mesh%midpos)//"_"//trim(mesh%resolution) 

  end subroutine namegrid

  subroutine initghostcells(mesh)
    !---------------------------------------------------------------------
    !
    ! INITGHOSTCELLS
    !
    ! Defines the number of ghost cells, as well the panel grid interior
    ! indexes
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Number of ghost cells
    mesh%nbgl = 3
    mesh%nbgr = 3
    mesh%nbg  = mesh%nbgl + mesh%nbgr
 
    ! Interior indexes
    mesh%i0   = mesh%nbgl + 1 
    mesh%iend = mesh%i0 + mesh%n - 1 
    mesh%j0   = mesh%i0 
    mesh%jend = mesh%iend

    ! Total number of cells along a coordinate axis including ghost cells
    mesh%ntotal = mesh%n + mesh%nbg

  end subroutine


  subroutine equiangular_cubedsphere_generation(mesh)
    !---------------------------------------------------------------------
    !
    ! EQUIANGULAR_CUBEDSPHERE_GENERATION
    !
    ! This routine generates the equiangular cubed-sphere vertices
    !
    ! Based on "C. Ronchi, R. Iacono, P.S. Paolucci, The “Cubed Sphere”:
    ! A New Method for the Solution of Partial Differential Equations
    ! in Spherical Geometry."
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Real aux vars
    real(r8), allocatable :: x(:) ! Local coordinates (angular)
    real(r8), allocatable :: y(:) ! Local coordinates (angular)
    real(r8), allocatable :: tanx(:)
    real(r8), allocatable :: tany(:)

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter
    integer(i4) :: i0, j0, iend, jend ! 2D grid interior indexes

    ! Number of cells along a coordinate axis
    n = mesh%n
 
    ! Number of cells along a coordinate axis + number of ghost cells
    ntotal = mesh%ntotal
   
    print*, 'Computing cubed-sphere vertices...'

    ! Allocation
    call r8_1darray_allocation(x, 0, ntotal)
    call r8_1darray_allocation(y, 0, ntotal)
    call r8_1darray_allocation(tanx, 0, ntotal)
    call r8_1darray_allocation(tany, 0, ntotal)

    ! Local coordinates grid size
    mesh%dx = pio2/mesh%n
    mesh%dy = pio2/mesh%n

    ! Init local coordinates grid
    x(0) = -pio4 - mesh%nbgl*mesh%dx
    do i = 1, ntotal
      x(i) = x(i-1) + mesh%dx
    end do
    y = x

    ! Angular coordinates are mapped to the cube
    tanx = dtan(x)
    tany = tanx

    ! Compute the gnomonic mapping at the vertices
    do panel = 1, nbfaces
      do i = 0, ntotal
        do j = 0, ntotal
          call equidistant_gnomonic_map(acube*tanx(i), acube*tany(j), mesh%po(i,j,panel)%p, panel)
        end do
      end do
    end do

    ! Deallocation
    deallocate(x, y, tanx, tany)
  end subroutine equiangular_cubedsphere_generation


  subroutine cubedsphere_properties(mesh)
    !---------------------------------------------------------------------
    !
    ! CUBEDSPHERE_PROPERTIES
    !
    ! This routine generates cubed-sphere properties, such as areas, lenghts,
    ! metric tensor, tangent vectors...
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p, panel    ! Panel counter
    integer(i4) :: i0, j0, iend, jend ! 2D grid interior indexes

    ! Number of cells along a coordinate axis
    n = mesh%n
 
    ! Number of cells along a coordinate axis + number of ghost cells
    ntotal = mesh%ntotal

    ! Compute the midpoints
    call compute_midpoints(mesh)

    ! Compute the local coordinates
    call compute_localcoords(mesh)

    if(mesh%resolution=="unif")then
      panel = 1 ! In this case, we only need to compute areas, metric tensor and etc in a single panel
    else
      panel = nbfaces
      print*, 'ERROR on meshallocation: invalid mesh resolution', mesh%resolution
      stop
    end if

    print*, 'Computing cubed-sphere areas...'

    ! Compute the quadrilateral areas
    do p = 1, panel
      do i = 1, ntotal
        do j = 1, ntotal
          mesh%area(i,j,p) = sphquadarea(mesh%po(i-1,j-1,1)%p, mesh%po(i-1,j,1)%p, mesh%po(i,j-1,1)%p, mesh%po(i,j,1)%p)
        end do
      end do
    end do

    !i0 = mesh%i0
    !j0 = mesh%j0
    !iend = mesh%iend
    !jend = mesh%jend
    !print*, n ,(6._r8*sum(mesh%area(i0:iend,j0:jend,1))-unitspharea)/unitspharea

  end subroutine cubedsphere_properties


  subroutine compute_midpoints(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_MIDPOINTS
    !
    ! Given the cell vertices already computed in mesh, this routines
    ! compute the cell centers and edges midpoints.
    !
    !
    !  Cell points are labeled as below
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
    !
    ! In other words, given po, this routine compute pu, pv and pc
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter
    ! Real aux vars
    real(r8), allocatable :: x(:) ! Local coordinates (angular)
    real(r8), allocatable :: y(:) ! Local coordinates (angular)
    real(r8), allocatable :: tanx(:)
    real(r8), allocatable :: tany(:)
    real(r8) :: pc(1:3), pc1(1:3), pc2(1:3)

    ! Number of cells along a coordinate axis
    n = mesh%n
 
    ! Number of cells along a coordinate axis + number of ghost cells
    ntotal = mesh%ntotal
 
    print*, 'Computing cubed-sphere midpoints...'

    select case(trim(mesh%midpos))
      case("local")
        ! Allocation
        call r8_1darray_allocation(x, 0, ntotal)
        call r8_1darray_allocation(y, 1, ntotal)
        call r8_1darray_allocation(tanx, 0, ntotal)
        call r8_1darray_allocation(tany, 1, ntotal)

        ! Local coordinates grid size
        mesh%dx = pio2/mesh%n
        mesh%dy = pio2/mesh%n

        ! Init local coordinates grid
        x(0) = -pio4 - mesh%nbgl*mesh%dx
        do i = 1, ntotal
          x(i) = x(i-1) + mesh%dx
        end do
        y = (x(0:ntotal-1)+x(1:ntotal))*0.5d0

        ! Angular coordinates are mapped to the cube
        tanx = dtan(x)
        tany = dtan(y)

        ! Compute pu
        do panel = 1, nbfaces
          do i = 0, ntotal
            do j = 1, ntotal
              call equidistant_gnomonic_map(acube*tanx(i), acube*tany(j), mesh%pu(i,j,panel)%p, panel)
            end do
          end do
        end do

        ! Compute pv
        do panel = 1, nbfaces
          do i = 1, ntotal
            do j = 0, ntotal
              call equidistant_gnomonic_map(acube*tany(i), acube*tanx(j), mesh%pv(i,j,panel)%p, panel)
            end do
          end do
        end do

        ! Compute pc
        do panel = 1, nbfaces
          do i = 1, ntotal
            do j = 1, ntotal
              call equidistant_gnomonic_map(acube*tany(i), acube*tany(j), mesh%pc(i,j,panel)%p , panel)
            end do
          end do
        end do

      case("geo")
        ! Compute pu
        do panel = 1, nbfaces
          do i = 0, ntotal
            do j = 1, ntotal
              mesh%pu(i,j,panel)%p = midpoint(mesh%po(i,j,panel)%p, mesh%po(i,j-1,panel)%p)
            end do
          end do
        end do

        ! Compute pv
        do panel = 1, nbfaces
          do i = 1, ntotal
            do j = 0, ntotal
              mesh%pv(i,j,panel)%p = midpoint(mesh%po(i,j,panel)%p, mesh%po(i-1,j,panel)%p)
            end do
          end do
        end do

        ! Compute pc
        do panel = 1, nbfaces
          do i = 1, ntotal
            do j = 1, ntotal
              pc1 = midpoint(mesh%pu(i-1,j,panel)%p, mesh%pu(i,j,panel)%p)
              pc2 = midpoint(mesh%pv(i,j-1,panel)%p, mesh%pv(i,j,panel)%p)
              pc  = midpoint(pc1, pc2)
              mesh%pc(i,j,panel)%p = pc
              call cart2sph(mesh%pc(i,j,panel)%p(1), mesh%pc(i,j,panel)%p(2), mesh%pc(i,j,panel)%p(3), &
                            mesh%pc(i,j,panel)%lat ,mesh%pc(i,j,panel)%lon)
            end do
          end do
        end do

      case default
        print*, 'ERROR on compute_midpoints: invalid midpoint position: ', mesh%midpos
        stop
      end select


      ! Convert po to latlon
      do panel = 1, nbfaces
        do i = 0, ntotal
          do j = 0, ntotal
            call cart2sph(mesh%po(i,j,panel)%p(1), mesh%po(i,j,panel)%p(2), mesh%po(i,j,panel)%p(3), &
                          mesh%po(i,j,panel)%lon , mesh%po(i,j,panel)%lat)
          end do
        end do
      end do

      ! Convert pc to latlon
      do panel = 1, nbfaces
        do i = 1, ntotal
          do j = 1, ntotal
            call cart2sph(mesh%pc(i,j,panel)%p(1), mesh%pc(i,j,panel)%p(2), mesh%pc(i,j,panel)%p(3), &
                          mesh%pc(i,j,panel)%lon , mesh%pc(i,j,panel)%lat)
          end do
        end do
      end do

      ! Convert pu to latlon
      do panel = 1, nbfaces
        do i = 0, ntotal
          do j = 1, ntotal
            call cart2sph(mesh%pu(i,j,panel)%p(1), mesh%pu(i,j,panel)%p(2), mesh%pu(i,j,panel)%p(3), &
                          mesh%pu(i,j,panel)%lon , mesh%pu(i,j,panel)%lat)
          end do
        end do
      end do

      ! Convert pv to latlon
      do panel = 1, nbfaces
        do i = 1, ntotal
          do j = 0, ntotal
            call cart2sph(mesh%pv(i,j,panel)%p(1), mesh%pv(i,j,panel)%p(2), mesh%pv(i,j,panel)%p(3), &
                          mesh%pv(i,j,panel)%lon , mesh%pv(i,j,panel)%lat)
          end do
        end do
      end do


  end subroutine compute_midpoints 

  subroutine compute_localcoords(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_LOCALCOORDS
    !
    ! Given the cell points po, pc, pu and pv, this routines compute
    ! their local coordinates
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter

    ! Real aux vars
    real(r8) :: x ! Local coordinates
    real(r8) :: y ! Local coordinates
    real(r8) :: p(1:3) ! sphere point

    ntotal = mesh%ntotal

    print*, 'Computing cubed-sphere local coordinates...'

    ! Compute po
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 1, ntotal
         ! print*,panel
          p =  mesh%pc(i,j,panel)%p
        !  call inverse_equidistant_gnomonic_map(x, y, p, mesh)
        !  mesh%po(i,j,panel)%x = datan2(x)
        !  mesh%po(i,j,panel)%y = datan2(y)
        !  print*,x*rad2deg,y*rad2deg
        end do
      end do
    end do
!stop

    ! Compute pu
    do panel = 1, nbfaces
      do i = 0, ntotal
        do j = 1, ntotal
          !mesh%pu(i,j,panel)%x
        end do
      end do
    end do

    ! Compute pv
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 0, ntotal
          !mesh%pv(i,j,panel)%x
        end do
      end do
    end do

    ! Compute pc
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 1, ntotal
          !mesh%pc(i,j,panel)%p
        end do
      end do
    end do


  end subroutine compute_localcoords 

  subroutine latlon_grid(mesh)
    !---------------------------------------------------------------------
    !
    ! LATLON_GRID
    !
    ! Creates a latlon grid and its nearest neighbors indexes
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    integer(i4) :: i, j !counters
    integer(i4) :: ix, jy, panel !latlon point index on the cubed-sphere

    ! Latlon grid spacing
    mesh%dlon = 2._r8*pi/mesh%nlon
    mesh%dlat = pi/mesh%nlat

    print*, 'Generating latlon grid...'
    ! Creates the ll grid
    do i = 0, mesh%nlon
      do j = 0, mesh%nlat
        mesh%ll(i,j)%lon = -pi   + i*mesh%dlon
        mesh%ll(i,j)%lat = -pio2 + j*mesh%dlat
        call sph2cart(mesh%ll(i,j)%lon, mesh%ll(i,j)%lat, mesh%ll(i,j)%p(1), mesh%ll(i,j)%p(2), mesh%ll(i,j)%p(3))
      end do
    end do

    print*, 'Perfoming binary search: latlon to cubed sphere...'
    ! Get nearest neighbor indexes
    do i = 0, mesh%nlon
      do j = 0, mesh%nlat
        call binary_search(mesh%ll(i,j)%p, mesh, ix, jy, panel) 
        mesh%panels_ll(i,j) = panel
        mesh%ix_ll(i,j) = ix
        mesh%jy_ll(i,j) = jy
      end do
    end do

    ! Deallocation
    deallocate(mesh%ll)
  end subroutine

end module cubed_sphere 

