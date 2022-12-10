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
      vector, &
      matrix, &
      cubedsphere

  !Input routines
  use input, only: &
      meshload, &
      meshread

  !Output routines
  use output, only: &
      meshstore, &
      printmesh
 
  !Data allocation
  use allocation, only: &
      meshallocation, &
      r8_1darray_allocation

  !Data deallocation
  use deallocation, only: &
      tgvectors_deallocation 

  !Spherical geometry
  use sphgeo, only: &
      equidistant_gnomonic_map, &
      norm, &
      midpoint, &
      sphquadarea, &
      inverse_equidistant_gnomonic_map, &
      sph2cart, &
      cart2sph, &
      binary_search, &
      arcintersec, &
      arclen, &
      proj_vec_sphere, &
      tangent_ll_lat, &
      tangent_ll_lon, &
      derivative_ydir_equiangular_gnomonic_map, &
      derivative_xdir_equiangular_gnomonic_map, &
      cross_product

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

      case("read") !Read vertices from file
         call meshread(mesh)

      case default
         print*, "MESH BUILD ERROR: Invalid mesh kind : ", mesh%kind
         stop
         return
      end select

      ! Generate the grid properties
      call cubedsphere_properties(mesh)

      ! Deallocate tg vectors and local coordinates
      !call tgvectors_deallocation(mesh)

      ! Create lat/lon grid
      call latlon_grid(mesh)

      ! Store the grid
      call meshstore(mesh, header)
  end if

    ! Print grid features on screen
    call printmesh(mesh)

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

    !---------------------------------------------------------------------
    ! Panel indexes distribution
    !      +---+
    !      | 5 |
    !  +---+---+---+---+
    !  | 4 | 1 | 2 | 3 |
    !  +---+---+---+---+
    !      | 6 |
    !      +---+
    !---------------------------------------------------------------------
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

    !Total number of cells (ignoring ghost cells)
    mesh%nbcells = 6*n*n
 
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
    integer(i4) :: p! Panel counter
    integer(i4) :: i0, j0, iend, jend ! 2D grid interior indexes

    ! Number of cells along a coordinate axis
    n = mesh%n
 
    ! Number of cells along a coordinate axis + number of ghost cells
    ntotal = mesh%ntotal

    ! Compute the midpoints
    call compute_midpoints(mesh)

    ! Compute the local coordinates
    call compute_localcoords(mesh)

    ! Compute tangent vectors
    call compute_tgvectors(mesh)

    ! Compute angles
    call compute_angles(mesh)

    ! Compute conversions between contravariant and latlon representation of vectors
    call compute_ll2contra(mesh)


    print*, 'Computing cubed-sphere areas, lenghts...'

    do p = 1, nbfaces
      do i = 1, ntotal
        do j = 1, ntotal
          ! Compute the quadrilateral areas
          mesh%area(i,j,p) = sphquadarea(mesh%po(i-1,j-1,1)%p, mesh%po(i-1,j,1)%p, mesh%po(i,j-1,1)%p, mesh%po(i,j,1)%p)

          ! Compute lenghts in x direction
          mesh%lx(i,j,p) = arclen(mesh%pu(i-1,j,1)%p, mesh%pu(i,j,1)%p)

          ! Compute lenghts in y direction
          mesh%ly(i,j,p) = arclen(mesh%pv(i,j-1,1)%p, mesh%pv(i,j,1)%p)
 
        end do
      end do
    end do


    ! save some stats
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    mesh%minarea  = minval(mesh%area(i0:iend,j0:jend,:))
    mesh%maxarea  = maxval(mesh%area(i0:iend,j0:jend,:))
    mesh%meanarea = sum(mesh%area(i0:iend,j0:jend,:))/mesh%nbcells

    mesh%mindist  = minval(mesh%lx(i0-1:iend,j0:jend,:))
    mesh%maxdist  = maxval(mesh%lx(i0-1:iend,j0:jend,:))
    mesh%meandist = sum(mesh%lx(i0-1:iend,j0:jend,:))/mesh%nbcells
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
    real(r8) :: pc(1:3), p1(1:3), p2(1:3), q1(1:3), q2(1:3)

    !logical 
    logical :: intersection_existence

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
              p1 = mesh%pu(i-1,j  , panel)%p
              p2 = mesh%pu(i  ,j  , panel)%p
              q1 = mesh%pv(i  ,j-1, panel)%p
              q2 = mesh%pv(i  ,j  , panel)%p
              intersection_existence = arcintersec(p1, p2, q1, q2, pc)
              if (intersection_existence)then
                mesh%pc(i,j,panel)%p = pc
              else
                print*, 'ERROR on compute_midpoints: geodesics do not intersect.'
                stop
              end if
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
    ! Given the vertices pv, this routines compute
    ! their local coordinates using the inverse of the gnomonic mapping
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p ! Panel counter

    ! Real aux vars
    real(r8) :: x ! Local coordinates
    real(r8) :: y ! Local coordinates
    real(r8) :: q(1:3) ! sphere point

    ntotal = mesh%ntotal

    print*, 'Computing cubed-sphere local coordinates...'

    ! Compute po
    do p = 1, nbfaces
      do i = 0, ntotal
        do j = 0, ntotal
          q =  mesh%po(i,j,p)%p
          !call inverse_equidistant_gnomonic_map(x, y, q, mesh, p)
          !mesh%x_po(i,j,p) = datan2(x)
          !mesh%y_po(i,j,p) = datan2(y)
          !print*,mesh%x_po(i,j,p)*rad2deg, mesh%y_po(i,j,p)*rad2deg, mesh%dy*rad2deg
        end do
        !stop
      end do
    end do

  end subroutine compute_localcoords 

  subroutine compute_tgvectors(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_TGVECTORS
    !
    ! This routine computes the tangent vector in x and y direction
    ! at pu and pv points
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter

    ! Real aux vars
    real(r8) :: p(1:3) ! sphere point
    real(r8) :: v(1:3) ! vector point
    real(r8) :: v1(1:3) ! vector point
    real(r8) :: v2(1:3) ! vector point
    real(r8) :: projv_p(1:3) ! vector point

    real(r8), allocatable :: x(:), y(:)
    real(r8) :: error, error1

    ntotal = mesh%ntotal

    ! Allocation
    call r8_1darray_allocation(x, 0, ntotal)
    call r8_1darray_allocation(y, 1, ntotal)

    ! Init local coordinates grid
    x(0) = -pio4 - mesh%nbgl*mesh%dx
    do i = 1, ntotal
      x(i) = x(i-1) + mesh%dx
    end do
    y = (x(0:ntotal-1)+x(1:ntotal))*0.5d0

    print*, 'Computing cubed-sphere tg vectors...'

    error = 0._r8

    ! Compute tgx at pu
    do panel = 1, nbfaces
      do i = 0, ntotal-1
        do j = 1, ntotal
          v = mesh%pu(i+1,j,panel)%p - mesh%pu(i,j,panel)%p
          p = mesh%pu(i,j,panel)%p
          projv_p = proj_vec_sphere(v, p)
          v1 = projv_p/norm(projv_p)

          call derivative_xdir_equiangular_gnomonic_map(x(i), y(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgx_pu(i,j,panel)%v = v2 
        end do
      end do
    end do

    do panel = 1, nbfaces
      do j = 1, ntotal
        v = mesh%pu(ntotal-1,j,panel)%p - mesh%pu(ntotal,j,panel)%p
        p = mesh%pu(ntotal,j,panel)%p
        projv_p = proj_vec_sphere(v, p)
        v1 = -projv_p/norm(projv_p)

        call derivative_xdir_equiangular_gnomonic_map(x(i), y(j), v2, panel)
        !v2 = v2/norm(v2)
        !error1 = norm(v1-v2)
        !error = max(error, error1)

        mesh%tgx_pu(ntotal,j,panel)%v = v2
     end do
    end do

    ! Compute tgy at pu
     do panel = 1, nbfaces
      do i = 0, ntotal
        do j = 1, ntotal
          v = mesh%po(i,j,panel)%p - mesh%pu(i,j,panel)%p
          p = mesh%pu(i,j,panel)%p
          projv_p = proj_vec_sphere(v, p) 
          v1 = projv_p/norm(projv_p) 

          call derivative_ydir_equiangular_gnomonic_map(x(i), y(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgy_pu(i,j,panel)%v = v2
       end do
      end do
    end do

   ! Compute tgy at pv
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 0, ntotal-1
          v = mesh%pv(i,j+1,panel)%p - mesh%pv(i,j,panel)%p
          p = mesh%pv(i,j,panel)%p 
          projv_p = proj_vec_sphere(v, p)
          v1 = projv_p/norm(projv_p)

          call derivative_ydir_equiangular_gnomonic_map(y(i), x(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgy_pv(i,j,panel)%v = v2
        end do
      end do
    end do

    do panel = 1, nbfaces
      do i = 1, ntotal
        v = mesh%pv(i,ntotal-1,panel)%p - mesh%pv(i,ntotal,panel)%p
        p = mesh%pv(i,ntotal,panel)%p 
        projv_p = proj_vec_sphere(v, p)
        v1 = -projv_p/norm(projv_p)

        call derivative_ydir_equiangular_gnomonic_map(y(i), x(j), v2, panel)
        !v2 = v2/norm(v2)
        !error1 = norm(v1-v2)
        !error = max(error, error1)

        mesh%tgy_pv(i,j,panel)%v = v2
 
     end do
    end do



   ! Compute tgx at pv
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 0, ntotal
          v = mesh%po(i,j,panel)%p - mesh%pv(i,j,panel)%p
          p = mesh%pv(i,j,panel)%p 
          projv_p = proj_vec_sphere(v, p)
          v1 = projv_p/norm(projv_p)

          call derivative_xdir_equiangular_gnomonic_map(y(i), x(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgx_pv(i,j,panel)%v = v2
       end do
      end do
    end do
    
    ! Compute tgx at pc
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 1, ntotal
          v = mesh%pu(i,j,panel)%p - mesh%pc(i,j,panel)%p
          p = mesh%pc(i,j,panel)%p 
          projv_p = proj_vec_sphere(v, p)
          v1 = projv_p/norm(projv_p)

          call derivative_xdir_equiangular_gnomonic_map(y(i), y(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgx_pc(i,j,panel)%v = v2
       end do
      end do
    end do


   ! Compute tgy at pc
    do panel = 1, nbfaces
      do i = 1, ntotal
        do j = 1, ntotal
          v = mesh%pv(i,j,panel)%p - mesh%pc(i,j,panel)%p
          p = mesh%pc(i,j,panel)%p 
          projv_p = proj_vec_sphere(v, p)
          v1 = projv_p/norm(projv_p)

          call derivative_ydir_equiangular_gnomonic_map(y(i), y(j), v2, panel)
          !v2 = v2/norm(v2)
          !error1 = norm(v1-v2)
          !error = max(error, error1)

          mesh%tgy_pc(i,j,panel)%v = v2
        end do
      end do
    end do

  end subroutine compute_tgvectors 

 

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

  subroutine compute_angles(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_ANGLES
    !
    ! This routine computes the angles at pc, pu, pv points
    ! by calculating the tangent vector cross products
    ! and dot produts
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: ntotal, n
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p! Panel counter

    ! Real aux vars
    real(r8) :: v1(1:3), v2(1:3), v3(1:3)  ! vectors

    ntotal = mesh%ntotal

    print*, 'Computing cubed-sphere angles...'

    ! Compute angles at pc
    do p = 1, nbfaces 
      do i = 1, ntotal
        do j = 1, ntotal
          ! Cross product
          v1 = mesh%tgx_pc(i,j,p)%v
          v2 = mesh%tgy_pc(i,j,p)%v
          v3 = cross_product(v1,v2)
          mesh%sinc(i,j,p) = norm(v3)
       end do
      end do
    end do

    ! Compute angles a pu
    do p = 1, nbfaces
      do i = 0, ntotal
        do j = 1, ntotal
          ! Cross product
          v1 = mesh%tgx_pu(i,j,p)%v
          v2 = mesh%tgy_pu(i,j,p)%v
          v3 = cross_product(v1,v2)
          mesh%sinu(i,j,p) = norm(v3)

          ! Dot product
          mesh%cosu(i,j,p) = dot_product(v1,v2)
       end do
      end do
    end do

   ! Compute angles at pv
    do p = 1, nbfaces
      do i = 1, ntotal
        do j = 0, ntotal
          ! Cross product
          v1 = mesh%tgx_pv(i,j,p)%v
          v2 = mesh%tgy_pv(i,j,p)%v
          v3 = cross_product(v1,v2)
          mesh%sinv(i,j,p) = norm(v3)

          ! Dot product
          mesh%cosv(i,j,p) = dot_product(v1,v2)
       end do
      end do
    end do

    end subroutine compute_angles

    subroutine compute_ll2contra(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_ll2contra
    !
    ! This routine computes the matrices that performs
    ! the conversions between contravariant and latlon representation
    ! of tangent vectors on the sphere
    !---------------------------------------------------------------------
      type(cubedsphere), intent(inout) :: mesh

      ! Integer auxs
      integer(i4) :: ntotal
      integer(i4) :: i, j ! 2D grid counters
      integer(i4) :: p ! Panel counter

      ! Real aux vars
      real(r8) :: elon(1:3), elat(1:3), ex(1:3), ey(1:3) ! vectors
      real(r8) :: lat, lon
      real(r8) :: a11, a12, a21, a22, det

      ntotal = mesh%ntotal

      print*, 'Computing conversion matrices...'

      ! Compute at pu
      do p = 1, nbfaces
        do i = 0, ntotal
          do j = 1, ntotal
            ex = mesh%tgx_pu(i,j,p)%v
            ey = mesh%tgy_pu(i,j,p)%v

            lat = mesh%pu(i,j,p)%lat
            lon = mesh%pu(i,j,p)%lon

            call tangent_ll_lon(lon, elon)
            call tangent_ll_lat(lon, lat, elat)

            a11 = dot_product(ex, elon)  
            a12 = dot_product(ey, elon)  
            a21 = dot_product(ex, elat)  
            a22 = dot_product(ey, elat) 
            det = a11*a22 - a21*a12

            ! Contra to latlon matrix
            mesh%contra2ll_pu(i,j,p)%M(1,1) = a11
            mesh%contra2ll_pu(i,j,p)%M(1,2) = a12
            mesh%contra2ll_pu(i,j,p)%M(2,1) = a21
            mesh%contra2ll_pu(i,j,p)%M(2,2) = a22

            ! latlon to contra matrix
            mesh%ll2contra_pu(i,j,p)%M(1,1) =  a22
            mesh%ll2contra_pu(i,j,p)%M(1,2) = -a12
            mesh%ll2contra_pu(i,j,p)%M(2,1) = -a21
            mesh%ll2contra_pu(i,j,p)%M(2,2) =  a11
            mesh%ll2contra_pu(i,j,p)%M(:,:) = mesh%ll2contra_pu(i,j,p)%M(:,:)/det
 
         end do
        end do
      end do

      ! Compute at pv
      do p = 1, nbfaces
        do i = 1, ntotal
          do j = 0, ntotal
            ex = mesh%tgx_pv(i,j,p)%v
            ey = mesh%tgy_pv(i,j,p)%v

            lat = mesh%pv(i,j,p)%lat
            lon = mesh%pv(i,j,p)%lon

            call tangent_ll_lon(lon, elon)
            call tangent_ll_lat(lon, lat, elat)

            a11 = dot_product(ex, elon)  
            a12 = dot_product(ey, elon)  
            a21 = dot_product(ex, elat)  
            a22 = dot_product(ey, elat) 
            det = a11*a22 - a21*a12

            ! Contra to latlon matrix
            mesh%contra2ll_pv(i,j,p)%M(1,1) = a11
            mesh%contra2ll_pv(i,j,p)%M(1,2) = a12
            mesh%contra2ll_pv(i,j,p)%M(2,1) = a21
            mesh%contra2ll_pv(i,j,p)%M(2,2) = a22

            ! latlon to contra matrix
            mesh%ll2contra_pv(i,j,p)%M(1,1) =  a22
            mesh%ll2contra_pv(i,j,p)%M(1,2) = -a12
            mesh%ll2contra_pv(i,j,p)%M(2,1) = -a21
            mesh%ll2contra_pv(i,j,p)%M(2,2) =  a11
            mesh%ll2contra_pv(i,j,p)%M(:,:) = mesh%ll2contra_pv(i,j,p)%M(:,:)/det
         end do
        end do
      end do

      end subroutine compute_ll2contra


end module cubed_sphere 

