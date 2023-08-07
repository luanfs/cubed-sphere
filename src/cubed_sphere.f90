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
    showonscreen, &
    i0, iend, &
    j0, jend, &
    n0, nend, &
    nghost, &
    hs

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
    inverse_equiangular_gnomonic_map, &
    sph2cart, &
    cart2sph, &
    tangent_ll_lat, &
    tangent_ll_lon, &
    derivative_ydir_equiangular_gnomonic_map, &
    derivative_xdir_equiangular_gnomonic_map

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
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    !Auxiliar variables
    character (len=256):: header
    logical:: ifile

    !Total number of cells (ignoring ghost cells)
    mesh%nbcells = 6*mesh%n**2
 
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
    header=trim(griddir)//trim(mesh%name)//"_pc.dat"
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

        ! Generate the grid properties
        call cubedsphere_properties(mesh)

    else
        !------------------------------------------------
        !Generate mesh
        !------------------------------------------------
        select case(trim(mesh%kind))
        case("equiangular") !Equiangular grid
            call equiangular_cubedsphere_generation(mesh)

            ! Compute tangent vectors
            call compute_tgvectors(mesh)

        case("read") !Read vertices from file
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

    ! Deallocate tg vectors and local coordinates
    call tgvectors_deallocation(mesh)

    ! save some stats
    mesh%minarea  = minval(mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx*mesh%dy
    mesh%maxarea  = maxval(mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx*mesh%dy
    mesh%meanarea = sum(mesh%mt_pc(i0:iend,j0:jend,:))/mesh%nbcells*mesh%dx*mesh%dy

    mesh%mindist  = minval(mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx
    mesh%maxdist  = maxval(mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx
    mesh%meandist = sum(mesh%mt_pc(i0:iend,j0:jend,:))/mesh%nbcells*mesh%dx

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

    !mesh%name = trim(mesh%name)

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

    ! Halo size
    mesh%halosize = 4
    hs = mesh%halosize 
    ! Interior indexes
    mesh%i0   = 1
    mesh%iend = mesh%n
    mesh%j0   = 1 
    mesh%jend = mesh%n

    ! Indexes for cells along a coordinate axis including ghost cells
    mesh%n0 = -mesh%halosize+1
    mesh%nend =  mesh%n + mesh%halosize

    ! global vars
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend
    nghost = mesh%halosize
    n0 = mesh%n0
    nend = mesh%nend
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
    real(kind=8), allocatable :: x(:) ! Local coordinates (angular)
    real(kind=8), allocatable :: y(:) ! Local coordinates (angular)
    real(kind=8), allocatable :: tanx(:)
    real(kind=8), allocatable :: tany(:)

    ! Integer auxs
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter

    print*, 'Computing cubed-sphere vertices...'

    ! Allocation
    call r8_1darray_allocation(x, n0, nend+1)
    call r8_1darray_allocation(y, n0, nend+1)
    call r8_1darray_allocation(tanx, n0, nend+1)
    call r8_1darray_allocation(tany, n0, nend+1)

    ! Local coordinates grid size
    mesh%dx = pio2/mesh%n
    mesh%dy = pio2/mesh%n

    ! Init local coordinates grid
    x(n0) = -pio4 - mesh%halosize*mesh%dx
    do i = n0+1, nend+1
        x(i) = x(i-1) + mesh%dx
    end do
    y = x

    ! Angular coordinates are mapped to the cube
    tanx = dtan(x)
    tany = tanx

    ! Compute the gnomonic mapping at the vertices
    do panel = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                call equidistant_gnomonic_map(acube*tanx(i), acube*tany(j), mesh%po(i,j,panel)%p, panel)
            end do
        end do
    end do

    ! Deallocation
    deallocate(x, y, tanx, tany)

    print*, 'Computing cubed-sphere midpoints...'

    ! Allocation
    call r8_1darray_allocation(x, n0, nend+1)
    call r8_1darray_allocation(y, n0, nend)
    call r8_1darray_allocation(tanx, n0, nend+1)
    call r8_1darray_allocation(tany, n0, nend)

    ! Init local coordinates grid
    x(n0) = -pio4 - mesh%halosize*mesh%dx
    do i = n0+1, nend+1
        x(i) = x(i-1) + mesh%dx
    end do
    y = (x(n0:nend)+x(n0+1:nend+1))*0.5d0

    ! Angular coordinates are mapped to the cube
    tanx = dtan(x)
    tany = dtan(y)

    ! Compute pu
    do panel = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                call equidistant_gnomonic_map(acube*tanx(i), acube*tany(j), mesh%pu(i,j,panel)%p, panel)
            end do
        end do
    end do

    ! Compute pv
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                call equidistant_gnomonic_map(acube*tany(i), acube*tanx(j), mesh%pv(i,j,panel)%p, panel)
            end do
        end do
    end do

    ! Compute pc
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                call equidistant_gnomonic_map(acube*tany(i), acube*tany(j), mesh%pc(i,j,panel)%p , panel)
            end do
        end do
    end do

    ! Convert po to latlon
    do panel = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                call cart2sph(mesh%po(i,j,panel)%p(1), mesh%po(i,j,panel)%p(2), mesh%po(i,j,panel)%p(3), &
                          mesh%po(i,j,panel)%lon , mesh%po(i,j,panel)%lat)
            end do
        end do
    end do

    ! Convert pc to latlon
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                call cart2sph(mesh%pc(i,j,panel)%p(1), mesh%pc(i,j,panel)%p(2), mesh%pc(i,j,panel)%p(3), &
                          mesh%pc(i,j,panel)%lon , mesh%pc(i,j,panel)%lat)
            end do
        end do
    end do

    ! Convert pu to latlon
    do panel = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                call cart2sph(mesh%pu(i,j,panel)%p(1), mesh%pu(i,j,panel)%p(2), mesh%pu(i,j,panel)%p(3), &
                          mesh%pu(i,j,panel)%lon , mesh%pu(i,j,panel)%lat)
            end do
        end do
    end do

    ! Convert pv to latlon
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                call cart2sph(mesh%pv(i,j,panel)%p(1), mesh%pv(i,j,panel)%p(2), mesh%pv(i,j,panel)%p(3), &
                          mesh%pv(i,j,panel)%lon , mesh%pv(i,j,panel)%lat)
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
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p ! Panel counter

    ! Compute metric tensor
    call compute_metric_tensor(mesh)

    ! Compute conversions between contravariant, contravariant and latlon representation of vectors
    call compute_conversion_matrices(mesh)
        

end subroutine cubedsphere_properties


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
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: panel    ! Panel counter

    ! Real aux vars
    real(kind=8) :: v(1:3) ! vector point
    real(kind=8), allocatable :: x(:), y(:)

    ! Allocation
    call r8_1darray_allocation(x, n0, nend+1)
    call r8_1darray_allocation(y, n0, nend)

    ! Init local coordinates grid
    x(n0) = -pio4 - mesh%halosize*mesh%dx
    do i = n0+1,nend+1
        x(i) = x(i-1) + mesh%dx
    end do
    y(n0:nend) = (x(n0:nend)+x(n0+1:nend+1))*0.5d0

    print*, 'Computing cubed-sphere tg vectors...'

    ! Compute tgx/tgy at pu
    do panel = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                call derivative_xdir_equiangular_gnomonic_map(x(i), y(j), v, panel)
                mesh%tgx_pu(i,j,panel)%v = v
                call derivative_ydir_equiangular_gnomonic_map(x(i), y(j), v, panel)
                mesh%tgy_pu(i,j,panel)%v = v
            end do
        end do
    end do

    ! Compute tgx/tgy at pv
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                call derivative_ydir_equiangular_gnomonic_map(y(i), x(j), v, panel)
                mesh%tgy_pv(i,j,panel)%v = v
                call derivative_xdir_equiangular_gnomonic_map(y(i), x(j), v, panel)
                mesh%tgx_pv(i,j,panel)%v = v
            end do
        end do
    end do

    ! Compute tgx/tgy at pc
    do panel = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                call derivative_xdir_equiangular_gnomonic_map(y(i), y(j), v, panel)
                mesh%tgx_pc(i,j,panel)%v = v
                call derivative_ydir_equiangular_gnomonic_map(y(i), y(j), v, panel)
                mesh%tgy_pc(i,j,panel)%v = v
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
    real(kind=8) :: x, y
    real(kind=8), allocatable :: xx(:), yy(:)

    ! Latlon grid spacing
    mesh%dlon = 2d0*pi/mesh%nlon
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

    print*, 'Search indexes: latlon to cubed sphere...'
    allocate(xx(1:mesh%n))
    allocate(yy(1:mesh%n))
    do i = 1, mesh%n
        xx(i) = -pio4 + (i-0.5)*mesh%dx
        yy(i) = -pio4 + (i-0.5)*mesh%dx
    end do
    ! Get nearest neighbor indexes
    do i = 0, mesh%nlon
        do j = 0, mesh%nlat
            call inverse_equiangular_gnomonic_map(x, y, panel, mesh%ll(i,j)%p, mesh)
            !ix = floor((x+pio4-mesh%dx*0.d0)/mesh%dx)
            !jy = floor((y+pio4-mesh%dx*0.d0)/mesh%dy)
            ix = minloc(abs(xx(:)-x),DIM=1)
            jy = minloc(abs(yy(:)-y),DIM=1)
            mesh%panels_ll(i,j) = panel
            mesh%ix_ll(i,j) = ix
            mesh%jy_ll(i,j) = jy

            !print*,ix,jy,x*rad2deg,y*rad2deg,panel
            !if(abs(x)*rad2deg>45d0 .or. abs(y)*rad2deg>45._kind=8)then
            !    print*, 'error1'
            !    stop
            !end if
            !if(ix<i0.or. ix>iend)then
            !    print*, 'error2'
            !    stop
            !end if
            !if(jy<j0 .or. jy>jend)then
            !    print*, 'error3'
            !    stop
            !end if
 
        end do
    end do

    ! Deallocation
    deallocate(mesh%ll)
    deallocate(xx, yy)
end subroutine

subroutine compute_metric_tensor(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_METRIC_TENSOR
    !
    ! This routine computes the metric tensor at pc, pu, pv points
    ! by calculating the tangent vector dot products
    !
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p! Panel counter

    ! Real aux vars
    real(kind=8) :: v1(1:3), v2(1:3) ! vectors

    print*, 'Computing cubed-sphere metric tensor...'

    ! Compute metric tensor at pc
    do p = 1, nbfaces 
        do i = n0, nend
            do j = n0, nend
                v1 = mesh%tgx_pc(i,j,p)%v
                v2 = mesh%tgy_pc(i,j,p)%v
                mesh%mt_pc(i,j,p) = (norm2(v1)**2)*(norm2(v2)**2) - dot_product(v1,v2)**2
                mesh%mt_pc(i,j,p) = dsqrt(mesh%mt_pc(i,j,p))
            end do
        end do
    end do

    ! Compute angles a pu
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                v1 = mesh%tgx_pu(i,j,p)%v
                v2 = mesh%tgy_pu(i,j,p)%v
                mesh%mt_pu(i,j,p) = (norm2(v1)**2)*(norm2(v2)**2) - dot_product(v1,v2)**2
                mesh%mt_pu(i,j,p) = dsqrt(mesh%mt_pu(i,j,p))
            end do
        end do
    end do

    ! Compute angles at pv
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                v1 = mesh%tgx_pv(i,j,p)%v
                v2 = mesh%tgy_pv(i,j,p)%v
                mesh%mt_pv(i,j,p) = (norm2(v1)**2)*(norm2(v2)**2) - dot_product(v1,v2)**2
                mesh%mt_pv(i,j,p) = dsqrt(mesh%mt_pv(i,j,p))
            end do
        end do
    end do

end subroutine compute_metric_tensor


subroutine compute_conversion_matrices(mesh)
    !---------------------------------------------------------------------
    !
    ! COMPUTE_CONVERSION_MATRICES
    !
    ! This routine computes the matrices that performs
    ! the conversions between contravariant, covariant and latlon representation
    ! of tangent vectors on the sphere
    !---------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Integer auxs
    integer(i4) :: i, j ! 2D grid counters
    integer(i4) :: p ! Panel counter

    ! Real aux vars
    real(kind=8) :: elon(1:3), elat(1:3), ex(1:3), ey(1:3) ! vectors
    real(kind=8) :: lat, lon
    real(kind=8) :: a11, a12, a21, a22, det

    print*, 'Computing conversion matrices...'

    ! Contravariant/latlon
    ! Compute at pu
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
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
        do i = n0, nend
            do j = n0, nend+1
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

    ! Compute at pc
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                ex = mesh%tgx_pc(i,j,p)%v
                ey = mesh%tgy_pc(i,j,p)%v

                lat = mesh%pc(i,j,p)%lat
                lon = mesh%pc(i,j,p)%lon

                call tangent_ll_lon(lon, elon)
                call tangent_ll_lat(lon, lat, elat)

                a11 = dot_product(ex, elon)
                a12 = dot_product(ey, elon)
                a21 = dot_product(ex, elat)
                a22 = dot_product(ey, elat)
                det = a11*a22 - a21*a12

                ! Contra to latlon matrix
                mesh%contra2ll_pc(i,j,p)%M(1,1) = a11
                mesh%contra2ll_pc(i,j,p)%M(1,2) = a12
                mesh%contra2ll_pc(i,j,p)%M(2,1) = a21
                mesh%contra2ll_pc(i,j,p)%M(2,2) = a22

                ! latlon to contra matrix
                mesh%ll2contra_pc(i,j,p)%M(1,1) =  a22
                mesh%ll2contra_pc(i,j,p)%M(1,2) = -a12
                mesh%ll2contra_pc(i,j,p)%M(2,1) = -a21
                mesh%ll2contra_pc(i,j,p)%M(2,2) =  a11
                mesh%ll2contra_pc(i,j,p)%M(:,:) = mesh%ll2contra_pc(i,j,p)%M(:,:)/det
            end do
        end do
    end do

    ! Contravariant/covariant
    ! Compute at pu
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                ex = mesh%tgx_pu(i,j,p)%v
                ey = mesh%tgy_pu(i,j,p)%v

                a11 = dot_product(ex, ex)
                a12 = dot_product(ey, ex)
                a21 = dot_product(ex, ey)
                a22 = dot_product(ey, ey)
                det = a11*a22 - a21*a12

                ! Contra to covari matrix
                mesh%contra2covari_pu(i,j,p)%M(1,1) = a11
                mesh%contra2covari_pu(i,j,p)%M(1,2) = a12
                mesh%contra2covari_pu(i,j,p)%M(2,1) = a21
                mesh%contra2covari_pu(i,j,p)%M(2,2) = a22

                ! Covari to contra matrix
                mesh%covari2contra_pu(i,j,p)%M(1,1) =  a22
                mesh%covari2contra_pu(i,j,p)%M(1,2) = -a12
                mesh%covari2contra_pu(i,j,p)%M(2,1) = -a21
                mesh%covari2contra_pu(i,j,p)%M(2,2) =  a11
                mesh%covari2contra_pu(i,j,p)%M(:,:) = mesh%covari2contra_pu(i,j,p)%M(:,:)/det
            end do
        end do
    end do

    ! Compute at pv
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                ex = mesh%tgx_pv(i,j,p)%v
                ey = mesh%tgy_pv(i,j,p)%v

                a11 = dot_product(ex, ex)
                a12 = dot_product(ey, ex)
                a21 = dot_product(ex, ey)
                a22 = dot_product(ey, ey)
                det = a11*a22 - a21*a12

                ! Contra to covari matrix
                mesh%contra2covari_pv(i,j,p)%M(1,1) = a11
                mesh%contra2covari_pv(i,j,p)%M(1,2) = a12
                mesh%contra2covari_pv(i,j,p)%M(2,1) = a21
                mesh%contra2covari_pv(i,j,p)%M(2,2) = a22

                ! Covari to contra matrix
                mesh%covari2contra_pv(i,j,p)%M(1,1) =  a22
                mesh%covari2contra_pv(i,j,p)%M(1,2) = -a12
                mesh%covari2contra_pv(i,j,p)%M(2,1) = -a21
                mesh%covari2contra_pv(i,j,p)%M(2,2) =  a11
                mesh%covari2contra_pv(i,j,p)%M(:,:) = mesh%covari2contra_pv(i,j,p)%M(:,:)/det
            end do
        end do
    end do

    ! Compute at pc
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                ex = mesh%tgx_pc(i,j,p)%v
                ey = mesh%tgy_pc(i,j,p)%v

                a11 = dot_product(ex, ex)
                a12 = dot_product(ey, ex)
                a21 = dot_product(ex, ey)
                a22 = dot_product(ey, ey)
                det = a11*a22 - a21*a12

                ! Contra to covari matrix
                mesh%contra2covari_pc(i,j,p)%M(1,1) = a11
                mesh%contra2covari_pc(i,j,p)%M(1,2) = a12
                mesh%contra2covari_pc(i,j,p)%M(2,1) = a21
                mesh%contra2covari_pc(i,j,p)%M(2,2) = a22

                ! Covari to contra matrix
                mesh%covari2contra_pc(i,j,p)%M(1,1) =  a22
                mesh%covari2contra_pc(i,j,p)%M(1,2) = -a12
                mesh%covari2contra_pc(i,j,p)%M(2,1) = -a21
                mesh%covari2contra_pc(i,j,p)%M(2,2) =  a11
                mesh%covari2contra_pc(i,j,p)%M(:,:) = mesh%covari2contra_pc(i,j,p)%M(:,:)/det
            end do
        end do
    end do



    ! latlon/covariant
    ! Compute at pu
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                ! latlon to covari matrix
                mesh%ll2covari_pu(i,j,p)%M = matmul(mesh%contra2covari_pu(i,j,p)%M, mesh%ll2contra_pu(i,j,p)%M )

                a11 = mesh%ll2covari_pu(i,j,p)%M(1,1)
                a12 = mesh%ll2covari_pu(i,j,p)%M(1,2)
                a21 = mesh%ll2covari_pu(i,j,p)%M(2,1)
                a22 = mesh%ll2covari_pu(i,j,p)%M(2,2)
                det = a11*a22 - a21*a12

                ! Covari to covari matrix
                mesh%covari2ll_pu(i,j,p)%M(1,1) =  a22
                mesh%covari2ll_pu(i,j,p)%M(1,2) = -a12
                mesh%covari2ll_pu(i,j,p)%M(2,1) = -a21
                mesh%covari2ll_pu(i,j,p)%M(2,2) =  a11
                mesh%covari2ll_pu(i,j,p)%M(:,:) = mesh%covari2ll_pu(i,j,p)%M(:,:)/det
            end do
        end do
    end do

    ! Compute at pv
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                ! latlon to covari matrix
                mesh%ll2covari_pv(i,j,p)%M = matmul(mesh%contra2covari_pv(i,j,p)%M, mesh%ll2contra_pv(i,j,p)%M )

                a11 = mesh%ll2covari_pv(i,j,p)%M(1,1)
                a12 = mesh%ll2covari_pv(i,j,p)%M(1,2)
                a21 = mesh%ll2covari_pv(i,j,p)%M(2,1)
                a22 = mesh%ll2covari_pv(i,j,p)%M(2,2)
                det = a11*a22 - a21*a12

                ! Covari to latlon matrix
                mesh%covari2ll_pv(i,j,p)%M(1,1) =  a22
                mesh%covari2ll_pv(i,j,p)%M(1,2) = -a12
                mesh%covari2ll_pv(i,j,p)%M(2,1) = -a21
                mesh%covari2ll_pv(i,j,p)%M(2,2) =  a11
                mesh%covari2ll_pv(i,j,p)%M(:,:) = mesh%covari2ll_pv(i,j,p)%M(:,:)/det
            end do
        end do
    end do

    ! Compute at pc
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                ! latlon to covari matrix
                mesh%ll2covari_pc(i,j,p)%M = matmul(mesh%contra2covari_pc(i,j,p)%M, mesh%ll2contra_pc(i,j,p)%M )

                a11 = mesh%ll2covari_pc(i,j,p)%M(1,1)
                a12 = mesh%ll2covari_pc(i,j,p)%M(1,2)
                a21 = mesh%ll2covari_pc(i,j,p)%M(2,1)
                a22 = mesh%ll2covari_pc(i,j,p)%M(2,2)
                det = a11*a22 - a21*a12

                ! Covari to latlon matrix
                mesh%covari2ll_pc(i,j,p)%M(1,1) =  a22
                mesh%covari2ll_pc(i,j,p)%M(1,2) = -a12
                mesh%covari2ll_pc(i,j,p)%M(2,1) = -a21
                mesh%covari2ll_pc(i,j,p)%M(2,2) =  a11
                mesh%covari2ll_pc(i,j,p)%M(:,:) = mesh%covari2ll_pc(i,j,p)%M(:,:)/det
            end do
        end do
    end do
end subroutine compute_conversion_matrices

end module cubed_sphere 
