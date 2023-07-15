module allocation
!========================================================================
!
! Module for memory allocation routines
!
! Routines based on iModel (https://github.com/pedrospeixoto/iModel)
!
! Luan Santos 2022
!========================================================================

!Global constants
use constants, only: &
  i4, &
  r8, &
  nbfaces, &
  n_lat, n_lon, &
  n0, nend

!Data structures
use datastruct, only: &
  point_structure, &
  vector, &
  matrix, &
  cubedsphere, &
  scalar_field, &
  vector_field, &
  simulation, &
  ppm_parabola

implicit none

contains 

subroutine r8_1darray_allocation(data, i0, iend)
    !---------------------------------------------------
    ! r8_1DARRAY_ALLOCATION
    ! allocates a real r8 1d array
    !
    ! allocates data(i0:iend)
    !--------------------------------------------------
    real(r8), allocatable:: data(:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend), stat=status)
        if(status>0)then
            print *, "ERROR on r8_1darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on r8_1darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine r8_1darray_allocation 

subroutine r8_2darray_allocation(data, i0, iend, j0, jend)
    !---------------------------------------------------
    ! r8_2DARRAY_ALLOCATION
    ! allocates a real r8 2d array
    !
    ! allocates data(i0:iend, j0:jend)
    !--------------------------------------------------
    real(r8), allocatable:: data(:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend), stat=status)
        if(status>0)then
            print *, "ERROR on r8_2darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on r8_2darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine r8_2darray_allocation 

subroutine i4_2darray_allocation(data, i0, iend, j0, jend)
    !---------------------------------------------------
    ! i4_2DARRAY_ALLOCATION
    ! allocates an integer 2d array
    !
    ! allocates data(i0:iend, j0:jend)
    !--------------------------------------------------
    integer(i4), allocatable:: data(:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend), stat=status)
        if(status>0)then
            print *, "ERROR on i4_2darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on i4_2darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine i4_2darray_allocation 



subroutine r8_3darray_allocation(data, i0, iend, j0, jend, k0, kend)
    !---------------------------------------------------
    ! r8_3DARRAY_ALLOCATION
    ! allocates a real r8 3d array
    !
    ! allocates r8(i0:iend,j0:jend,k0:kend)
    !--------------------------------------------------
    real(r8), allocatable:: data(:,:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend
    integer(i4), intent(in):: k0, kend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend,k0:kend), stat=status)
        if(status>0)then
            print *, "ERROR on r8_3darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on r8_3darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine r8_3darray_allocation 


subroutine point_3darray_allocation(data, i0, iend, j0, jend, k0, kend)
    !---------------------------------------------------
    ! POINT_3DARRAY_ALLOCATION
    ! allocates a point structure 3d array
    !
    ! allocates data(i0:iend,j0:jend,k0:kend)
    !--------------------------------------------------
    type(point_structure), allocatable:: data(:,:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend
    integer(i4), intent(in):: k0, kend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend,k0:kend), stat=status)
        if(status>0)then
            print *, "ERROR on point_3darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on point_3darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine point_3darray_allocation 

subroutine point_2darray_allocation(data, i0, iend, j0, jend)
    !---------------------------------------------------
    ! POINT_2DARRAY_ALLOCATION
    ! allocates a point structure 2d array
    !
    ! allocates data(i0:iend,j0:jend)
    !--------------------------------------------------
    type(point_structure), allocatable:: data(:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend), stat=status)
        if(status>0)then
            print *, "ERROR on point_2darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on point_2darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine point_2darray_allocation 

subroutine vector_3darray_allocation(data, i0, iend, j0, jend, k0, kend)
    !---------------------------------------------------
    ! VECTOR_3DARRAY_ALLOCATION
    ! allocates a vector structure 3d array
    !
    ! allocates data(i0:iend,j0:jend,k0:kend)
    !--------------------------------------------------
    type(vector), allocatable:: data(:,:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend
    integer(i4), intent(in):: k0, kend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend,k0:kend), stat=status)
        if(status>0)then
            print *, "ERROR on vector_3darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on vector_3darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine vector_3darray_allocation 

subroutine matrix_3darray_allocation(data, i0, iend, j0, jend, k0, kend)
    !---------------------------------------------------
    ! MATRIX_3DARRAY_ALLOCATION
    ! allocates a matrix structure 3d array
    !
    ! allocates data(i0:iend,j0:jend,k0:kend)
    !--------------------------------------------------
    type(matrix), allocatable:: data(:,:,:) 

    ! Dimension
    integer(i4), intent(in):: i0, iend
    integer(i4), intent(in):: j0, jend
    integer(i4), intent(in):: k0, kend

    ! aux vars
    integer (i4):: status

    if(.not.allocated(data))then
        allocate(data(i0:iend,j0:jend,k0:kend), stat=status)
        if(status>0)then
            print *, "ERROR on matrix_3darray_allocation: Allocation problem ", status
            stop
        end if
    else
        print *, "ERROR on matrix_3darray_allocation: Trying to allocate an already allocated variable "
        stop
    end if
    return
end subroutine matrix_3darray_allocation 



subroutine meshallocation(mesh)
    !---------------------------------------------------
    ! MESHALLOCATION
    ! allocate all the needed mesh atributtes
    !--------------------------------------------------
    type(cubedsphere):: mesh

    ! aux vars
    integer (i4):: n0, nf

    ! Number of cells along a coordinate axis (including ghost cells)
    n0 = mesh%n0
    nf = mesh%nend

    ! Allocate the cubed-sphere points
    call point_3darray_allocation(mesh%pc, n0, nf  , n0, nf  , 1, nbfaces) 
    call point_3darray_allocation(mesh%pu, n0, nf+1, n0, nf  , 1, nbfaces) 
    call point_3darray_allocation(mesh%pv, n0, nf  , n0, nf+1, 1, nbfaces) 
    call point_3darray_allocation(mesh%po, n0, nf+1, n0, nf+1, 1, nbfaces)

    ! Allocate the cubed-sphere tangent vectors
    call vector_3darray_allocation(mesh%tgx_pc, n0, nf  , n0, nf  , 1, nbfaces) 
    call vector_3darray_allocation(mesh%tgy_pc, n0, nf  , n0, nf  , 1, nbfaces) 
    call vector_3darray_allocation(mesh%tgx_pu, n0, nf+1, n0, nf  , 1, nbfaces) 
    call vector_3darray_allocation(mesh%tgy_pu, n0, nf+1, n0, nf  , 1, nbfaces) 
    call vector_3darray_allocation(mesh%tgx_pv, n0, nf  , n0, nf+1, 1, nbfaces) 
    call vector_3darray_allocation(mesh%tgy_pv, n0, nf  , n0, nf+1, 1, nbfaces) 

    ! Metric tensor
    call r8_3darray_allocation(mesh%mt_pc, n0, nf  , n0, nf  , 1, nbfaces)
    call r8_3darray_allocation(mesh%mt_pu, n0, nf+1, n0, nf  , 1, nbfaces)
    call r8_3darray_allocation(mesh%mt_pv, n0, nf  , n0, nf+1, 1, nbfaces)


    ! Allocate the cubed-sphere latlon/contravariant conversion matrix
    call matrix_3darray_allocation(mesh%contra2ll_pu, n0, nf+1, n0, nf  , 1, nbfaces) 
    call matrix_3darray_allocation(mesh%ll2contra_pu, n0, nf+1, n0, nf  , 1, nbfaces) 
    call matrix_3darray_allocation(mesh%contra2ll_pv, n0, nf  , n0, nf+1, 1, nbfaces) 
    call matrix_3darray_allocation(mesh%ll2contra_pv, n0, nf  , n0, nf+1, 1, nbfaces) 

    ! Latlon-grid allocation
    mesh%nlon = n_lon
    mesh%nlat = n_lat
    call point_2darray_allocation(mesh%ll, 0, mesh%nlon, 0, mesh%nlat) 
    call i4_2darray_allocation(mesh%ix_ll, 0, mesh%nlon, 0, mesh%nlat)
    call i4_2darray_allocation(mesh%jy_ll, 0, mesh%nlon, 0, mesh%nlat)
    call i4_2darray_allocation(mesh%panels_ll, 0, mesh%nlon, 0, mesh%nlat)

    return
end subroutine meshallocation 

subroutine scalar_field_allocation(Q, mesh, pos)
    !---------------------------------------------------
    ! SCALAR_FIELD_ALLOCATION
    !
    ! This routines allocates the scalar field Q
    ! at position 'pos', including ghost cells.
    !
    ! mesh must be already allocated.
    !
    ! Position of the values relative to a mesh
    !   0 - Centers (pc)
    !   1 - Vertices (po) 
    !   2 - Midpoint at u position (pu)
    !   3 - Midpoint at v position (pv)
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in):: mesh
    type(scalar_field), intent(inout) :: Q
   
    ! Position
    integer, intent(in) :: pos

    ! aux vars
    integer (i4):: n0, nend

    ! Number of cells along a coordinate axis (including ghost cells)
    n0 = mesh%n0
    nend = mesh%nend

    ! Allocate the cubed-sphere points
    select case(pos)
        case(0) ! centers
            call r8_3darray_allocation(Q%f, n0, nend, n0, nend, 1, nbfaces) 

        case(1) ! vertices
            call r8_3darray_allocation(Q%f, n0, nend+1, n0, nend+1, 1, nbfaces) 

        case(2) ! mid u
            call r8_3darray_allocation(Q%f, n0, nend+1, n0, nend, 1, nbfaces) 

        case(3) ! mid v
            call r8_3darray_allocation(Q%f, n0, nend, n0, nend+1, 1, nbfaces) 

        case default
            print*, 'ERROR on scalar_field_allocation: invalid position', pos
    end select
      
    ! Store the position
    Q%pos = pos
    return
end subroutine scalar_field_allocation 

subroutine vector_field_allocation(V, mesh, pos)
    !---------------------------------------------------
    ! VECTOR_FIELD_ALLOCATION
    !
    ! This routines allocates the vector field V
    ! at position 'pos', including ghost cells.

    ! mesh must be already allocated.
    !
    ! Position of the values relative to a mesh
    !   0 - Centers (pc)
    !   1 - Vertices (po) 
    !   2 - Midpoint at u position (pu)
    !   3 - Midpoint at v position (pv)
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in):: mesh
    type(vector_field), intent(inout) :: V
   
    ! Position
    integer, intent(in) :: pos
     
    ! Allocate all the vector representations
    call scalar_field_allocation(V%u, mesh, pos) 
    call scalar_field_allocation(V%v, mesh, pos) 
    call scalar_field_allocation(V%ucontra, mesh, pos) 
    call scalar_field_allocation(V%vcontra, mesh, pos) 
    call scalar_field_allocation(V%ucontra_av, mesh, pos) 
    call scalar_field_allocation(V%vcontra_av, mesh, pos) 
    call scalar_field_allocation(V%ucovari, mesh, pos) 
    call scalar_field_allocation(V%vcovari, mesh, pos) 
     
    ! Store the position
    V%pos = pos
    return
end subroutine vector_field_allocation 


subroutine ppm_parabola_allocation(p, mesh)
    !---------------------------------------------------
    !   PPM_PARABOLA_ALLOCATION
    !
    ! This routines allocates the scalar field Q
    ! at position 'pos', including ghost cells.
    !
    ! mesh must be already allocated.
    !
    ! Direction of the values relative to a mesh
    !   1 - x direction
    !   2 - y direction 
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in):: mesh
    type(ppm_parabola), intent(inout) :: p
   

    call r8_3darray_allocation(p%q_L, n0, nend, n0, nend, 1, nbfaces) 
    call r8_3darray_allocation(p%q_R, n0, nend, n0, nend, 1, nbfaces) 
    call r8_3darray_allocation(p%dq,  n0, nend, n0, nend, 1, nbfaces) 
    call r8_3darray_allocation(p%q6,  n0, nend, n0, nend, 1, nbfaces) 
    call r8_3darray_allocation(p%df,  n0, nend, n0, nend, 1, nbfaces) 
    call scalar_field_allocation(p%Q, mesh, 0)

    !---------------------------------------------------
    ! Allocate the cubed-sphere points
    select case(p%dir)
        case(1) ! x direction
            call r8_3darray_allocation(p%f_L,   n0, nend+1, n0, nend, 1, nbfaces) 
            call r8_3darray_allocation(p%f_R,   n0, nend+1, n0, nend, 1, nbfaces) 
            call r8_3darray_allocation(p%f_upw, n0, nend+1, n0, nend, 1, nbfaces) 

      case(2) ! y direction
            call r8_3darray_allocation(p%f_L,   n0, nend, n0, nend+1, 1, nbfaces) 
            call r8_3darray_allocation(p%f_R,   n0, nend, n0, nend+1, 1, nbfaces) 
            call r8_3darray_allocation(p%f_upw, n0, nend, n0, nend+1, 1, nbfaces) 

      case default
            print*, 'ERROR on ppm_parabola_allocation: invalid direction', p%dir

      end select
      return
end subroutine ppm_parabola_allocation 


subroutine allocate_adv_vars(mesh)
    use advection_vars
    !---------------------------------------------------
    ! ALLOCATE_ADV_VARS 
    !
    ! This routine allocated the variables of 
    ! the advection simulation testcase
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh

    ! Allocate variables
    call scalar_field_allocation(Q, mesh, 0)
    call scalar_field_allocation(Q_exact, mesh, 0)
    call vector_field_allocation(wind_pu, mesh, 2)
    call vector_field_allocation(wind_pv, mesh, 3)
    call scalar_field_allocation(cx_pu, mesh, 2)
    call scalar_field_allocation(cy_pv, mesh, 3)
    call scalar_field_allocation(div_ugq, mesh, 0)
    call scalar_field_allocation(div_ugq_exact, mesh, 0)
    call scalar_field_allocation(div_ugq_error, mesh, 0)

    !Operator splitting variables
    call scalar_field_allocation(Qx, mesh, 0)
    call scalar_field_allocation(Qy, mesh, 0)

    !PPM vars
    call ppm_parabola_allocation(px, mesh)
    call ppm_parabola_allocation(py, mesh)

end subroutine
end module allocation 

