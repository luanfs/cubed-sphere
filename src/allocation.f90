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
      n_lat, n_lon

  !Data structures
  use datastruct, only: &
      point_structure, &
      cubedsphere

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


  subroutine meshallocation(mesh)
    !---------------------------------------------------
    ! MESHALLOCATION
    ! allocate all the needed mesh atributtes
    !--------------------------------------------------
    type(cubedsphere):: mesh
    
    ! aux vars
    integer (i4):: status
    integer (i4):: n

    ! Number of cells along a coordinate axis (including ghost cells)
    n = mesh%ntotal

    ! Allocate the cubed-sphere points
    call point_3darray_allocation(mesh%pc, 1, n, 1, n, 1, nbfaces) 
    call point_3darray_allocation(mesh%pu, 0, n, 1, n, 1, nbfaces) 
    call point_3darray_allocation(mesh%pv, 1, n, 0, n, 1, nbfaces) 
    call point_3darray_allocation(mesh%po, 0, n, 0, n, 1, nbfaces)

    ! Latlon-grid allocation
    mesh%nlon = n_lon
    mesh%nlat = n_lat
    call point_2darray_allocation(mesh%ll, 0, mesh%nlon, 0, mesh%nlat) 
    call i4_2darray_allocation(mesh%ix_ll, 0, mesh%nlon, 0, mesh%nlat)
    call i4_2darray_allocation(mesh%jy_ll, 0, mesh%nlon, 0, mesh%nlat)
    call i4_2darray_allocation(mesh%panels_ll, 0, mesh%nlon, 0, mesh%nlat)


    if (mesh%kind=="equiangular")then
      ! Areas allocation
      call r8_3darray_allocation(mesh%area, 1, n, 1, n, 1, 1)

      ! Metric tensor allocation
      call r8_3darray_allocation(mesh%gc, 1, n, 1, n, 1, 1)

    else
      print*, 'ERROR on meshallocation: invalid mesh kind', mesh%kind
    end if
    return
  end subroutine meshallocation 

end module allocation 

