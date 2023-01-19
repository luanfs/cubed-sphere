module simulpack
  !========================================================================
  !
  ! Module for simulation package
  ! Includes: grid quality, divergence and advection simulations
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
      pi

  !Data structures
  use datastruct, only: &
      cubedsphere, &
      scalar_field, &
      vector_field, &
      simulation

  ! Allocation
  use allocation, only: &
      scalar_field_allocation

  ! Deallocation
  use deallocation, only: &
      adv_deallocation

  ! Output
  use output, only: &
      plot_scalarfield, &
      write_final_errors

  ! Input
  use input, only: &
      getadvparameters

  ! Initial conditions
  use advection_ic, only: &
      init_adv_vars, &
      compute_exact_div

  ! Discrete operators
  use discrete_operators, only: &
      cfl_x, cfl_y, &
      divergence

  ! Linear algebra
  use linear_algebra, only: &
      error_norm_max_rel, &
      error_norm_1_rel, &
      error_norm_2_rel


  implicit none

  contains 

  subroutine grid_quality(mesh)
    !---------------------------------------------------
    ! GRID_QUALITY
    !--------------------------------------------------
    type(cubedsphere), intent(in):: mesh

    ! Scalar fields
    type(scalar_field):: area
    type(scalar_field):: length
    type(scalar_field):: sinc

    ! aux vars
    real(r8) :: erad_km, erad2_km ! Earth radius in km
    integer(i4) :: p
    integer(i4) :: i   , j
    integer(i4) :: i0  , j0
    integer(i4) :: iend, jend

    print*, 'Grid quality testing...'
    print*

    ! Earth radius in km
    erad_km  = erad/1e3_r8
    erad2_km = erad_km*erad_km

    ! Scalar fields allocation
    call scalar_field_allocation(area, mesh, 0)
    call scalar_field_allocation(length, mesh, 0)
    call scalar_field_allocation(sinc, mesh, 0)

    ! Interior indexes
    i0   = mesh%i0
    iend = mesh%iend
    j0   = mesh%j0
    jend = mesh%jend
  
    if(mesh%resolution == 'unif') then ! use information from a single panel
      ! Get areas assuming earth radius
      area%f(i0:iend,j0:jend,1) = erad2_km*mesh%area(i0:iend,j0:jend,1)
      
      do p = 2, nbfaces
        area%f(i0:iend,j0:jend,p) =  area%f(i0:iend, j0:jend,1) 
      end do

      ! Get lenghts assuming earth radius
      length%f(i0:iend,j0:jend,1) = erad_km*(mesh%lx(i0:iend,j0:jend,1) + mesh%ly(i0:iend,j0:jend,1))*0.5_r8

      do p = 2, nbfaces
        length%f(i0:iend,j0:jend,p) =  length%f(i0:iend, j0:jend,1) 
      end do

      ! Get metric tensor (unit sphere)
      sinc%f(i0:iend,j0:jend,1) = mesh%sinc(i0:iend,j0:jend,1) 

      do p = 2, nbfaces
        sinc%f(i0:iend,j0:jend,p) =  sinc%f(i0:iend, j0:jend,1) 
      end do


    else 
      print*, 'ERROR on grid_quality: invalid mesh resolution: ', mesh%resolution
      stop
    end if

    ! Name scalar fields
    area%name = "area"
    length%name = "length"
    sinc%name = "sinc"

    ! Plot scalar fields
    call plot_scalarfield(area, mesh)
    call plot_scalarfield(length, mesh)
    call plot_scalarfield(sinc, mesh)

    deallocate(area%f, length%f, sinc%f)
  end subroutine grid_quality 

  subroutine div_test(mesh)
    use advection_vars
    !---------------------------------------------------
    ! DIV_TEST
    !--------------------------------------------------
    type(cubedsphere),intent(inout):: mesh

    ! aux integer
    integer(i4) :: i, j, p
    integer(i4) :: i0, iend
    integer(i4) :: j0, jend

    !File name for output
    character (len=256):: filename

    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    ! Get test parameter from par/advection.par
    call getadvparameters(advsimul)

    ! For divergence testing, we make Q = 1
    advsimul%ic = 1

    ! Initialize the variables (allocation, initial condition,...)
    call init_adv_vars(mesh)

    ! Test name
    advsimul%name = "div_"//"vf"//trim(advsimul%vf_name)//"_"//trim(advsimul%recon1d)
    !print*, advsimul%name

    ! CFL number
    call cfl_x(mesh, wind_pu, cx_pu, advsimul%dt)
    call cfl_y(mesh, wind_pv, cy_pv, advsimul%dt)

    ! Multiply Q by the metric tensor
    Q%f = 1._r8

    ! Compute the divergence
    call divergence(Q, wind_pu, wind_pv, mesh)

    ! Exact divergence
    call compute_exact_div(div_ugq_exact, mesh, advsimul)

    ! Compute the errors
    call compute_errors_field(div_ugq, div_ugq_exact, div_ugq_error, &
      advsimul%linf_error, advsimul%l1_error, advsimul%l2_error, mesh)

    ! Name scalar fields
    div_ugq_error%name = trim(advsimul%name)//"_error"
    div_ugq%name = trim(advsimul%name)

    ! Plot scalar fields
    call plot_scalarfield(div_ugq, mesh)
    call plot_scalarfield(div_ugq_error, mesh)

    ! Deallocate vars
    call adv_deallocation()

    ! Print errors on screen
    print*
    print '(a22, 3e16.8)','linf, l1, l2 errors:', advsimul%linf_error, advsimul%l1_error, advsimul%l2_error

    ! Write errors in a file
    filename = trim(advsimul%name)//"_"//trim(mesh%name)//"_errors"
    call  write_final_errors(advsimul, mesh, filename) 

  end subroutine div_test


  subroutine compute_errors_field(Q, Q_ref, Q_error, linf, l1, l2, mesh)
    !---------------------------------------------------
    ! COMPUTE_ERRORS_FIELD
    !
    ! Given the scalar field reference field Q_ref,
    ! this routine compute the l1, l2, linf errors
    ! of Q
    !--------------------------------------------------
    type(cubedsphere),intent(in) :: mesh
    type(scalar_field), intent(in) :: Q      ! numerical approximation
    type(scalar_field), intent(in) :: Q_ref   ! reference solution
    type(scalar_field), intent(out) :: Q_error ! error
    real(r8), intent(out) :: linf, l1, l2
    ! aux vars
    integer (i4) :: i0, iend
    integer (i4) :: j0, jend

    select case(Q%pos)
      case(0) ! centers
        i0 = mesh%i0
        j0 = mesh%j0
        iend = mesh%iend
        jend = mesh%jend

      case(1) ! vertices
        i0 = mesh%i0-1
        j0 = mesh%j0-1
        iend = mesh%iend
        jend = mesh%jend

      case(2) ! mid u
        i0 = mesh%i0-1
        j0 = mesh%j0
        iend = mesh%iend
        jend = mesh%jend
 
      case(3) ! mid v
        i0 = mesh%i0
        j0 = mesh%j0-1
        iend = mesh%iend
        jend = mesh%jend
 
      case default
        print*, 'ERROR on compute_errors_field: invalid position', Q%pos

    end select
     
    ! Compute the errors
    Q_error%f = Q_ref%f - Q%f
    linf = error_norm_max_rel(Q%f(i0:iend,j0:jend,:), Q_ref%f(i0:iend,j0:jend,:)) 
    l1   = error_norm_1_rel  (Q%f(i0:iend,j0:jend,:), Q_ref%f(i0:iend,j0:jend,:)) 
    l2   = error_norm_2_rel  (Q%f(i0:iend,j0:jend,:), Q_ref%f(i0:iend,j0:jend,:)) 

 
  end subroutine compute_errors_field

end module simulpack 
