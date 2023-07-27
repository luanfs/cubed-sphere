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
    getadvparameters, &
    getinterp_parameters

! Initial conditions
use advection_ic, only: &
    init_adv_vars, &
    compute_exact_div, &
    compute_ic_adv

! Advection timestep
use advection_timestep, only: &
    adv_timestep

! Linear algebra
use linear_algebra, only: &
    error_norm_max_rel, &
    error_norm_1_rel, &
    error_norm_2_rel

! Diagnostics 
use diagnostics, only: &
    mass_computation

! Duogrid interpolation
use duogrid_interpolation, only: &
    dg_interp

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
    type(scalar_field):: mt_pc

    ! aux vars
    real(r8) :: erad_km, erad2_km ! Earth radius in km
    integer(i4) :: p
    integer(i4) :: i   , j
    print*, 'Grid quality testing...'
    print*

    ! Earth radius in km
    erad_km  = erad/1e3_r8
    erad2_km = erad_km*erad_km

    ! Scalar fields allocation
    call scalar_field_allocation(area, mesh, 0)
    call scalar_field_allocation(length, mesh, 0)
    call scalar_field_allocation(mt_pc, mesh, 0)

    ! Get metric tensor
    mt_pc%f(i0:iend,j0:jend,:) = mesh%mt_pc(i0:iend,j0:jend,:)

    ! Get areas assuming earth radius
    area%f(i0:iend,j0:jend,:) = erad2_km*mesh%mt_pc(i0:iend,j0:jend,:)*mesh%dx*mesh%dy

    ! Get lenghts assuming earth radius
    length%f(i0:iend,j0:jend,:) = dsqrt(area%f(i0:iend,j0:jend,:))

    ! Name scalar fields
    area%name = "area"
    length%name = "length"
    mt_pc%name = "metric_tensor"

    ! Plot scalar fields
    call plot_scalarfield(area, mesh)
    call plot_scalarfield(length, mesh)
    call plot_scalarfield(mt_pc, mesh)

    deallocate(area%f, length%f, mt_pc%f)
end subroutine grid_quality 


subroutine interpolation_test(mesh)
    use advection_vars
    !---------------------------------------------------
    ! INTEPORLATION_TEST
    ! routine to test the lagrange interpolation
    ! at ghost cell positions
    !--------------------------------------------------
    type(cubedsphere),intent(inout):: mesh

    !File name for output
    character (len=256):: filename

    !Errors
    real(r8) :: error_q

    ! Get test parameter from par/interpolation.par
    call getinterp_parameters(advsimul)

    ! Initialize the variables (allocation, initial condition,...)
    call init_adv_vars(mesh)
    !advsimul%name = "div_"//trim(advsimul%name)

    ! Duogrid nterpolation of the scalar field Q 
    call dg_interp(Q, L_pc)

    ! Compute the error
    error_q = maxval(abs(Q_exact%f(:,:,:)-Q%f(:,:,:)))

    ! Print errors on screen
    print*
    print '(a22, 2e16.8)','(q, u) errors:', error_q, error_q
    stop
    ! Write errors in a file
    !filename = trim(advsimul%name)//"_"//trim(mesh%name)//"_errors"
    !call  write_final_errors(advsimul, mesh, filename) 

end subroutine interpolation_test


subroutine div_test(mesh)
    use advection_vars
    !---------------------------------------------------
    ! DIV_TEST
    !--------------------------------------------------
    type(cubedsphere),intent(inout):: mesh

    ! aux integer
    integer(i4) :: i, j, p

    !File name for output
    character (len=256):: filename

    ! Get test parameter from par/advection.par
    call getadvparameters(advsimul)

    ! For divergence testing, we make Q = 1
    advsimul%ic = 1

    ! Initialize the variables (allocation, initial condition,...)
    call init_adv_vars(mesh)

    advsimul%name = "div_"//trim(advsimul%name)

    ! Multiply Q by the metric tensor
    Q%f = 1._r8

    ! Compute the divergence obtained in one timestep
    call adv_timestep(mesh)

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
    advsimul%mass_variation = mass_computation(div_ugq, mesh)

    ! Deallocate vars
    call adv_deallocation()

    ! Print errors on screen
    print*
    print '(a22, 3e16.8)','linf, l1, l2 errors:', advsimul%linf_error, advsimul%l1_error, advsimul%l2_error
    print '(a22, 1e16.8)','div mass:', advsimul%mass_variation
    ! Write errors in a file
    filename = trim(advsimul%name)//"_"//trim(mesh%name)//"_errors"
    call write_final_errors(advsimul, mesh, filename) 

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
    type(scalar_field), intent(inout) :: Q      ! numerical approximation
    type(scalar_field), intent(inout) :: Q_ref   ! reference solution
    type(scalar_field), intent(inout) :: Q_error ! error
    real(r8), intent(out) :: linf, l1, l2
    ! aux vars
    integer (i4) :: x0, xend
    integer (i4) :: y0, yend

    select case(Q%pos)
        case(0) ! pc
            x0 = i0
            y0 = j0
            xend = iend
            yend = jend

        case(1) ! po
            x0 = i0+1
            y0 = j0+1
            xend = iend+1
            yend = jend+1

        case(2) ! pu
            x0 = i0+1
            y0 = j0
            xend = iend+1
            yend = jend

        case(3) ! pv
            x0 = i0
            y0 = j0+1
            xend = iend
            yend = jend+1

        case default
            print*, 'ERROR on compute_errors_field: invalid position', Q%pos

    end select

    ! Compute the errors
    Q_error%f = Q_ref%f - Q%f
    linf = error_norm_max_rel(Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 
    l1   = error_norm_1_rel  (Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 
    l2   = error_norm_2_rel  (Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 


end subroutine compute_errors_field

end module simulpack 
