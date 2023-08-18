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
    adv_deallocation, &
    swm_deallocation

! Output
use output, only: &
    plot_scalarfield, &
    write_final_errors_adv, &
    write_final_errors_swm, &
    write_final_errors_interp, &
    output_adv, &
    output_swm, &
    compute_errors_field, &
    print_advparameters, &
    print_swmparameters

! Input
use input, only: &
    getadvparameters, &
    getswmparameters, &
    getinterp_parameters

! Initial conditions
use advection_ic, only: &
    init_adv_vars, &
    compute_exact_div, &
    compute_ic_adv

use swm_ic, only: &
    init_swm_vars
 
! Advection timestep
use advection_timestep, only: &
    adv_timestep, &
    adv_update


! Shallow water timestep
use swm_timestep, only: &
    sw_timestep, &
    sw_timestep_Cgrid

! Diagnostics 
use diagnostics, only: &
    mass_computation

! Duogrid interpolation
use duogrid_interpolation, only: &
    dg_interp, &
    interp_C2Aduogrid, &
    interp_C2Agrid, &
    interp_A2Cduogrid, &
    interp_A2Cgrid, &
    interp_D2Aduogrid

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
    real(kind=8) :: erad_km, erad2_km ! Earth radius in km
    integer(i4) :: p
    integer(i4) :: i   , j
    print*, 'Grid quality testing...'
    print*

    ! Earth radius in km
    erad_km  = erad/1000.d0
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
    real(kind=8) :: error_q
    real(kind=8) :: error_ucontra, error_vcontra
    real(kind=8) :: error_ucovari, error_vcovari

    ! Get test parameter from par/interpolation.par
    call getadvparameters(advsimul)
    call getinterp_parameters(advsimul)

    ! Initialize the variables (allocation, initial condition,...)
    call init_adv_vars(mesh)
    advsimul%avd = 3
    !advsimul%name = "div_"//trim(advsimul%name)

    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the scalar field Q 
    call dg_interp(Q%f, L_pc)

    ! Compute the error
    error_q = maxval(abs(Q_exact%f(:,:,:)-Q%f(:,:,:)))


    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the vector field on a C grid to its ghost cell values
    ! first we interpolate from C grid to the A grid ghost cells
    call interp_C2Aduogrid(wind_pu%ucontra%f, wind_pv%vcontra%f,&
    wind_pc%u%f, wind_pc%v%f,wind_pc%ucontra%f, wind_pc%vcontra%f,&
    L_pc, mesh%contra2ll_pc, mesh%ll2contra_pc)

    ! then we interpolate from C grid to the A grid inner cells
    call interp_C2Agrid(wind_pu%ucontra%f, wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f,  advsimul%avd)

    ! now we fill the ghost cell C grid
    call interp_A2Cduogrid(wind_pu%ucontra%f, wind_pu%vcontra%f, wind_pv%ucontra%f, &
    wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f)
 
    ! then we interpolate from A grid to the C grid inner cells
    call interp_A2Cgrid(wind_pu%ucontra%f, wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f, advsimul%avd)

    error_ucontra = maxval(abs(wind_pu%ucontra%f(i0-1:iend+2,n0:nend,:)-wind_pu%ucontra_old%f(i0-1:iend+2,n0:nend,:)))
    error_vcontra = maxval(abs(wind_pv%vcontra%f(n0:nend,j0-1:jend+2,:)-wind_pv%vcontra_old%f(n0:nend,j0-1:jend+2,:)))
    error_ucontra = max(error_ucontra, error_vcontra)

    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the vector field on a D grid
    ! first we interpolate to the A grid ghost cells)
    !call interp_D2Aduogrid(wind_pu%vcovari%f, wind_pv%ucovari%f, &
    !wind_pc%v%f, wind_pc%u%f, wind_pc%vcovari%f, wind_pc%ucovari%f, &
    !L_pc, mesh%covari2ll_pc, mesh%ll2covari_pc)
    call interp_D2Aduogrid(wind_pu, wind_pv, wind_pc, L_pc, mesh)

    ! then we interpolate from d grid to the a grid inner cells
    call interp_c2agrid(wind_pu%vcovari%f, wind_pv%ucovari%f, wind_pc%vcovari%f, wind_pc%ucovari%f,  advsimul%avd)

    ! now we fill the ghost cell D grid
    call interp_A2Cduogrid(wind_pu%ucovari%f, wind_pu%vcovari%f, wind_pv%ucovari%f, &
    wind_pv%vcovari%f, wind_pc%ucovari%f, wind_pc%vcovari%f)
 
    ! then we interpolate from A grid to the D grid inner cell
    call interp_A2Cgrid(wind_pu%ucovari%f, wind_pv%vcovari%f, wind_pc%ucovari%f, wind_pc%vcovari%f, advsimul%avd)

    error_ucovari = maxval(abs(wind_pu%ucovari%f(i0-1:iend+2,n0:nend,:)-wind_pu%ucovari_old%f(i0-1:iend+2,n0:nend,:)))
    error_vcovari = maxval(abs(wind_pv%vcovari%f(n0:nend,j0-1:jend+2,:)-wind_pv%vcovari_old%f(n0:nend,j0-1:jend+2,:)))
    error_ucovari = max(error_ucovari, error_vcovari)

    print*
    print '(a22, 3e16.8)','(q, u, v) errors:', error_q, error_ucontra, error_ucovari

    ! Write errors in a file
    filename = "interp_ic"//trim(advsimul%ic_name)//"_vf"//trim(advsimul%vf_name)&
    //"_id"//trim(advsimul%id_name)//"_"//trim(mesh%name)//"_errors"
    call  write_final_errors_interp(filename, error_q, error_ucontra, error_ucovari)

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

    ! Initialize the variables (allocation, initial condition,...)
    advsimul%ic = 1
    advsimul%n = 1
    call init_adv_vars(mesh)

    ! print parameters
    call print_advparameters(advsimul)

    ! Compute the divergence obtained in one timestep
    call adv_timestep(mesh)

    ! Exact divergence
    call compute_exact_div(div_ugq_exact, mesh, advsimul)

    ! Compute the errors
    call compute_errors_field(div_ugq, div_ugq_exact, div_ugq_error, &
      advsimul%linf_error_h, advsimul%l1_error_h, advsimul%l2_error_h, mesh)

    ! Plot scalar fields
    call plot_scalarfield(div_ugq, mesh)
    call plot_scalarfield(div_ugq_error, mesh)

    advsimul%mass_variation = mass_computation(div_ugq, mesh)

    ! Print errors on screen
    print*
    print '(a22, 3e16.8)','linf, l1, l2 errors:', advsimul%linf_error_h, advsimul%l1_error_h, advsimul%l2_error_h
    print '(a22, 1e16.8)','div mass:', advsimul%mass_variation

    ! Write errors in a file
    filename = "div_"//trim(advsimul%name)//"_"//trim(mesh%name)//"_errors"
    call write_final_errors_adv(advsimul, mesh, filename) 

    ! Deallocate vars
    call adv_deallocation()

end subroutine div_test


subroutine adv_test(mesh)
    use advection_vars
    !---------------------------------------------------
    ! ADV_TEST
    !--------------------------------------------------
    type(cubedsphere),intent(inout):: mesh

    ! aux integer
    integer(i4) :: n

    !File name for output
    character (len=256):: filename

    ! Get test parameter from par/advection.par
    call getadvparameters(advsimul)

    ! Initialize the variables (allocation, initial condition,...)
    call init_adv_vars(mesh)

    ! print parameters
    call print_advparameters(advsimul)

    ! Initial output
    call output_adv(mesh)

    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the vector field on a C grid to its ghost cell values
    ! first we interpolate from C grid to the A grid ghost cells
    call interp_C2Aduogrid(wind_pu%ucontra%f, wind_pv%vcontra%f,&
    wind_pc%u%f, wind_pc%v%f,wind_pc%ucontra%f, wind_pc%vcontra%f,&
    L_pc, mesh%contra2ll_pc, mesh%ll2contra_pc)

    ! now we fill the ghost cell C grid
    call interp_A2Cduogrid(wind_pu%ucontra%f, wind_pu%vcontra%f, wind_pv%ucontra%f, &
    wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f)
    wind_pu%ucontra_old%f(:,:,:) = wind_pu%ucontra%f(:,:,:)
    wind_pv%vcontra_old%f(:,:,:) = wind_pv%vcontra%f(:,:,:)

    ! Temporal loop
    advsimul%t = 0.d0
    do n = 1, advsimul%nsteps
        ! Update time
        advsimul%t = advsimul%t + advsimul%dt
        advsimul%n = n

        ! Update the solution
        call adv_timestep(mesh)

        ! Output
        call output_adv(mesh)

        ! Update the velocity field for the next time step - only for variable velocity
        if(advsimul%vf>=2)then
            call adv_update(wind_pu, wind_pv, mesh, advsimul%vf, advsimul%t)
        end if
    end do

    ! Write errors in a file
    filename = "adv_"//trim(advsimul%name)//"_"//trim(mesh%name)//"_errors"
    call write_final_errors_adv(advsimul, mesh, filename) 

    ! Deallocate vars
    call adv_deallocation()

end subroutine adv_test

subroutine swm_test(mesh)
    use swm_vars
    !---------------------------------------------------
    ! SWM_TEST
    !--------------------------------------------------
    type(cubedsphere),intent(inout):: mesh

    ! aux integer
    integer(i4) :: n

    !File name for output
    character (len=256):: filename

    ! Get test parameter from par/swm.par
    call getswmparameters(swm_simul)

    ! Initialize the variables (allocation, initial condition,...)
    call init_swm_vars(mesh)

    ! Initialize velocity at ghost cells
    call sw_timestep_Cgrid(mesh)
    !print*, maxval(abs(wind_pc%ucovari%f(:,:,:)-wind_pc%ucovari_old%f(:,:,:)))
    !print*, maxval(abs(wind_pc%vcovari%f(:,:,:)-wind_pc%vcovari_old%f(:,:,:)))
    !print*, maxval(abs(wind_pc%ucontra%f(:,:,:)-wind_pc%ucontra_old%f(:,:,:)))
    !print*, maxval(abs(wind_pc%vcontra%f(:,:,:)-wind_pc%vcontra_old%f(:,:,:)))
    !print*, maxval(abs(wind_pu%ucontra%f(i0-1:iend+2,:,:)-wind_pu%ucontra_old%f(i0-1:iend+2,:,:)))
    !print*, maxval(abs(wind_pv%vcontra%f(:,j0-1:jend+2,:)-wind_pv%vcontra_old%f(:,j0-1:jend+2,:)))

    ! CFL number
    swm_simul%cfl = maxval(abs(wind_pu%ucontra%f))
    swm_simul%cfl = max(swm_simul%cfl, maxval(abs(wind_pv%vcontra%f)))
    swm_simul%cfl = swm_simul%cfl*swm_simul%dt/mesh%dx
 
    ! print parameters
    call print_swmparameters(swm_simul)

    ! Initial output
    call output_swm(mesh)

    wind_pu%ucontra_old%f(:,:,:) = wind_pu%ucontra%f(:,:,:)
    wind_pv%vcontra_old%f(:,:,:) = wind_pv%vcontra%f(:,:,:)

    ! Temporal loop
    swm_simul%t = 0.d0
    do n = 1, swm_simul%nsteps
        ! Update time
        swm_simul%t = swm_simul%t + swm_simul%dt
        swm_simul%n = n

        ! Update the solution
        call sw_timestep(mesh)

        ! Output
        call output_swm(mesh)
    end do

    ! Write errors in a file
    filename = "swm_"//trim(swm_simul%name)//"_"//trim(mesh%name)//"_errors"
    call write_final_errors_swm(swm_simul, mesh, filename) 

    ! Deallocate vars
    call swm_deallocation()

end subroutine swm_test


end module simulpack 
