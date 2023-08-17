module swm_timestep
!===============================================================================================
! Module for routines needed in a shallow water model timestep
!
! Luan da Fonseca Santos - 2023
! (luan.santos@usp.br)
!===============================================================================================

!Global constants
use constants, only: &
    i4, &
    pi, &
    nbfaces, &
    i0, iend, &
    j0, jend, &
    n0, nend

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  deg2rad

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field, &
  simulation

! Discrete operators 
use discrete_operators, only: &
    divergence, &
    vorticity_fluxes, &
    cfl_x, cfl_y

! duo grid interpolation
use duogrid_interpolation, only: &
    interp_D2Aduogrid, &
    interp_A2Cduogrid, &
    interp_C2Agrid, &
    interp_A2Cgrid, &
    dg_interp

! Model variables
use swm_vars

! Departure point
use departure_point, only: &
    swm_time_averaged_wind

implicit none

contains 

subroutine sw_timestep(mesh)
    !--------------------------------------------------
    ! Compute one time step for the shallow water
    ! problem on the sphere
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    call sw_timestep_Cgrid(mesh)
    call sw_timestep_Dgrid(mesh)

end subroutine sw_timestep


subroutine sw_timestep_Dgrid(mesh)
    !--------------------------------------------------
    ! Compute one time step for the shallow water
    ! problem on the sphere on the D grid
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Discrete divergence
    dginterp = .false. 
    call divergence(div_ugH, H, wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, swm_simul, mesh, L_pc, dginterp)

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(H, swm_simul, div_ugH)
    ! Update the fluid depth
    H%f = H%f - swm_simul%dt*div_ugH%f
    !$OMP END PARALLEL WORKSHARE

    ! interpolate depth to ghost cells
    if(swm_simul%et=='duogrid') then
        call dg_interp(H%f, L_pc)
    end if

    ! compute the vorticity fluxes
    if(swm_simul%n==-1)then
        !----------------------------------------------------------------------------------
        ! compute relative vorticity using the winds from D grid !formula 19 from LR1997
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(rel_vort, wind_pu, wind_pv, mesh, i0, iend, j0, jend)
        rel_vort%f(i0:iend,j0:jend,:) = &
        (wind_pu%vcovari%f(i0+1:iend+1,j0:jend,:) - wind_pu%vcovari%f(i0:iend,j0:jend,:))/mesh%dx - &
        (wind_pv%ucovari%f(i0:iend,j0+1:jend+1,:) - wind_pv%ucovari%f(i0:iend,j0:jend,:))/mesh%dy
        rel_vort%f(i0:iend,j0:jend,:) = rel_vort%f(i0:iend,j0:jend,:)/mesh%mt_pc(i0:iend,j0:jend,:)
        !$OMP END PARALLEL WORKSHARE
        !----------------------------------------------------------------------------------

        ! interpolate relative vorticity to ghost cells
        call dg_interp(rel_vort%f, L_pc)

        !----------------------------------------------------------------------------------
        ! absolute vorticity
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(abs_vort, rel_vort, fcoriolis_pc, i0, iend, j0, jend)
        abs_vort%f(:,:,:) = rel_vort%f(:,:,:) + fcoriolis_pc%f(:,:,:)
        !$OMP END PARALLEL WORKSHARE
        !----------------------------------------------------------------------------------

    end if

    call vorticity_fluxes(div_abs_vort, abs_vort_flux_pu, abs_vort_flux_pv, &
                          rel_vort, abs_vort, fcoriolis_pc,&
                          wind_pu, wind_pv, cx_pu, cy_pv, &
                          px, py, Qx, Qy, swm_simul, mesh, L_pc)

    !print*, maxval(abs(div_abs_vort%f(i0:iend,j0:jend,:)))

    !div_abs_vort%f(i0:iend,j0:jend,:) = &
    !((abs_vort_flux_pu%f(i0+1:iend+1,j0:jend,:)  - &
    !  abs_vort_flux_pu%f(i0:iend    ,j0:jend,:)) +  &
    !( abs_vort_flux_pv%f(i0:iend,j0+1:jend+1,:) - &
    !  abs_vort_flux_pv%f(i0:iend,j0:jend,:) ))/mesh%mt_pc(i0:iend,j0:jend,:)/mesh%dx
    !print*, maxval(abs(div_abs_vort%f(i0:iend,j0:jend,:)))

    !print*, maxval(abs(abs_vort_flux_pu%f(i0+1:iend+1,j0:jend,:)))
    !print*, maxval(abs(abs_vort_flux_exact_pu%f(i0:iend+1,j0:jend,:)))
    !print*, maxval(abs(abs_vort_flux_pu%f(i0:iend+1,j0:jend,:)))
    !print*, maxval(abs(abs_vort_flux_exact_pu%f(i0:iend+1,j0:jend,:) - &
    !                   abs_vort_flux_pu%f(i0:iend+1,j0:jend,:)))/maxval(abs(abs_vort_flux_exact_pu%f(i0:iend+1,j0:jend,:)))
    !print*
    !div_abs_vort%f(i0:iend,j0:jend,:) = &
    !((abs_vort_flux_exact_pu%f(i0+1:iend+1,j0:jend,:)  - &
    !  abs_vort_flux_exact_pu%f(i0:iend    ,j0:jend,:)) +  &
    !( abs_vort_flux_exact_pv%f(i0:iend,j0+1:jend+1,:) - &
    !  abs_vort_flux_exact_pv%f(i0:iend,j0:jend,:) ))/mesh%mt_pc(i0:iend,j0:jend,:)/mesh%dx

    !print*, maxval(abs( &
    !abs_vort_flux_exact_pu%f(i0+1:iend+1,j0:jend,:) - &
    !abs_vort_flux_pu%f(i0+1:iend+1,j0:jend,:) ))/maxval(abs(abs_vort_flux_exact_pu%f(i0+1:iend+1,j0:jend,:)))

    !print*, maxval(abs(abs_vort_flux_exact_pu%f(i0+1:iend+1,j0:jend,:)))
    !print*, maxval(abs(abs_vort_flux_pu%f(i0+1:iend+1,j0:jend,:)))
    !print*, maxval(abs(div_abs_vort%f(i0:iend,j0:jend,:)))
    !stop
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(abs_vort, swm_simul, div_abs_vort)
    ! Update the absolute vorticity
    !abs_vort%f = abs_vort%f - swm_simul%dt*div_abs_vort%f
    !abs_vort%f = div_abs_vort%f
    !$OMP END PARALLEL WORKSHARE


end subroutine sw_timestep_Dgrid



subroutine sw_timestep_Cgrid(mesh)
    !--------------------------------------------------
    ! Compute one time step for the shallow water
    ! problem on the sphere on the C grid
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the vector field on a D grid
    ! first we interpolate to the A grid ghost cells
    call interp_D2Aduogrid(wind_pu, wind_pv, wind_pc, L_pc, mesh)

    ! then we interpolate from D grid to the a grid inner cells
    call interp_C2Agrid(wind_pu%vcovari%f, wind_pv%ucovari%f, wind_pc%vcovari%f, wind_pc%ucovari%f, swm_simul%avd)

    ! Convert from covariant to contravariant
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(n0, nend) &
    !$OMP SHARED(wind_pc, wind_pu, wind_pv, mesh)
    wind_pc%ucontra%f(n0:nend,n0:nend,:) = &
    wind_pc%ucovari%f(n0:nend,n0:nend,:)*mesh%covari2contra_pc(n0:nend,n0:nend,:)%M(1,1)+&
    wind_pc%vcovari%f(n0:nend,n0:nend,:)*mesh%covari2contra_pc(n0:nend,n0:nend,:)%M(1,2) 

    wind_pc%vcontra%f(n0:nend,n0:nend,:) =&
    wind_pc%ucovari%f(n0:nend,n0:nend,:)*mesh%covari2contra_pc(n0:nend,n0:nend,:)%M(2,1)+&
    wind_pc%vcovari%f(n0:nend,n0:nend,:)*mesh%covari2contra_pc(n0:nend,n0:nend,:)%M(2,2) 
    !$OMP END PARALLEL WORKSHARE

    ! now we fill the ghost cell C grid
    call interp_A2Cduogrid(wind_pu%ucontra%f, wind_pu%vcontra%f, wind_pv%ucontra%f, &
    wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f)

    ! then we interpolate from A grid to the C grid inner cells
    call interp_A2Cgrid(wind_pu%ucontra%f, wind_pv%vcontra%f, wind_pc%ucontra%f, wind_pc%vcontra%f, swm_simul%avd)


    ! Compute time-averaged wind
    call swm_time_averaged_wind(wind_pu, wind_pv, wind_pc, swm_simul%dp, &
        swm_simul%dto2, mesh%dx, mesh, L_pc)

    ! CFL number
    call cfl_x(mesh, wind_pu%ucontra_time_av, cx_pu, swm_simul%dt)
    call cfl_y(mesh, wind_pv%vcontra_time_av, cy_pv, swm_simul%dt)

end subroutine sw_timestep_Cgrid



end module swm_timestep
