module discrete_operators
!===============================================================================================
!   This module contains all the routines that compute discrete operators
!===============================================================================================

! Constants
use constants, only: &
  i4, &
  nbfaces, &
  i0, iend, &
  j0, jend, &
  n0, nend, &
  erad

!Data structures
use datastruct, only: &
  scalar_field, &
  velocity_field, &
  cubedsphere, &
  ppm_parabola, &
  simulation, &
  lagrange_poly_cs

! 1d fluxes 
use ppm_flux, only: &
    ppm_flux_pu, &
    ppm_flux_pv, &
    ppm_fluxes_PL07, &
    numerical_flux_ppm_pu, &
    numerical_flux_ppm_pv

! reconstruction 
use ppm_reconstruction, only: &
    ppm_reconstruction_x, &
    ppm_reconstruction_y

! Mass fixer
use mass_fixer

! Interpolation
use duogrid_interpolation, only: &
    dg_interp


implicit none

contains 

subroutine F_operator(Q, wind_pu, cx_pu, px, mesh, dt)
    !---------------------------------------------------
    !
    ! Dimension splitting operators implementation
    !
    ! References:
    ! Lin, S., & Rood, R. B. (1996). Multidimensional Flux-Form Semi-Lagrangian
    ! Transport Schemes, Monthly Weather Review, 124(9), 2046-2070, from
    ! https://journals.ametsoc.org/view/journals/mwre/124/9/1520-0493_1996_124_2046_mffslt_2_0_co_2.xml
    !
    ! Flux operator in x direction
    ! Formula 2.7 from Lin and Rood 1996
    !---------------------------------------------------
    type(velocity_field), intent(inout) :: wind_pu
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(cubedsphere), intent(inout) :: mesh
    real(kind=8), intent(in) :: dt

    ! Compute fluxes
    if(px%et=='duogrid') then
        call ppm_flux_pu(Q, px, wind_pu%ucontra_time_av, cx_pu, mesh)
    end if

    ! F operator
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(px, mesh, i0, iend, dt)
    px%df(i0:iend,:,:) = -(dt/mesh%dx)*(px%f_upw(i0+1:iend+1,:,:)-px%f_upw(i0:iend,:,:))
    !$OMP END PARALLEL WORKSHARE

end subroutine F_operator

subroutine G_operator(Q, wind_pv, cy_pv, py, mesh, dt)
    !---------------------------------------------------
    !
    ! Dimension splitting operators implementation
    !
    ! References:
    ! Lin, S., & Rood, R. B. (1996). Multidimensional Flux-Form Semi-Lagrangian
    ! Transport Schemes, Monthly Weather Review, 124(9), 2046-2070, from
    ! https://journals.ametsoc.org/view/journals/mwre/124/9/1520-0493_1996_124_2046_mffslt_2_0_co_2.xml
    !
    ! Flux operator in y direction
    ! Formula 2.8 from Lin and Rood 1996
    !---------------------------------------------------
    type(velocity_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(cubedsphere), intent(inout) :: mesh
    real(kind=8), intent(in) :: dt

    ! Compute fluxes
    if(py%et=='duogrid') then
        call ppm_flux_pv(Q, py, wind_pv%vcontra_time_av, cy_pv, mesh)
    end if

    ! G operator
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(py, mesh, j0, jend, dt)
     py%df(:,j0:jend,:) = -(dt/mesh%dy)*(py%f_upw(:,j0+1:jend+1,:)-py%f_upw(:,j0:jend,:))
    !$OMP END PARALLEL WORKSHARE

end subroutine G_operator

subroutine inner_f_operator(Q, wind_pu, cx_pu, px, mesh, dt, sp)
    !---------------------------------------------------
    ! Inner flux operator in x direction
    !---------------------------------------------------
    type(velocity_field), intent(inout) :: wind_pu
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: sp ! splitting method
    real(kind=8), intent(in) :: dt

    ! Compute fluxes
    if(px%et=='duogrid') then
        call ppm_flux_pu(Q, px, wind_pu%ucontra_time_av, cx_pu, mesh)
    end if

    ! F operator
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(px, mesh, i0, iend, dt)
    px%df(i0:iend,:,:) = -(dt/mesh%dx)*(px%f_upw(i0+1:iend+1,:,:)-px%f_upw(i0:iend,:,:))
    !$OMP END PARALLEL WORKSHARE

    ! Inner operator
    select case (sp)
        case ('avlt')
            ! nothing to do here
            !px%df = px%df

        case ('pl07')
            ! PL07 - equation 17 and 18
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(px, mesh, n0, nend, cx_pu, dt)
            px%df = (-px%Q%f + (px%Q%f + px%df)/&
            (1.d0-(cx_pu%f(n0+1:,:,:)*mesh%mt_pu(n0+1:,:,:)-cx_pu%f(:nend,:,:)*mesh%mt_pu(:nend,:,:))))
            !$OMP END PARALLEL WORKSHARE

        case default
            print*, 'ERROR in inner_f_operator: invalid operator splitting,  ', sp 
            stop
    end select
 
end subroutine inner_f_operator

subroutine inner_g_operator(Q, wind_pv, cy_pv, py, mesh, dt, sp)
    !---------------------------------------------------
    ! Inner flux operator in y direction
    !---------------------------------------------------
    type(velocity_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: sp ! splitting method
    real(kind=8), intent(in) :: dt

    ! Compute fluxes
    if(py%et=='duogrid') then
        call ppm_flux_pv(Q, py, wind_pv%vcontra_time_av, cy_pv, mesh)
    end if

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(py, mesh, j0, jend, dt)
    ! G operator
    py%df(:,j0:jend,:) = -(dt/mesh%dy)*(py%f_upw(:,j0+1:jend+1,:)-py%f_upw(:,j0:jend,:))
    !$OMP END PARALLEL WORKSHARE

    ! Inner operator
    select case (sp)
        case ('avlt')
            ! nothing to do here
            !py%df = py%df

        case ('pl07')
            ! PL07 - equation 17 and 18
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(py, mesh, n0, nend, cy_pv, dt)
            py%df = (-py%Q%f + (py%Q%f + py%df)/&
            (1.d0-(cy_pv%f(:,n0+1:,:)*mesh%mt_pv(:,n0+1:,:)-cy_pv%f(:,:nend,:)*mesh%mt_pv(:,:nend,:))))
            !$OMP END PARALLEL WORKSHARE

        case default
            print*, 'ERROR in inner_g_operator: invalid operator splitting,  ', sp 
            stop
    end select
 
end subroutine inner_g_operator


subroutine divergence(div_ugq, Q, wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, advsimul, mesh, L_pc)
    !---------------------------------------------------
    !
    ! Computes the divergence of u*q using the
    ! dimension splliting method
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(inout) :: advsimul
    type(velocity_field), intent(inout) :: wind_pu
    type(velocity_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(scalar_field), intent(inout) :: div_ugq
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(scalar_field), intent(inout) :: Qx ! variable to advect in x direction
    type(scalar_field), intent(inout) :: Qy ! variable to advect in y direction
    type(lagrange_poly_cs), intent(inout) :: L_pc ! lagrange polynomial

    if(advsimul%et=='duogrid') then
        ! Interpolate scalar field to ghost cells
        call dg_interp(Q%f, L_pc)

        ! Dimension splitting operators
        call inner_f_operator(Q, wind_pu, cx_pu, px, mesh, advsimul%dt, advsimul%opsplit)
        call inner_g_operator(Q, wind_pv, cy_pv, py, mesh, advsimul%dt, advsimul%opsplit)

        ! Compute next splitting input
        ! Metric tensor scheme
        select case (advsimul%mt)
        case ('mt0')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(Qx, Qy, px, py, mesh)
            Qx%f = px%Q%f + 0.5d0*px%df
            Qy%f = py%Q%f + 0.5d0*py%df
            Qx%f = Qx%f/mesh%mt_pc
            Qy%f = Qy%f/mesh%mt_pc
            !$OMP END PARALLEL WORKSHARE
        case ('pl07')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(Qx, Qy, px, py, mesh)
            Qx%f = px%Q%f + 0.5d0*px%df
            Qy%f = py%Q%f + 0.5d0*py%df
            !$OMP END PARALLEL WORKSHARE
        case default
            print*, 'ERROR in divergence: invalid metric tensor formulation,  ', advsimul%mt
            stop
        end select

        ! Compute fluxes
        call F_operator(Qy, wind_pu, cx_pu, px, mesh, advsimul%dt)
        call G_operator(Qx, wind_pv, cy_pv, py, mesh, advsimul%dt)

        ! Applies mass fixer (average at cube interfaces)
        if (advsimul%mf=='af') then
            call average_flux_at_cube_intefaces(px, py, mesh%dx, mesh%dy, advsimul%dt)
        end if
     
        ! Compute the divergence
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(div_ugq, px, py, advsimul, mesh)
        div_ugq%f = -(px%df + py%df)/advsimul%dt/mesh%mt_pc
        !$OMP END PARALLEL WORKSHARE

        ! Applies mass fixer (project divergence in nullspace)
        if (advsimul%mf=='lpr' .or. advsimul%mf=='gpr') then
            call divergence_projection(div_ugq, advsimul, mesh)
        end if

    else if(advsimul%et == 'pl07')then
        ! Compute fluxes
        call ppm_fluxes_PL07(Q, Q, px, py, wind_pu%ucontra_time_av, &
        wind_pv%vcontra_time_av, cx_pu, cy_pv, mesh)

        ! Dimension splitting operators
        call inner_f_operator(Q, wind_pu, cx_pu, px, mesh, advsimul%dt, advsimul%opsplit)
        call inner_g_operator(Q, wind_pv, cy_pv, py, mesh, advsimul%dt, advsimul%opsplit)

        ! Compute next splitting input
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(Qx, Qy, px, py)
        Qx%f = px%Q%f+0.5d0*px%df
        Qy%f = py%Q%f+0.5d0*py%df
        !$OMP END PARALLEL WORKSHARE

        ! Metric tensor scheme
        select case (advsimul%mt)
        case ('mt0')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(Qx, Qy, mesh)
            Qx%f = Qx%f/mesh%mt_pc
            Qy%f = Qy%f/mesh%mt_pc
            !$OMP END PARALLEL WORKSHARE
        case ('pl07')
            ! Nothing to do here
            !Qx%f = Qx%f
            !Qy%f = Qy%f

        case default
            print*, 'ERROR in divergence: invalid metric tensor forumalion,  ', advsimul%mt
            stop
        end select

        ! Compute fluxes
        call ppm_fluxes_PL07(Qy, Qx, px, py, wind_pu%ucontra_time_av, &
        wind_pv%vcontra_time_av, cx_pu, cy_pv, mesh)

        call F_operator(Qy, wind_pu, cx_pu, px, mesh, advsimul%dt)
        call G_operator(Qx, wind_pv, cy_pv, py, mesh, advsimul%dt)

        ! Compute the divergence
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(div_ugq, px, py, advsimul, mesh)
        div_ugq%f = -(px%df + py%df)/advsimul%dt/mesh%mt_pc
        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in divergence: invalid edge treatment ', advsimul%et
        stop
    end if
end subroutine divergence


subroutine vorticity_fluxes(div_abs_vort,abs_vort_flux_pu, abs_vort_flux_pv, &
                      rel_vort, abs_vort, fcoriolis_pc,&
                      wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, swmsimul, mesh, L_pc)
    !---------------------------------------------------
    ! Computes the absolute vorticity fluxes at pu and pv
    ! dimension splliting method
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(inout) :: swmsimul
    type(scalar_field), intent(inout) :: abs_vort_flux_pu
    type(scalar_field), intent(inout) :: abs_vort_flux_pv
    type(scalar_field), intent(inout) :: rel_vort
    type(scalar_field), intent(inout) :: abs_vort
    type(scalar_field), intent(inout) :: div_abs_vort
    type(scalar_field), intent(inout) :: fcoriolis_pc
    type(velocity_field), intent(inout) :: wind_pu
    type(velocity_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: cy_pv
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(scalar_field), intent(inout) :: Qx ! variable to advect in x direction
    type(scalar_field), intent(inout) :: Qy ! variable to advect in y direction
    type(lagrange_poly_cs), intent(inout) :: L_pc ! lagrange polynomial

    if(swmsimul%et=='duogrid') then
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

        ! Interpolate scalar field to ghost cells
        ! Dimension splitting operators
        call inner_f_operator(abs_vort, wind_pu, cx_pu, px, mesh, swmsimul%dt, swmsimul%opsplit)
        call inner_g_operator(abs_vort, wind_pv, cy_pv, py, mesh, swmsimul%dt, swmsimul%opsplit)

        ! Compute next splitting input
        ! Metric tensor scheme
        select case (swmsimul%mt)
        case ('mt0')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(Qx, Qy, px, py, mesh)
            Qx%f = px%Q%f + 0.5d0*px%df
            Qy%f = py%Q%f + 0.5d0*py%df
            Qx%f = Qx%f/mesh%mt_pc
            Qy%f = Qy%f/mesh%mt_pc
            !$OMP END PARALLEL WORKSHARE
        case ('pl07')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(Qx, Qy, px, py, mesh)
            Qx%f = px%Q%f + 0.5d0*px%df
            Qy%f = py%Q%f + 0.5d0*py%df
            !$OMP END PARALLEL WORKSHARE
        case default
            print*, 'ERROR in vorticity_fluxes: invalid metric tensor formulation,  ', swmsimul%mt
            stop
        end select

        ! Compute fluxes
        call F_operator(Qy, wind_pu, cx_pu, px, mesh, swmsimul%dt)
        call G_operator(Qx, wind_pv, cy_pv, py, mesh, swmsimul%dt)

        ! Compute the divergence
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(div_abs_vort, px, py, swmsimul, mesh)
        div_abs_vort%f = -(px%df + py%df)/swmsimul%dt/mesh%mt_pc
        !$OMP END PARALLEL WORKSHARE

        abs_vort_flux_pu%f = px%f_upw
        !abs_vort_flux_pu%f = abs_vort_flux_pu%f/swmsimul%dt
        !abs_vort_flux_pu%f = abs_vort_flux_pu%f/mesh%mt_pc

        abs_vort_flux_pv%f = py%f_upw
        !abs_vort_flux_pv%f = abs_vort_flux_pv%f/swmsimul%dt
        !abs_vort_flux_pv%f = abs_vort_flux_pv%f/mesh%mt_pc

    else
        print*, 'ERROR in vorticity_fluxes: invalid edge treatment ', swmsimul%et
        stop
    end if
end subroutine vorticity_fluxes



subroutine ke_fluxes(wind_pu, wind_pv, wind_po, cx_po, cy_po, &
                      Ku_px, Kv_py, K_po, swmsimul, mesh)
    !---------------------------------------------------
    ! Computes the kinetic energy fluxes at po
    ! time averaged winds must be already computed
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(inout) :: swmsimul
    type(velocity_field), intent(inout) :: wind_pu
    type(velocity_field), intent(inout) :: wind_pv
    type(velocity_field), intent(inout) :: wind_po
    type(scalar_field), intent(inout) :: K_po
    type(scalar_field), intent(inout) :: cx_po
    type(scalar_field), intent(inout) :: cy_po
    type(ppm_parabola), intent(inout) :: Ku_px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: Kv_py ! ppm in y direction

    ! PPM reconstructions
    call ppm_reconstruction_x(wind_pv%ucovari, Ku_px)
    call ppm_reconstruction_y(wind_pu%vcovari, Kv_py)

    ! Compute the fluxes
    call numerical_flux_ppm_pu(wind_pv%ucovari, Ku_px, wind_po%ucontra_time_av, cx_po, mesh)
    call numerical_flux_ppm_pv(wind_pu%vcovari, Kv_py, wind_po%vcontra_time_av, cy_po, mesh)
    ! now we can compute the kinetic energy
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(K_po, Ku_px, Kv_py)
    K_po%f = (Ku_px%f_upw + Kv_py%f_upw)*0.5d0
    !$OMP END PARALLEL WORKSHARE
 
end subroutine ke_fluxes


subroutine cfl_x(mesh, ucontra_pu, cx_pu, dt)
    !-----------------------------------------------------
    ! Routine for computing the CFL number at pu in x 
    ! direction
    !----------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(in)  :: ucontra_pu
    type(scalar_field), intent(inout) :: cx_pu
    real(kind=8), intent(in)::dt

    ! Compute CFL
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(cx_pu, ucontra_pu, dt, mesh)
    cx_pu%f = ucontra_pu%f*(dt/mesh%dx)
    !$OMP END PARALLEL WORKSHARE

end subroutine cfl_x


subroutine cfl_y(mesh, vcontra_pv, cy_pv, dt)
    !-----------------------------------------------------
    ! Routine for computing the CFL number at pv in y 
    ! direction
    !----------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(in)  :: vcontra_pv
    type(scalar_field), intent(inout) :: cy_pv
    real(kind=8), intent(in)::dt

    ! Compute CFL
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(cy_pv, vcontra_pv, dt, mesh)
    cy_pv%f = vcontra_pv%f*(dt/mesh%dy)
    !$OMP END PARALLEL WORKSHARE

end subroutine cfl_y

end module discrete_operators 
