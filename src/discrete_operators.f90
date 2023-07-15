module discrete_operators
!===============================================================================================
!   This module contains all the routines that compute discrete operators
!===============================================================================================

! Constants
use constants, only: &
  i4, &
  r8, &
  nbfaces, &
  i0, iend, &
  j0, jend, &
  n0, nend

!Data structures
use datastruct, only: &
  scalar_field, &
  vector_field, &
  cubedsphere, &
  ppm_parabola, &
  simulation

! 1d fluxes 
use ppm_flux, only: &
  ppm_flux_pu, &
  ppm_flux_pv

implicit none

contains 

subroutine F_operator(Q, wind_pu, cx_pu, px, mesh, dt, mt)
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
    type(vector_field), intent(inout) :: wind_pu
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: mt ! metric tensor formulation
    real(r8), intent(in) :: dt

    ! Compute fluxes
    call ppm_flux_pu(Q, px, wind_pu%ucontra_av, cx_pu, mesh, mt)

    ! F operator
    px%df(i0:iend,:,:) = -(dt/mesh%dx)*(px%f_upw(i0+1:iend+1,:,:)-px%f_upw(i0:iend,:,:))

end subroutine F_operator

subroutine G_operator(Q, wind_pv, cy_pv, py, mesh, dt, mt)
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
    type(vector_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: mt ! metric tensor formulation
    real(r8), intent(in) :: dt

    ! Compute fluxes
    call ppm_flux_pv(Q, py, wind_pv%vcontra_av, cy_pv, mesh, mt)

    ! G operator
    py%df(:,j0:jend,:) = -(dt/mesh%dy)*(py%f_upw(:,j0+1:jend+1,:)-py%f_upw(:,j0:jend,:))

end subroutine G_operator

subroutine inner_f_operator(Q, wind_pu, cx_pu, px, mesh, dt, mt, sp)
    !---------------------------------------------------
    ! Inner flux operator in x direction
    !---------------------------------------------------
    type(vector_field), intent(inout) :: wind_pu
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: mt ! metric tensor formulation
    character(len=16) :: sp ! splitting method
    real(r8), intent(in) :: dt

    ! Compute fluxes
    call ppm_flux_pu(Q, px, wind_pu%ucontra_av, cx_pu, mesh, mt)

    ! F operator
    px%df(i0:iend,:,:) = -(dt/mesh%dx)*(px%f_upw(i0+1:iend+1,:,:)-px%f_upw(i0:iend,:,:))

    ! Inner operator
    select case (sp)
        case ('avlt')
            ! nothing to do here
            !px%df = px%df

        case ('pl07')
            ! PL07 - equation 17 and 18
            px%df = (-px%Q%f + (px%Q%f + px%df)/&
            (1._r8-(cx_pu%f(n0+1:,:,:)*mesh%mt_pu(n0+1:,:,:)-cx_pu%f(:nend,:,:)*mesh%mt_pu(:nend,:,:))))

        case default
            print*, 'ERROR inner_f_operator: invalid operator splitting,  ', sp 
            stop
    end select
 
end subroutine inner_f_operator

subroutine inner_g_operator(Q, wind_pv, cy_pv, py, mesh, dt, mt, sp)
    !---------------------------------------------------
    ! Inner flux operator in y direction
    !---------------------------------------------------
    type(vector_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(cubedsphere), intent(inout) :: mesh
    character(len=16) :: mt ! metric tensor formulation
    character(len=16) :: sp ! splitting method
    real(r8), intent(in) :: dt

    ! Compute fluxes
    call ppm_flux_pv(Q, py, wind_pv%vcontra_av, cy_pv, mesh, mt)

    ! G operator
    py%df(:,j0:jend,:) = -(dt/mesh%dy)*(py%f_upw(:,j0+1:jend+1,:)-py%f_upw(:,j0:jend,:))

    ! Inner operator
    select case (sp)
        case ('avlt')
            ! nothing to do here
            !py%df = py%df

        case ('pl07')
            ! PL07 - equation 17 and 18
            py%df = (-py%Q%f + (py%Q%f + py%df)/&
            (1._r8-(cy_pv%f(:,n0+1:,:)*mesh%mt_pv(:,n0+1:,:)-cy_pv%f(:,:nend,:)*mesh%mt_pv(:,:nend,:))))

        case default
            print*, 'ERROR inner_f_operator: invalid operator splitting,  ', sp 
            stop
    end select
 
end subroutine inner_g_operator


subroutine divergence(div_ugq, Q, wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, advsimul, mesh)
    !---------------------------------------------------
    !
    ! Computes the divergence of u*q using the
    ! dimension splliting method
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(in) :: advsimul
    type(vector_field), intent(inout) :: wind_pu
    type(vector_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(inout) :: Q
    type(scalar_field), intent(inout) :: div_ugq
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    type(scalar_field), intent(inout) :: Qx ! variable to advect in x direction
    type(scalar_field), intent(inout) :: Qy ! variable to advect in y direction

    ! initialize some vars
    wind_pu%ucontra_av%f = wind_pu%ucontra%f
    wind_pv%vcontra_av%f = wind_pv%vcontra%f
 
    ! Dimension splitting operators
    call inner_f_operator(Q, wind_pu, cx_pu, px, mesh, advsimul%dt, advsimul%mt, advsimul%opsplit)
    call inner_g_operator(Q, wind_pv, cy_pv, py, mesh, advsimul%dt, advsimul%mt, advsimul%opsplit)

    ! Compute next splitting input
    Qx%f = px%Q%f+0.5_r8*px%df
    Qy%f = py%Q%f+0.5_r8*py%df

    ! Metric tensor scheme
    select case (advsimul%mt)
    case ('mt0')
        Qx%f = Qx%f/mesh%mt_pc
        Qy%f = Qy%f/mesh%mt_pc

    case ('pl07')
        ! Nothing to do here
        !Qx%f = Qx%f
        !Qy%f = Qy%f

    case default
        print*, 'ERROR in divergence: invalid operator splitting,  ', advsimul%opsplit
        stop
    end select
 
    ! Compute fluxes
    call F_operator(Qy, wind_pu, cx_pu, px, mesh, advsimul%dt, advsimul%mt)
    call G_operator(Qx, wind_pv, cy_pv, py, mesh, advsimul%dt, advsimul%mt)

    ! Compute the divergence
    div_ugq%f = -(px%df + py%df)/advsimul%dt/mesh%mt_pc

end subroutine divergence


subroutine cfl_x(mesh, wind_pu, cx_pu, dt)
    !-----------------------------------------------------
    ! Routine for computing the CFL number at pu in x 
    ! direction
    !----------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(vector_field), intent(in)  :: wind_pu
    type(scalar_field), intent(inout) :: cx_pu
    real(r8), intent(in)::dt

    ! Compute CFL
    cx_pu%f = wind_pu%ucontra%f*(dt/mesh%dx)

end subroutine cfl_x


subroutine cfl_y(mesh, wind_pv, cy_pv, dt)
    !-----------------------------------------------------
    ! Routine for computing the CFL number at pv in y 
    ! direction
    !----------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(vector_field), intent(in)  :: wind_pv
    type(scalar_field), intent(inout) :: cy_pv
    real(r8), intent(in)::dt

    ! Compute CFL
    cy_pv%f = wind_pv%vcontra%f*(dt/mesh%dy)

end subroutine cfl_y

end module discrete_operators 
