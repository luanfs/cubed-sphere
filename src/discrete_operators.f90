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

subroutine F_operator(F_gQ, f_pu, mesh, dt)
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
    real(r8), intent(inout), allocatable :: F_gQ(:,:,:)
    real(r8), intent(in), allocatable :: f_pu(:,:,:)
    type(cubedsphere), intent(in) :: mesh
    real(r8), intent(in) :: dt

    F_gQ(i0:iend,:,:) = -(dt/mesh%dx)*(f_pu(i0+1:iend+1,:,:) - f_pu(i0:iend,:,:))

end subroutine F_operator

subroutine G_operator(G_gQ, f_pv, mesh, dt)
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
    real(r8), intent(inout), allocatable :: G_gQ(:,:,:)
    real(r8), intent(in), allocatable :: f_pv(:,:,:)
    type(cubedsphere), intent(in) ::mesh
    real(r8), intent(in) :: dt

    G_gQ(:,j0:jend,:) = -(dt/mesh%dy)*(f_pv(:,j0+1:jend+1,:) - f_pv(:,j0:jend,:))

end subroutine G_operator



subroutine divergence(gQ, wind_pu, wind_pv, cx_pu, cy_pv, &
                        px, py, advsimul, mesh)
    !---------------------------------------------------
    !
    ! Computes the divergence of u*q using the
    ! dimension splliting method
    ! The divergence is given by px%df + py%df
    !---------------------------------------------------
    use advection_vars, only: &
        Qx, Qy

    type(cubedsphere), intent(in) :: mesh
    type(simulation), intent(in) :: advsimul
    type(vector_field), intent(inout) :: wind_pu
    type(vector_field), intent(inout) :: wind_pv
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(in) :: gQ
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
 
    ! Compute fluxes
    wind_pu%ucontra_av%f = wind_pu%ucontra%f
    wind_pv%vcontra_av%f = wind_pv%vcontra%f
    call ppm_flux_pu(gQ, px, wind_pu%ucontra_av, cx_pu)
    call ppm_flux_pv(gQ, py, wind_pv%vcontra_av, cy_pv)

    ! Dimension splitting operators
    call F_operator(px%df, px%f_upw, mesh, advsimul%dt)
    call G_operator(py%df, py%f_upw, mesh, advsimul%dt)

    ! Splitting scheme
    select case (advsimul%opsplit)
    case ('avlt')
        Qx%f = gQ%f+0.5_r8*px%df
        Qy%f = gQ%f+0.5_r8*py%df
        !Qx%f = Qx%f/mesh%mt_pc
        !Qy%f = Qy%f/mesh%mt_pc
    Case ('pl07')
        ! PL07 - equation 17 and 18
        Qx%f = 0.5_r8*(gQ%f + (gQ%f + px%df)/(1._r8-(cx_pu%f(n0+1:,:,:)-cx_pu%f(:nend,:,:))))
        Qy%f = 0.5_r8*(gQ%f + (gQ%f + py%df)/(1._r8-(cy_pv%f(:,n0+1:,:)-cy_pv%f(:,:nend,:))))

    case default
        print*, 'ERROR in divergence: invalid operator splitting,  ', advsimul%opsplit
        stop
    end select

    ! Compute fluxes
    call ppm_flux_pu(Qy, px, wind_pu%ucontra_av, cx_pu)
    call ppm_flux_pv(Qx, py, wind_pv%vcontra_av, cy_pv)

    ! Dimension splitting operators
    call F_operator(px%df, px%f_upw, mesh, advsimul%dt)
    call G_operator(py%df, py%f_upw, mesh, advsimul%dt)
 
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
