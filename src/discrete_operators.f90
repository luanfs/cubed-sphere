module discrete_operators
  !===============================================================================================
  !   This module contains all the miscellaneous routines
  !===============================================================================================

  ! Constants
  use constants, only: &
      i4, &
      r8, &
      nbfaces, &
      i0, iend, &
      j0, jend

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

  subroutine F_operator(F_gQ, wind_pu, f_pu, mesh, dt)
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
    type(vector_field), intent(in) :: wind_pu
    type(cubedsphere), intent(in) :: mesh

    !aux
    real(r8) :: dt

    F_gQ(i0:iend,:,:) = -(dt/mesh%area(i0:iend,:,:))*mesh%dy* &
                         (wind_pu%ucontra%f(i0:iend,:,:)*f_pu(i0:iend,:,:) - &
                          wind_pu%ucontra%f(i0-1:iend-1,:,:)*f_pu(i0-1:iend-1,:,:))

  end subroutine F_operator

  subroutine G_operator(G_gQ, wind_pv, f_pv, mesh, dt)
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
    type(vector_field), intent(in) :: wind_pv
    type(cubedsphere), intent(in) ::mesh

    !aux
    real(r8) :: dt

    G_gQ(:,j0:jend,:) = -(dt/mesh%area(:,j0:jend,:))*mesh%dx* &
                         (wind_pv%vcontra%f(:,j0:jend,:)*f_pv(:,j0:jend,:) - &
                          wind_pv%vcontra%f(:,j0-1:jend-1,:)*f_pv(:,j0-1:jend-1,:))

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
    type(vector_field), intent(in) :: wind_pu
    type(vector_field), intent(in) :: wind_pv
    type(scalar_field), intent(inout) :: cx_pu
    type(scalar_field), intent(inout) :: cy_pv
    type(scalar_field), intent(in) :: gQ
    type(ppm_parabola), intent(inout) :: px ! ppm in x direction
    type(ppm_parabola), intent(inout) :: py ! ppm in y direction
    integer(i4) :: n, ng

    n = mesh%n
    ng = mesh%nbg
 
    ! Compute fluxes
    call ppm_flux_pu(gQ, px, cx_pu)
    call ppm_flux_pv(gQ, py, cy_pv)

    ! Dimension splitting operators
    call F_operator(px%df, wind_pu, px%f_upw, mesh, advsimul%dt)
    call G_operator(py%df, wind_pv, py%f_upw, mesh, advsimul%dt)

    ! Splitting scheme
    select case (advsimul%opsplit)
    case ('lr96')
        Qx%f = gQ%f+0.5_r8*px%df
        Qy%f = gQ%f+0.5_r8*py%df

    case ('l04')
        ! L04 equation 7 and 8
        px%dF = px%df + (cx_pu%f(1:,:,:)-cx_pu%f(:N+ng,:,:))*gQ%f
        py%dF = py%df + (cy_pv%f(:,1:,:)-cy_pv%f(:,:N+ng,:))*gQ%f
        Qx%f = gQ%f+0.5_r8*px%df
        Qy%f = gQ%f+0.5_r8*py%df
    case ('pl07')
        ! PL07 - equation 17 and 18
        Qx%f = 0.5_r8*(gQ%f + (gQ%f + px%df)/(1._r8-(cx_pu%f(1:,:,:)-cx_pu%f(:N+ng,:,:))))
        Qy%f = 0.5_r8*(gQ%f + (gQ%f + py%df)/(1._r8-(cy_pv%f(:,1:,:)-cy_pv%f(:,:N+ng,:))))

    case default
        print*, 'ERROR in divergence: invalid operator splitting,  ', advsimul%opsplit
    end select

    ! Compute fluxes
    call ppm_flux_pu(Qy, px, cx_pu)
    call ppm_flux_pv(Qx, py, cy_pv)

    ! Dimension splitting operators
    call F_operator(px%df, wind_pu, px%f_upw, mesh, advsimul%dt)
    call G_operator(py%df, wind_pv, py%f_upw, mesh, advsimul%dt)
 
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
