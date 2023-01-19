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
      cubedsphere

  ! 1d fluxes 
  use ppm_flux, only: &
      ppm_flux_pu, &
      ppm_flux_pv

 implicit none

contains 

  subroutine F_operator(F_gQ, gQ, wind_pu, f_pu, mesh, dt)
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
    type(scalar_field), intent(inout) :: F_gQ
    type(scalar_field), intent(in) :: gQ
    type(vector_field), intent(in) :: wind_pu
    real(r8), intent(in), allocatable :: f_pu(:,:,:)
    type(cubedsphere), intent(in) :: mesh

    !aux
    real(r8) :: dt

    F_gQ%f(i0:iend,:,:) = -(dt/mesh%area(i0:iend,:,:))*mesh%dy* &
                         (wind_pu%ucontra%f(i0:iend,:,:)*f_pu(i0:iend,:,:) - &
                          wind_pu%ucontra%f(i0-1:iend-1,:,:)*f_pu(i0-1:iend-1,:,:))

  end subroutine F_operator

  subroutine G_operator(G_gQ, gQ, wind_pv, f_pv, mesh, dt)
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
    type(scalar_field), intent(inout) :: G_gQ
    type(scalar_field), intent(in) :: gQ
    type(vector_field), intent(in) :: wind_pv
    real(r8), intent(in), allocatable :: f_pv(:,:,:)
    type(cubedsphere), intent(in) ::mesh

    !aux
    real(r8) :: dt

    G_gQ%f(:,j0:jend,:) = -(dt/mesh%area(:,j0:jend,:))*mesh%dx* &
                         (wind_pv%vcontra%f(:,j0:jend,:)*f_pv(:,j0:jend,:) - &
                          wind_pv%vcontra%f(:,j0-1:jend-1,:)*f_pv(:,j0-1:jend-1,:))

  end subroutine G_operator



  subroutine divergence(Q, wind_pu, wind_pv, mesh)
    !---------------------------------------------------
    !
    ! Computes the divergence of u*q using the
    ! dimension splliting method
    !
    !---------------------------------------------------
    use advection_vars, only: &
        G_gQ,   F_gQ, &
        FG_gQ, GF_gQ, &
        cx_pu, cy_pv, &
        px, py, gQ, &
        Qx, Qy, &
        div_ugq, &
        advsimul

    type(cubedsphere), intent(in) :: mesh
    type(vector_field), intent(in) :: wind_pu
    type(vector_field), intent(in) :: wind_pv
    type(scalar_field), intent(in) :: Q
    integer(i4) :: n

    n = mesh%n

    ! Multiply Q by the metric tensor
    gQ%f = Q%f*mesh%sinc
 
    ! Compute fluxes
    call ppm_flux_pu(gQ, px, cx_pu, wind_pu%ucontra)
    call ppm_flux_pv(gQ, py, cy_pv, wind_pv%vcontra)

    ! Dimension splliting operators
    call F_operator(F_gQ, gQ, wind_pu, px%f_upw, mesh, advsimul%dt)
    call G_operator(G_gQ, gQ, wind_pv, py%f_upw, mesh, advsimul%dt)

    ! Advect in x and y directions
    Qx%f = gQ%f + 0.5_r8*F_gQ%f 
    Qy%f = gQ%f + 0.5_r8*G_gQ%f

    ! Compute fluxes
    call ppm_flux_pu(Qy, px, cx_pu, wind_pu%ucontra)
    call ppm_flux_pv(Qx, py, cy_pv, wind_pv%vcontra)

    ! Dimension spllinting operators
    call F_operator(FG_gQ, Qy, wind_pu, px%f_upw, mesh, advsimul%dt)
    call G_operator(GF_gQ, Qx, wind_pv, py%f_upw, mesh, advsimul%dt)
 
    ! Compute the divergence
    div_ugq%f(i0:iend,j0:jend,:) = (FG_gQ%f(i0:iend,j0:jend,:) + GF_gQ%f(i0:iend,j0:jend,:))/advsimul%dt

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
