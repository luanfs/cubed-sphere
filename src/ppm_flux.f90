module ppm_flux
  !===============================================================================================
  !
  ! Piecewise Parabolic Method (PPM) flux evaluation module
  !
  ! This module contains routines that performs the PPM flux computation 
  ! on the cubed-sphere at pu and pv points
  !
  ! The PPM reconstruction is performed at ppm_reconstruction module
  !
  ! Luan da Fonseca Santos - 2022
  !
  ! References:
  ! -  Phillip Colella, Paul R Woodward, The Piecewise Parabolic Method (PPM) for gas-dynamical simulations,
  ! Journal of Computational Physics, Volume 54, Issue 1, 1984, Pages 174-201, ISSN 0021-9991,
  ! https://doi.org/10.1016/0021-9991(84)90143-8.
  !
  ! -  Carpenter , R. L., Jr., Droegemeier, K. K., Woodward, P. R., & Hane, C. E. (1990).
  ! Application of the Piecewise Parabolic Method (PPM) to Meteorological Modeling, Monthly Weather Review, 118(3),
  ! 586-612.
  !
  !===============================================================================================

  ! Constants
  use constants, only: &
      i4, &
      r8, &
      nbfaces, &
      i0, iend, &
      j0, jend, &
      nghost

  !Data structures
  use datastruct, only: &
      scalar_field, &
      ppm_parabola

  !PPM reconstruction
  use ppm_reconstruction, only: &
        ppm_reconstruction_x, &
        ppm_reconstruction_y

 implicit none

contains 

  subroutine ppm_flux_pu(Q, px, cx_pu)
    !---------------------------------------------------------------------------------
    ! PPM_FLUX_PU
    !
    ! Given the average values of a scalar field Q, this routine reconstructs
    ! the PPM aproximation using the reconstruction routine
    ! from module ppm_reconstruction and evaluates the flux
    ! at pu points
    !--------------------------------------------------------------------------------
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: px    ! parabola
    type(scalar_field), intent(inout) :: cx_pu       ! CFL in x direction at pu points 

    integer(i4)::N
    select case(px%recon)
      case('ppm', 'hyppm')
        ! Reconstructs the values of Q using a piecewise parabolic polynomial
        call ppm_reconstruction_x(Q, px)

        ! Compute the fluxes
        call numerical_flux_ppm_pu(px, cx_pu)

        ! Get upwind flux
        call upwind_flux_pu(px, cx_pu)
 
     case default
      print*, 'ERROR on ppm_flux_pu: invalid 1D flux method: ', px%recon 
      stop
    end select
    return 

  end subroutine ppm_flux_pu


  subroutine ppm_flux_pv(Q, py, cy_pv)
    !---------------------------------------------------------------------------------
    ! PPM_FLUX_PV
    !
    ! Given the average values of a scalar field Q, this routine reconstructs
    ! the PPM aproximation using the reconstruction routine
    ! from module ppm_reconstruction and evaluates the flux
    ! at pv points
    !--------------------------------------------------------------------------------
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: py       ! parabola
    type(scalar_field), intent(inout) :: cy_pv    ! CFL in y direction at pv points 

    select case(py%recon)
      case('ppm', 'hyppm')
        ! Reconstructs the values of Q using a piecewise parabolic polynomial
        call ppm_reconstruction_y(Q, py)

        ! Compute the fluxes
        call numerical_flux_ppm_pv(py, cy_pv)

        ! Get upwind flux
        call upwind_flux_pv(py, cy_pv)

     case default
      print*, 'ERROR on ppm_flux_pv: invalid 1D flux method: ', py%recon
      stop
    end select
    return 

  end subroutine ppm_flux_pv

  subroutine numerical_flux_ppm_pu(px, cx_pu)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PU
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pu points
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: px ! parabola
    type(scalar_field), intent(inout) :: cx_pu    ! CFL in x direction at pu points 
    integer(i4) :: N

    ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
    ! Flux at left edges
    px%f_L(i0-1:iend,:,:) = px%q_R(i0-1:iend,:,:) - cx_pu%f(i0-1:iend,:,:)*0.5*&
                       (px%dq(i0-1:iend,:,:) - (1.0-(2.0/3.0)*cx_pu%f(i0-1:iend,:,:))*px%q6(i0-1:iend,:,:))

    ! Flux at right edges
    px%f_R(i0-1:iend,:,:) = px%q_L(i0:iend+1,:,:) - cx_pu%f(i0-1:iend,:,:)*0.5*&
                       (px%dq(i0:iend+1,:,:) + (1.0+(2.0/3.0)*cx_pu%f(i0-1:iend,:,:) )*px%q6(i0:iend+1,:,:))

 end subroutine numerical_flux_ppm_pu

  subroutine numerical_flux_ppm_pv(py, cy_pv)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PV
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pv points
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: py ! parabola
    type(scalar_field), intent(inout) :: cy_pv    ! CFL in y direction at pv points 

    ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
    ! Flux at left edges
    py%f_L(:,i0-1:iend,:) = py%q_R(:,i0-1:iend,:) - cy_pv%f(:,i0-1:iend,:)*0.5*(py%dq(:,i0-1:iend,:) - &
                       (1.0-(2.0/3.0)*cy_pv%f(:,i0-1:iend,:))*py%q6(:,i0-1:iend,:))

    ! Flux at right edges
    py%f_R(:,i0-1:iend,:) = py%q_L(:,i0:iend+1,:) - cy_pv%f(:,i0-1:iend,:)*0.5*(py%dq(:,i0:iend+1,:) + &
                       (1.0+(2.0/3.0)*cy_pv%f(:,i0-1:iend,:))*py%q6(:,i0:iend+1,:))

 end subroutine numerical_flux_ppm_pv


  subroutine upwind_flux_pu(px, cx_pu)
    !---------------------------------------------------------------------------------
    ! UPWIND_FLUX_PPM_PU
    !
    ! Given the fluxes at right and left, this routine computes the upwind flux at pu
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: px    ! flux values at edges (pu points)
    type(scalar_field), intent(in) :: cx_pu  ! cfl of u (contravariant velocity) at pu points 

    ! aux
    integer(i4) :: i, j, p, N

    N = px%n

    do p = 1, nbfaces
      do i = i0-1, iend
        do j = 0, N+nghost
          ! Upwind flux
          if(cx_pu%f(i,j,p) > 0._r8)then
            px%f_upw(i,j,p) = px%f_L(i,j,p)
          else
            px%f_upw(i,j,p) = px%f_R(i,j,p)
          end if
        end do
      end do
    end do
  end subroutine upwind_flux_pu

  subroutine upwind_flux_pv(py, cy_pv)
    !---------------------------------------------------------------------------------
    ! UPWIND_FLUX_PPM_PV
    !
    ! Given the fluxes at right and left, this routine computes the upwind flux at pv
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: py
    type(scalar_field), intent(in) :: cy_pv ! cfl of v (contravariant velocity) at pu points 

    ! aux
    integer(i4) :: i, j, p, N

    N = py%n
    do p = 1, nbfaces
      do i = 0, N+nghost
        do j = j0-1, jend
          ! Upwind flux
          if(cy_pv%f(i,j,p) >= 0._r8)then
            py%f_upw(i,j,p) = py%f_L(i,j,p)
          else
            py%f_upw(i,j,p) = py%f_R(i,j,p)
          end if
        end do
      end do
    end do
  end subroutine upwind_flux_pv


end module ppm_flux
