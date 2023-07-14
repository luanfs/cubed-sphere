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
  n0, nend

!Data structures
use datastruct, only: &
  scalar_field, &
  ppm_parabola, &
  cubedsphere

!PPM reconstruction
use ppm_reconstruction, only: &
    ppm_reconstruction_x, &
    ppm_reconstruction_y

implicit none

contains 

subroutine ppm_flux_pu(Q, px, V_pu_av, cx_pu, mesh, mt)
    !---------------------------------------------------------------------------------
    ! PPM_FLUX_PU
    !
    ! Given the average values of a scalar field Q, this routine reconstructs
    ! the PPM aproximation using the reconstruction routine
    ! from module ppm_reconstruction and evaluates the flux
    ! at pu points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px      ! parabola
    type(scalar_field), intent(inout) :: V_pu_av ! time averaged wind at pu points 
    type(scalar_field), intent(inout) :: cx_pu   ! CFL number of V_pu_av
    character(len=16) :: mt ! metric tensor formulation

    select case(px%recon)
        case('ppm', 'hyppm')
            ! Reconstructs the values of Q using a piecewise parabolic polynomial
            call ppm_reconstruction_x(Q, px)

            ! Compute the fluxes
            call numerical_flux_ppm_pu(Q, px, cx_pu, mesh)

            ! Get upwind flux
            call upwind_flux_pu(px, V_pu_av, cx_pu, mesh, mt)
        case default
            print*, 'ERROR on ppm_flux_pu: invalid 1D flux method: ', px%recon 
            stop
    end select
    return 

end subroutine ppm_flux_pu


subroutine ppm_flux_pv(Q, py, V_pv_av, cy_pv, mesh, mt)
    !---------------------------------------------------------------------------------
    ! PPM_FLUX_PV
    !
    ! Given the average values of a scalar field Q, this routine reconstructs
    ! the PPM aproximation using the reconstruction routine
    ! from module ppm_reconstruction and evaluates the flux
    ! at pv points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py      ! parabola
    type(scalar_field), intent(inout) :: V_pv_av ! time averaged wind at pv points 
    type(scalar_field), intent(inout) :: cy_pv   ! CFL number of V_pv_av
    character(len=16) :: mt ! metric tensor formulation

    select case(py%recon)
        case('ppm', 'hyppm')
            ! Reconstructs the values of Q using a piecewise parabolic polynomial
            call ppm_reconstruction_y(Q, py)

            ! Compute the fluxes
            call numerical_flux_ppm_pv(Q, py, cy_pv, mesh)

            ! Get upwind flux
            call upwind_flux_pv(py, V_pv_av, cy_pv, mesh, mt)

        case default
            print*, 'ERROR on ppm_flux_pv: invalid 1D flux method: ', py%recon
            stop
    end select
    return 

end subroutine ppm_flux_pv

subroutine numerical_flux_ppm_pu(Q, px, cx_pu, mesh)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PU
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pu points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: px ! parabola
    type(scalar_field), intent(inout) :: cx_pu   ! CFL number of V_pu_av
    integer(i4) :: N

    select case(px%mt)
        case('mt0')
            px%q_L(i0-1:iend+1,:,:) = px%q_L(i0-1:iend+1,:,:)*mesh%mt_pu(i0-1:iend+1,:,:)
            px%q_R(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:)*mesh%mt_pu(i0:iend+2,:,:)
            px%Q%f(i0-1:iend+1,:,:) = Q%f(i0-1:iend+1,:,:)*mesh%mt_pc(i0-1:iend+1,:,:)

        case('pl07')
            px%Q%f(i0-1:iend+1,:,:) = Q%f(i0-1:iend+1,:,:)

        case default
            print*, 'ERROR on numerical_flux_ppm_pu: invalid 1D metric tensor method: ', px%mt
            stop
    end select

    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    px%dq(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:) - px%q_L(i0-1:iend+1,:,:)
    px%q6(i0-1:iend+1,:,:) = 6._r8*px%Q%f(i0-1:iend+1,:,:) - 3._r8*(px%q_R(i0-1:iend+1,:,:) + px%q_L(i0-1:iend+1,:,:))

    ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
    ! Flux at left edges
    px%f_L(i0:iend+1,:,:) = px%q_R(i0-1:iend,:,:) - cx_pu%f(i0:iend+1,:,:)*0.5*&
    (px%dq(i0-1:iend,:,:) - (1.0-(2.0/3.0)*cx_pu%f(i0:iend+1,:,:))*px%q6(i0-1:iend,:,:))

    ! Flux at right edges
    px%f_R(i0:iend+1,:,:) = px%q_L(i0:iend+1,:,:) - cx_pu%f(i0:iend+1,:,:)*0.5*&
    (px%dq(i0:iend+1,:,:) + (1.0+(2.0/3.0)*cx_pu%f(i0:iend+1,:,:) )*px%q6(i0:iend+1,:,:))

    return 
end subroutine numerical_flux_ppm_pu

subroutine numerical_flux_ppm_pv(Q, py, cy_pv, mesh)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PV
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pv points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(inout) :: Q
    type(ppm_parabola), intent(inout) :: py ! parabola
    type(scalar_field), intent(inout) :: cy_pv   ! CFL number of V_pv_av

    select case(py%mt)
        case('mt0')
            py%q_L(:,j0-1:jend+1,:) = py%q_L(:,j0-1:jend+1,:)*mesh%mt_pv(:,j0-1:jend+1,:)
            py%q_R(:,j0-1:jend+1,:) = py%q_R(:,j0-1:jend+1,:)*mesh%mt_pv(:,j0:jend+2,:)
            py%Q%f(:,j0-1:jend+1,:) = Q%f(:,j0-1:jend+1,:)*mesh%mt_pc(:,j0-1:jend+1,:)

        case('pl07')
            py%Q%f(:,j0-1:jend+1,:) = Q%f(:,j0-1:jend+1,:)

        case default
            print*, 'ERROR on numerical_flux_ppm_pv: invalid 1D metric tensor method: ', py%mt
            stop
    end select

    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    py%dq(:,j0-1:jend+1,:) = py%q_R(:,j0-1:jend+1,:) - py%q_L(:,j0-1:jend+1,:)
    py%q6(:,j0-1:jend+1,:) = 6._r8*py%Q%f(:,j0-1:jend+1,:) - 3._r8*(py%q_R(:,j0-1:jend+1,:) + py%q_L(:,j0-1:jend+1,:))


    ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
    ! Flux at left edges
    py%f_L(:,j0:jend+1,:) = py%q_R(:,j0-1:jend,:) - cy_pv%f(:,j0:jend+1,:)*0.5*(py%dq(:,i0-1:iend,:) - &
    (1.0-(2.0/3.0)*cy_pv%f(:,j0-1:jend,:))*py%q6(:,j0-1:jend,:))

    ! Flux at right edges
    py%f_R(:,j0:jend+1,:) = py%q_L(:,j0:jend+1,:) - cy_pv%f(:,j0:jend+1,:)*0.5*(py%dq(:,j0:jend+1,:) + &
    (1.0+(2.0/3.0)*cy_pv%f(:,j0:jend+1,:))*py%q6(:,j0:jend+1,:))

end subroutine numerical_flux_ppm_pv


subroutine upwind_flux_pu(px, V_pu_av, cx_pu, mesh, mt)
    !---------------------------------------------------------------------------------
    ! UPWIND_FLUX_PPM_PU
    !
    ! Given the fluxes at right and left, this routine computes the upwind flux at pu
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(ppm_parabola), intent(inout) :: px    ! flux values at edges (pu points)
    type(scalar_field), intent(in) :: V_pu_av  ! time averaged u (contravariant velocity) at pu points 
    type(scalar_field), intent(in) :: cx_pu  ! cfl of u (contravariant velocity) at pu points 
    character(len=16) :: mt ! metric tensor formulation
    ! aux
    integer(i4) :: i, j, p

    do p = 1, nbfaces
        do i = i0, iend+1
            do j = n0, nend
                ! Upwind flux
                if(cx_pu%f(i,j,p) >= 0._r8)then
                    px%f_upw(i,j,p) = V_pu_av%f(i,j,p)*px%f_L(i,j,p)
                else
                    px%f_upw(i,j,p) = V_pu_av%f(i,j,p)*px%f_R(i,j,p)
                end if
            end do
        end do
    end do

    if (mt == 'pl07') then
        px%f_upw(:,:,:) = px%f_upw(:,:,:)*mesh%mt_pu(:,:,:)
    end if

end subroutine upwind_flux_pu

subroutine upwind_flux_pv(py, V_pv_av, cy_pv, mesh, mt)
    !---------------------------------------------------------------------------------
    ! UPWIND_FLUX_PPM_PV
    !
    ! Given the fluxes at right and left, this routine computes the upwind flux at pv
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(ppm_parabola), intent(inout) :: py
    type(scalar_field), intent(in) :: V_pv_av  ! time averaged v (contravariant velocity) at pu points 
    type(scalar_field), intent(in) :: cy_pv ! cfl of v (contravariant velocity) at pv points 
    character(len=16) :: mt ! metric tensor formulation
    ! aux
    integer(i4) :: i, j, p

    do p = 1, nbfaces
        do i = n0, nend
            do j = j0, jend+1
                ! Upwind flux
                if(cy_pv%f(i,j,p) >= 0._r8)then
                    py%f_upw(i,j,p) = V_pv_av%f(i,j,p)*py%f_L(i,j,p)
                else
                    py%f_upw(i,j,p) = V_pv_av%f(i,j,p)*py%f_R(i,j,p)
                end if
            end do
        end do
    end do

    if (mt == 'pl07') then
        py%f_upw(:,:,:) = py%f_upw(:,:,:)*mesh%mt_pv(:,:,:)
    end if

end subroutine upwind_flux_pv

end module ppm_flux
