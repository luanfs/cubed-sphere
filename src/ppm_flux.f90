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
    ppm_reconstruction_y, &
    edges_extrapolation

!edge treatment for pl07
use duogrid_interpolation, only: &
    gethalodata_PL07

implicit none

contains 

subroutine ppm_flux_pu(Q, px, V_pu_av, cx_pu, mesh)
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

    select case(px%recon)
        case('ppm', 'hyppm')
            ! Reconstructs the values of Q using a piecewise parabolic polynomial
            call ppm_reconstruction_x(Q, px)

            ! Compute the fluxes
            call numerical_flux_ppm_pu(Q, px, V_pu_av, cx_pu, mesh)

        case default
            print*, 'ERROR on ppm_flux_pu: invalid 1D reconstruction method: ', px%recon 
            stop
    end select
    return 

end subroutine ppm_flux_pu


subroutine ppm_flux_pv(Q, py, V_pv_av, cy_pv, mesh)
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

    select case(py%recon)
        case('ppm', 'hyppm')
            ! Reconstructs the values of Q using a piecewise parabolic polynomial
            call ppm_reconstruction_y(Q, py)

            ! Compute the fluxes
            call numerical_flux_ppm_pv(Q, py, V_pv_av, cy_pv, mesh)

        case default
            print*, 'ERROR on ppm_flux_pv: invalid 1D reconstruction method: ', py%recon
            stop
    end select
    return 

end subroutine ppm_flux_pv


subroutine ppm_fluxes_PL07(Qx, Qy, px, py, V_pu_av, V_pv_av, cx_pu, cy_pv, mesh)
    !---------------------------------------------------------------------------------
    ! PPM_FLUX_PL07
    !
    ! Given the average values of  scalar field Qx and Qy, this routine reconstructs
    ! the PPM aproximation of Qx and Qy using the reconstruction routine
    ! from module ppm_reconstruction and evaluates the flux
    ! at pu and pv points
    ! Uses the formalation from PL07
    !
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: Qx, Qy
    type(ppm_parabola), intent(inout) :: px, py  ! parabola
    type(scalar_field), intent(inout) :: V_pu_av, V_pv_av !time averaged wind at pu points 
    type(scalar_field), intent(inout) :: cx_pu, cy_pv   ! CFL number

    if(px%et=='pl07') then
        select case(px%recon)
            case('ppm', 'hyppm')
                ! Reconstructs the values of Qx and Qy using a piecewise parabolic polynomial
                call gethalodata_PL07(Qx, Qy)
                call ppm_reconstruction_x(Qx, px)
                call ppm_reconstruction_y(Qy, py)
                !call edges_extrapolation(Qx, Qy, px, py)

                ! Compute the fluxes
                call numerical_flux_ppm_pu(Qx, px, V_pu_av, cx_pu, mesh)
                call numerical_flux_ppm_pv(Qy, py, V_pv_av, cy_pv, mesh)

            case default
                print*, 'ERROR on ppm_fluxes_pl07: invalid 1D reconstruction method: ', px%recon 
                stop
        end select
    end if
    return 

end subroutine ppm_fluxes_PL07


subroutine numerical_flux_ppm_pu(Q, px, V_pu_av, cx_pu, mesh)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PU
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pu points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: px ! parabola
    type(scalar_field), intent(in) :: V_pu_av ! time averaged wind at pu points 
    type(scalar_field), intent(in) :: cx_pu   ! CFL number of V_pu_av
    integer(i4) :: i, j, p

    select case(px%mt)
        case('mt0')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(px, Q, mesh, i0, iend)
            px%q_L(i0-1:iend+1,:,:) = px%q_L(i0-1:iend+1,:,:)*mesh%mt_pu(i0-1:iend+1,:,:)
            px%q_R(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:)*mesh%mt_pu(i0:iend+2,:,:)
            px%Q%f(i0-1:iend+1,:,:) = Q%f(i0-1:iend+1,:,:)*mesh%mt_pc(i0-1:iend+1,:,:)
            !$OMP END PARALLEL WORKSHARE

        case('pl07')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(px, Q, i0, iend)
            px%Q%f(i0-1:iend+1,:,:) = Q%f(i0-1:iend+1,:,:)
            !$OMP END PARALLEL WORKSHARE

        case default
            print*, 'ERROR on numerical_flux_ppm_pu: invalid 1D metric tensor method: ', px%mt
            stop
    end select

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(px, i0, iend)
    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    px%dq(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:) - px%q_L(i0-1:iend+1,:,:)
    px%q6(i0-1:iend+1,:,:) = 6.d0*px%Q%f(i0-1:iend+1,:,:) - 3.d0*(px%q_R(i0-1:iend+1,:,:) + px%q_L(i0-1:iend+1,:,:))
    !$OMP END PARALLEL WORKSHARE

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(V_pu_av, cx_pu, px) & 
    !$OMP SHARED(n0, nend, i0, iend, nbfaces) &
    !$OMP PRIVATE(i, j, p) &
    !$OMP SCHEDULE(static)
    ! Upwind fluxes
    do i = i0, iend+1
        do j = n0, nend
            do p = 1, nbfaces
                ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
                if(cx_pu%f(i,j,p) >= 0.d0)then
                    ! Flux at left edges
                    px%f_upw(i,j,p) = V_pu_av%f(i,j,p)*(px%q_R(i-1,j,p) + &
                    0.5d0*cx_pu%f(i,j,p)*(px%q6(i-1,j,p) - px%dq(i-1,j,p)) - &
                    (cx_pu%f(i,j,p)*cx_pu%f(i,j,p)/3.d0)*px%q6(i-1,j,p))

                else
                    ! Flux at right edges
                    px%f_upw(i,j,p) = V_pu_av%f(i,j,p)*(px%q_L(i,j,p) - &
                    0.5d0*cx_pu%f(i,j,p)*(px%q6(i,j,p) + px%dq(i,j,p)) - &
                    (cx_pu%f(i,j,p)*cx_pu%f(i,j,p)/3.d0)*px%q6(i,j,p))

                end if
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! metric tensor multiplication
    if (px%mt == 'pl07') then
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(px, mesh)
        px%f_upw(:,:,:) = px%f_upw(:,:,:)*mesh%mt_pu(:,:,:)
        !$OMP END PARALLEL WORKSHARE
    end if

    return 
end subroutine numerical_flux_ppm_pu

subroutine numerical_flux_ppm_pv(Q, py, V_pv_av, cy_pv, mesh)
    !---------------------------------------------------------------------------------
    ! NUMERICAL_FLUX_PPM_PV
    !
    ! Given the ppm coefficients, this routine computes the fluxes
    ! at pv points
    !--------------------------------------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: py ! parabola
    type(scalar_field), intent(in) :: V_pv_av ! time averaged wind at pv points 
    type(scalar_field), intent(in) :: cy_pv   ! CFL number of V_pv_av
    integer(i4) :: i, j, p

    select case(py%mt)
        case('mt0')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(py, Q, mesh, j0, jend)
            py%q_L(:,j0-1:jend+1,:) = py%q_L(:,j0-1:jend+1,:)*mesh%mt_pv(:,j0-1:jend+1,:)
            py%q_R(:,j0-1:jend+1,:) = py%q_R(:,j0-1:jend+1,:)*mesh%mt_pv(:,j0:jend+2,:)
            py%Q%f(:,j0-1:jend+1,:) = Q%f(:,j0-1:jend+1,:)*mesh%mt_pc(:,j0-1:jend+1,:)
            !$OMP END PARALLEL WORKSHARE

        case('pl07')
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(py, Q, j0, jend)
            py%Q%f(:,j0-1:jend+1,:) = Q%f(:,j0-1:jend+1,:)
            !$OMP END PARALLEL WORKSHARE

        case default
            print*, 'ERROR on numerical_flux_ppm_pv: invalid 1D metric tensor method: ', py%mt
            stop
    end select

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(py, j0, jend)
    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    py%dq(:,j0-1:jend+1,:) = py%q_R(:,j0-1:jend+1,:) - py%q_L(:,j0-1:jend+1,:)
    py%q6(:,j0-1:jend+1,:) = 6.d0*py%Q%f(:,j0-1:jend+1,:) - 3.d0*(py%q_R(:,j0-1:jend+1,:) + py%q_L(:,j0-1:jend+1,:))
    !$OMP END PARALLEL WORKSHARE

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(V_pv_av, cy_pv, py) & 
    !$OMP SHARED(n0, nend, j0, jend, nbfaces) &
    !$OMP PRIVATE(i, j, p) &
    !$OMP SCHEDULE(static)
    do i = n0, nend
        do j = j0, jend+1
            do p = 1, nbfaces
                ! Compute the fluxes (formula 1.12 from Collela and Woodward 1984)
                if(cy_pv%f(i,j,p) >= 0.d0)then
                    ! Flux at left edges
                    py%f_upw(i,j,p) = V_pv_av%f(i,j,p)*(py%q_R(i,j-1,p) + &
                    0.5d0*cy_pv%f(i,j,p)*(py%q6(i,j-1,p) - py%dq(i,j-1,p)) - &
                    (cy_pv%f(i,j,p)*cy_pv%f(i,j,p)/3.d0)*py%q6(i,j-1,p))

                else
                    ! Flux at right edges
                    py%f_upw(i,j,p) = V_pv_av%f(i,j,p)*(py%q_L(i,j,p) - &
                    0.5d0*cy_pv%f(i,j,p)*(py%q6(i,j,p) + py%dq(i,j,p)) - &
                    (cy_pv%f(i,j,p)*cy_pv%f(i,j,p)/3.d0)*py%q6(i,j,p))

                end if
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! metric tensor multiplication
    if (py%mt == 'pl07') then
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(py, mesh) 
        py%f_upw(:,:,:) = py%f_upw(:,:,:)*mesh%mt_pv(:,:,:)
        !$OMP END PARALLEL WORKSHARE
    end if

end subroutine numerical_flux_ppm_pv

end module ppm_flux
