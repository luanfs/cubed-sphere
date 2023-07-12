module ppm_reconstruction
!===============================================================================================
!
! Piecewise Parabolic Method (PPM) polynomial reconstruction module
!
! This module contains routines that performs the PPM reconstruction
! on the cubed-sphere in x and y direction.
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

!Global constants
use constants, only: &
    i4, &
    r8, &
    i0, iend, &
    j0, jend

!Data structures
use datastruct, only: &
  scalar_field, &
  ppm_parabola

implicit none

contains 

subroutine ppm_reconstruction_x(Q, px)
    !---------------------------------------------------------------------------------
    ! PPM_RECONSTRUCTION_X
    !
    ! Given the average values of a scalar field Q, this routine constructs
    ! a piecewise parabolic aproximation of Q using its average value.
    ! The reconstruction is perfomed in x direction of each cubedsphere panel
    !--------------------------------------------------------------------------------
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: px

    !aux
    real(r8) :: a1, a2, a3, a4, a5

    select case(px%recon)
        case('ppm') ! PPM from CW84 paper
            ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
            a1 = 7._r8/12._r8
            a2 = a1
            a3 = -1._r8/12._r8
            a4 = a3

            ! Assign values of q_R and q_L
            px%q_R(i0-2:iend+2,:,:) = a1*Q%f(i0-2:iend+2,:,:) + a2*Q%f(i0-3:iend+1,:,:)&
                                    + a3*Q%f(i0-1:iend+3,:,:) + a4*Q%f(i0-4:iend,:,:)
            px%q_L(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:)
            px%q_R(i0-1:iend+1,:,:) = px%q_R(i0:iend+2,:,:)

        case('hyppm') !Hybrid PPM from PL07
            ! coeffs from equations 41 and 42 from Putman and Lin 2007
            a1 =   2._r8/60._r8
            a2 = -13._r8/60._r8
            a3 =  47._r8/60._r8
            a4 =  27._r8/60._r8
            a5 =  -3._r8/60._r8

            ! Assign values of Q_R and Q_L
            px%q_R(i0-1:iend+1,:,:) = a1*Q%f(i0-3:iend-1,:,:) + a2*Q%f(i0-2:iend,:,:) + a3*Q%f(i0-1:iend+1,:,:)&
                                    + a4*Q%f(i0:iend+2,:,:) + a5*Q%f(i0+1:iend+3,:,:)
            px%q_L(i0-1:iend+1,:,:) = a5*Q%f(i0-3:iend-1,:,:) + a4*Q%f(i0-2:iend,:,:) + a3*Q%f(i0-1:iend+1,:,:)&
                                    + a2*Q%f(i0:iend+2,:,:) + a1*Q%f(i0+1:iend+3,:,:)

        case default
            print*, 'ERROR on ppm_reconstruction_x: invalid 1D reconstruction method: ', px%recon 
            stop
        end select
 
        ! Compute the polynomial coefs
        ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
        px%dq(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:) - px%q_L(i0-1:iend+1,:,:)
        px%q6(i0-1:iend+1,:,:) = 6._r8*Q%f(i0-1:iend+1,:,:) - 3._r8*(px%q_R(i0-1:iend+1,:,:) + px%q_L(i0-1:iend+1,:,:))

       return 

end subroutine ppm_reconstruction_x

subroutine ppm_reconstruction_y(Q, py)
    !---------------------------------------------------------------------------------
    ! PPM_RECONSTRUCTION_Y
    !
    ! Given the average values of a scalar field Q, this routine constructs
    ! a piecewise parabolic aproximation of Q using its average value.
    ! The reconstruction is perfomed in y direction of each cubedsphere panel
    !--------------------------------------------------------------------------------
    type(scalar_field), intent(in) :: Q
    type(ppm_parabola), intent(inout) :: py

    !aux
    real(r8) :: a1, a2, a3, a4, a5

    select case(py%recon)
        case('ppm') ! PPM from CW84 paper
            ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
            a1 = 7._r8/12._r8
            a2 = a1
            a3 = -1._r8/12._r8
            a4 = a3

            ! Assign values of q_R and q_L
            py%q_R(:,i0-2:iend+2,:) = a1*Q%f(:,i0-2:iend+2,:) + a2*Q%f(:,i0-3:iend+1,:)&
                                    + a3*Q%f(:,i0-1:iend+3,:) + a4*Q%f(:,i0-4:iend,:)
            py%q_L(:,i0-1:iend+1,:) = py%q_R(:,i0-1:iend+1,:)
            py%q_R(:,i0-1:iend+1,:) = py%q_R(:,i0:iend+2,:)

        case('hyppm') !Hybrid PPM from PL07
            ! coeffs from equations 41 and 42 from Putman and Lin 2007
            a1 =   2._r8/60._r8
            a2 = -13._r8/60._r8
            a3 =  47._r8/60._r8
            a4 =  27._r8/60._r8
            a5 =  -3._r8/60._r8

            ! Assign values of Q_R and Q_L
            py%q_R(:,i0-1:iend+1,:) = a1*Q%f(:,i0-3:iend-1,:) + a2*Q%f(:,i0-2:iend,:) + a3*Q%f(:,i0-1:iend+1,:)&
                                    + a4*Q%f(:,i0:iend+2,:) + a5*Q%f(:,i0+1:iend+3,:)
            py%q_L(:,i0-1:iend+1,:) = a5*Q%f(:,i0-3:iend-1,:) + a4*Q%f(:,i0-2:iend,:) + a3*Q%f(:,i0-1:iend+1,:)& 
                                    + a2*Q%f(:,i0:iend+2,:) + a1*Q%f(:,i0+1:iend+3,:)

        case default
            print*, 'ERROR on ppm_reconstruction_y: invalid 1D reconstruction method: ', py%recon
            stop
    end select
 
    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    py%dq(:,i0-1:iend+1,:) = py%q_R(:,i0-1:iend+1,:) - py%q_L(:,i0-1:iend+1,:)
    py%q6(:,i0-1:iend+1,:) = 6._r8*Q%f(:,i0-1:iend+1,:) - 3._r8*(py%q_R(:,i0-1:iend+1,:) + py%q_L(:,i0-1:iend+1,:))
    return 

end subroutine ppm_reconstruction_y


end module ppm_reconstruction 
