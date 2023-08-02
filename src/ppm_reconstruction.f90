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
    real(kind=8) :: a1, a2, a3, a4, a5

    select case(px%recon)
        case('ppm') ! PPM from CW84 paper
            ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
            a1 = 7.d0/12.d0
            a2 = a1
            a3 = -1.d0/12.d0
            a4 = a3

            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(px, Q, a1, a2, a3, a4, i0, iend)
            ! Assign values of q_R and q_L
            px%q_R(i0-2:iend+2,:,:) = a1*Q%f(i0-2:iend+2,:,:) + a2*Q%f(i0-3:iend+1,:,:)&
                                    + a3*Q%f(i0-1:iend+3,:,:) + a4*Q%f(i0-4:iend,:,:)
            px%q_L(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:)
            px%q_R(i0-1:iend+1,:,:) = px%q_R(i0:iend+2,:,:)
            !$OMP END PARALLEL WORKSHARE

        case('hyppm') !Hybrid PPM from PL07
            ! coeffs from equations 41 and 42 from Putman and Lin 2007
            a1 =   2.d0/60.d0
            a2 = -13.d0/60.d0
            a3 =  47.d0/60.d0
            a4 =  27.d0/60.d0
            a5 =  -3.d0/60.d0

            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(px, Q, a1, a2, a3, a4, a5, i0, iend)
            ! Assign values of Q_R and Q_L
            px%q_R(i0-1:iend+1,:,:) = a1*Q%f(i0-3:iend-1,:,:) + a2*Q%f(i0-2:iend,:,:) + a3*Q%f(i0-1:iend+1,:,:)&
                                    + a4*Q%f(i0:iend+2,:,:) + a5*Q%f(i0+1:iend+3,:,:)
            px%q_L(i0-1:iend+1,:,:) = a5*Q%f(i0-3:iend-1,:,:) + a4*Q%f(i0-2:iend,:,:) + a3*Q%f(i0-1:iend+1,:,:)&
                                    + a2*Q%f(i0:iend+2,:,:) + a1*Q%f(i0+1:iend+3,:,:)
            !$OMP END PARALLEL WORKSHARE

        case default
            print*, 'ERROR on ppm_reconstruction_x: invalid 1D reconstruction method: ', px%recon 
            stop
        end select
 
        ! Compute the polynomial coefs
        ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
        !px%dq(i0-1:iend+1,:,:) = px%q_R(i0-1:iend+1,:,:) - px%q_L(i0-1:iend+1,:,:)
        !px%q6(i0-1:iend+1,:,:) = 6.d0*Q%f(i0-1:iend+1,:,:) - 3.d0*(px%q_R(i0-1:iend+1,:,:) + px%q_L(i0-1:iend+1,:,:))

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
    real(kind=8) :: a1, a2, a3, a4, a5

    select case(py%recon)
        case('ppm') ! PPM from CW84 paper
            ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
            a1 = 7.d0/12.d0
            a2 = a1
            a3 = -1.d0/12.d0
            a4 = a3

            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(py, Q, a1, a2, a3, a4, j0, jend)
            ! Assign values of q_R and q_L
            py%q_R(:,j0-2:jend+2,:) = a1*Q%f(:,j0-2:jend+2,:) + a2*Q%f(:,j0-3:jend+1,:)&
                                    + a3*Q%f(:,j0-1:jend+3,:) + a4*Q%f(:,j0-4:jend,:)
            py%q_L(:,j0-1:jend+1,:) = py%q_R(:,j0-1:jend+1,:)
            py%q_R(:,j0-1:jend+1,:) = py%q_R(:,j0:jend+2,:)
            !$OMP END PARALLEL WORKSHARE

        case('hyppm') !Hybrid PPM from PL07
            ! coeffs from equations 41 and 42 from Putman and Lin 2007
            a1 =   2.d0/60.d0
            a2 = -13.d0/60.d0
            a3 =  47.d0/60.d0
            a4 =  27.d0/60.d0
            a5 =  -3.d0/60.d0

            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(py, Q, a1, a2, a3, a4, a5, j0, jend)
            ! Assign values of Q_R and Q_L
            py%q_R(:,j0-1:jend+1,:) = a1*Q%f(:,j0-3:jend-1,:) + a2*Q%f(:,j0-2:jend,:) + a3*Q%f(:,j0-1:jend+1,:)&
                                    + a4*Q%f(:,j0:jend+2,:) + a5*Q%f(:,j0+1:jend+3,:)
            py%q_L(:,j0-1:jend+1,:) = a5*Q%f(:,j0-3:jend-1,:) + a4*Q%f(:,j0-2:jend,:) + a3*Q%f(:,j0-1:jend+1,:)& 
                                    + a2*Q%f(:,j0:jend+2,:) + a1*Q%f(:,j0+1:jend+3,:)

            !$OMP END PARALLEL WORKSHARE

         case default
            print*, 'ERROR on ppm_reconstruction_y: invalid 1D reconstruction method: ', py%recon
            stop
    end select

    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    !py%dq(:,i0-1:iend+1,:) = py%q_R(:,i0-1:iend+1,:) - py%q_L(:,i0-1:iend+1,:)
    !py%q6(:,i0-1:iend+1,:) = 6.d0*Q%f(:,i0-1:iend+1,:) - 3.d0*(py%q_R(:,i0-1:iend+1,:) + py%q_L(:,i0-1:iend+1,:))
    return 

end subroutine ppm_reconstruction_y

subroutine average_parabolas_at_cube_intefaces(px, py)
    !---------------------------------------------------------------------------------
    ! AVERAGE_PARABOLAS_AT_CUBE_INTERFACES
    !
    ! Given the PPM parabolas px and py, this routine
    ! average the values of the edge reconstruction at the
    ! cube edge points
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: px, py
    real(kind=8) :: a, b

    a = 0.5d0
    b = 1.d0 - a

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(px, py, i0, iend, j0, jend, a, b)
    ! Average panels 1-2,2-3,3-4,4-1
    px%q_R(iend+1,j0:jend,1:3) = a*px%q_R(iend+1,j0:jend,1:3) + b*px%q_L(i0,j0:jend,2:4)
    px%q_L(i0,j0:jend,2:4)     = px%q_R(iend+1,j0:jend,1:3)

    px%q_R(iend+1,j0:jend,4) = a*px%q_R(iend+1,j0:jend,4) + b*px%q_L(i0,j0:jend,1)
    px%q_L(i0,j0:jend,1)  = px%q_R(iend+1,j0:jend,4)

    ! Average panels 1-5
    py%q_L(i0:iend,j0,5)   = a*py%q_L(i0:iend,j0,5) + b*py%q_R(i0:iend,jend+1,1)
    py%q_R(i0:iend,jend+1,1) = py%q_L(i0:iend,j0,5)

    ! Average panels 2-5
    px%q_R(iend+1,j0:jend,5) = a*px%q_R(iend+1,j0:jend,5) + b*py%q_R(i0:iend,jend+1,2)
    py%q_R(i0:iend,jend+1,2) = px%q_R(iend+1,j0:jend,5)

    ! Average panels 3-5
    py%q_R(i0:iend,jend+1,5) = a*py%q_R(i0:iend,jend+1,5) + b*py%q_R(iend:i0:-1,jend+1,3)
    py%q_R(i0:iend,jend+1,3) = py%q_R(iend:i0:-1,jend+1,5)

    ! Average panels 4-5
    px%q_L(i0,j0:jend,5)    = a*px%q_L(i0,j0:jend,5) + b*(py%q_R(iend:i0:-1,jend+1,4))
    py%q_R(i0:iend,jend+1,4) = (px%q_L(i0,jend:j0:-1,5))

    ! Average panels 1-6
    py%q_R(i0:iend,jend+1,6) = a*py%q_R(i0:iend,jend+1,6) + b*py%q_L(i0:iend,j0,1)
    py%q_L(i0:iend,j0,1)   = py%q_R(i0:iend,jend+1,6)

    ! Average panels 2-6
    py%q_L(i0:iend,j0,2)     = a*py%q_L(i0:iend,j0,2) + b*px%q_R(iend+1,jend:j0:-1,6)
    px%q_R(iend+1,j0:jend,6) = py%q_L(iend:i0:-1,j0,2)

    ! Average panels 3-6
    py%q_L(i0:iend,j0,3) = a*py%q_L(i0:iend,j0,3) + b*py%q_L(iend:i0:-1,j0,6)
    py%q_L(i0:iend,j0,6) = py%q_L(iend:i0:-1,j0,3)

    ! Average panels 4-6
    py%q_L(i0:iend,j0,4) = a*px%q_L(i0,j0:jend,6) + b*py%q_L(i0:iend,j0,4)
    px%q_L(i0,j0:jend,6) = py%q_L(i0:iend,j0,4)

    !$OMP END PARALLEL WORKSHARE

end subroutine average_parabolas_at_cube_intefaces

subroutine edges_extrapolation(Qx, Qy, px, py)
    !---------------------------------------------------------------------------------
    ! AVERAGE_PARABOLAS_AT_CUBE_INTERFACES
    !
    ! Given the PPM parabolas px and py, this routine
    ! average the values of the edge reconstruction at the
    ! cube edge points
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: px, py
    type(scalar_field), intent(inout) :: Qx, Qy

    call average_parabolas_at_cube_intefaces(px, py)


end subroutine edges_extrapolation






end module ppm_reconstruction 
