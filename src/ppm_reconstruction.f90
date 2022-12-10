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
        r8

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
    integer(i4) :: N

    N = px%n
    select case(px%recon)
      case('ppm') ! PPM from CW84 paper
        ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
        a1 = 7._r8/12._r8
        a2 = a1
        a3 = -1._r8/12._r8
        a4 = a3

        ! Assign values of q_R and q_L
        px%q_R(2:N+5,:,:) = a1*Q%f(2:N+5,:,:) + a2*Q%f(1:N+4,:,:) + a3*Q%f(3:N+6,:,:) + a4*Q%f(0:N+3,:,:)
        px%q_L(3:N+4,:,:) = px%q_R(3:N+4,:,:)
        px%q_R(3:N+4,:,:) = px%q_R(4:N+5,:,:)

     case('hyppm') !Hybrid PPM from PL07
        ! coeffs from equations 41 and 42 from Putman and Lin 2007
        a1 =   2._r8/60._r8
        a2 = -13._r8/60._r8
        a3 =  47._r8/60._r8
        a4 =  27._r8/60._r8
        a5 =  -3._r8/60._r8

        ! Assign values of Q_R and Q_L
        px%q_R(3:N+4,:,:) = a1*Q%f(1:N+2,:,:) + a2*Q%f(2:N+3,:,:) + a3*Q%f(3:N+4,:,:) + a4*Q%f(4:N+5,:,:) + a5*Q%f(5:N+6,:,:)
        px%q_L(3:N+4,:,:) = a5*Q%f(1:N+2,:,:) + a4*Q%f(2:N+3,:,:) + a3*Q%f(3:N+4,:,:) + a2*Q%f(4:N+5,:,:) + a1*Q%f(5:N+6,:,:)

    case default
      print*, 'ERROR on ppm_reconstruction_x: invalid 1D reconstruction method: ', px%recon 
      stop
    end select
 
    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    px%dq(3:N+4,:,:) = px%q_R(3:N+4,:,:) - px%q_L(3:N+4,:,:)
    px%q6(3:N+4,:,:) = 6._r8*Q%f(3:N+4,:,:) - 3._r8*(px%q_R(3:N+4,:,:) + px%q_L(3:N+4,:,:))

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
    integer(i4) :: N

    N = py%n
    select case(py%recon)
      case('ppm') ! PPM from CW84 paper
        ! Values of Q at right edges (q_(j+1/2)) - Formula 1.9 from Collela and Woodward 1984
        a1 = 7._r8/12._r8
        a2 = a1
        a3 = -1._r8/12._r8
        a4 = a3

        ! Assign values of q_R and q_L
        py%q_R(:,2:N+5,:) = a1*Q%f(:,2:N+5,:) + a2*Q%f(:,1:N+4,:) + a3*Q%f(:,3:N+6,:) + a4*Q%f(:,0:N+3,:)
        py%q_L(:,3:N+4,:) = py%q_R(:,3:N+4,:)
        py%q_R(:,3:N+4,:) = py%q_R(:,4:N+5,:)

     case('hyppm') !Hybrid PPM from PL07
        ! coeffs from equations 41 and 42 from Putman and Lin 2007
        a1 =   2._r8/60._r8
        a2 = -13._r8/60._r8
        a3 =  47._r8/60._r8
        a4 =  27._r8/60._r8
        a5 =  -3._r8/60._r8

        ! Assign values of Q_R and Q_L
        !py%q_R(2:N+3,:,:) = a1*Q%f(0:N+1,:,:) + a2*Q%f(1:N+2,:,:) + a3*Q%f(2:N+3,:,:) + a4*Q%f(3:N+4,:,:) + a5*Q%f(4:N+5,:,:)
        !py%q_L(2:N+3,:,:) = a5*Q%f(0:N+1,:,:) + a4*Q%f(1:N+2,:,:) + a3*Q%f(2:N+3,:,:) + a2*Q%f(3:N+4,:,:) + a1*Q%f(4:N+5,:,:)
        py%q_R(:,3:N+4,:) = a1*Q%f(:,1:N+2,:) + a2*Q%f(:,2:N+3,:) + a3*Q%f(:,3:N+4,:) + a4*Q%f(:,4:N+5,:) + a5*Q%f(:,5:N+6,:)
        py%q_L(:,3:N+4,:) = a5*Q%f(:,1:N+2,:) + a4*Q%f(:,2:N+3,:) + a3*Q%f(:,3:N+4,:) + a2*Q%f(:,4:N+5,:) + a1*Q%f(:,5:N+6,:)


    case default
      print*, 'ERROR on ppm_reconstruction_y: invalid 1D reconstruction method: ', py%recon
      stop
    end select
 
    ! Compute the polynomial coefs
    ! q(x) = q_L + z*(dq + q6*(1-z)) z in [0,1]
    py%dq(:,3:N+4,:) = py%q_R(:,3:N+4,:) - py%q_L(:,3:N+4,:)
    py%q6(:,3:N+4,:) = 6._r8*Q%f(:,3:N+4,:) - 3._r8*(py%q_R(:,3:N+4,:) + py%q_L(:,3:N+4,:))

    return 

  end subroutine ppm_reconstruction_y


end module ppm_reconstruction 
