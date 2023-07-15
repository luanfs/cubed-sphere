module departure_point
!===============================================================================================
! Module for routines that are used to compute the departure point
!
! Luan da Fonseca Santos - 2022
! (luan.santos@usp.br)
!===============================================================================================


!Global constants
use constants, only: &
    i4, &
    r8, &
    nbfaces, &
    i0, iend, &
    j0, jend

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field
implicit none

contains

subroutine adv_time_averaged_wind(V_pu, V_pv, mesh, dp)
    !--------------------------------------------------
    ! Compute the time average velocity needed
    ! for the departure point scheme
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(vector_field), intent(out) :: V_pu, V_pv
    character(len=16) :: dp ! departute point method
    ! aux
    integer(i4) :: i, j, p

    ! Inner operator
    select case (dp)
        case ('rk1')
            ! nothing to do here
            !px%df = px%df

        case ('rk2')
            ! PL07 - equation 17 and 18
            !px%df = (-px%Q%f + (px%Q%f + px%df)/&

        case default
            print*, 'ERROR in adv_time_averaged_wind: invalid operator splitting,  ', dp
            stop
    end select
 
end subroutine adv_time_averaged_wind


end module departure_point
