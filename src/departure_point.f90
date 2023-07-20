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
    n0, nend, &
    i0, iend, &
    j0, jend

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field

implicit none

contains

subroutine adv_time_averaged_wind(wind_pu, wind_pv, dp, dto2, dx)
    !--------------------------------------------------
    ! Compute the time average velocity needed
    ! for the departure point scheme
    !
    !--------------------------------------------------
    type(vector_field), intent(inout) :: wind_pu, wind_pv
    character(len=16) :: dp ! departute point method
    real(r8) :: dto2 ! time step over two
    real(r8) :: dx
    ! aux
    real(r8) :: a, a1, a2, u1, u2
    integer(i4) :: i, j, p

    ! Inner operator
    select case (dp)
        case ('rk1')
            ! RK1
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(wind_pu, wind_pv)
            wind_pu%ucontra_time_av%f = wind_pu%ucontra%f
            wind_pv%vcontra_time_av%f = wind_pv%vcontra%f
            !wind_pu%ucontra_time_av%f = 1.5_r8*wind_pu%ucontra%f-0.5_r8*wind_pu%ucontra_old%f
            !wind_pv%vcontra_time_av%f = 1.5_r8*wind_pv%vcontra%f-0.5_r8*wind_pv%vcontra_old%f
            !$OMP END PARALLEL WORKSHARE
        case ('rk2')
            ! RK2
            ! time extrapolation to obtaind the wind centererd at time

            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(wind_pu, wind_pv)
            wind_pu%ucontra_time_centered%f = 1.5_r8*wind_pu%ucontra%f-0.5_r8*wind_pu%ucontra_old%f
            wind_pv%vcontra_time_centered%f = 1.5_r8*wind_pv%vcontra%f-0.5_r8*wind_pv%vcontra_old%f
            !$OMP END PARALLEL WORKSHARE

            ! wind for dp in x direction
            !$OMP PARALLEL DO &
            !$OMP DEFAULT(NONE) & 
            !$OMP SHARED(wind_pu) & 
            !$OMP SHARED(a, a1, a2, u1, u2, dto2, dx) &
            !$OMP SHARED(n0, nend, i0, iend, nbfaces) &
            !$OMP PRIVATE(i, j, p) &
            !$OMP SCHEDULE(static)
            do j = n0, nend
                do i = i0, iend+1
                    do p = 1, nbfaces
                        ! Linear interpolation weight
                        a = wind_pu%ucontra%f(i,j,p)*dto2/dx

                        ! Upwind linear interpolation
                        if (wind_pu%ucontra%f(i,j,p)>0._r8) then
                            u1 = wind_pu%ucontra_time_centered%f(i-1,j,p)
                            u2 = wind_pu%ucontra_time_centered%f(i,j,p)
                            a1 = a
                            a2 = 1._r8-a
                        else
                            u1 = wind_pu%ucontra_time_centered%f(i,j,p)
                            u2 = wind_pu%ucontra_time_centered%f(i+1,j,p)
                            a1 = 1._r8+a
                            a2 = -a
                        end if
                        wind_pu%ucontra_time_av%f(i,j,p) = a1*u1 + a2*u2
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

            ! wind for dp in y direction
            !$OMP PARALLEL DO &
            !$OMP DEFAULT(NONE) & 
            !$OMP SHARED(wind_pv) & 
            !$OMP SHARED(a, a1, a2, u1, u2, dto2, dx) &
            !$OMP SHARED(n0, nend, j0, jend, nbfaces) &
            !$OMP PRIVATE(i, j, p) &
            !$OMP SCHEDULE(static)
            do i = n0, nend
                do j = j0, jend+1
                    do p = 1, nbfaces
                        ! Linear interpolation weight
                        a = wind_pv%vcontra%f(i,j,p)*dto2/dx

                        ! Upwind linear interpolation
                        if (wind_pv%vcontra%f(i,j,p)>0._r8) then
                            u1 = wind_pv%vcontra_time_centered%f(i,j-1,p)
                            u2 = wind_pv%vcontra_time_centered%f(i,j,p)
                            a1 = a
                            a2 = 1._r8-a
                        else
                            u1 = wind_pv%vcontra_time_centered%f(i,j,p)
                            u2 = wind_pv%vcontra_time_centered%f(i,j+1,p)
                            a1 = 1._r8+a
                            a2 = -a
                        end if
                        wind_pv%vcontra_time_av%f(i,j,p) = a1*u1 + a2*u2
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

        case default
            print*, 'ERROR in adv_time_averaged_wind: invalid operator splitting,  ', dp
            stop
    end select
 
end subroutine adv_time_averaged_wind


end module departure_point
