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
    nbfaces, &
    n0, nend, &
    i0, iend, &
    j0, jend

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field, &
  lagrange_poly_cs

! Interpolation
use duogrid_interpolation, only: &
    interp_A2Cduogrid, &
    interp_C2Aduogrid

implicit none

contains

subroutine adv_time_averaged_wind(wind_pu, wind_pv, wind_pc, dp, dto2, dx, mesh, L_pc)
    !--------------------------------------------------
    ! Compute the time average velocity needed
    ! for the departure point scheme
    !
    !--------------------------------------------------
    type(vector_field), intent(inout) :: wind_pu, wind_pv, wind_pc
    type(cubedsphere), intent(inout) :: mesh
    type(lagrange_poly_cs), intent(inout) :: L_pc
    character(len=16) :: dp ! departute point method
    real(kind=8) :: dto2 ! time step over two
    real(kind=8) :: dx
    ! aux
    real(kind=8) :: a, a1, a2, u1, u2
    integer(i4) :: i, j, p

    ! Interpolation of the wind at ghost cells
    ! first we interpolate to the A grid ghost cells
    call interp_C2Aduogrid(wind_pu%ucontra%f, wind_pv%vcontra%f, wind_pc%u%f, wind_pc%v%f, &
    wind_pc%ucontra%f, wind_pc%vcontra%f, L_pc, mesh%contra2ll_pc, mesh%ll2contra_pc)

    ! now we fill the ghost cell C grid
    call interp_A2Cduogrid(wind_pu%ucontra%f, wind_pu%vcontra%f, wind_pv%ucontra%f, wind_pv%vcontra%f, &
    wind_pc%ucontra%f, wind_pc%vcontra%f)

    select case (dp)
        case ('rk1')
            ! RK1
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(wind_pu, wind_pv)
            wind_pu%ucontra_time_av%f = wind_pu%ucontra_old%f
            wind_pv%vcontra_time_av%f = wind_pv%vcontra_old%f
            !$OMP END PARALLEL WORKSHARE
        case ('rk2')
            ! RK2

            ! time extrapolation to obtaind the wind centered at time
            !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
            !$OMP SHARED(wind_pu, wind_pv)
            wind_pu%ucontra_time_centered%f = 1.5d0*wind_pu%ucontra%f-0.5d0*wind_pu%ucontra_old%f
            wind_pv%vcontra_time_centered%f = 1.5d0*wind_pv%vcontra%f-0.5d0*wind_pv%vcontra_old%f
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
                        if (wind_pu%ucontra%f(i,j,p)>0.d0) then
                            u1 = wind_pu%ucontra_time_centered%f(i-1,j,p)
                            u2 = wind_pu%ucontra_time_centered%f(i,j,p)
                            a1 = a
                            a2 = 1.d0-a
                        else
                            u1 = wind_pu%ucontra_time_centered%f(i,j,p)
                            u2 = wind_pu%ucontra_time_centered%f(i+1,j,p)
                            a1 = 1.d0+a
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
                        if (wind_pv%vcontra%f(i,j,p)>0.d0) then
                            u1 = wind_pv%vcontra_time_centered%f(i,j-1,p)
                            u2 = wind_pv%vcontra_time_centered%f(i,j,p)
                            a1 = a
                            a2 = 1.d0-a
                        else
                            u1 = wind_pv%vcontra_time_centered%f(i,j,p)
                            u2 = wind_pv%vcontra_time_centered%f(i,j+1,p)
                            a1 = 1.d0+a
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
