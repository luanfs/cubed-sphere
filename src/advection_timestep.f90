module advection_timestep
!===============================================================================================
! Module for routines needed in a advection timestep
!
! Luan da Fonseca Santos - 2022
! (luan.santos@usp.br)
!===============================================================================================

!Global constants
use constants, only: &
    i4, &
    r8, &
    pi, &
    nbfaces

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  deg2rad, &
  ll2contra, &
  contra2ll

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field

! Discrete operators 
use discrete_operators, only: &
    divergence, &
    cfl_x, cfl_y

! Advection initial condition
use advection_ic, only: &
  velocity_adv

! Model variables
use advection_vars

implicit none

contains 

subroutine adv_timestep(mesh)
    use advection_vars
    !--------------------------------------------------
    ! Compute the velocity update needed in a timestep 
    ! for the advection problem on the sphere
    !
    ! ic - initial conditions
    ! vf - velocity field
    !
    ! Q - scalar field average values on cells
    ! V_pu - velocity at pu
    ! V_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    
    ! CFL number
    call cfl_x(mesh, wind_pu, cx_pu, advsimul%dt)
    call cfl_y(mesh, wind_pv, cy_pv, advsimul%dt)

    ! Discrete divergence
    call divergence(div_ugq, Q, wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, advsimul, mesh)

    ! Update the solution
    Q%f = Q%f - advsimul%dt*div_ugq%f
end subroutine adv_timestep


subroutine adv_update(Q, V_pu, V_pv, mesh, ic, vf)
    !--------------------------------------------------
    ! Compute the velocity update needed in a timestep 
    ! for the advection problem on the sphere
    !
    ! ic - initial conditions
    ! vf - velocity field
    !
    ! Q - scalar field average values on cells
    ! V_pu - velocity at pu
    ! V_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(scalar_field), intent(out) :: Q
    type(vector_field), intent(out) :: V_pu, V_pv
    integer(i4), intent(in) :: ic, vf

    ! aux
    integer(i4) :: i0, iend, j0, jend
    integer(i4) :: i, j, p

    !aux
    real(r8) :: lat, lon
    real(r8) :: ulon, vlat, ucontra, vcontra
    

    ! interior grid indexes
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    ! Vector field at pu
    do p = 1 , nbfaces
        do i = i0-1, iend
            do j = j0, jend
                lat  = mesh%pu(i,j,p)%lat
                lon  = mesh%pu(i,j,p)%lon

                call velocity_adv(ulon, vlat, lat, lon, 0._r8, vf)
                call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pu(i,j,p)%M)
                V_pu%ucontra%f(i,j,p) = ucontra
                V_pu%vcontra%f(i,j,p) = vcontra
            end do
        end do
    end do    

    ! Vector field at pv
    do p = 1 , nbfaces
        do i = i0, iend
            do j = j0-1, jend
                lat  = mesh%pv(i,j,p)%lat
                lon  = mesh%pv(i,j,p)%lon

                call velocity_adv(ulon, vlat, lat, lon, 0._r8, vf)
                call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pv(i,j,p)%M)
                V_pv%ucontra%f(i,j,p) = ucontra
                V_pv%vcontra%f(i,j,p) = vcontra
            end do
        end do
    end do
 
end subroutine adv_update

end module advection_timestep
