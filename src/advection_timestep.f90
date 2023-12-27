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
    pi, &
    nbfaces, &
    i0, iend, &
    j0, jend, &
    n0, nend

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  deg2rad

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  velocity_field, &
  simulation

! Discrete operators 
use discrete_operators, only: &
    divergence, &
    cfl_x, cfl_y

! Advection initial condition
use advection_ic, only: &
  velocity_adv

! Model variables
use advection_vars

! Departure point
use departure_point, only: &
    adv_time_averaged_wind

! Interpolation
use duogrid_interpolation, only: &
    dg_interp


implicit none

contains 

subroutine adv_timestep(mesh)
    !--------------------------------------------------
    ! Compute one time step for the advection
    ! problem on the sphere
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! Time averaged wind - when vf==1, we only need to compute for n=1
    if(advsimul%vf>=2 .or. advsimul%n==1)then
        ! Compute time centered wind
        call adv_centered_wind(wind_pu, wind_pv, mesh)
        
        ! Compute time-averaged wind
        call adv_time_averaged_wind(wind_pu, wind_pv, wind_pc, advsimul%dp, &
        advsimul%dto2, mesh%dx, mesh, L_pc)

        ! CFL number
        call cfl_x(mesh, wind_pu%ucontra_time_av, cx_pu, advsimul%dt)
        call cfl_y(mesh, wind_pv%vcontra_time_av, cy_pv, advsimul%dt)
    end if

    ! Discrete divergence
    call divergence(div_ugq, Q, gQ, wind_pu, wind_pv, cx_pu, cy_pv, &
                      px, py, Qx, Qy, advsimul, mesh, L_pc)

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, advsimul, div_ugq)
    ! Update the solution
    Q%f = Q%f - advsimul%dt*div_ugq%f
    !$OMP END PARALLEL WORKSHARE
end subroutine adv_timestep


subroutine adv_update(V_pu, V_pv, mesh, vf, t)
    !--------------------------------------------------
    ! Compute the velocity update needed in a timestep 
    ! for the advection problem on the sphere
    !
    !
    ! Q - scalar field average values on cells
    ! V_pu - velocity at pu
    ! V_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(velocity_field), intent(inout) :: V_pu, V_pv
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: t

    ! aux
    integer(i4) :: i, j, p

    !aux
    real(kind=8) :: lat, lon
    real(kind=8) :: ulon, vlat, ucontra, vcontra

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(V_pu, mesh) & 
    !$OMP SHARED(i0, iend, j0, jend, nbfaces, vf, t) &
    !$OMP SHARED(n0, nend) &
    !$OMP PRIVATE(i, j, p, ulon, vlat, ucontra, vcontra, lat, lon) &
    !$OMP SCHEDULE(static) 
    ! Vector field at pu
    do i = n0, nend+1
        do j = n0, nend
            do p = 1, nbfaces
                lat  = mesh%pu(i,j,p)%lat
                lon  = mesh%pu(i,j,p)%lon

                ! Update velocity
                call velocity_adv(ulon, vlat, lat, lon, t, vf)

                ! LL2contra
                ucontra = mesh%ll2contra_pu(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(1,2)*vlat
                vcontra = mesh%ll2contra_pu(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(2,2)*vlat
                V_pu%ucontra%f(i,j,p) = ucontra
                V_pu%vcontra%f(i,j,p) = vcontra
            end do
        end do
    end do    
    !$OMP END PARALLEL DO

    ! Vector field at pv
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(V_pv, mesh) & 
    !$OMP SHARED(i0, iend, j0, jend, nbfaces, vf, t) &
    !$OMP SHARED(n0, nend) &
    !$OMP PRIVATE(i, j, p, ulon, vlat, ucontra, vcontra, lat, lon) &
    !$OMP SCHEDULE(static) 
    do i = n0, nend
        do j = n0, nend+1
            do p = 1, nbfaces
                lat  = mesh%pv(i,j,p)%lat
                lon  = mesh%pv(i,j,p)%lon

                ! Update velocity
                call velocity_adv(ulon, vlat, lat, lon, t, vf)

                ! LL2contra
                ucontra = mesh%ll2contra_pv(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(1,2)*vlat
                vcontra = mesh%ll2contra_pv(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(2,2)*vlat

                V_pv%ucontra%f(i,j,p) = ucontra
                V_pv%vcontra%f(i,j,p) = vcontra
            end do
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine adv_update


subroutine adv_centered_wind(V_pu, V_pv, mesh)
    !--------------------------------------------------
    ! Compute the velocity update needed in a timestep 
    ! for the advection problem on the sphere
    !
    !
    ! Q - scalar field average values on cells
    ! V_pu - velocity at pu
    ! V_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(velocity_field), intent(inout) :: V_pu, V_pv

    ! aux
    integer(i4) :: i, j, p

    !aux
    real(kind=8) :: lat, lon
    real(kind=8) :: ulon, vlat, ucontra, vcontra

    if(advsimul%vf==1)then
       !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
       !$OMP SHARED(V_pu, V_pv)
       V_pu%ucontra_time_centered%f(:,:,:) = V_pu%ucontra%f(:,:,:)
       V_pv%vcontra_time_centered%f(:,:,:) = V_pv%vcontra%f(:,:,:)
       !$OMP END PARALLEL WORKSHARE
    else
       !$OMP PARALLEL DO &
       !$OMP DEFAULT(NONE) & 
       !$OMP SHARED(V_pu, mesh, advsimul) & 
       !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
       !$OMP SHARED(n0, nend) &
       !$OMP PRIVATE(i, j, p, ulon, vlat, ucontra, vcontra, lat, lon) &
       !$OMP SCHEDULE(static) 
       ! Vector field at pu
       do i = n0, nend+1
           do j = n0, nend
               do p = 1, nbfaces
                   lat  = mesh%pu(i,j,p)%lat
                   lon  = mesh%pu(i,j,p)%lon

                   ! Update velocity
                   call velocity_adv(ulon, vlat, lat, lon, advsimul%t-advsimul%dto2, advsimul%vf)

                   ! LL2contra
                   ucontra = mesh%ll2contra_pu(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(1,2)*vlat
                   vcontra = mesh%ll2contra_pu(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(2,2)*vlat
                   V_pu%ucontra_time_centered%f(i,j,p) = ucontra
                   V_pu%vcontra_time_centered%f(i,j,p) = vcontra
               end do
           end do
       end do    
       !$OMP END PARALLEL DO

       ! Vector field at pv
       !$OMP PARALLEL DO &
       !$OMP DEFAULT(NONE) & 
       !$OMP SHARED(V_pv, mesh, advsimul) & 
       !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
       !$OMP SHARED(n0, nend) &
       !$OMP PRIVATE(i, j, p, ulon, vlat, ucontra, vcontra, lat, lon) &
       !$OMP SCHEDULE(static) 
       do i = n0, nend
           do j = n0, nend+1
               do p = 1, nbfaces
                   lat  = mesh%pv(i,j,p)%lat
                   lon  = mesh%pv(i,j,p)%lon

                   ! Update velocity
                   call velocity_adv(ulon, vlat, lat, lon, advsimul%t-advsimul%dto2, advsimul%vf)

                   ! LL2contra
                   ucontra = mesh%ll2contra_pv(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(1,2)*vlat
                   vcontra = mesh%ll2contra_pv(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(2,2)*vlat

                   V_pv%ucontra_time_centered%f(i,j,p) = ucontra
                   V_pv%vcontra_time_centered%f(i,j,p) = vcontra
               end do
           end do
       end do
       !$OMP END PARALLEL DO
    endif
end subroutine adv_centered_wind


end module advection_timestep
