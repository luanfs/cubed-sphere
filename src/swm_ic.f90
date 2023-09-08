module swm_ic
!===============================================================================================
! Module for shallow water test case set up (initial condition, exact solution and etc)
!
! Luan da Fonseca Santos - October 2023
! (luan.santos@usp.br)
!===============================================================================================

!Global constants
use constants, only: &
    i4, &
    pi, &
    nbfaces, &
    i0, iend, &
    j0, jend, &
    n0, nend, &
    grav, gravi, &
    omega

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  deg2rad

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  velocity_field, &
  simulation, &
  lagrange_poly_cs

! Data allocation
use allocation, only: &
  scalar_field_allocation, &
  velocity_field_allocation, &
  allocate_swm_vars

use diagnostics, only: &
  mass_computation

! duo grid interpolation
use duogrid_interpolation, only: &
    compute_lagrange_cs, &
    interp_D2Aduogrid, &
    interp_A2Cduogrid, &
    interp_C2Agrid, &
    interp_A2Cgrid

use swm_timestep, only: &
    sw_timestep_Cgrid

! Constants
use constants,  only: &
    erad, &
    day2sec, &
    sec2day
implicit none

contains 


subroutine init_swm_vars(mesh)
    use swm_vars
    !--------------------------------------------------
    ! Initialize the swm simulation variables
    ! This routine also allocate the fields
    !
    ! ic - initial conditions
    !
    ! H - scalar field average values on cells
    ! wind_pu - velocity at pu
    ! wind_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    ! aux
    integer(i4) :: i, j, p

    ! ppm direction
    px%dir = 1 ! x direction
    py%dir = 2 ! y direction
    Ku_px%dir = 1 ! x direction
    Kv_py%dir = 2 ! y direction

    ! ppm reference point
    px%point = 1 ! pc
    py%point = 1 ! pc
    Ku_px%point = 2 ! po
    Kv_py%point = 2 ! po

    ! Reconstruction scheme
    px%recon = swm_simul%recon1d
    py%recon = swm_simul%recon1d
    Ku_px%recon = swm_simul%recon1d
    Kv_py%recon = swm_simul%recon1d


    ! Metric tensor formulation
    px%mt = swm_simul%mt
    py%mt = swm_simul%mt
    px%et = swm_simul%et
    py%et = swm_simul%et
    Ku_px%mt = 'plane'
    Kv_py%mt = 'plane'
    Ku_px%et = swm_simul%et
    Kv_py%et = swm_simul%et


    ! N
    px%n = mesh%n
    py%n = mesh%n
    Ku_px%n = mesh%n
    Kv_py%n = mesh%n


    ! Lagrange polynomial at centers
    L_pc%degree =  swm_simul%id
    L_pc%order =  L_pc%degree+1
    L_pc%pos = 1

    ! Allocate the variables
    call allocate_swm_vars(mesh)

    ! Compute lagrange polynomials
    call compute_lagrange_cs(L_pc, mesh)
    !swm_simul%avd = 3

    ! Time step over 2
    swm_simul%dto2 = swm_simul%dt*0.5d0
    swm_simul%dt2 = swm_simul%dt*swm_simul%dt

    ! Final time step converted to seconds
    swm_simul%tf = swm_simul%tf*day2sec

    ! Number of time steps
    swm_simul%nsteps = int(swm_simul%tf/swm_simul%dt)

    ! adjust time step
    swm_simul%dt = swm_simul%Tf/swm_simul%nsteps
    swm_simul%nplot = swm_simul%nsteps/(swm_simul%nplot-1)
    swm_simul%plotcounter = 0

    if (swm_simul%ic == 0) then
        swm_simul%nsteps = 1
    end if

    ! Compute the initial conditions
    call compute_ic_swm(H_exact, wind_pu, wind_pv, wind_pc, mesh, swm_simul, L_pc)
    H%f(i0:iend,j0:jend,:) = H_exact%f(i0:iend,j0:jend,:)

    ! Compute initial mass
    swm_simul%mass0 = mass_computation(H, mesh)

    ! var used in pr mass fixer
    if(swm_simul%mf == 'gpr')then
        swm_simul%a2 = sum(mesh%mt_pc(i0:iend,j0:jend,:)*mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx*mesh%dy
    else if(swm_simul%mf == 'lpr')then
        swm_simul%a2 = sum(mesh%mt_pc(i0  ,j0:jend,:)*mesh%mt_pc(i0  ,j0:jend,:)) + &
        sum(mesh%mt_pc(iend,j0:jend,:)*mesh%mt_pc(iend,j0:jend,:)) + &
        sum(mesh%mt_pc(i0+1:iend-1,j0  ,:)*mesh%mt_pc(i0+1:iend-1,j0  ,:))+&
        sum(mesh%mt_pc(i0+1:iend-1,jend,:)*mesh%mt_pc(i0+1:iend-1,jend,:))
        swm_simul%a2 = swm_simul%a2*mesh%dx*mesh%dy
    else
        if(swm_simul%mf .ne. 'none' .and. swm_simul%mf .ne. 'af')then
            print*, 'ERROR in  swm_ic: invalid mass fixer: ', swm_simul%mf
            stop
        end if
    end if

    ! Define wheter exact solution at all time steps is available or not
    if(swm_simul%ic <= 2)then
        swm_simul%exactsolution = .true.
    else
        swm_simul%exactsolution = .false.
    end if

    ! Filename (for outputs)
    write (swm_simul%ic_name, *) swm_simul%ic
    swm_simul%ic_name = adjustl(swm_simul%ic_name)

    write (swm_simul%id_name, *) swm_simul%id
    swm_simul%id_name = adjustl(swm_simul%id_name)

    write (swm_simul%avd_name, *) swm_simul%avd
    swm_simul%avd_name = adjustl(swm_simul%avd_name)

    swm_simul%name = "ic"//trim(swm_simul%ic_name)//"_"//trim(swm_simul%opsplit) &
    //"_"//trim(swm_simul%recon1d)//"_mt"//trim(swm_simul%mt)//"_"//trim(swm_simul%dp) &
    //"_mf"//trim(swm_simul%mf)//"_et"//trim(swm_simul%et)//"_id"//&
    trim(swm_simul%id_name)//"_av"//trim(swm_simul%avd_name)

    ! basename
    !swm_simul%name = trim(swm_simul%name)

    ! Name scalar fields
    div_ugH_error%name = "swm_"//trim(swm_simul%name)//"_div_error"
    div_ugH%name = "swm_"//trim(swm_simul%name)//"_div"
    rel_vort_error%name = "swm_"//trim(swm_simul%name)//"_rel_vort_error"
    rel_vort%name = "swm_"//trim(swm_simul%name)//"_rel_vort"
 
    H%name = "swm_"//trim(swm_simul%name)//"_H"
    H_error%name = "swm_"//trim(swm_simul%name)//"_H_error"
end subroutine init_swm_vars



subroutine compute_ic_swm(H, V_pu, V_pv, V_pc, mesh, swm_simul, L_pc)
    use swm_vars, only: div_ugH_exact, rel_vort_exact, fcoriolis_pc, &
    abs_vort_exact, abs_vort_flux_exact_pu, abs_vort_flux_exact_pv, abs_vort, &
    H_po_exact, H_pu_exact, H_pv_exact, Ku_po_exact, Kv_po_exact, K_po_exact, &
    wind_po, vcovari_pu_exact, ucovari_pv_exact
    !--------------------------------------------------
    ! Compute the initial conditions for the shallow water
    ! problem on the sphere
    ! This routine also allocate the fields
    !
    ! ic - initial conditions
    !
    ! H -  average values of fluid depth
    ! V_pu - velocity at pu
    ! V_pv - velocity at pv
    !
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(inout) :: swm_simul
    type(lagrange_poly_cs), intent(inout) :: L_pc
    type(scalar_field), intent(inout) :: H
    type(velocity_field), intent(inout) :: V_pu, V_pv, V_pc

    ! aux
    integer(i4) :: i, j, p

    !aux
    real(kind=8) :: lat, lon
    real(kind=8) :: ulon, vlat, ucontra, vcontra, ucovari, vcovari
    real(kind=8) :: rv_pu, rv_pv, av_pu, av_pv, f_pu, f_pv, gradh_lat, gradh_lon

    !debug - check if wind conversion is correct
    real(kind=8) :: ull, vll, error1, error
    error = 0.d0


    ! Scalar field at pc
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                lat  = mesh%pc(i,j,p)%lat
                lon  = mesh%pc(i,j,p)%lon

                ! Fluid depth
                H%f(i,j,p) = h0_swm(lat, lon, swm_simul%ic)

                ! Coriolis force
                fcoriolis_pc%f(i,j,p) = fcoriolis(lat, lon, swm_simul%ic)

                ! Compute velocity
                call velocity_swm(ulon, vlat, lat, lon, 0.d0, swm_simul%ic)

                ! LL2contra
                ucontra = mesh%ll2contra_pc(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pc(i,j,p)%M(1,2)*vlat
                vcontra = mesh%ll2contra_pc(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pc(i,j,p)%M(2,2)*vlat

                ! LL2covari
                ucovari = mesh%ll2covari_pc(i,j,p)%M(1,1)*ulon + mesh%ll2covari_pc(i,j,p)%M(1,2)*vlat
                vcovari = mesh%ll2covari_pc(i,j,p)%M(2,1)*ulon + mesh%ll2covari_pc(i,j,p)%M(2,2)*vlat

                V_pc%ucontra_old%f(i,j,p) = ucontra
                V_pc%vcontra_old%f(i,j,p) = vcontra
                V_pc%ucovari_old%f(i,j,p) = ucovari
                V_pc%vcovari_old%f(i,j,p) = vcovari

            end do
        end do
    end do

    ! Vector field at pu
    do p = 1 , nbfaces
        do i = n0, nend+1
            do j = n0, nend
                lat  = mesh%pu(i,j,p)%lat
                lon  = mesh%pu(i,j,p)%lon

                ! Compute velocity
                call velocity_swm(ulon, vlat, lat, lon, 0.d0, swm_simul%ic)

                ! LL2contra
                ucontra = mesh%ll2contra_pu(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(1,2)*vlat
                vcontra = mesh%ll2contra_pu(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pu(i,j,p)%M(2,2)*vlat


                ! LL2covari
                ucovari = mesh%ll2covari_pu(i,j,p)%M(1,1)*ulon + mesh%ll2covari_pu(i,j,p)%M(1,2)*vlat
                vcovari = mesh%ll2covari_pu(i,j,p)%M(2,1)*ulon + mesh%ll2covari_pu(i,j,p)%M(2,2)*vlat

                V_pu%ucontra_old%f(i,j,p) = ucontra
                V_pu%vcontra_old%f(i,j,p) = vcontra
                V_pu%ucovari_old%f(i,j,p) = ucovari
                V_pu%vcovari_old%f(i,j,p) = vcovari
                vcovari_pu_exact%f(i,j,p) = vcovari
            end do
        end do
    end do

    V_pu%vcovari%f(i0:iend+1,j0:jend,:) = V_pu%vcovari_old%f(i0:iend+1,j0:jend,:)
    V_pu%ucovari%f(i0:iend+1,j0:jend,:) = V_pu%ucovari_old%f(i0:iend+1,j0:jend,:)
 
    ! Vector field at pv
    do p = 1 , nbfaces
        do i = n0, nend
            do j = n0, nend+1
                lat  = mesh%pv(i,j,p)%lat
                lon  = mesh%pv(i,j,p)%lon

                ! Compute velocity
                call velocity_swm(ulon, vlat, lat, lon, 0.d0, swm_simul%ic)

                ! LL2contra
                ucontra = mesh%ll2contra_pv(i,j,p)%M(1,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(1,2)*vlat
                vcontra = mesh%ll2contra_pv(i,j,p)%M(2,1)*ulon + mesh%ll2contra_pv(i,j,p)%M(2,2)*vlat


                ! LL2covari
                ucovari = mesh%ll2covari_pv(i,j,p)%M(1,1)*ulon + mesh%ll2covari_pv(i,j,p)%M(1,2)*vlat
                vcovari = mesh%ll2covari_pv(i,j,p)%M(2,1)*ulon + mesh%ll2covari_pv(i,j,p)%M(2,2)*vlat

                V_pv%ucontra_old%f(i,j,p) = ucontra
                V_pv%vcontra_old%f(i,j,p) = vcontra
                V_pv%ucovari_old%f(i,j,p) = ucovari
                V_pv%vcovari_old%f(i,j,p) = vcovari
                ucovari_pv_exact%f(i,j,p) = ucovari
            end do
        end do
    end do

    V_pv%vcovari%f(i0:iend,j0:jend+1,:) = V_pv%vcovari_old%f(i0:iend,j0:jend+1,:)
    V_pv%ucovari%f(i0:iend,j0:jend+1,:) = V_pv%ucovari_old%f(i0:iend,j0:jend+1,:)

    ! vars to check consistency
    if (swm_simul%ic==0)then
        ! pc fields
        do p = 1, nbfaces
            do i = n0, nend
                do j = n0, nend
                    lat  = mesh%pc(i,j,p)%lat
                    lon  = mesh%pc(i,j,p)%lon

                    ! divergence
                    call div_swm(div_ugh_exact%f(i,j,p), lat, lon, swm_simul%ic)

                    ! relative vorticity
                    call rel_vort_swm(rel_vort_exact%f(i,j,p), lat, lon, swm_simul%ic)

                    ! absolute vorticity
                    abs_vort_exact%f(i,j,p) = rel_vort_exact%f(i,j,p) + fcoriolis_pc%f(i,j,p)
                end do
            end do
        end do
        !abs_vort%f(i0:iend,j0:jend,:) = abs_vort_exact%f(i0:iend,j0:jend,:)

        ! po fields
        do p = 1, nbfaces
            do i = n0, nend+1
                do j = n0, nend+1
                    lat  = mesh%po(i,j,p)%lat
                    lon  = mesh%po(i,j,p)%lon

                    ! depth
                    H_po_exact%f(i,j,p) = h0_swm(lat, lon, swm_simul%ic)
                end do
            end do
        end do
 
        ! Fields at pu
        do p = 1, nbfaces
            do i = n0, nend+1
                do j = n0, nend
                    lat  = mesh%pu(i,j,p)%lat
                    lon  = mesh%pu(i,j,p)%lon

                    ! Fluid depth
                    H_pu_exact%f(i,j,p) = h0_swm(lat, lon, swm_simul%ic)

                    ! relative vorticity at pu
                    call rel_vort_swm(rv_pu, lat, lon, swm_simul%ic)

                    ! coriolis at pu
                    f_pu =  fcoriolis(lat, lon, swm_simul%ic)

                    ! absolute vorticity at pu
                    av_pu = rv_pu + f_pu 

                    ! absolute vorticity flux at pu
                    abs_vort_flux_exact_pu%f(i,j,p) = av_pu*V_pu%ucontra_old%f(i,j,p)
                    abs_vort_flux_exact_pu%f(i,j,p) = abs_vort_flux_exact_pu%f(i,j,p)*mesh%mt_pu(i,j,p)

                end do
            end do
        end do

        ! Fields at pv
        do p = 1, nbfaces
            do i = n0, nend
                do j = n0, nend+1
                    lat  = mesh%pv(i,j,p)%lat
                    lon  = mesh%pv(i,j,p)%lon

                    ! Fluid depth
                    H_pv_exact%f(i,j,p) = h0_swm(lat, lon, swm_simul%ic)

                    ! relative vorticity at pv
                    call rel_vort_swm(rv_pv, lat, lon, swm_simul%ic)

                    ! coriolis at pv
                    f_pv =  fcoriolis(lat, lon, swm_simul%ic)

                    ! absolute vorticity at pv
                    av_pv = rv_pv + f_pv

                    ! absolute vorticity flux at pv
                    abs_vort_flux_exact_pv%f(i,j,p) = av_pv*V_pv%vcontra_old%f(i,j,p)
                    abs_vort_flux_exact_pv%f(i,j,p) = abs_vort_flux_exact_pv%f(i,j,p)*mesh%mt_pv(i,j,p)

                end do
            end do
        end do

        ! Vector field at po
        do p = 1 , nbfaces
            do i = n0, nend+1
                do j = n0, nend+1
                    lat  = mesh%po(i,j,p)%lat
                    lon  = mesh%po(i,j,p)%lon

                    ! Compute velocity
                    call velocity_swm(ulon, vlat, lat, lon, 0.d0, swm_simul%vf)

                    ! LL2contra
                    ucontra = mesh%ll2contra_po(i,j,p)%M(1,1)*ulon + mesh%ll2contra_po(i,j,p)%M(1,2)*vlat
                    vcontra = mesh%ll2contra_po(i,j,p)%M(2,1)*ulon + mesh%ll2contra_po(i,j,p)%M(2,2)*vlat

                    ! LL2covari
                    ucovari = mesh%ll2covari_po(i,j,p)%M(1,1)*ulon + mesh%ll2covari_po(i,j,p)%M(1,2)*vlat
                    vcovari = mesh%ll2covari_po(i,j,p)%M(2,1)*ulon + mesh%ll2covari_po(i,j,p)%M(2,2)*vlat


                    wind_po%ucontra_old%f(i,j,p) = ucontra
                    wind_po%vcontra_old%f(i,j,p) = vcontra

                    wind_po%ucovari_old%f(i,j,p) = ucovari
                    wind_po%vcovari_old%f(i,j,p) = vcovari

                    Ku_po_exact%f(i,j,p) = ucontra*ucovari
                    Kv_po_exact%f(i,j,p) = vcontra*vcovari
                    K_po_exact%f(i,j,p) = (Ku_po_exact%f(i,j,p)+Kv_po_exact%f(i,j,p))*0.5d0
                    !error = max(error, abs(Ku_po%f(i,j,p)-ulon*ulon))
                    !error = max(error, abs(Kv_po%f(i,j,p)-vlat*vlat))
                    !error = max(error, abs(Ku_po%f(i,j,p)+Kv_po%f(i,j,p)-ulon**2-vlat**2))
                end do
            end do
        end do

        !print*, error
        !stop

    end if
   !stop
end subroutine compute_ic_swm


function h0_swm(lat, lon, ic)
    !--------------------------------------------------
    ! Compute the initial condition of the swm
    ! problem on the sphere
    ! 
    ! Possible initial conditions (ic)
    ! 1 - one gaussian hill
    ! 2 - steady state from Will et al 1992 paper 
    !--------------------------------------------------
    integer(i4), intent(in) :: ic
    real(kind=8), intent(in) :: lat, lon
    real(kind=8) :: h0_swm

    ! aux vars
    real(kind=8) :: alpha ! rotation angle
    real(kind=8) :: x, y, z ! r3 coordinates
    real(kind=8) :: x0, y0, z0 ! r3 coordinates of center point
    real(kind=8) :: x1, y1, z1 ! r3 coordinates of center point
    real(kind=8) :: lat0, lon0 ! latlon coordinates of center point
    real(kind=8) :: lat1, lon1 ! latlon coordinates of center point
    real(kind=8) :: b0 ! Gaussian width
    real(kind=8) :: u0, h0
    integer(i4):: m, n

    select case(ic)
        case(1) ! one Gaussian hill
            call sph2cart(lon, lat, x, y, z)
            ! Gaussian center
            lon0 = pi*0.25d0
            lat0 = pi/6.d0
            call sph2cart(lon0, lat0, x0, y0, z0)
            b0 = 10.d0
            h0 = 1000.d0
            h0_swm = h0*0.5d0*(1.d0 + dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2)))

        case(0, 2) ! steady state from will92
            alpha = -45.d0*deg2rad ! Rotation angle
            u0    =  erad*2.d0*pi/12.d0*sec2day ! Wind speed
            h0 = 2.94d0*10000.d0*gravi
            h0_swm = h0 - gravi*(erad*omega*u0 + u0*u0*0.5d0) &
            *(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))**2

        case default
            print*, "ERROR on h0_swm: invalid initial condition."
            stop
    end select
    return
end function h0_swm



subroutine velocity_swm(ulon, vlat, lat, lon, time, vf)
    !--------------------------------------------------
    ! Compute the velocity field of the shallow water
    ! problem on the sphere
    ! 
    ! Possible velocity fields (vf)

    ! 2 - rotated zonal wind
    ! 3 - non divergent
    ! 4 - non divergent
    ! 5 - divergent
    ! 6 - trinometric field
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: lat, lon, time
    real(kind=8), intent(inout) :: ulon, vlat

    ! aux vars
    real(kind=8) :: alpha ! rotation angle
    real(kind=8) :: u0, T, k, lonp
    integer(i4) :: n, m

    select case(vf)
        case(0, 1, 2) ! rotated zonal wind
            alpha = -45.d0*deg2rad ! Rotation angle
            u0    =  erad*2.d0*pi/12.d0*sec2day ! Wind speed
            ulon  =  u0*(dcos(lat)*dcos(alpha) + dsin(lat)*dcos(lon)*dsin(alpha))
            vlat  = -u0*dsin(lon)*dsin(alpha)
            ! divide by earth radius to map the winds to unit sphere
            ulon = ulon/erad
            vlat = vlat/erad
        case default
            print*, "ERROR on velocity_swm: invalid vector field"
            stop
    end select
    return
end subroutine velocity_swm

subroutine div_swm(div, lat, lon, vf)
    !--------------------------------------------------
    ! Compute the exact divergence of the velocity field * depth
    ! 
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: lat, lon
    real(kind=8), intent(out) :: div

    select case(vf)
        case(0, 1, 2)
            div = 0.d0
        case default
            print*, "ERROR on div_swm: invalid vector field"
            stop
    end select
    return
end subroutine div_swm


subroutine rel_vort_swm(rel_vort, lat, lon, vf)
    !--------------------------------------------------
    ! Compute the exact relative vorcity of the velocity fields
    ! 
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: lat, lon
    real(kind=8), intent(out) :: rel_vort
    real(kind=8) :: u0, alpha

    select case(vf)
        case(0, 1, 2)
            alpha = -45.d0*deg2rad ! Rotation angle
            u0    =  erad*2.d0*pi/12.d0*sec2day ! Wind speed
            rel_vort = 2.d0*u0*(-dcos(lat)*dcos(lon)*dsin(alpha) + dsin(lat)*dcos(alpha))
            rel_vort = rel_vort/erad
        case default
            print*, "ERROR on rel_vort_swm: invalid vector field"
            stop
    end select
    return
end subroutine rel_vort_swm

function fcoriolis(lat, lon, ic)
    !--------------------------------------------------
    ! Compute the Coriolis term
    !--------------------------------------------------
    integer(i4), intent(in) :: ic
    real(kind=8), intent(in) :: lat, lon
    real(kind=8) :: fcoriolis

    ! aux vars
    real(kind=8) :: alpha ! rotation angle

    select case(ic)
        case(0, 1, 2)
            alpha = -45.d0*deg2rad ! Rotation angle
            fcoriolis = 2.d0*omega*(-dcos(lat)*dcos(lon)*dsin(alpha) + dsin(lat)*dcos(alpha))

        case(3)
            fcoriolis  =  2.d0*omega*dsin(lat)
        case default
            print*, "ERROR on fcoriolis: invalid initial condition."
            stop
    end select
    return
end function fcoriolis


subroutine gradh_swm(gradh_lon, gradh_lat, lat, lon, vf)
    !--------------------------------------------------
    ! Compute the exact divergence of the depth field
    ! 
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: lat, lon
    real(kind=8), intent(out) :: gradh_lon, gradh_lat
    real(kind=8) :: alpha ! rotation angle
    real(kind=8) :: u0

    select case(vf)
        case(0, 1, 2)
            alpha = -45.d0*deg2rad ! Rotation angle
            u0    =  erad*2.d0*pi/12.d0*sec2day ! Wind speed

            gradh_lon = 2.d0*dcos(lon)*dcos(lat)*dsin(alpha)*dsin(lon)*dcos(lat)*dsin(alpha) + &
            2.d0*dsin(lat)*dcos(alpha)*dsin(lon)*dcos(lat)*dsin(alpha)
            gradh_lon = gradh_lon/dcos(lat)

            gradh_lat = 2.d0*(-dcos(lon)*dcos(lat)*dsin(alpha)+&
                        dsin(lat)*dcos(alpha))*(dcos(lon)*dsin(lat)*dsin(alpha)-&
                        dcos(lat)*dcos(alpha))


            gradh_lat = -gravi*(erad*omega*u0 + u0*u0*0.5d0)*gradh_lat
            gradh_lon = -gravi*(erad*omega*u0 + u0*u0*0.5d0)*gradh_lon
        case default
            print*, "ERROR on gradh_swm: invalid vector field"
            stop
    end select
    return
end subroutine gradh_swm


end module swm_ic
