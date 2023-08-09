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
    n0, nend

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  deg2rad

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field, &
  simulation, &
  lagrange_poly_cs

! Data allocation
use allocation, only: &
  scalar_field_allocation, &
  vector_field_allocation, &
  allocate_swm_vars

use diagnostics, only: &
  mass_computation

! duo grid interpolation
use duogrid_interpolation, only: &
    compute_lagrange_cs, &
    interp_D2Aduogrid, &
    interp_C2Agrid, &
    interp_A2Cgrid


implicit none

contains 

function h0_swm(lat, lon, ic)
    !--------------------------------------------------
    ! Compute the initial condition of the swm
    ! problem on the sphere
    ! 
    ! Possible initial conditions (ic)
    ! 1 - constant field
    ! 2 - one gaussian hill
    ! 3 - two gaussian hills
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
    real(kind=8) :: u0, f
    integer(i4):: m, n

    select case(ic)
        case(1) ! constant scalar field
            h0_swm = 1.d0

        case(2) ! one Gaussian hill
            call sph2cart(lon, lat, x, y, z)
            ! Gaussian center
            lon0 = pi*0.25d0
            lat0 = pi/6.d0
            call sph2cart(lon0, lat0, x0, y0, z0)
            b0 = 10.d0
            h0_swm = dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))
   
        case(3) ! two Gaussian hills
            call sph2cart(lon, lat, x, y, z)
            ! Gaussian hill centers
            lon0 = -pi/6.d0
            lat0 = 0.d0
            lon1 = pi/6.d0
            lat1 = 0.d0
            call sph2cart(lon0, lat0, x0, y0, z0)
            call sph2cart(lon1, lat1, x1, y1, z1)
            b0 = 5.0
            h0_swm = dexp(-b0*((x-x1)**2+ (y-y1)**2 + (z-z1)**2)) + &
                     dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))

        case(4) ! steady state from will92
            alpha = -45.d0*deg2rad ! Rotation angle
            u0 = 2.d0*pi/5.d0     ! Wind speed
            f = (-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))
            h0_swm = 1.0 - f*f

        case(5) ! Trigonometric field
            m = 1
            n = 1
            h0_swm = (-dcos(lon) * dsin(m * lon) * m * dcos(n * lat) ** 4 / dcos(lat) - &
            dsin(lon) * dcos(m * lon) * m ** 2 * dcos(n * lat) ** 4 / dcos(lat) + &
            12.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 2 * dsin(n *lat) ** 2 * n ** 2 * dcos(lat) - &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 4 * n ** 2 * dcos(lat) + &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 3 * dsin(n * lat) * n * dsin(lat)) / dcos(lat)

        case default
            print*, "ERROR on q0_swm: invalid initial condition."
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
        case(1)! rotated zonal wind
            alpha = -45.d0*deg2rad ! Rotation angle
            u0    =  2.d0*pi/5.d0 ! Wind speed
            ulon  =  u0*(dcos(lat)*dcos(alpha) + dsin(lat)*dcos(lon)*dsin(alpha))
            vlat  = -u0*dsin(lon)*dsin(alpha)

        case(2) ! Non divergent field 4 from Nair and Lauritzen 2010
            T = 5.d0 ! Period
            k = 2.d0
            lonp = lon-2.d0*pi*time/T
            ulon = k*(dsin((lonp+pi))**2)*(dsin(2.*lat))*(dcos(pi*time/T))+2.*pi*dcos(lat)/T
            vlat = k*(dsin(2*(lonp+pi)))*(dcos(lat))*(dcos(pi*time/T))

        case(3)! Divergent field 3 from Nair and Lauritzen 2010
            T = 5.d0 ! Period
            k = 1.d0
            ulon = -k*(dsin((lon+pi)/2.d0)**2)*(dsin(2.d0*lat))*(dcos(lat)**2)*(dcos(pi*time/T))
            vlat = (k/2.d0)*(dsin((lon+pi)))*(dcos(lat)**3)*(dcos(pi*time/T))

        case(4) ! trigonometric field
            m = 1
            n = 1
            ulon = m*(dsin(lon)*dsin(m*lon)*dcos(n*lat)**3)!/np.cos(lat)
            vlat = 4*n*(dcos(n*lat)**3)*dsin(n*lat)*dcos(m*lon)*dsin(lon)
        case default
            print*, "ERROR on velocity_swm: invalid vector field"
            stop
    end select
    return
end subroutine velocity_swm

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

    ! Reconstruction scheme
    px%recon = swm_simul%recon1d
    py%recon = swm_simul%recon1d

    ! Metric tensor formulation
    px%mt = swm_simul%mt
    py%mt = swm_simul%mt
    px%et = swm_simul%et
    py%et = swm_simul%et

    ! N
    px%n = mesh%n
    py%n = mesh%n

    ! Lagrange polynomial at centers
    L_pc%degree =  swm_simul%id
    L_pc%order =  L_pc%degree+1
    L_pc%pos = 1

    ! Allocate the variables
    call allocate_swm_vars(mesh)

    ! Compute lagrange polynomials
    call compute_lagrange_cs(L_pc, mesh)
    swm_simul%id_d2a = 3

    ! Time step over 2
    swm_simul%dto2 = swm_simul%dt*0.5d0

    ! Final time step - we adopt 5 units for all simulations
    swm_simul%tf = 5.d0

    ! Number of time steps
    swm_simul%nsteps = int(swm_simul%tf/swm_simul%dt)
    swm_simul%nplot = swm_simul%nsteps/(swm_simul%nplot-1)
    swm_simul%plotcounter = 0

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
    if(swm_simul%ic == 1 .and. swm_simul%vf <= 2)then
        swm_simul%exactsolution = .true.
    else if(swm_simul%ic .ne. 3 .and. swm_simul%vf == 1)then
        swm_simul%exactsolution = .true.
     else
        swm_simul%exactsolution = .false.
    end if

    ! Filename (for outputs)
    write (swm_simul%ic_name, *) swm_simul%ic
    swm_simul%ic_name = adjustl(swm_simul%ic_name)

    write (swm_simul%id_name, *) swm_simul%id
    swm_simul%id_name = adjustl(swm_simul%id_name)

    swm_simul%name = "ic"//trim(swm_simul%ic_name)//"_"//trim(swm_simul%opsplit) &
    //"_"//trim(swm_simul%recon1d)//"_mt"//trim(swm_simul%mt)//"_"//trim(swm_simul%dp) &
    //"_mf"//trim(swm_simul%mf)//"_et"//trim(swm_simul%et)//"_id"//trim(swm_simul%id_name)

    ! basename
    !swm_simul%name = trim(swm_simul%name)

    ! Name scalar fields
    div_ugH_error%name = "swm_"//trim(swm_simul%name)//"_div_error"
    div_ugH%name = "swm_"//trim(swm_simul%name)//"_div"
    H%name = "swm_"//trim(swm_simul%name)//"_Q"
    H_error%name = "swm_"//trim(swm_simul%name)//"_Q_error"

end subroutine init_swm_vars

subroutine compute_ic_swm(H, V_pu, V_pv, V_pc, mesh, swm_simul, L_pc)
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
    type(vector_field), intent(inout) :: V_pu, V_pv, V_pc

    ! aux
    integer(i4) :: i, j, p

    !aux
    real(kind=8) :: lat, lon
    real(kind=8) :: ulon, vlat, ucontra, vcontra, ucovari, vcovari
    
    !debug - check if wind conversion is correct
    real(kind=8) :: ull, vll, error1, error
    error = 0.d0


    ! Scalar field at pc
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                lat  = mesh%pc(i,j,p)%lat
                lon  = mesh%pc(i,j,p)%lon
                H%f(i,j,p) = h0_swm(lat, lon, swm_simul%ic)

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
            end do
        end do
    end do

    V_pv%vcovari%f(i0:iend,j0:jend+1,:) = V_pv%vcovari_old%f(i0:iend,j0:jend+1,:)
    V_pv%ucovari%f(i0:iend,j0:jend+1,:) = V_pv%ucovari_old%f(i0:iend,j0:jend+1,:)

    !-----------------------------------------------------------------------------
    ! Duogrid interpolation of the vector field on a D grid
    ! first we interpolate to the A grid (including A grid ghost cells)
    call interp_D2Aduogrid(V_pu, V_pv, V_pc, L_pc, mesh)

    ! then we interpolate from D grid to the A grid inner cells
    call interp_C2Agrid(V_pu%vcovari%f, V_pv%ucovari%f, V_pc%vcovari%f, V_pc%ucovari%f,  swm_simul%id_d2a)

    ! now we fill the ghost cell D grid
    !call interp_A2Dduogrid(V_pu, V_pv, V_pc)

    ! then we interpolate from A grid to the D grid inner cells
    call interp_A2Cgrid(V_pu%ucovari%f, V_pv%vcovari%f, V_pc%ucovari%f, V_pc%vcovari%f, swm_simul%id_d2a)



    ! Convert from covariant to contravariant
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(n0, nend) &
    !$OMP SHARED(V_pc, V_pu, V_pv, mesh)
    V_pu%ucontra%f(n0:nend+1,n0:nend,:) = &
    V_pu%ucovari%f(n0:nend+1,n0:nend,:)*mesh%covari2contra_pu(n0:nend+1,n0:nend,:)%M(1,1)+&
    V_pu%vcovari%f(n0:nend+1,n0:nend,:)*mesh%covari2contra_pu(n0:nend+1,n0:nend,:)%M(1,2) 

    V_pv%vcontra%f(n0:nend,n0:nend+1,:) =&
    V_pv%ucovari%f(n0:nend,n0:nend+1,:)*mesh%covari2contra_pv(n0:nend,n0:nend+1,:)%M(2,1)+&
    V_pv%vcovari%f(n0:nend,n0:nend+1,:)*mesh%covari2contra_pv(n0:nend,n0:nend+1,:)%M(2,2) 
    !$OMP END PARALLEL WORKSHARE

    print*, maxval(abs(V_pc%ucovari%f(:,:,:)-V_pc%ucovari_old%f(:,:,:)))
    print*, maxval(abs(V_pc%vcovari%f(:,:,:)-V_pc%vcovari_old%f(:,:,:)))

    print*, maxval(abs(V_pu%ucovari%f(i0-1:iend+2,:,:)-V_pu%ucovari_old%f(i0-1:iend+2,:,:)))
    print*, maxval(abs(V_pu%vcovari%f(i0-1:iend+2,:,:)-V_pu%vcovari_old%f(i0-1:iend+2,:,:)))

    print*, maxval(abs(V_pv%ucovari%f(:,j0-1:jend+2,:)-V_pv%ucovari_old%f(:,j0-1:jend+2,:)))
    print*, maxval(abs(V_pv%vcovari%f(:,j0-1:jend+2,:)-V_pv%vcovari_old%f(:,j0-1:jend+2,:)))
    print*, maxval(abs(V_pu%ucontra%f(i0-1:iend+2,:,:)-V_pu%ucontra_old%f(i0-1:iend+2,:,:)))
    print*, maxval(abs(V_pv%vcontra%f(:,j0-1:jend+2,:)-V_pv%vcontra_old%f(:,j0-1:jend+2,:)))



    ! CFL number
    swm_simul%cfl = maxval(abs(V_pu%ucontra%f))
    swm_simul%cfl = max(swm_simul%cfl, maxval(abs(V_pv%vcontra%f)))
    swm_simul%cfl = swm_simul%cfl*swm_simul%dt/mesh%dx
    !print*,error
    !stop
end subroutine compute_ic_swm

subroutine div_swm(div, lat, lon, vf)
    !--------------------------------------------------
    ! Compute the exact divergence of the velocity fields
    ! 
    ! Possible velocity fields (vf)
    ! 1 - zonal wind
    ! 2 - rotated zonal wind
    ! 3 - non divergent
    ! 4 - non divergent
    ! 5 - divergent
    ! 6 - trigometric field
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(kind=8), intent(in) :: lat, lon
    real(kind=8), intent(out) :: div
    integer(i4) :: m, n

    select case(vf)
        case(1, 2)
            div = 0.d0
        case(3)
            print*, 'error on div_swm: div is not implemented for this vector field: ', vf
            stop
        case(4)
            m = 1
            n = 1
            div = -(-dcos(lon) * dsin(m * lon) * m * dcos(n * lat) ** 4 / dcos(lat) - &
            dsin(lon) * dcos(m * lon) * m ** 2 * dcos(n * lat) ** 4 / dcos(lat) + &
            12.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 2 * dsin(n *lat) ** 2 * n ** 2 * dcos(lat) - &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 4 * n ** 2 * dcos(lat) + &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 3 * dsin(n * lat) * n * dsin(lat)) / dcos(lat)
        case default
            print*, "ERROR on div_swm: invalid vector field"
            stop
      end select
    return
end subroutine div_swm

 
subroutine swm_compute_exact_div(div, mesh, swm_simul)
    !--------------------------------------------------
    ! Compute the exact divergence at cell centers
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(simulation), intent(in) :: swm_simul
    type(scalar_field), intent(inout) :: div

    ! aux
    integer(i4) :: i0, iend, j0, jend, nt
    integer(i4) :: i, j, p

    !aux
    real(kind=8) :: lat, lon

    ! interior grid indexes
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    ! Scalar field at pc
    do p = 1, nbfaces
        do i = i0, iend
            do j = j0, jend
                lat  = mesh%pc(i,j,p)%lat
                lon  = mesh%pc(i,j,p)%lon
                call div_swm(div%f(i,j,p), lat, lon, swm_simul%vf)
            end do
        end do
    end do
end subroutine swm_compute_exact_div

end module swm_ic
