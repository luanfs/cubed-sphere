module advection_ic
!===============================================================================================
! Module for advection test case set up (initial condition, exact solution and etc)
!
! Luan da Fonseca Santos - October 2022
! (luan.santos@usp.br)
! Test cases are based in the paper "A class of deformational ﬂow test cases for linear
! transport problems  on the sphere", 2010, Ramachandran D. Nair and Peter H. Lauritzen
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
  deg2rad, &
  ll2contra, &
  contra2ll, &
  covari2contra, &
  contra2covari

! Data structures
use datastruct, only: &
  cubedsphere, &
  scalar_field, &
  vector_field, &
  simulation

! Data allocation
use allocation, only: &
  scalar_field_allocation, &
  vector_field_allocation, &
  allocate_adv_vars

! Diagnostics
use diagnostics, only: &
  mass_computation, &
  adv_diagnostics

! duo grid interpolation
use duogrid_interpolation, only: &
    compute_lagrange_cs

implicit none

contains 

function q0_adv(lat, lon, ic)
    !--------------------------------------------------
    ! Compute the initial condition of the advection
    ! problem on the sphere
    ! 
    ! Possible initial conditions (ic)
    ! 1 - constant field
    ! 2 - one gaussian hill
    ! 3 - two gaussian hills
    !--------------------------------------------------
    integer(i4), intent(in) :: ic
    real(kind=8), intent(in) :: lat, lon
    real(kind=8) :: q0_adv

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
            q0_adv = 1.d0

        case(2) ! one Gaussian hill
            call sph2cart(lon, lat, x, y, z)
            ! Gaussian center
            lon0 = pi*0.25d0
            lat0 = pi/6.d0
            call sph2cart(lon0, lat0, x0, y0, z0)
            b0 = 10.d0
            q0_adv = dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))
   
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
            q0_adv = dexp(-b0*((x-x1)**2+ (y-y1)**2 + (z-z1)**2)) + &
                     dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))

        case(4) ! steady state from will92
            alpha = -45.d0*deg2rad ! Rotation angle
            u0 = 2.d0*pi/5.d0     ! Wind speed
            f = (-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))
            q0_adv = 1.0 - f*f

        case(5) ! Trigonometric field
            m = 1
            n = 1
            q0_adv = (-dcos(lon) * dsin(m * lon) * m * dcos(n * lat) ** 4 / dcos(lat) - &
            dsin(lon) * dcos(m * lon) * m ** 2 * dcos(n * lat) ** 4 / dcos(lat) + &
            12.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 2 * dsin(n *lat) ** 2 * n ** 2 * dcos(lat) - &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 4 * n ** 2 * dcos(lat) + &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 3 * dsin(n * lat) * n * dsin(lat)) / dcos(lat)

        case default
            print*, "ERROR on q0_adv: invalid initial condition."
            stop
    end select
    return
end function q0_adv



subroutine velocity_adv(ulon, vlat, lat, lon, time, vf)
    !--------------------------------------------------
    ! Compute the velocity field of the advection
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
            print*, "ERROR on velocity_adv: invalid vector field"
            stop
    end select
    return
end subroutine velocity_adv

subroutine init_adv_vars(mesh)
    use advection_vars
    !--------------------------------------------------
    ! Initialize the advection simulation variables
    ! This routine also allocate the fields
    !
    ! ic - initial conditions
    ! vf - velocity field
    !
    ! Q - scalar field average values on cells
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
    px%recon = advsimul%recon1d
    py%recon = advsimul%recon1d

    ! Metric tensor formulation
    px%mt = advsimul%mt
    py%mt = advsimul%mt
    px%et = advsimul%et
    py%et = advsimul%et

    ! N
    px%n = mesh%n
    py%n = mesh%n

    ! Lagrange polynomial at centers
    L_pc%degree =  advsimul%id
    L_pc%order =  L_pc%degree+1
    L_pc%pos = 1

    ! Allocate the variables
    call allocate_adv_vars(mesh)

    ! Compute lagrange polynomials
    call compute_lagrange_cs(L_pc, mesh)
    advsimul%id_d2a = 3

    ! Time step over 2
    advsimul%dto2 = advsimul%dt*0.5d0

    ! Final time step - we adopt 5 units for all simulations
    advsimul%tf = 5.d0

    ! Number of time steps
    advsimul%nsteps = int(advsimul%tf/advsimul%dt)
    advsimul%nplot = advsimul%nsteps/(advsimul%nplot-1)
    advsimul%plotcounter = 0

    ! Compute the initial conditions
    call compute_ic_adv(Q_exact, wind_pu, wind_pv, wind_pc, mesh, advsimul)
    Q%f(i0:iend,j0:jend,:) = Q_exact%f(i0:iend,j0:jend,:)

    ! Compute initial mass
    advsimul%mass0 = mass_computation(Q, mesh)

    ! var used in pr mass fixer
    if(advsimul%mf == 'gpr')then
        advsimul%a2 = sum(mesh%mt_pc(i0:iend,j0:jend,:)*mesh%mt_pc(i0:iend,j0:jend,:))*mesh%dx*mesh%dy
    else if(advsimul%mf == 'lpr')then
        advsimul%a2 = sum(mesh%mt_pc(i0  ,j0:jend,:)*mesh%mt_pc(i0  ,j0:jend,:)) + &
        sum(mesh%mt_pc(iend,j0:jend,:)*mesh%mt_pc(iend,j0:jend,:)) + &
        sum(mesh%mt_pc(i0+1:iend-1,j0  ,:)*mesh%mt_pc(i0+1:iend-1,j0  ,:))+&
        sum(mesh%mt_pc(i0+1:iend-1,jend,:)*mesh%mt_pc(i0+1:iend-1,jend,:))
        advsimul%a2 = advsimul%a2*mesh%dx*mesh%dy
    else
        if(advsimul%mf .ne. 'none' .and. advsimul%mf .ne. 'af')then
            print*, 'ERROR in  advection_ic: invalid mass fixer: ', advsimul%mf
            stop
        end if
    end if

    ! Define wheter exact solution at all time steps is available or not
    if(advsimul%ic == 1 .and. advsimul%vf <= 2)then
        advsimul%exactsolution = .true.
    else if(advsimul%ic .ne. 3 .and. advsimul%vf == 1)then
        advsimul%exactsolution = .true.
     else
        advsimul%exactsolution = .false.
    end if

    ! Filename (for outputs)
    write (advsimul%ic_name, *) advsimul%ic
    advsimul%ic_name = adjustl(advsimul%ic_name)

    write (advsimul%vf_name, *) advsimul%vf
    advsimul%vf_name = adjustl(advsimul%vf_name)

    write (advsimul%id_name, *) advsimul%id
    advsimul%id_name = adjustl(advsimul%id_name)

    advsimul%name = "ic"//trim(advsimul%ic_name)//"_vf"//trim(advsimul%vf_name)//"_"//trim(advsimul%opsplit) &
    //"_"//trim(advsimul%recon1d)//"_mt"//trim(advsimul%mt)//"_"//trim(advsimul%dp) &
    //"_mf"//trim(advsimul%mf)//"_et"//trim(advsimul%et)//"_id"//trim(advsimul%id_name)

    ! basename
    !advsimul%name = trim(advsimul%name)

    ! Name scalar fields
    div_ugq_error%name = "div_"//trim(advsimul%name)//"_div_error"
    div_ugq%name = "div_"//trim(advsimul%name)//"_div"
    Q%name = "adv_"//trim(advsimul%name)//"_Q"
    Q_error%name = "adv_"//trim(advsimul%name)//"_Q_error"

end subroutine init_adv_vars

subroutine compute_ic_adv(Q, V_pu, V_pv, V_pc, mesh, advsimul)
    !--------------------------------------------------
    ! Compute the initial conditions for the advection
    ! problem on the sphere
    ! This routine also allocate the fields
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
    type(simulation), intent(inout) :: advsimul
    type(scalar_field), intent(inout) :: Q
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
                Q%f(i,j,p) = q0_adv(lat, lon, advsimul%ic)

                call velocity_adv(ulon, vlat, lat, lon, 0.d0, advsimul%vf)
                call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pc(i,j,p)%M)
                call contra2covari(ucovari, vcovari, ucontra, vcontra, mesh%contra2covari_pc(i,j,p)%M)
                V_pc%ucontra_old%f(i,j,p) = ucovari
                V_pc%vcontra_old%f(i,j,p) = vcovari

            end do
        end do
    end do
    ! Vector field at pu
    do p = 1 , nbfaces
        do i = n0, nend+1
            do j = n0, nend
                lat  = mesh%pu(i,j,p)%lat
                lon  = mesh%pu(i,j,p)%lon

                call velocity_adv(ulon, vlat, lat, lon, 0.d0, advsimul%vf)
                call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pu(i,j,p)%M)

                V_pu%ucontra_old%f(i,j,p) = ucontra
                V_pu%vcontra_old%f(i,j,p) = vcontra

                call contra2covari(ucovari, vcovari, ucontra, vcontra, mesh%contra2covari_pu(i,j,p)%M)
                V_pu%ucovari%f(i,j,p) = ucovari
                V_pu%vcovari%f(i,j,p) = vcovari
                ! debug 
                !error1 = abs( (ulon**2+vlat**2) -(ucontra*ucovari + vcontra*vcovari) ) 
                !call contra2ll(ull, vll, ucontra, vcontra, mesh%contra2ll_pu(i,j,p)%M)
                !error1 =  abs(ulon-ull)
                !error1 = max(error1, abs(vll-vlat))
                !error = max(error, error1)
            end do
        end do
    end do

    V_pu%ucontra%f(i0:iend+1,j0:jend,:) = V_pu%ucontra_old%f(i0:iend+1,j0:jend,:)
    V_pu%vcontra%f(i0:iend+1,j0:jend,:) = V_pu%vcontra_old%f(i0:iend+1,j0:jend,:)
 
    ! Vector field at pv
    do p = 1 , nbfaces
        do i = n0, nend
            do j = n0, nend+1
                lat  = mesh%pv(i,j,p)%lat
                lon  = mesh%pv(i,j,p)%lon

                call velocity_adv(ulon, vlat, lat, lon, 0.d0, advsimul%vf)
                call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pv(i,j,p)%M)

                V_pv%ucontra_old%f(i,j,p) = ucontra
                V_pv%vcontra_old%f(i,j,p) = vcontra

                call contra2covari(ucovari, vcovari, ucontra, vcontra, mesh%contra2covari_pv(i,j,p)%M)
                V_pv%ucovari%f(i,j,p) = ucovari
                V_pv%vcovari%f(i,j,p) = vcovari


                ! debug 
                !error1 = abs( (ulon**2+vlat**2) -(ucontra*ucovari + vcontra*vcovari) ) 
                !call contra2ll(ull, vll, ucontra, vcontra, mesh%contra2ll_pv(i,j,p)%M)
                !error1 =  abs(ulon-ull)
                !error1 = max(error1, abs(vll-vlat))
                !error = max(error, error1)
            end do
        end do
    end do

    V_pv%ucontra%f(i0:iend,j0:jend+1,:) = V_pv%ucontra_old%f(i0:iend,j0:jend+1,:)
    V_pv%vcontra%f(i0:iend,j0:jend+1,:) = V_pv%vcontra_old%f(i0:iend,j0:jend+1,:)
 
    ! CFL number
    advsimul%cfl = maxval(abs(V_pu%ucontra%f))
    advsimul%cfl = max(advsimul%cfl, maxval(abs(V_pv%vcontra%f)))
    advsimul%cfl = advsimul%cfl*advsimul%dt/mesh%dx
    !print*,error
    !stop
end subroutine compute_ic_adv

subroutine div_adv(div, lat, lon, vf)
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
            print*, 'error on div_adv: div is not implemented for this vector field: ', vf
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
            print*, "ERROR on div_adv: invalid vector field"
            stop
      end select
    return
end subroutine div_adv

 
subroutine compute_exact_div(div, mesh, advsimul)
    !--------------------------------------------------
    ! Compute the exact divergence at cell centers
    !--------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    type(simulation), intent(in) :: advsimul
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
                call div_adv(div%f(i,j,p), lat, lon, advsimul%vf)
            end do
        end do
    end do
end subroutine compute_exact_div

end module advection_ic
