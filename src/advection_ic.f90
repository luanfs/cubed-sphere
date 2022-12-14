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
    real(r8), intent(in) :: lat, lon
    real(r8) :: q0_adv

    ! aux vars
    real(r8) :: alpha ! rotation angle
    real(r8) :: x, y, z ! r3 coordinates
    real(r8) :: x0, y0, z0 ! r3 coordinates of center point
    real(r8) :: x1, y1, z1 ! r3 coordinates of center point
    real(r8) :: lat0, lon0 ! latlon coordinates of center point
    real(r8) :: lat1, lon1 ! latlon coordinates of center point
    real(r8) :: b0 ! Gaussian width 

    select case(ic)

      case(1) ! constant scalar field
        q0_adv = 1._r8

      case(2) ! one Gaussian hill
          call sph2cart(lon, lat, x, y, z)
          ! Gaussian center
          lon0 = 0._r8
          lat0 = 0._r8
          call sph2cart(lon0, lat0, x0, y0, z0)
          b0 = 5._r8
          q0_adv = dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))
   
      case(3) ! two Gaussian hills
        call sph2cart(lon, lat, x, y, z)
        ! Gaussian hill centers
        lon0 = -pi/6._r8
        lat0 = 0._r8
        lon1 = pi/6._r8
        lat1 = 0._r8
        call sph2cart(lon0, lat0, x0, y0, z0)
        call sph2cart(lon1, lat1, x1, y1, z1)
        b0 = 5.0
        q0_adv = dexp(-b0*((x-x1)**2+ (y-y1)**2 + (z-z1)**2)) + &
                 dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))

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
    ! 1 - zonal wind
    ! 2 - rotated zonal wind
    ! 3 - non divergent
    ! 4 - non divergent
    ! 5 - divergent
    ! 6 - trinometric field
    !--------------------------------------------------
    integer(i4), intent(in) :: vf
    real(r8), intent(in) :: lat, lon, time
    real(r8), intent(inout) :: ulon, vlat

    ! aux vars
    real(r8) :: alpha ! rotation angle
    real(r8) :: u0, T, k, lonp
    integer(i4) :: n, m

    select case(vf)
      case(1, 2)! zonal wind
        if(vf == 1) then
           alpha = 0._r8*deg2rad   ! Rotation angle
        else 
           alpha = -45._r8*deg2rad ! Rotation angle
        end if
        u0    =  2.0*pi/5.0 ! Wind speed
        ulon  =  u0*(dcos(lat)*dcos(alpha) + dsin(lat)*dcos(lon)*dsin(alpha))
        vlat  = -u0*dsin(lon)*dsin(alpha)

      case(3) !  Non divergent field 2 from Nair and Lauritzen 2010
        T = 5._r8 ! Period
        k = 2._r8
        ulon = k*dsin(lon+pi)**2 * dsin(2.0*lat) * dcos(pi*time/T)
        vlat = k*dsin(2*(lon+pi)) * dcos(lat) * dcos(pi*time/T)

      case(4) ! Non divergent field 4 from Nair and Lauritzen 2010
        T = 5.0 ! Period
        k = 2.0
        lonp = lon-2*pi*time/T
        ulon = k*(dsin((lonp+pi))**2)*(dsin(2.*lat))*(dcos(pi*time/T))+2.*pi*dcos(lat)/T
        vlat = k*(dsin(2*(lonp+pi)))*(dcos(lat))*(dcos(pi*time/T))

      case(5)! Divergent field 3 from Nair and Lauritzen 2010
        T = 5.0 ! Period
        k = 1.0
        ulon = -k*(dsin((lon+pi)/2.0)**2)*(dsin(2.0*lat))*(dcos(lat)**2)*(dcos(pi*time/T))
        vlat = (k/2.0)*(dsin((lon+pi)))*(dcos(lat)**3)*(dcos(pi*time/T))

      case(6) ! trigonometric field
        m = 1
        n = 1
        ulon = -m*(dsin(lon)*dsin(m*lon)*dcos(n*lat)**3)!/np.cos(lat)
        vlat = -4*n*(dcos(n*lat)**3)*dsin(n*lat)*dcos(m*lon)*dsin(lon)
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
    integer(i4) :: i0, iend, j0, jend
    integer(i4) :: i, j, p

    ! ppm direction
    px%dir = 1 ! x direction
    py%dir = 2 ! y direction

    ! Reconstruction scheme
    px%recon = advsimul%recon1d
    py%recon = advsimul%recon1d

    ! N
    px%n = mesh%n
    py%n = mesh%n

    ! Allocate the variables
    call allocate_adv_vars(mesh)

    ! Compute the initial conditions
    call compute_ic_adv(Q, wind_pu, wind_pv, mesh, advsimul)

    ! Time step
    advsimul%dt = 0.5_r8*mesh%dx/maxval(abs(wind_pu%ucontra%f))

    ! Compute initial mass
    advsimul%mass0 = mass_computation(Q, mesh)

   ! Define wheter exact solution is available or not
    if(advsimul%vf == 1 .or. advsimul%vf == 2)then
      advsimul%exactsolution = .true.
    else if(advsimul%ic == 1 .and. advsimul%vf <= 4)then
      advsimul%exactsolution = .true.
    else
      advsimul%exactsolution = .false.
    end if

    ! Filename (for outputs)
    write (advsimul%ic_name, *) advsimul%ic
    advsimul%ic_name = adjustl(advsimul%ic_name)

    write (advsimul%vf_name, *) advsimul%vf
    advsimul%vf_name = adjustl(advsimul%vf_name)

    advsimul%name = "adv_"//"ic"//trim(advsimul%ic_name)//"_vf"//trim(advsimul%vf_name)//"_"//trim(advsimul%recon1d)
    !print*, trim(advsimul%name)
  end subroutine init_adv_vars

  subroutine compute_ic_adv(Q, V_pu, V_pv, mesh, advsimul)
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
    type(simulation), intent(in) :: advsimul
    type(scalar_field), intent(inout) :: Q
    type(vector_field), intent(inout) :: V_pu, V_pv

    ! aux
    integer(i4) :: i0, iend, j0, jend, nt
    integer(i4) :: i, j, p, N

    !aux
    real(r8) :: lat, lon
    real(r8) :: ulon, vlat, ucontra, vcontra
    
    !debug - check if wind conversion is correct
    real(r8) :: ull, vll, error1, error
    error = 0._r8

    ! interior grid indexes
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend
    nt = mesh%ntotal
    N = mesh%n

    ! Scalar field at pc
    do p = 1, nbfaces
      do i = i0, iend
        do j = j0, jend
          lat  = mesh%pc(i,j,p)%lat
          lon  = mesh%pc(i,j,p)%lon
          Q%f(i,j,p) = q0_adv(lat, lon, advsimul%ic)
        end do
      end do
    end do

    ! Vector field at pu
    do p = 1 , nbfaces
      do i = i0-1, iend
        do j = 1, nt
          lat  = mesh%pu(i,j,p)%lat
          lon  = mesh%pu(i,j,p)%lon

          call  velocity_adv(ulon, vlat, lat, lon, 0._r8, advsimul%vf)
          call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pu(i,j,p)%M)
          V_pu%ucontra%f(i,j,p) = ucontra
          V_pu%vcontra%f(i,j,p) = vcontra

          ! debug 
          !call contra2ll(ull, vll, ucontra, vcontra, mesh%contra2ll_pu(i,j,p)%M)
          !error1 =  abs(ulon-ull)
          !error1 = max(error1, abs(vll-vlat))
          !error = max(error, error1)
        end do
      end do
    end do    

    ! Vector field at pv
    do p = 1 , nbfaces
      do i = 1, nt
        do j = j0-1, jend
          lat  = mesh%pv(i,j,p)%lat
          lon  = mesh%pv(i,j,p)%lon

          call  velocity_adv(ulon, vlat, lat, lon, 0._r8, advsimul%vf)
          call ll2contra(ulon, vlat, ucontra, vcontra, mesh%ll2contra_pv(i,j,p)%M)

          V_pv%ucontra%f(i,j,p) = ucontra
          V_pv%vcontra%f(i,j,p) = vcontra

          ! debug 
          !call contra2ll(ull, vll, ucontra, vcontra, mesh%contra2ll_pv(i,j,p)%M)
          !error1 =  abs(ulon-ull)
          !error1 = max(error1, abs(vll-vlat))
          !error = max(error, error1)
        end do
      end do
    end do
 
    !print*, N
    !print*, V_pu%ucontra%f(3:N+3,6,6)
    !print*, shape(V_pu%ucontra%f)
    !stop   !print*, error
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
    real(r8), intent(in) :: lat, lon
    real(r8), intent(out) :: div
    integer(i4) :: m, n

    select case(vf)
      case(1, 2, 3, 4)
        div = 0._r8
      case(5)
        print*, 'error on div_adv: div is not implemented for this vector field: ', vf
        stop
      case(6)
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
    real(r8) :: lat, lon

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
