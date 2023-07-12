module sphgeo
!========================================================================
    !
! This module contains all the spherical geometry routines 
!
! Based on module 'smeshpack' from iModel (https://github.com/pedrospeixoto/iModel)
!
! Luan Santos
!========================================================================

!Global constants
use constants, only: &
     deg2rad, &
     eps, &
     eps2, &
     i4, &
     pardir, &
     pi, &
     pi2, &
     pio2, &
     r8, &
     r4, &
     r16, &
     rad2deg, &
     acube, &
     nbfaces

!Data structures
use datastruct, only: &
  cubedsphere

!Numerical linear algebra
use linear_algebra, only: &
  norm, &
  cross_product, &
  det

implicit none

contains 

!---------------------------------------------------------------------
! panel indexes distribution
!      +---+
!      | 5 |
!  +---+---+---+---+
!  | 4 | 1 | 2 | 3 |
!  +---+---+---+---+
!      | 6 |
!      +---+
!---------------------------------------------------------------------

subroutine equidistant_gnomonic_map(x, y, p, panel)
!---------------------------------------------------------------
! this routine computes the gnomonic mapping based on the equidistant projection
! defined by rancic et al (96) for each panel
! - x, y are the variables defined in [-a, a].
! - the projection is applied on the points (x,y)
! - returns the cartesian coordinates of the
! projected point p.
!
! references: 
! - rancic, m., purser, r.j. and mesinger, f. (1996), a global shallow-water model using an expanded
!  spherical cube: gnomonic versus conformal coordinates. q.j.r. meteorol. soc., 122: 959-982. 
!  https://doi.org/10.1002/qj.49712253209
! - nair, r. d., thomas, s. j., & loft, r. d. (2005). a discontinuous galerkin transport scheme on the
! cubed sphere, monthly weather review, 133(4), 814-828. retrieved feb 7, 2022, 
! from https://journals.ametsoc.org/view/journals/mwre/133/4/mwr2890.1.xml
!
!---------------------------------------------------------------
    real(r8), intent(in) :: x ! local coordinates
    real(r8), intent(in) :: y ! local coordinates
    real(r8), intent(out) :: p(1:3) ! projected point

    ! aux vars
    real(r8) :: r, r2, invr,  xor, yor, aor

    ! panel
    integer(i4), intent(in) :: panel

    ! auxiliary variables
    r2   = acube**2 + x**2 + y**2
    r    = dsqrt(r2)
    invr = 1._r8/r
    xor  = invr*x
    yor  = invr*y
    aor  = invr*acube

    ! compute the cartesian coordinates for each panel
    ! with the aid of the auxiliary variables  
    select case(panel)
        case(1)
            p(1) =  aor
            p(2) =  xor
            p(3) =  yor

        case(2)
            p(1) = -xor
            p(2) =  aor
            p(3) =  yor

        case(3)
            p(1) = -aor
            p(2) = -xor
            p(3) =  yor

        case(4)
            p(1) =  xor
            p(2) = -aor
            p(3) =  yor      

        case(5)
            p(1) = -yor
            p(2) =  xor
            p(3) =  aor

        case(6)
            p(1) =  yor
            p(2) =  xor
            p(3) = -aor

        case default
            print*, 'error on equidistant_gnomonic_map: invalid panel, ', panel
    end select
    return
end subroutine equidistant_gnomonic_map

subroutine derivative_xdir_equidistant_gnomonic_map(x, y, p, panel)
    !---------------------------------------------------------------
    ! this routine computes the derivative in x direction of the
    ! gnomonic mapping based on the equidistant projection
    !---------------------------------------------------------------
    real(r8), intent(in) :: x ! local coordinates
    real(r8), intent(in) :: y ! local coordinates
    real(r8), intent(out) :: p(1:3) ! projected point

    ! aux vars
    real(r8) :: r2, r32, invr, ax, xy, y2, a, a2

    ! panel
    integer(i4), intent(in) :: panel

    a = acube

    !auxiliary variables
    r2   = a**2 + x**2 + y**2
    r32  = dsqrt(r2)**3
    invr = 1._r8/r32
    ax = a*x
    xy = x*y
    y2 = y*y
    a2 = a*a

    ! compute the cartesian coordinates for each panel
    ! with the aid of the auxiliary variables  
    select case(panel)
        case(1)
            p(1) = -ax*invr
            p(2) =  (a2 + y2)*invr
            p(3) = -xy*invr

        case(2)
            p(1) = -(a2 + y2)*invr
            p(2) = -ax*invr
            p(3) = -xy*invr

        case(3)
            p(1) =  ax*invr
            p(2) = -(a2 + y2)*invr
            p(3) = -xy*invr

        case(4)
            p(1) =  (a2 + y2)*invr
            p(2) =  ax*invr
            p(3) = -xy*invr 

        case(5)
            p(1) =  xy*invr
            p(2) =  (a2 + y2)*invr
            p(3) = -ax*invr

        case(6)
            p(1) = -xy*invr
            p(2) =  (a2 + y2)*invr
            p(3) =  ax*invr

        case default
        print*, 'error on derivative_xdir_equidistant_gnomonic_map: invalid panel, ', panel
        end select
    return
end subroutine derivative_xdir_equidistant_gnomonic_map


subroutine derivative_ydir_equidistant_gnomonic_map(x, y, p, panel)
    !---------------------------------------------------------------
    ! this routine computes the derivative in y direction of the
    ! gnomonic mapping based on the equidistant projection
    !---------------------------------------------------------------
    real(r8), intent(in) :: x ! local coordinates
    real(r8), intent(in) :: y ! local coordinates
    real(r8), intent(out) :: p(1:3) ! projected point

    ! aux vars
    real(r8) :: r2, r32, invr, ay, xy, x2, a, a2

    ! panel
    integer(i4), intent(in) :: panel

    a = acube

    !auxiliary variables
    r2   = a**2 + x**2 + y**2
    r32  = dsqrt(r2)**3
    invr = 1._r8/r32
    ay = a*y
    xy = x*y
    x2 = x*x
    a2 = a*a

    ! compute the cartesian coordinates for each panel
    ! with the aid of the auxiliary variables  
    select case(panel)
        case(1)
            p(1) = -ay*invr
            p(2) = -xy*invr
            p(3) = (a2 + x2)*invr

        case(2)
            p(1) =  xy*invr
            p(2) = -ay*invr
            p(3) =  (a2 + x2)*invr

        case(3)
            p(1) =  ay*invr
            p(2) =  xy*invr
            p(3) = (a2 + x2)*invr
   
        case(4)
            p(1) = -xy*invr
            p(2) =  ay*invr
            p(3) = (a2 + x2)*invr 

        case(5)
            p(1) = -(a2 + x2)*invr
            p(2) = -xy*invr
            p(3) = -ay*invr

        case(6)
            p(1) =  (a2 + x2)*invr
            p(2) = -xy*invr
            p(3) =  ay*invr

        case default
            print*, 'error on derivative_ydir_equidistant_gnomonic_map: invalid panel, ', panel
    end select
    return
end subroutine derivative_ydir_equidistant_gnomonic_map

subroutine derivative_xdir_equiangular_gnomonic_map(x, y, p, panel)
    !---------------------------------------------------------------
    ! this routine computes the derivative in x direction of the
    ! gnomonic mapping based on the equiangular projection
    !---------------------------------------------------------------
    real(r8), intent(in) :: x ! angular coordinates
    real(r8), intent(in) :: y ! angular coordinates
    real(r8), intent(out) :: p(1:3) ! tangent vector
    real(r8) :: a, cos2x, cos2y, tanx, tany

    ! panel
    integer(i4), intent(in) :: panel

    a = acube
    tanx = a*dtan(x)
    tany = a*dtan(y)
    cos2x = dcos(x)**2
    cos2y = dcos(y)**2

    call derivative_xdir_equidistant_gnomonic_map(tanx, tany, p, panel)

    p = p/cos2x
    !p = p/cos2y
    p = a*p
end subroutine derivative_xdir_equiangular_gnomonic_map


subroutine derivative_ydir_equiangular_gnomonic_map(x, y, p, panel)
    !---------------------------------------------------------------
    ! this routine computes the derivative in y direction of the
    ! gnomonic mapping based on the equiangular projection
    !---------------------------------------------------------------
    real(r8), intent(in) :: x ! angular coordinates
    real(r8), intent(in) :: y ! angular coordinates
    real(r8), intent(out) :: p(1:3) ! tangent vector
    real(r8) :: a, cos2x, cos2y, tanx, tany

    ! panel
    integer(i4), intent(in) :: panel

    a = acube
    tanx = a*dtan(x)
    tany = a*dtan(y)
    cos2x = dcos(x)**2
    cos2y = dcos(y)**2

    call derivative_ydir_equidistant_gnomonic_map(tanx, tany, p, panel)

    !p = p/cos2x
    p = p/cos2y
    p = a*p
    return
end subroutine derivative_ydir_equiangular_gnomonic_map


subroutine inverse_equiangular_gnomonic_map(x, y, panel, p, mesh)
    !---------------------------------------------------------------
    ! this routine computes the inverse of the gnomonic mapping based on the equidistant projection
    ! defined by rancic et al (96) for each panel
    ! - returns the cube coordinates of the
    ! projected points.
    !
    !---------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    
    real(r8), intent(out) :: x ! local coordinates
    real(r8), intent(out) :: y ! local coordinates
    real(r8), intent(in) :: p(1:3) ! projected point

    ! aux vars
    real(r8) :: tanx, tany

    ! panel
    integer(i4), intent(inout) :: panel

    call get_cubedsphere_panel(p, panel, mesh)

    ! compute the local coordinates for each panel
    select case(panel)
        case(1)
            tanx = p(2)/p(1)
            tany = p(3)/p(1)
        case(2)
            tanx = -p(1)/p(2)
            tany =  p(3)/p(2)
        case(3)
            tanx =  p(2)/p(1)
            tany = -p(3)/p(1)
        case(4)
            tanx = -p(1)/p(2)
            tany = -p(3)/p(2)
        case(5)
            tanx =  p(2)/p(3)
            tany = -p(1)/p(3)
        case(6)
            tanx = -p(2)/p(3)
            tany = -p(1)/p(3)

        case default
            print*, 'error on inverse_equidistant_gnomonic_map: invalid panel, ', panel
            stop
        end select
        x = datan(tanx)
        y = datan(tany)
    return
end subroutine inverse_equiangular_gnomonic_map


subroutine get_cubedsphere_panel(p, panel, mesh)
!---------------------------------------------------------------
! get_cubedsphere_panel
!
! given a point on the sphere, this routine finds the panel
! where p is located
!
!---------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    real(r8), intent(in) :: p(1:3) ! given point
    real(r8)  :: p1(1:3), p2(1:3), p3(1:3), p4(1:3) ! aux vars

    integer(i4), intent(inout):: panel
    integer(i4) :: i0, j0, iend, jend, indexes, panel2

    logical :: insidepanel

    !----------------------------------------------------------------------------------
    ! reference: lauritzen, p. h., bacmeister, j. t., callaghan, p. f., and taylor,
    !   m. a.: ncar_topo (v1.0): ncar global model topography generation software
    !   for unstructured grids, geosci. model dev., 8, 3975–3986,
    !   https://doi.org/10.5194/gmd-8-3975-2015, 2015.
    !----------------------------------------------------------------------------------
    panel = -1

    if(abs(p(1))>=abs(p(2)) .and. abs(p(1))>=abs(p(3)) .and. p(1)>0._r8)then
        panel = 1 

    else if(abs(p(2))>=abs(p(1)) .and. abs(p(2))>=abs(p(3)) .and. p(2)>0._r8)then
        panel = 2

    else if(abs(p(1))>=abs(p(2)) .and. abs(p(1))>=abs(p(3)) .and. p(1)<0._r8)then
        panel = 3

    else if(abs(p(2))>=abs(p(1)) .and. abs(p(2))>=abs(p(3)) .and. p(2)<0._r8)then
        panel = 4

    else if(abs(p(3))>=abs(p(1)) .and. abs(p(3))>=abs(p(2)) .and. p(3)>0._r8)then
        panel = 5

    else if(abs(p(3))>=abs(p(1)) .and. abs(p(3))>=abs(p(2)) .and. p(3)<=0._r8)then
        panel = 6
    end if

    if(panel == -1)then
        print*, 'error on get_cubedsphere_panel: couldnt find the panel.'
        stop
    end if

end subroutine get_cubedsphere_panel


subroutine ll2contra(u, v, ucontra, vcontra, M)
!---------------------------------------------------------------
! LL2CONTRA
!
! Given the components (u,v) (u-east, v-north) a vector in the
! geographical coordinates system and conversion
! matrix between latlon and the cubedsphere coordinates system,
! this routine computes the contravariant components of the vector
!---------------------------------------------------------------
    real(r8), intent(in):: M(1:2,1:2)
    real(r8), intent(in):: u, v
    real(r8), intent(out):: ucontra, vcontra

    ucontra = M(1,1)*u + M(1,2)*v
    vcontra = M(2,1)*u + M(2,2)*v

    return
end subroutine

subroutine contra2ll(u, v, ucontra, vcontra, M)
!---------------------------------------------------------------
! CONTRA2LL
!
! Given the contravariant components of  a vector in the
! cubedsphere coordinates system and the conversion
! matrix between the cubedsphere coordinates and lalot coordinates,
! this routine computes the geographical components of the vector
! (u,v) (u-east, v-north)
!---------------------------------------------------------------
    real(r8), intent(in):: M(1:2,1:2)
    real(r8), intent(in):: ucontra, vcontra
    real(r8), intent(out):: u, v

    u = M(1,1)*ucontra + M(1,2)*vcontra
    v = M(2,1)*ucontra + M(2,2)*vcontra

    return

end subroutine

subroutine tangent_ll_lon (lon, e_lon)
    !------------------------------------------------------------------------------------
    !
    ! generate the unit tangent (r^3) vector in the geographic coordinates in longitude direction
    !
    !------------------------------------------------------------------------------------
    real (r8), intent(in)  :: lon
    real (r8), intent(out) :: e_lon(1:3)

    e_lon(1) = -dsin(lon)
    e_lon(2) =  dcos(lon)
    e_lon(3) =  0._r8
    return
end subroutine tangent_ll_lon

subroutine tangent_ll_lat (lon, lat, e_lat)
    !------------------------------------------------------------------------------------
    !
    ! generate the unit tangent (r^3) vector in the geographic coordinates in latitude direction
    !
    !------------------------------------------------------------------------------------
    real (r8), intent(in)  :: lat, lon
    real (r8), intent(out) :: e_lat(1:3)
    e_lat(1) = -dsin(lat)*dcos(lon)
    e_lat(2) = -dsin(lat)*dsin(lon)
    e_lat(3) =  dcos(lat)
    return
end subroutine tangent_ll_lat


!===============================================================================================
!   the following routines were taken from imodel  https://github.com/pedrospeixoto/imodel
!===============================================================================================

subroutine sph2cart (lon, lat, x, y, z )
  !------------------------------------------------------------------------------------
  ! sph2cart 
  !
  !     transforms geographical coordinates (lat,lon) to cartesian coordinates.
  !     similar to stripack's trans
  !
  !    input: lat, latitudes of the node in radians [-pi/2,pi/2]
  !           lon, longitudes of the nodes in radians [-pi,pi]
  !
  !    output:  x, y, z, the coordinates in the range -1 to 1. 
  !                    x**2 + y**2 + z**2 = 1 
  !---------------------------------------------------------------------
  real (r8), intent(in) :: lon
  real (r8), intent(in) :: lat
  real (r8), intent(out) :: x
  real (r8), intent(out) :: y
  real (r8), intent(out) :: z
  real (r8):: coslat

  coslat = dcos (lat)
  x = coslat * dcos (lon)
  y = coslat * dsin (lon)
  z = dsin (lat)

  return
end subroutine sph2cart

subroutine cart2sph ( x, y, z, lon, lat )
  !---------------------------------------------------------------------
  ! cart2sph 
  !     transforms  cartesian coordinates to geographical (lat,lon) coordinates .
  !     similar to stripack's scoord
  !
  !    input:  x, y, z, the coordinates in the range -1 to 1. 
  !
  !    output, lon, longitude of the node in radians [-pi,pi].
  !                       lon=0 if point lies on the z-axis.  
  !    output, lat, latitude of the node in radians [-pi/2,pi/2].
  !                       lat=0 if   x**2 + y**2 + z**2 = 0.
  !------------------------------------------------------------------------------------
  real    (r8), intent(in) :: x
  real    (r8), intent(in) :: y
  real    (r8), intent(in) :: z
  real    (r8), intent(out) :: lat
  real    (r8), intent(out) :: lon
  real    (r8):: pnrm

  pnrm = dsqrt ( x **2 + y**2 + z**2 )
  if ( x /= 0.0d+00 .or. y /= 0.0d+00 ) then
     lon = datan2 ( y, x )
  else
     lon = 0.0d+00
  end if

  if ( pnrm /= 0.0d+00 ) then
     lat = dasin ( z / pnrm )
  else
     print*, "cart2sph warning: point not in the unit sphere. norm= ", pnrm
     lat = 0.0d+00
  end if

  return
end subroutine cart2sph

subroutine convert_vec_cart2sph(p, v , vlon, vlat)
    !---------------------------------------------------------------------
    !	convert_vec_cart2sph
    !
    !   recieves a point p=(x,y,z) and a vector (v) at this point
    !   returns the vector at this point in (lon,lat) reference
    !      vlon=west-east direction
    !      vlat=south-north direction
    !
    !   the vector must be tangent to the sphere, that is, v.(x,y,z)=0
    !   this is done in order to plot the vectorfield on gmt
    !---------------------------------------------------------------------
    !point cartesian coords
    real(r8), intent(in) :: p(1:3)
    !cartesian vector on point
    real(r8), intent(in) :: v(1:3)
    !spherical coord vector on point
    real(r8), intent(out) :: vlat
    real(r8), intent(out) :: vlon
    !auxiliar variables
    real(r8):: r
    real(r8):: rho
    real(r8):: rvec(1:3)
    real(r8):: latvec(1:3)
    real(r8):: lonvec(1:3)
    real(r8):: zero
    real(r8):: test

    zero=0
    r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
    rho=dsqrt(p(1)**2+p(2)**2)

    !case where the point is in the north or south pole
    if(rho<10*eps)then
       !print*, "pole:", v
       vlon=v(2)
       vlat=v(1)
       return
    end if

    rvec=(/p(1),p(2),p(3)/)
    rvec=rvec/r

    latvec=(/-p(1)*p(3),-p(2)*p(3),rho**2/)
    latvec=latvec/rho

    lonvec=(/-p(2), p(1), zero /)
    lonvec=lonvec/rho

    test=dot_product(v,rvec)
    if(abs(test)>10e-5)then
       print *,"CONVERT_VEC_CART2SPH Warning: Vector not tangent to sphere."
       print '(a,3f10.6)',"Vector:",v(1:3)
       print '(a,3f10.6)',"Point:",rvec(1:3)
       print '(a,f10.6)',"Dot Product:", test
       stop
    end if
    vlat=dot_product(v,latvec)
    vlon=dot_product(v,lonvec)

    return    
end subroutine convert_vec_cart2sph

subroutine convert_vec_sph2cart(vlon, vlat, p, v)
    !---------------------------------------------------------------------
    !	CONVERT_VEC_SPH2CART
    !
    !   Recieves a point p=(x,y,z) and a vector at this point in lon, lat system
    !   in radians (vlon, vlat), ie, 
    !      vlon=West-East direction
    !      vlat=South-North direction
    !   Returns the vector at this point in cartesian reference (v)
    !---------------------------------------------------------------------
    !Point cartesian coords
    real(r8), intent(in) :: p(1:3)
    real(r8), intent(in) :: vlon
    real(r8), intent(in) :: vlat
    !Cartesian vector on point
    real(r8), intent(out) :: v(1:3)
    !Auxiliar variables
    real(r8):: r
    real(r8):: rho

    r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
    rho=dsqrt(p(1)**2+p(2)**2)

    !Case where the point is in the north or south pole
    if(rho==0)then
       v(1)=vlat
       v(2)=vlon
       v(3)=0
       return
    else    
       !The minus sign in vlat is due to the diference between spherical coords and
       !   geographical coords
       v(1)=-vlon*(p(2)/rho) - (vlat*p(1)*p(3))/(rho)
       v(2)=vlon*(p(1)/rho) - (vlat*p(2)*p(3))/(rho)
       v(3)=vlat*rho
    end if

    return    
end subroutine convert_vec_sph2cart

end module sphgeo
