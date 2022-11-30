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


  subroutine equidistant_gnomonic_map(x, y, p, panel)
      !---------------------------------------------------------------
      ! This routine computes the Gnomonic mapping based on the equidistant projection
      ! defined by Rancic et al (96) for each panel
      ! - x, y are the variables defined in [-a, a].
      ! - The projection is applied on the points (x,y)
      ! - Returns the Cartesian coordinates of the
      ! projected point p.
      !
      ! References: 
      ! - Rancic, M., Purser, R.J. and Mesinger, F. (1996), A global shallow-water model using an expanded
      !  spherical cube: Gnomonic versus conformal coordinates. Q.J.R. Meteorol. Soc., 122: 959-982. 
      !  https://doi.org/10.1002/qj.49712253209
      ! - Nair, R. D., Thomas, S. J., & Loft, R. D. (2005). A Discontinuous Galerkin Transport Scheme on the
      ! Cubed Sphere, Monthly Weather Review, 133(4), 814-828. Retrieved Feb 7, 2022, 
      ! from https://journals.ametsoc.org/view/journals/mwre/133/4/mwr2890.1.xml
      !
      !---------------------------------------------------------------
      real(r8), intent(in) :: x ! Local coordinates
      real(r8), intent(in) :: y ! Local coordinates
      real(r8), intent(out) :: p(1:3) ! Projected point

      ! aux vars
      real(r8) :: r, r2, invr,  xor, yor, aor

      ! panel
      integer(i4), intent(in) :: panel


      ! Auxiliary variables
      r2   = acube**2 + x**2 + y**2
      r    = dsqrt(r2)
      invr = 1._r8/r
      xor  = invr*x
      yor  = invr*y
      aor  = invr*acube

      ! Compute the Cartesian coordinates for each panel
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
            print*, 'ERROR on equidistant_gnomonic_map: invalid panel, ', panel
         end select
         return
   end subroutine equidistant_gnomonic_map

   subroutine inverse_equidistant_gnomonic_map(x, y, p, mesh)
      !---------------------------------------------------------------
      ! This routine computes the inverse of the Gnomonic mapping based on the equidistant projection
      ! defined by Rancic et al (96) for each panel
      ! - Returns the cube coordinates of the
      ! projected points.
      !
      !---------------------------------------------------------------
      type(cubedsphere), intent(inout) :: mesh

      real(r8), intent(out) :: x ! Local coordinates
      real(r8), intent(out) :: y ! Local coordinates
      real(r8), intent(inout) :: p(1:3) ! Projected point

      ! aux vars
      real(r8) :: r

      ! panel
      integer(i4) :: panel

      ! Get the panel
      call get_cubedsphere_panel(p, panel, mesh)
      print*, panel
      print*
      ! Compute the local coordinates for each panel
      select case(panel)
         case(1)
            r = acube/p(1)
            x = p(2)*r
            y = p(3)*r

         case(2)
            r =  acube/p(2)
            x = -p(1)*r
            y =  p(3)*r

         case(3)
            r = -acube/p(1)
            x = -p(2)*r
            y =  p(3)*r
     
         case(4)
            r = -acube/p(2)
            x =  p(1)*r
            y =  p(3)*r

         case(5)
            r =  acube/p(3)
            x =  p(2)*r
            y = -p(1)*r

         case(6)
            r = -acube/p(3)
            x =  p(2)*r
            y =  p(1)*r

         case default
            print*, 'ERROR on inverse_equidistant_gnomonic_map: invalid panel, ', panel
            stop
         end select
         return
   end subroutine inverse_equidistant_gnomonic_map


   subroutine get_cubedsphere_panel(p, panel, mesh)
      !---------------------------------------------------------------
      ! GET_CUBEDSPHERE_PANEL
      !
      ! Given a point on the sphere, this routine finds the panel
      ! where p is located
      !
      !---------------------------------------------------------------
      type(cubedsphere), intent(inout) :: mesh

      real(r8), intent(in) :: p(1:3) ! Given point
      real(r8)  :: p1(1:3), p2(1:3), p3(1:3), p4(1:3) ! aux vars

      integer(i4), intent(inout):: panel
      integer(i4) :: i0, j0, iend, jend, indexes, panel2

      logical :: insidepanel

      !----------------------------------------------------------------------------------
      ! Reference: Lauritzen, P. H., Bacmeister, J. T., Callaghan, P. F., and Taylor,
      !   M. A.: NCAR_Topo (v1.0): NCAR global model topography generation software
      !   for unstructured grids, Geosci. Model Dev., 8, 3975–3986,
      !   https://doi.org/10.5194/gmd-8-3975-2015, 2015.
      !----------------------------------------------------------------------------------
      !panel = -1

      !if(abs(p(1))>=abs(p(2)) .and. abs(p(1))>=abs(p(3)) .and. p(1)>0._r8)then
      !   panel = 1 

      !else if(abs(p(2))>=abs(p(1)) .and. abs(p(2))>=abs(p(3)) .and. p(2)>0._r8)then
      !   panel = 2

      !else if(abs(p(1))>=abs(p(2)) .and. abs(p(1))>=abs(p(3)) .and. p(1)<0._r8)then
      !   panel = 3

      !else if(abs(p(2))>=abs(p(1)) .and. abs(p(2))>=abs(p(3)) .and. p(2)<0._r8)then
      !   panel = 4

      !else if(abs(p(3))>=abs(p(1)) .and. abs(p(3))>=abs(p(2)) .and. p(3)>0._r8)then
      !   panel = 5

      !else if(abs(p(3))>=abs(p(1)) .and. abs(p(3))>=abs(p(2)) .and. p(3)<=0._r8)then
      !   panel = 6
      !end if

      !if(panel == -1)then
      !  print*, 'ERROR on get_cubedsphere_panel: couldnt find the panel.'
      !  stop
      !end if

      ! Check if p in inside the panel quadrilateral
      i0   = mesh%i0
      iend = mesh%iend
      j0   = mesh%j0
      jend = mesh%jend

      do panel2 = 1, nbfaces
          p1 = mesh%po(i0-1, j0-1, panel2)%p
          p2 = mesh%po(iend, j0-1, panel2)%p
          p3 = mesh%po(iend, jend, panel2)%p
          p4 = mesh%po(i0-1, jend, panel2)%p

          ! check if p is inside the panel
          insidepanel = insidequad(p, p1, p2, p3, p4) 

          if (insidepanel)then
               panel = panel2
               !print*,'found ', panel2- panel
               !if(panel2 .ne. panel )then
               !  print*, 'error, ', panel2, panel
                !read(*,*)
               !  stop
               !endif
             return
          end if
      end do

      if(panel2 == 7)then
         print*, 'ERROR on get_cubedsphere_panel: couldnt find the panel.'
         stop
      end if

  end subroutine get_cubedsphere_panel


  subroutine binary_search(p, mesh, ix, jy, panel)
      !---------------------------------------------------------------
      ! BINARY_SEARCH 
      !
      ! Given a point p on the sphere and the panel where the point is located, 
      ! this routine finds quadrilateral that contains p perfoming
      ! a binary search
      !
      !---------------------------------------------------------------
      type(cubedsphere), intent(inout) :: mesh

      real(r8), intent(in) :: p(1:3) ! Given point
      real(r8) :: p1(1:3), p2(1:3), p3(1:3), p4(1:3), p5(1:3)  ! aux vars
      real(r8) :: p6(1:3), p7(1:3), p8(1:3), p9(1:3)

      integer(i4), intent(inout) :: ix, jy, panel
      
      ! Aux vars
      integer(i4) :: i, j, i0, j0, iend, jend, imid, jmid, panel2

      ! Logical vars
      logical :: insidequad1
      logical :: insidequad2
      logical :: insidequad3
      logical :: insidequad4

      ! Get the panel
      call get_cubedsphere_panel(p, panel, mesh)

      i0   = mesh%i0-1
      iend = mesh%iend
      j0   = mesh%j0-1
      jend = mesh%jend
      panel2=0

      do while(iend-i0>1 .or. jend-j0>1)
         imid = (i0+iend)/2
         jmid = (j0+jend)/2

         ! Get the points
         p1 = mesh%po(i0  , j0  , panel)%p
         p2 = mesh%po(imid, j0  , panel)%p
         p3 = mesh%po(iend, j0  , panel)%p

         p4 = mesh%po(i0  , jmid, panel)%p
         p5 = mesh%po(imid, jmid, panel)%p
         p6 = mesh%po(iend, jmid, panel)%p

         p7 = mesh%po(i0  , jend, panel)%p
         p8 = mesh%po(imid, jend, panel)%p
         p9 = mesh%po(iend, jend, panel)%p

         ! check if p is inside the quadrilateral formed by p1, p2, p5, p4
         insidequad1 = insidequad(p, p1, p2, p5, p4) 
 
         if (.not. insidequad1) then
            ! check if p is inside the quadrilateral formed by p2, p3, p6, p5
            insidequad2 = insidequad(p, p2, p3, p6, p5)

            if (.not. insidequad2)then
               ! check if p is inside the quadrilateral formed by p5, p6, p9, p8
               insidequad3 = insidequad(p, p5, p6, p9, p8)    

               if (.not. insidequad3)then
                  ! check if p is inside the quadrilateral formed by p4, p5, p8, p7
                  insidequad4 = insidequad(p, p4, p5, p8, p7)

                  if (.not. insidequad4)then
                     print*, 'ERROR on binary_search: point is not contained in any quadrilateral. '
                     stop

                  else

                     ! p is inside the quadrilateral formed by p4, p5, p8, p7
                     iend = imid
                     j0   = jmid

                  end if
               else

                  ! p is inside the quadrilateral formed by p5, p6, p9, p8
                  i0 = imid
                  j0 = jmid

               end if

            else

               ! p is inside the quadrilateral formed by p2, p3, p6, p5
               i0   = imid
               jend = jmid

            end if
         else

            ! p is inside the quadrilateral formed by p1, p2, p5, p4
            iend = imid
            jend = jmid

         end if
      end do

      ! save indexes
      ix = iend
      jy = jend
  end subroutine




  function midpoint(p1, p2) 
    !--------------------------------------------------------------
    ! MIDPOINT 
    !  Calculates the midpoint of a geodesic arc
    !  Returns a 'point_structure' kind
    !-------------------------------------------------------------
    real(r8), intent(in) :: p1(1:3), p2(1:3)
    real(r8):: midpoint(1:3)

    ! Aux var
    real(r8):: p(1:3)

    !Position
    p = (p1 + p2)*0.5_r8
    p = p/norm(p)
    midpoint = p
    return
  end function midpoint



  !===============================================================================================
  !   The following routines were taken from iModel  https://github.com/pedrospeixoto/iModel
  !===============================================================================================

   subroutine sph2cart (lon, lat, x, y, z )
      !------------------------------------------------------------------------------------
      ! SPH2CART 
      !
      !     Transforms geographical coordinates (lat,lon) to Cartesian coordinates.
      !     Similar to stripack's TRANS
      !
      !    Input: LAT, latitudes of the node in radians [-pi/2,pi/2]
      !           LON, longitudes of the nodes in radians [-pi,pi]
      !
      !    Output:  X, Y, Z, the coordinates in the range -1 to 1. 
      !                    X**2 + Y**2 + Z**2 = 1 
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
      ! CART2SPH 
      !     Transforms  Cartesian coordinates to geographical (lat,lon) coordinates .
      !     Similar to stripack's SCOORD
      !
      !    Input:  X, Y, Z, the coordinates in the range -1 to 1. 
      !
      !    Output, LON, longitude of the node in radians [-pi,pi].
      !                       LON=0 if point lies on the Z-axis.  
      !    Output, LAT, latitude of the node in radians [-pi/2,pi/2].
      !                       LAT=0 if   X**2 + Y**2 + Z**2 = 0.
      !------------------------------------------------------------------------------------
      real    (r8), intent(in) :: x
      real    (r8), intent(in) :: y
      real    (r8), intent(in) :: z
      real    (r8), intent(out) :: lat
      real    (r8), intent(out) :: lon
      real    (r8):: pnrm

      pnrm = dsqrt ( x **2 + y**2 + z**2 )
      if ( x /= 0.0D+00 .or. y /= 0.0D+00 ) then
         lon = datan2 ( y, x )
      else
         lon = 0.0D+00
      end if

      if ( pnrm /= 0.0D+00 ) then
         lat = dasin ( z / pnrm )
      else
         print*, "CART2SPH Warning: Point not in the unit sphere. Norm= ", pnrm
         lat = 0.0D+00
      end if

      return
    end subroutine cart2sph

    subroutine convert_vec_cart2sph(p, v , vlon, vlat)
      !---------------------------------------------------------------------
      !	CONVERT_VEC_CART2SPH
      !
      !   Recieves a point p=(x,y,z) and a vector (v) at this point
      !   Returns the vector at this point in (lon,lat) reference
      !      vlon=West-East direction
      !      vlat=South-North direction
    !
    !   The vector must be tangent to the sphere, that is, v.(x,y,z)=0
    !   This is done in order to plot the vectorfield on GMT
    !---------------------------------------------------------------------
    !Point cartesian coords
    real(r8), intent(in) :: p(1:3)
    !Cartesian vector on point
    real(r8), intent(in) :: v(1:3)
    !Spherical coord vector on point
    real(r8), intent(out) :: vlat
    real(r8), intent(out) :: vlon
    !Auxiliar variables
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

    !Case where the point is in the north or south pole
    if(rho<10*eps)then
       !print*, "Pole:", v
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

  function proj_vec_sphere(v, p)
    !-----------------------------------------------------------
    !  Projects a vector 'v' on the plane tangent to a sphere
    !   Uses the the tangent plane relative to the unit sphere's
    !   point 'p', in cartesian coords
    !-----------------------------------------------------------
    real (r8), intent(in), dimension(1:3) :: v
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), dimension(1:3)  :: proj_vec_sphere

    proj_vec_sphere(1:3)=&
         v(1:3)-dot_product(v,p)*p(1:3)/norm(p)

    return
  end function proj_vec_sphere

  function arcdistll(lon1, lat1, lon2 , lat2)
    !----------------------------------------------------
    ! ARCDISTLL
    !   Arc length between p1 and p2 points on the sphere
    !     Receives latitude longitude coordinates in radians 
    !	  Returns the angle between points in radians
    !     Lat [-pi/2, pi/2] , lon [-pi, pi[   
    !
    !	Uses Vincenty formula
    !-----------------------------------------------------
    real (r8), intent(in) :: lat1
    real (r8), intent(in) :: lon1
    real (r8), intent(in) :: lat2
    real (r8), intent(in) :: lon2
    real (r8):: arcdistll
    real (r8):: dlat
    real (r8):: dlon

    dlat=lat1-lat2
    dlon=lon1-lon2

    arcdistll=0

    arcdistll= datan2(dsqrt( (dcos(lat2)*dsin(dlon))**2 + &
         (dcos(lat1)*dsin(lat2)-dsin(lat1)*dcos(lat2)*dcos(dlon))**2 ),&
         (dsin(lat1)*dsin(lat2)+dcos(lat1)*dcos(lat2)*dcos(dlon)) )
    !See wikipedia - Vincenty formula for the sphere
    arcdistll=abs(arcdistll)
    return

  end function arcdistll

  function arcdistxyz(p1, p2)
    !----------------------------------------------------
    ! ARCDISTXYZ - Prefer to use ARCLEN
    !    Arc length between p1 and p2 points on the sphere
    !     Receives cartesian coords. of points
    !     Returns the angle between points in radians
    !     It if preferable to use the routine ARCLEN
    !      to avoid non standard inverse cossine 
    !-----------------------------------------------------
    real (r8), dimension(3), intent(in) :: p1
    real (r8), dimension(3), intent(in) :: p2
    real (r8):: arcdistxyz
    real (r8):: p1norm
    real (r8):: p2norm

    arcdistxyz=0
    p1norm=dsqrt(dot_product(p1,p1))
    p2norm=dsqrt(dot_product(p2,p2))

    if(p1norm==0.or.p2norm==0) then
       print*
       print*,"ARCDIST Warning: Distance calculation error"
       print*,"p1:", p1, " p2:", p2
       return
    end if
    arcdistxyz=dacos(dot_product(p1,p2)/(p1norm*p2norm))

    if((p1norm >= p2norm+100*eps).or. &
         (p1norm <= p2norm-100*eps)) then
       print*, "ARCDIST WARNING: Vectors with different norms (dif):",p1norm-p2norm
       print*, "ARCDIST WARNING: Vectors are not on the same sphere! "
       print*
       return
    end if

    return
  end function arcdistxyz

  function arclen(p, q)
    !-----------------------------------------------------------
    ! ARCLEN from ssrfpack by r. renka            
    !
    !   This function computes the arc-length (angle in radians)            
    !    between a pair of points on the unit sphere. It is similar to
    !    arcdistxyz, but uses a calculation to avoid inverse cosine function
    !    p,q = arrays of length 3 containing the x, y, and z             
    !             coordinates (in that order) of points on the              
    !             unit sphere.                                                 
    !   Returns the angle in radians between the unit vectors              
    !       p and q.  0 <= arclen <= pi.                       
    !---------------------------------------------------------
    real (r8), intent(in) :: p(3)
    real (r8), intent(in) :: q(3)
    real (r8):: arclen

    !Dot product
    real (r8):: d

    d = dot_product(p+q, p+q)

    if (d==0.) then 
       ! p and q are separated by 180 degrees.
       arclen = pi
    elseif (d>=4.) then 
       ! p and q coincide.
       arclen = 0._r8
    else 
       arclen = 2._r8 * datan (dsqrt ( (4._r8 - d) / d) )
    endif

    return 
  end function arclen


  function sphtriarea(p1, p2, p3)
    !------------------------------------------
    !  SPHTRIAREA
    ! Calculates the area of 
    !    a spherical triangle on the unit sphere
    !    given the cartesian coords of nodes
    !--------------------------------------------

    !Cartesian coordinates of the nodes of spherical triangle
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Variables for spherical triangle area calculation
    ! a, b, c are lengths of the the 3 sides
    ! s is the semi-perimiter s=(a+b+c)/2
    ! e is the spherical excess of triangle
    !   e = a + b + c - 180, but we will useing L'Huilier's Theorem
    real (r8):: a
    real (r8):: b
    real (r8):: c
    real (r8):: s
    real (r8):: e
    real (r8):: tmp

    !Spherical triangle area
    real (r8):: sphtriarea

    !Check for degeneration
    if(norm(p1-p2)<eps/100 .or. norm(p2-p3)<eps/100 .or. norm(p1-p3)<eps/100)then
       sphtriarea=0._r8
       return
    end if

    !Calculate the sides length's and the semiperimiter
    s=0.
    a=arclen(p1,p2) !arcdistll(lon(1),lat(1), lon(2),lat(2))
    b=arclen(p2,p3)  !arcdistll(lon(2),lat(2), lon(3),lat(3))
    c=arclen(p3,p1)  !arcdistll(lon(3),lat(3), lon(1),lat(1))
    s=(a+b+c)/2

    !Calculate spherical triangle excess using L'Huilier's Theorem
    tmp=dtan(s/2)*dtan((s-a)/2)*dtan((s-b)/2)*dtan((s-c)/2)
    !Round off error might give almost zero negative numbers => assume zero
    if(tmp<0)then
       e=0.
    else
       e = 4*datan(dsqrt(tmp))
    end if
    sphtriarea=e

    return
  end function sphtriarea

  function sphquadarea(p1, p2, p3, p4)
    !------------------------------------------
    !  SPHQUADAREA
    ! Calculates the area of 
    !    a spherical quadrilateral on the unit sphere
    !    given the cartesian coords of nodes
    ! Compute the area splitting the quadrilateral
    ! as two triangles
    !--------------------------------------------

    !Cartesian coordinates of the nodes of spherical quadrilateral
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)
    real (r8),intent(in) :: p4(1:3)

    !Spherical triangle areas
    real (r8):: sphtriarea1
    real (r8):: sphtriarea2

    !Spherical quadrilateral area
    real (r8):: sphquadarea
 
    sphtriarea1 = sphtriarea(p1, p2, p3)
    sphtriarea2 = sphtriarea(p2, p3, p4)
    sphquadarea = sphtriarea1 + sphtriarea2
    return
  end function sphquadarea


  function sphtriangles(p1, p2, p3)
    !------------------------------------------
    !  SPHTRIANGLES
    ! Calculates the angles of 
    !    a spherical triangle on the unit sphere
    !--------------------------------------------

    !Latitudes and longitudes of nodes on the sphere
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Variables for spherical triangle area calculation
    ! a, b, c are lengths of the the 3 sides
    real (r8):: a
    real (r8):: b
    real (r8):: c

    !Spherical triangle area
    real (r8):: sphtriangles(1:3)

    !Calculate the sides length's and the semiperimiter
    a=arclen(p1,p2) !arcdistll(lon(1),lat(1), lon(2),lat(2))
    b=arclen(p2,p3)  !arcdistll(lon(2),lat(2), lon(3),lat(3))
    c=arclen(p3,p1)  !arcdistll(lon(3),lat(3), lon(1),lat(1))
    if((pi/2-a)<eps.or.(pi/2-b)<eps.or.(pi/2-c)<eps)then
       print*, "SPHTRIANGLES WARNING: Triangle with internal angle too large."
       print*, "Angles:", a, b, c
    end if
    !Using spherical law of cosine 
    sphtriangles(1)=dacos( (dcos(b)-dcos(a)*dcos(c)) / (dsin(a)*dsin(c)) )
    sphtriangles(2)=dacos( (dcos(c)-dcos(a)*dcos(b)) / (dsin(a)*dsin(b)) )
    sphtriangles(3)=dacos( (dcos(a)-dcos(b)*dcos(c)) / (dsin(b)*dsin(c)) )

    return
  end function sphtriangles

  function insidetr(p, p1, p2, p3)
    !----------------------------------------------------------
    ! INSIDETR
    !
    ! Checks if 'p' is inside the geodesical triangle formed bu p1, p2, p3
    !  The vertices of the triangle must be given ccwisely
    ! The algorithm checks if the point left of each edge
    !--------------------------------------------------------------
    !Point
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), intent(in), dimension(1:3) :: p1
    real (r8), intent(in), dimension(1:3) :: p2
    real (r8), intent(in), dimension(1:3) :: p3

    !Return true if point in triangle
    logical:: insidetr

    !Aux
    real (r8):: left

    insidetr=.false.

    !For every edge
    !Check if point is at left of edge
    left=det(p, p1, p2)
    if(left< -eps2 )then
       !Point not in triangle
       return
    end if

    left=det(p, p2, p3)
    if(left< -eps2 )then
       !Point not in triangle
       return
    end if

    left=det(p, p3, p1)
    if(left< -eps2 )then
       !Point not in triangle
       return
    end if

    !If left <=0  to all edge, then the point
    !  is inside the triangle, or on the edge
    insidetr=.true.

    return
  end function insidetr

  function insidequad(p, p1, p2, p3, p4)
    !----------------------------------------------------------
    ! INSIDEQUAD
    !
    ! Checks if 'p' is inside the geodesical quadrilateral formed
    ! by p1, p2, p3, p4
    ! The vertices of the quadrilateral must be given ccwisely
    ! Uses the routine insidetr
    !--------------------------------------------------------------
    !Point
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), intent(in), dimension(1:3) :: p1
    real (r8), intent(in), dimension(1:3) :: p2
    real (r8), intent(in), dimension(1:3) :: p3
    real (r8), intent(in), dimension(1:3) :: p4

    !Return true if point in quadrilateral
    logical:: insidequad

    insidequad = .false.

    insidequad = insidetr(p, p1, p2, p3)

    if(.not. insidequad) then
        insidequad = insidetr(p, p1, p3, p4)
    end if

    !if(.not. insidequad) then
    !    insidequad = insidetr(p, p2, p3, p4)
    !end if

    !if(.not. insidequad) then
    !    insidequad = insidetr(p, p2, p1, p4)
    !end if



    return
  end function insidequad


  function arcintersec(p1, p2, q1, q2, intpt)
    !-------------------------------------------------------------
    !  ARCINTERSEC
    !   
    !  Tests for existence of intersection between 2 geodesic arcs
    !  given by p1->p2 and q1->q2. Returns a logical .true. in case of
    !  intersection. An optional returning argument is the intersection point 'intpt'
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), dimension(1:3), intent(in) :: q1
    real(r8), dimension(1:3), intent(in) :: q2

    !intpt -> intersectin point
    real(r8), dimension(1:3), optional :: intpt
    real(r8), dimension(1:3) :: r

    !np -> normal to (p1, p2, zero) plane
    !nq -> normal to (q1, q2, zero) plane
    !nl -> np x nq
    real(r8), dimension(1:3) :: np
    real(r8), dimension(1:3) :: nq
    real(r8), dimension(1:3) :: nl

    !arcpr
    real (r8):: arcpr
    real (r8):: arcqr
    real (r8):: arcpq
    real (r8):: normnl
    real (r8), parameter :: eps=10e-6
    logical:: arcintersec
    logical:: arcintp

    !Initate variables
    arcintersec=.false.
    if(present(intpt))then
       intpt=0._r8
    endif
    r=0._r8

    !Calculate normal to planes formed by
    !  the arc points and the origin
    np=cross_product(p1,p2)
    np=np/norm(np)
    nq=cross_product(q1,q2)
    nq=nq/norm(nq)

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes 
    ! defined by np and nq 
    nl=cross_product(np,nq)

    !print*
    !print '(6f8.2)',p1,p2
    !print '(6f8.2)', q1, q2
    !print*,norm(nl)

    normnl=norm(nl)
    if(normnl<eps)then
       !Arcs are on the same great circle

       !Check if q1 is in arc (p1,p2)
       arcpq=arclen(p1,q1)+arclen(p2,q1)
       arcpr=arclen(p1,p2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if q2 is in arc (p1,p2)
       arcpq=arclen(p1,q2)+arclen(p2,q2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if p1 is in arc (q1,q2)
       arcpq=arclen(q1,p1)+arclen(q2,p1)
       arcpr=arclen(q1,q2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if p2 is in arc (q1,q2)
       arcpq=arclen(q1,p2)+arclen(q2,p2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

    else
       !Arcs are on diferent great circles

       !The intersection point must be in the intersection of
       ! the planes defined by np and nq and also belong to the sphere
       !So this point is +/- nl/norm(nl)
       r=nl/norm(nl)

       arcintp=.false.
       !Check if r is in arc (p1,p2)
       arcpr=arclen(p1,r)+arclen(p2,r)
       if(arcpr<=arclen(p1,p2)+eps)then
          arcintp=.true.
       end if
       if(.not.arcintp)then
          !Check the other point (-r)
          arcpr=arclen(p1,-r)+arclen(p2,-r)
          if(arcpr<=arclen(p1,p2)+eps)then
             arcintp=.true.
             r=-r
          else
             !Neither points (+/-r) are in (p1,p2) arc
             !  There cannot be intersection
             return
          end if
       endif

       !Check if the intersec point of (p1,p2) is in arc (q1,q2)
       arcqr=arclen(q1,r)+arclen(q2,r)
       if(arcintp .and. arcqr<=arclen(q1,q2)+eps)then
          arcintersec=.true.
       end if

    end if

    if(present(intpt).and. arcintersec)then
       intpt=r
    endif

    return
  end function arcintersec

  function gcircarcintersec(n, p1, p2)
    !-------------------------------------------------------------
    !  GCIRCARCINTERSEC
    !   
    !  Calculates the intersection between a great circle, 
    !  given by its normal component (n), and a geodesic arc
    !  given by p1->p2. Returns the intersection point.
    !
    !  If no point of intersection exists, r=(0, 0, 0)
    !  If arc on great circle, r=(9,0,0)
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: n
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2

    !intpt -> intersectin point
    real(r8), dimension(1:3) :: gcircarcintersec
    real(r8), dimension(1:3) :: r

    !np -> normal to (p1, p2, zero) plane
    !nl -> np x nq
    real(r8), dimension(1:3) :: np
    real(r8), dimension(1:3) :: nl

    !arcpr
    real (r8):: arcpr
    real (r8):: arcp1r
    real (r8):: arcp2r
    real (r8):: arcpq
    real (r8):: normnl
    !real (r8), parameter :: eps=10e-6
    logical:: lcircarcintersec

    !Initate variables
    lcircarcintersec=.false.
    gcircarcintersec=0._r8
    r=0._r8

    !Calculate normal to plane formed by
    !  the arc points and the origin
    np=cross_product(p1,p2)

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes 
    ! defined by np and n
    nl=cross_product(np,n)
    normnl=norm(nl)

    if(normnl<eps)then
       !Arc is on the great circle

       !There are infinite the intersection points
       ! intersection point is returned as (9, 0, 0)
       lcircarcintersec=.true.
       r(1)=9

    else
       !Arc not on the great circles

       !The intersection point must be in the intersection of
       ! the planes defined by np and nq and also belong to the sphere
       !So this point is +/- nl/norm(nl)
       r=nl/norm(nl)

       !Check if r is in arc (p1,p2)
       arcp1r=arclen(p1,r)
       arcp2r=arclen(p2,r)
       arcpr=arclen(p1,r)+arclen(p2,r)
       arcpq=arclen(p1,p2)
       if(arcp1r<=arcpq+eps .and. arcp2r<=arcpq+eps )then
          lcircarcintersec=.true.
       end if
       if(.not.lcircarcintersec)then
          !Check the other point (-r)
          arcp1r=arclen(p1,-r)
          arcp2r=arclen(p2,-r)
          arcpr=arclen(p1,-r)+arclen(p2,-r)
          if(arcp1r<=arcpq+eps .and. arcp2r<=arcpq+eps )then
             lcircarcintersec=.true.
             r=-r
          else
             !Neither points (+/-r) are in (p1,p2) arc
             !  There cannot be intersection
             return
          end if
       endif
    end if

    if(lcircarcintersec)then
       gcircarcintersec=r
    endif

    return
  end function gcircarcintersec

  function gcircintersec(p, np, nq)
    !-------------------------------------------------------------
    !  GCIRCINTERSEC
    !
    !  Calculates the intersection between 2 great circle,
    !  given by its normal component (np, nq) and
    !  returns the nearest intersection to the point p
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p
    real(r8), dimension(1:3), intent(in) :: np
    real(r8), dimension(1:3), intent(in) :: nq

    !intpt -> intersectin point
    real(r8), dimension(1:3) :: gcircintersec
    real(r8), dimension(1:3) :: r1
    real(r8), dimension(1:3) :: r2

    !nl -> np x nq
    real(r8), dimension(1:3) ::  nl

    !Initate variables
    gcircintersec=0._r8

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes
    ! defined by np and n
    nl=cross_product(np,nq)

    !The intersection point must be in the intersection of
    ! the planes defined by np and nq and also belong to the sphere
    !So this point is +/- nl/norm(nl)
    r1=nl/norm(nl)
    r2=-nl/norm(nl)

    if(norm(r1-p) < norm(r2-p) ) then
       gcircintersec=r1
    else
       gcircintersec=r2
    endif

    return
  end function gcircintersec

  function ptingarc(p, p1, p2, limit)
    !-------------------------------------------------------------
    !  ptingarc
    !
    !  Evaluates if a a point p is in a geodesic arc formed
    !    by p1-p2
    !
    !  Returns true if p in arc, or false if not
    !
    !  Limit : maximum angle (rad) to be considered still in the arc
    !      ex: 10e-7
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), intent(in) :: limit
    !real(r8) :: limitcos, limitsin

    !Logical output
    logical:: ptingarc

    !Normal - v1 x v2
    real(r8), dimension(1:3) :: nl
    real(r8):: normnl

    !Dot product between plane normal and point
    real(r8):: dp

    !Projection of point on plane
    !real(r8), dimension(1:3) :: proj

    !Angle
    real(r8):: d
    real(r8):: d1
    real(r8):: d2

    !Default return value
    ptingarc=.false.

    !Normal to plane formed by p1, p2
    nl=cross_product(p1,p2)
    normnl=norm(nl)
    !limitsin=dsin(limit)

    if(normnl<eps)then
       print*, "ptingarc warning: Points are too close or are the same"
       print*, "P1:", p1
       print*, "P2:", p2
       print*, "Norm dif:", normnl
       return
    end if

    !Calculate plane p1, 0, p2, normal
    nl=nl/normnl
    !Project the point on the plane
    dp=dot_product(nl, p)
    !proj=p-dp*nl
    !Put the projection on the sphere
    !proj=proj/norm(proj)

    !print*, abs(dp), limit, limitcos

    !Check if point is in great circle
    if(abs(dp)>eps)then
       !print*, "Point not in great circle"
       return
    end if

    !Calculate distances to points
    d1=arclen(p,p1)
    d2=arclen(p,p2)
    d=arclen(p1,p2)

    !print*, dp
    !print*, d1
    !print*, d2
    !print*, d

    !Test if point in arc
    if(d1+d2<d+limit)then
       ptingarc=.true.
    end if

    return
  end function ptingarc


  function dist2lines(p1, v1, p2, v2)
    !-------------------------------------------------------------
    !  DIST2LINES
    !
    !  Calculates the distânce between 2 striaght lines in R3
    !
    !  Line 1 passes in point p1 and is parallel to v1
    !  Line 2 passes in point p2 and is parallel to v2
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), dimension(1:3), intent(in) :: v1
    real(r8), dimension(1:3), intent(in) :: v2

    !distance
    real(r8):: dist2lines

    !Normal - v1 x v2
    real(r8), dimension(1:3) :: nl
    real(r8):: normnl

    !p1p2
    real(r8), dimension(1:3) :: p1p2

    !Initialize variable
    dist2lines=0._r8

    !Normal to plane
    nl=cross_product(v1,v2)
    normnl=norm(nl)

    !P2-P1
    p1p2=p2-p1

    !Calculate:
    !  Based on http://pages.pacificcoast.net/~cazelais/251/distance.pdf
    ! or equivalent: http://www.easycalculation.com/analytical/shortest-distance-between-lines.php
    if( normnl < eps ) then
       dist2lines=norm(cross_product(p1p2, v1))/norm(v1)
    else
       dist2lines=abs(dot_product(p1p2, nl))/normnl
    endif

    return
  end function dist2lines

  function rot_point (p, theta)
    !-------------------------------------------------------
    !  ROT_POINT
    !
    !   This subroutine applies a rotation 'gimble like'
    ! around x, y, z-axis respectively using angles theta(x, y, z),
    ! of the point p=(x, y, z)
    !
    ! On input:
    !       p=(x,y,z) = coordinates of a point on the unit sphere.
    !       theta_? = angles of rotation in radians
    !
    ! On output:
    !       pr=(xr,yr,zr) = rotated point
    !---------------------------------------------------------
    real (r8):: p(1:3)
    real (r8):: theta(1:3)
    real (r8):: pr(1:3, 1)
    real (r8):: rot_point(1:3)
    real(r8):: R(1:3, 1:3)

    !Save rotated point
    pr(1:3,1)=p(1:3)

    !Rotate around x axis
    R(1,1)=1; R(1,2)=0;             R(1,3) =0
    R(2,1)=0; R(2,2)=dcos(theta(1)); R(2,3) = -dsin(theta(1))
    R(3,1)=0; R(3,2)=dsin(theta(1)); R(3,3) = dcos(theta(1))
    pr=matmul(R, pr)

    !Rotate around y axis
    R(1,1)=dcos(theta(2));  R(1,2)=0; R(1,3) = dsin(theta(2))
    R(2,1)=0;              R(2,2)=1; R(2,3) =0
    R(3,1)=-dsin(theta(2)); R(3,2)=0; R(3,3) = dcos(theta(2))
    pr=matmul(R, pr)

    !Rotate around z axis
    R(1,1)=dcos(theta(3));  R(1,2)=-dsin(theta(3)); R(1,3) = 0
    R(2,1)=dsin(theta(3)); R(2,2)=dcos(theta(3)); R(2,3) = 0
    R(3,1)=0;              R(3,2)=0; R(3,3) = 1
    pr=matmul(R, pr)

    rot_point(1:3)=pr(1:3, 1)
    return
  end function rot_point

  subroutine ortogonalarc(p1, p2, p, n)
    !-------------------------------------------------------------
    !  ORTOGONALARC
    !
    !  Given an arc (p1,p2), computes the bisection point p and the
    !  normal n relative to the arc that crosses perpendicularly
    !  (p1,p2) at p.
    !-------------------------------------------------------------
    !Arc points
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    !Midpoint, normal relative to arc ortogonal to (p1,p2)
    ! Obviously n is tangent to (p1,p2)
    real(r8), dimension(1:3), intent(out) :: p
    real(r8), dimension(1:3), intent(out) :: n

    p=(p1+p2)/2._r8
    p=p/norm(p)

    n=proj_vec_sphere(p-p1, p)
    n=n/norm(n)

    return
  end subroutine ortogonalarc

end module sphgeo
