module duogrid_interpolation
!========================================================================
!
! Module for Lagrange interpolation at duogrid
!
! Luan Santos 2022
! Routines based on iModel (https://github.com/pedrospeixoto/iModel)
!========================================================================

!Global constants
use constants, only: &
    griddir, &
    i4, &
    pardir, &
    r8, &
    showonscreen, &
    simulcase, &
    nbfaces, &
    erad, &
    pi, pio2, pio4, deg2rad,&
    i0, iend, &
    j0, jend, &
    n0, nend, &
    nghost, &
    hs

!Data structures
use datastruct, only: &
    cubedsphere, &
    scalar_field, &
    velocity_field, &
    matrix, &
    simulation, &
    lagrange_poly_cs

! CS mapping
use sphgeo, only: &
    inverse_equiangular_gnomonic_map2

implicit none

contains 

subroutine gethalodata(Q, L)
    !---------------------------------------------------
    ! Fill the center ghost cell values of Q that are needed
    ! for duogrid interpolation
    !--------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: Q(:,:,:)
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4) :: p, east, north, south, west


    ! --------------------- Panel 1 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 1
    north = 5
    south = 6
    east  = 2
    west  = 4
 
    ! Data of panel 1 from east
    L%halodata_east(1:nghost,n0:nend,p) = Q(i0:i0+nghost-1, n0:nend, east) ! Panel 2

    ! Data of panel 1 from  west
    L%halodata_west(1:nghost,n0:nend,p) = Q(iend-nghost+1:iend, n0:nend, west) ! Panel 4

    ! Data of panel 1 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q(n0:nend, j0:j0+nghost-1, north) ! Panel 5

    ! Data of panel 1 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q(n0:nend, jend-nghost+1:jend, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE

    ! --------------------- Panel 2 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 2
    north = 5
    south = 6
    east  = 3
    west  = 1

    ! Data of panel 2 from east
    L%halodata_east(1:nghost,n0:nend,p) = Q(i0:i0+nghost-1, n0:nend, east) ! Panel 3

    ! Data of panel 2 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q(iend-nghost+1:iend, n0:nend, west) ! Panel 1

    ! Data of panel 2 from north
    L%halodata_north(n0:nend,1:nghost,p) = transpose(Q(jend:jend-nghost+1:-1, n0:nend, north)) ! Panel 5

    ! Data of panel 2 from south
    L%halodata_south(n0:nend,1:nghost,p) = transpose(Q(jend-nghost+1:jend, nend:n0:-1,south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 3 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 3
    north = 5
    south = 6
    east  = 4
    west  = 2

    ! Data of panel 3 from east
    L%halodata_east(1:nghost,n0:nend,p) = Q(i0:i0+nghost-1, n0:nend, east) ! Panel 4

    ! Data of panel 3 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q(iend-nghost+1:iend, n0:nend, west) ! Panel 2

    ! Data of panel 3 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 5

    ! Data of panel 3 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 4 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 4
    north = 5
    south = 6
    east  = 1
    west  = 3

    ! Data of panel 4 from east
    L%halodata_east(1:nghost,n0:nend,p) = Q(i0:i0+nghost-1, n0:nend, east) ! Panel 1

    ! Data of panel 4 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q(iend-nghost+1:iend, n0:nend, west) ! Panel 3

    ! Data of panel 4 from north
    L%halodata_north(n0:nend,1:nghost,p) = transpose(Q(i0:i0+nghost-1, nend:n0:-1, north)) ! Panel 5

    ! Data of panel 4 from south
    L%halodata_south(n0:nend,1:nghost,p) = transpose(Q(i0+nghost-1:i0:-1, n0:nend, south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 5 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 5
    north = 3
    south = 1
    east  = 2
    west  = 4

    ! Data of panel 5 from east
    L%halodata_east(1:nghost,n0:nend,p) = transpose(Q(n0:nend, jend:jend-nghost+1:-1, east)) ! Panel 2

    ! Data of panel 5 from west
    L%halodata_west(1:nghost,n0:nend,p) = transpose(Q(nend:n0:-1,jend-nghost+1:jend, west)) ! Panel 4

    ! Data of panel 5 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 3

    ! Data of panel 5 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q(n0:nend, jend-nghost+1:jend, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE

    ! --------------------- Panel 6 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q, L, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 6
    north = 1
    south = 3
    east  = 2
    west  = 4

    ! Data of panel 6 from east
    L%halodata_east(1:nghost,n0:nend,p) = transpose(Q(nend:n0:-1, j0:j0+nghost-1, east)) ! Panel 2

    ! Data of panel 6 from west
    L%halodata_west(1:nghost,n0:nend,p) = transpose(Q(n0:nend,j0+nghost-1:j0:-1, west)) ! Panel 4

    ! Data of panel 6 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q(n0:nend, j0:j0+nghost-1, north) ! Panel 3

    ! Data of panel 6 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE


end subroutine gethalodata



subroutine dg_interp(Q, L)
    !---------------------------------------------------
    !   duogrid interpolation of scalar field Q
    ! (ghost cells are defined at cell centers)
    !--------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: Q(:,:,:)
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4) :: i, j, p, g, d, g2, h

    ! Fill ghost cell centers
    call gethalodata(Q, L)
    do p = 1, nbfaces
        !--------------------------------------------------------------------------
        ! East panel interpolation
        do g = 1, nghost
            h = g-1
            do j = j0-h, jend+h
                ! Store in f the support points used in Lagrange interpolation
                L%f_nodes(j,g,:) = L%halodata_east(g, L%k0(j,g):L%kend(j,g), p)
                !print*, L%f_nodes(j,g,:)
                ! Does the interpolation
                Q(iend+g,j,p) = 0.d0
                do d = 1, L%order
                    Q(iend+g,j,p) = Q(iend+g,j,p) + L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
                end do
                !print*,Q%f(iend+g,j,p)
            end do
            !read(*,*)
        end do
        !stop
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! West panel interpolation
        ! Does the interpolation
        do g = 1, nghost
            g2 = nghost-g+1
            h = g2-1
            do j = j0-h, jend+h
                ! Store in f the support points used in Lagrange interpolation
                L%f_nodes(j,g2,:)= L%halodata_west(g, L%k0(j,g2):L%kend(j,g2), p)
                Q(i0-g2,j,p) = 0.d0
                do d = 1, L%order
                    Q(i0-g2,j,p) = Q(i0-g2,j,p) + L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
                end do
            end do
        end do
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! North panel interpolation
        do g = 1, nghost
            h = g-1
            do i = i0-h, iend+h
                ! Store in f the support points used in Lagrange interpolation
                L%f_nodes(i,g,:) = L%halodata_north(L%k0(i,g):L%kend(i,g), g, p)
                ! Does the interpolation
                Q(i,jend+g,p) = 0.d0
                do d = 1, L%order
                    Q(i,jend+g,p) = Q(i,jend+g,p) + L%f_nodes(i,g,d)*L%p_nodes(i,g,d)
                end do
            end do
        end do
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! South panel interpolation
        ! Does the interpolation
        do g = 1, nghost
            g2 = nghost-g+1
            h = g2-1
            do i = i0-h, iend+h
                ! Store in f the support points used in Lagrange interpolation
                L%f_nodes(i,g2,:)= L%halodata_south(L%k0(i,g2):L%kend(i,g2), g, p)
                Q(i,j0-g2,p) = 0.d0
                do d = 1, L%order
                    Q(i,j0-g2,p) = Q(i,j0-g2,p) + L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
                end do
            end do
        end do
 
        !--------------------------------------------------------------------------
    end do

    ! Get corners values
    call gethalodata(Q, L)

    ! Interpolate reiaming values at corner
    do p = 1, nbfaces
        !--------------------------------------------------------------------------
        ! East/North panel corner interpolation
        ! East
        do g = 1, nghost
            j = jend+g
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(j,g,:) = L%halodata_east(g, L%k0(j,g):L%kend(j,g), p)
            ! Does the interpolation
            Q(iend+g,j,p) = 0.d0
            do d = 1, L%order
                Q(iend+g,j,p) = Q(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
            end do
        end do

        ! North
        do g = 1, nghost
            j = jend+g
            ! North
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(j,g,:) = L%halodata_north(L%k0(j,g):L%kend(j,g), g, p)
            ! Does the interpolation
            do d = 1, L%order
                Q(iend+g,j,p) = Q(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
            end do
        end do
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! East/South panel corner interpolation
        ! East
        do g = 1, nghost
            j = j0-g
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(j,g,:) = L%halodata_east(g, L%k0(j,g):L%kend(j,g), p)
            ! Does the interpolation
            Q(iend+g,j,p) = 0.d0
            do d = 1, L%order
                Q(iend+g,j,p) = Q(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
            end do
        end do

        ! North
        do g = 1, nghost
            g2 = nghost-g+1
            i = iend+g2
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g2,:)= L%halodata_south(L%k0(i,g2):L%kend(i,g2), g, p)
            do d = 1, L%order
                Q(i,j0-g2,p) = Q(i,j0-g2,p) + 0.5d0*L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
            end do
        end do

        !--------------------------------------------------------------------------
        ! West/North panel corner interpolation
        ! West
        do g = 1, nghost
            g2 = nghost-g+1
            j = jend+g2
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(j,g2,:)= L%halodata_west(g, L%k0(j,g2):L%kend(j,g2), p)
            Q(i0-g2,j,p) = 0.d0
            do d = 1, L%order
                Q(i0-g2,j,p) = Q(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
            end do
        end do

        ! North
        do g = 1, nghost
            i = i0-g
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g,:) = L%halodata_north(L%k0(i,g):L%kend(i,g), g, p)
            ! Does the interpolation
            do d = 1, L%order
                Q(i,jend+g,p) = Q(i,jend+g,p) + 0.5d0*L%f_nodes(i,g,d)*L%p_nodes(i,g,d)
            end do
        end do
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! West/South panel corner interpolation
        ! West
        do g = 1, nghost
            g2 = nghost-g+1
            j = i0-g2
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(j,g2,:)= L%halodata_west(g, L%k0(j,g2):L%kend(j,g2), p)
            Q(i0-g2,j,p) = 0.d0
            do d = 1, L%order
                Q(i0-g2,j,p) = Q(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
            end do
        end do

        ! South
        ! Does the interpolation
        do g = 1, nghost
            g2 = nghost-g+1
            i = i0-g2
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g2,:)= L%halodata_south(L%k0(i,g2):L%kend(i,g2), g, p)
            do d = 1, L%order
                Q(i,j0-g2,p) = Q(i,j0-g2,p) + 0.5d0*L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
            end do
        end do
        !--------------------------------------------------------------------------
    end do

end subroutine dg_interp


subroutine interp_C2Aduogrid(u_pu, v_pv, ulon_pc, vlat_pc, u_pc, v_pc, L, A, Ainv)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at C grid (contravariant/covariant)
    ! to the A grid (contravariant/covariant) at ghost cell values at centers
    !---------------------------------------------------
    type(matrix), allocatable, intent(in) :: A(:,:,:), Ainv(:,:,:) ! conversion matrix
    real(kind=8), allocatable, intent(inout) :: u_pu(:,:,:), v_pv(:,:,:)
    real(kind=8), allocatable, intent(inout) :: u_pc(:,:,:), v_pc(:,:,:)
    real(kind=8), allocatable, intent(inout) :: ulon_pc(:,:,:), vlat_pc(:,:,:)
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4):: i, j, p, h
    real(kind=8) :: a1, a2, a3, a4
    real(kind=8) :: b1, b2, b3, b4
    real(kind=8) :: c1, c2

    h = hs-1
    ! cubic interpolation coeffs
    a1 =  5.d0/16.d0
    a2 = 15.d0/16.d0
    a3 = -5.d0/16.d0
    a4 =  1.d0/16.d0

    b1 = -1.d0/16.d0
    b2 =  9.d0/16.d0
    b3 =  9.d0/16.d0
    b4 = -1.d0/16.d0
 
    c1 =  9.d0/16.d0
    c2 = -1.d0/16.d0

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend) &
    !$OMP SHARED(a1, a2, a3, a4) &
    !$OMP SHARED(b1, b2, b3, b4) &
    !$OMP SHARED(h, hs) &
    !$OMP SHARED(u_pc, v_pc, ulon_pc, vlat_pc, u_pu, v_pv, A)
    ! west boundary
    u_pc(i0,j0:jend,:) = a1*u_pu(i0  ,j0:jend,:) &
                       + a2*u_pu(i0+1,j0:jend,:) &
                       + a3*u_pu(i0+2,j0:jend,:) &
                       + a4*u_pu(i0+3,j0:jend,:)

    u_pc(i0+1:i0+h,j0:jend,:) = b1*u_pu(i0:i0-1+h  ,j0:jend,:) &
                              + b2*u_pu(i0+1:i0+h  ,j0:jend,:) &
                              + b3*u_pu(i0+2:i0+1+h,j0:jend,:) &
                              + b4*u_pu(i0+3:i0+2+h,j0:jend,:)

    v_pc(i0:i0+h,j0+hs:jend-h,:) = b1*v_pv(i0:j0+h,j0+hs-1:jend-h-1,:) &
                                 + b2*v_pv(i0:i0+h,j0+hs:jend-h    ,:) &
                                 + b3*v_pv(i0:i0+h,j0+hs+1:jend-h+1,:) &
                                 + b4*v_pv(i0:i0+h,j0+hs+2:jend-h+2,:)

    ! east boundary
    u_pc(iend,j0:jend,:) = a4*u_pu(iend-2,j0:jend,:) &
                         + a3*u_pu(iend-1,j0:jend,:) &
                         + a2*u_pu(iend,j0:jend,:) &
                         + a1*u_pu(iend+1,j0:jend,:)

    u_pc(iend-hs:iend-1,j0:jend,:) = b4*u_pu(iend-hs-1:iend-2,j0:jend,:) &
                                   + b3*u_pu(iend-hs+0:iend-1,j0:jend,:) &
                                   + b2*u_pu(iend-hs+1:iend+0,j0:jend,:) &
                                   + b1*u_pu(iend-hs+2:iend+1,j0:jend,:)

    v_pc(iend-hs:iend,j0+hs:jend-h,:) = b1*v_pv(iend-hs:iend,j0+hs-1:jend-h-1,:) &
                                      + b2*v_pv(iend-hs:iend,j0+hs:jend-h    ,:) &
                                      + b3*v_pv(iend-hs:iend,j0+hs+1:jend-h+1,:) &
                                      + b4*v_pv(iend-hs:iend,j0+hs+2:jend-h+2,:)

    ! south boundary
    v_pc(i0:iend,j0,:) = a1*v_pv(i0:iend,j0  ,:) &
                       + a2*v_pv(i0:iend,j0+1,:) &
                       + a3*v_pv(i0:iend,j0+2,:) &
                       + a4*v_pv(i0:iend,j0+3,:)

    v_pc(i0:iend,j0+1:j0+h,:) = b1*v_pv(i0:iend,j0:j0-1+h  ,:) &
                              + b2*v_pv(i0:iend,j0+1:j0+h  ,:) &
                              + b3*v_pv(i0:iend,j0+2:j0+1+h,:) &
                              + b4*v_pv(i0:iend,j0+3:j0+2+h,:)

    u_pc(i0+hs:iend-h,j0:j0+h,:) = b1*u_pu(i0+hs-1:iend-h-1,j0:j0+h,:) &
                                 + b2*u_pu(i0+hs:iend-h    ,j0:j0+h,:) &
                                 + b3*u_pu(i0+hs+1:iend-h+1,j0:j0+h,:) &
                                 + b4*u_pu(i0+hs+2:iend-h+2,j0:j0+h,:) 

    ! north boundary
    v_pc(i0:iend,jend,:) = a4*v_pv(i0:iend,jend-2,:) &
                         + a3*v_pv(i0:iend,jend-1,:) &
                         + a2*v_pv(i0:iend,jend-0,:) &
                         + a1*v_pv(i0:iend,jend+1,:)

    v_pc(i0:iend,jend-hs:jend-1,:) = b4*v_pv(i0:iend,jend-hs-1:jend-2,:) &
                                   + b3*v_pv(i0:iend,jend-hs+0:jend-1,:) &
                                   + b2*v_pv(i0:iend,jend-hs+1:jend,:) &
                                   + b1*v_pv(i0:iend,jend-hs+2:jend+1,:)

    u_pc(i0+hs:iend-h,jend-hs:jend,:) = b1*u_pu(i0+hs-1:iend-h-1,jend-hs:jend,:) &
                                      + b2*u_pu(i0+hs:iend-h    ,jend-hs:jend,:) &
                                      + b3*u_pu(i0+hs+1:iend-h+1,jend-hs:jend,:) &
                                      + b4*u_pu(i0+hs+2:iend-h+2,jend-hs:jend,:)

    ! Convert from contravariant to latlon
    ulon_pc(i0:i0+h,j0:jend,:) = &
    u_pc(i0:i0+h,j0:jend,:)*A(i0:i0+h,j0:jend,:)%M(1,1)  + &
    v_pc(i0:i0+h,j0:jend,:)*A(i0:i0+h,j0:jend,:)%M(1,2) 
    vlat_pc(i0:i0+h,j0:jend,:) = &
    u_pc(i0:i0+h,j0:jend,:)*A(i0:i0+h,j0:jend,:)%M(2,1) + &
    v_pc(i0:i0+h,j0:jend,:)*A(i0:i0+h,j0:jend,:)%M(2,2) 

    ulon_pc(i0:iend,j0:j0+h,:) = &
    u_pc(i0:iend,j0:j0+h,:)*A(i0:iend,j0:j0+h,:)%M(1,1)  + &
    v_pc(i0:iend,j0:j0+h,:)*A(i0:iend,j0:j0+h,:)%M(1,2) 
    vlat_pc(i0:iend,j0:j0+h,:) = &
    u_pc(i0:iend,j0:j0+h,:)*A(i0:iend,j0:j0+h,:)%M(2,1) + &
    v_pc(i0:iend,j0:j0+h,:)*A(i0:iend,j0:j0+h,:)%M(2,2) 

    ulon_pc(iend-h:iend,j0:jend,:) = &
    u_pc(iend-h:iend,j0:jend,:)*A(iend-h:iend,j0:jend,:)%M(1,1)+&
    v_pc(iend-h:iend,j0:jend,:)*A(iend-h:iend,j0:jend,:)%M(1,2) 
    vlat_pc(iend-h:iend,j0:jend,:) = &
    u_pc(iend-h:iend,j0:jend,:)*A(iend-h:iend,j0:jend,:)%M(2,1)+&
    v_pc(iend-h:iend,j0:jend,:)*A(iend-h:iend,j0:jend,:)%M(2,2) 

    ulon_pc(i0:iend,jend-h:jend,:) = &
    u_pc(i0:iend,jend-h:jend,:)*A(i0:iend,jend-h:jend,:)%M(1,1)+&
    v_pc(i0:iend,jend-h:jend,:)*A(i0:iend,jend-h:jend,:)%M(1,2) 
    vlat_pc(i0:iend,jend-h:jend,:) = &
    u_pc(i0:iend,jend-h:jend,:)*A(i0:iend,jend-h:jend,:)%M(2,1)+&
    v_pc(i0:iend,jend-h:jend,:)*A(i0:iend,jend-h:jend,:)%M(2,2) 


    !$OMP END PARALLEL WORKSHARE

    ! Interpolate latlon to the ghost cell centers
    call dg_interp(ulon_pc, L)
    call dg_interp(vlat_pc, L)

    ! Convert from latlon to contravariant
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
    !$OMP SHARED(h, hs) &
    !$OMP SHARED(u_pc, v_pc, ulon_pc, vlat_pc, u_pu, v_pv, Ainv)
    u_pc(i0-hs:i0-1,n0:nend,:) = &
    ulon_pc(i0-hs:i0-1,n0:nend,:)*Ainv(i0-hs:i0-1,n0:nend,:)%M(1,1) + &
    vlat_pc(i0-hs:i0-1,n0:nend,:)*Ainv(i0-hs:i0-1,n0:nend,:)%M(1,2) 
    v_pc(i0-hs:i0-1,n0:nend,:) =&
    ulon_pc(i0-hs:i0-1,n0:nend,:)*Ainv(i0-hs:i0-1,n0:nend,:)%M(2,1) + &
    vlat_pc(i0-hs:i0-1,n0:nend,:)*Ainv(i0-hs:i0-1,n0:nend,:)%M(2,2) 



    u_pc(iend+1:iend+hs,n0:nend,:) = &
    ulon_pc(iend+1:iend+hs,n0:nend,:)*Ainv(iend+1:iend+hs,n0:nend,:)%M(1,1) + &
    vlat_pc(iend+1:iend+hs,n0:nend,:)*Ainv(iend+1:iend+hs,n0:nend,:)%M(1,2) 
    v_pc(iend+1:iend+hs,n0:nend,:) =&
    ulon_pc(iend+1:iend+hs,n0:nend,:)*Ainv(iend+1:iend+hs,n0:nend,:)%M(2,1) + &
    vlat_pc(iend+1:iend+hs,n0:nend,:)*Ainv(iend+1:iend+hs,n0:nend,:)%M(2,2) 



    u_pc(i0:iend,i0-hs:i0-1,:) = &
    ulon_pc(i0:iend,i0-hs:i0-1,:)*Ainv(i0:iend,i0-hs:i0-1,:)%M(1,1) + &
    vlat_pc(i0:iend,i0-hs:i0-1,:)*Ainv(i0:iend,i0-hs:i0-1,:)%M(1,2) 
    v_pc(i0:iend,i0-hs:i0-1,:) = &
    ulon_pc(i0:iend,i0-hs:i0-1,:)*Ainv(i0:iend,i0-hs:i0-1,:)%M(2,1) + &
    vlat_pc(i0:iend,i0-hs:i0-1,:)*Ainv(i0:iend,i0-hs:i0-1,:)%M(2,2) 



    u_pc(i0:iend,iend+1:iend+hs,:) = &
    ulon_pc(i0:iend,iend+1:iend+hs,:)*Ainv(i0:iend,iend+1:iend+hs,:)%M(1,1) + &
    vlat_pc(i0:iend,iend+1:iend+hs,:)*Ainv(i0:iend,iend+1:iend+hs,:)%M(1,2) 
    v_pc(i0:iend,iend+1:iend+hs,:) = &
    ulon_pc(i0:iend,iend+1:iend+hs,:)*Ainv(i0:iend,iend+1:iend+hs,:)%M(2,1) + &
    vlat_pc(i0:iend,iend+1:iend+hs,:)*Ainv(i0:iend,iend+1:iend+hs,:)%M(2,2) 
    !$OMP END PARALLEL WORKSHARE


end subroutine interp_C2Aduogrid


subroutine interp_C2Agrid(u_pu, v_pv, u_pc, v_pc, id)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at C grid (contra/covar/latlon)
    ! to the A grid (contra/covar/latlon),
    ! but not including its ghost cell values at centers
    !---------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: u_pu(:,:,:), v_pv(:,:,:)
    real(kind=8), allocatable, intent(inout) :: u_pc(:,:,:), v_pc(:,:,:)
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2, a3, a4
    real(kind=8) :: b1, b2, b3, b4

    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(u_pu, v_pv, u_pc, v_pc)
        u_pc(i0:iend,j0:jend,:) = &
        (u_pu(i0:iend,j0:jend,:) + u_pu(i0+1:iend+1,j0:jend,:))*0.5d0
        v_pc(i0:iend,j0:jend,:) = &
        (v_pv(i0:iend,j0:jend,:) + v_pv(i0:iend,j0+1:jend+1,:))*0.5d0
        !$OMP END PARALLEL WORKSHARE

    elseif(id == 3) then
        ! cubic interpolation coeffs
        a1 =  5.d0/16.d0
        a2 = 15.d0/16.d0
        a3 = -5.d0/16.d0
        a4 =  1.d0/16.d0

        b1 = -1.d0/16.d0
        b2 =  9.d0/16.d0
        b3 =  9.d0/16.d0
        b4 = -1.d0/16.d0
     
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(a1, a2, a3, a4) &
        !$OMP SHARED(b1, b2, b3, b4) &
        !$OMP SHARED(u_pu, v_pv, u_pc, v_pc)
        ! west boundary
        u_pc(i0,j0:jend,:) = a1*u_pu(i0,j0:jend,:) &
                           + a2*u_pu(i0+1,j0:jend,:) &
                           + a3*u_pu(i0+2,j0:jend,:) &
                           + a4*u_pu(i0+3,j0:jend,:)

        ! east boundary
        u_pc(iend,j0:jend,:) = a4*u_pu(iend-2,j0:jend,:) &
                             + a3*u_pu(iend-1,j0:jend,:) &
                             + a2*u_pu(iend,j0:jend,:) &
                             + a1*u_pu(iend+1,j0:jend,:)

        ! south boundary
        v_pc(i0:iend,j0,:) = a1*v_pv(i0:iend,j0,:) &
                           + a2*v_pv(i0:iend,j0+1,:) &
                           + a3*v_pv(i0:iend,j0+2,:) &
                           + a4*v_pv(i0:iend,j0+3,:)

        ! north boundary
        v_pc(i0:iend,jend,:) = a4*v_pv(i0:iend,jend-2,:) &
                             + a3*v_pv(i0:iend,jend-1,:) &
                             + a2*v_pv(i0:iend,jend-0,:) &
                             + a1*v_pv(i0:iend,jend+1,:)

        ! remaing cells
        u_pc(i0+1:iend-1,j0:jend,:) = b1*u_pu(i0:iend-2,j0:jend,:) &
                                    + b2*u_pu(i0+1:iend-1,j0:jend,:) &
                                    + b3*u_pu(i0+2:iend  ,j0:jend,:) &
                                    + b4*u_pu(i0+3:iend+1,j0:jend,:)

        v_pc(i0:iend,j0+1:jend-1,:) = b1*v_pv(i0:iend,j0:jend-2,:) &
                                    + b2*v_pv(i0:iend,j0+1:jend-1,:) &
                                    + b3*v_pv(i0:iend,j0+2:jend  ,:) &
                                    + b4*v_pv(i0:iend,j0+3:jend+1,:)

        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in interp_C2Agrid: invalid id, ', id
        stop
    end if
end subroutine interp_C2Agrid


subroutine interp_A2Cduogrid(u_pu, u_pv, v_pu, v_pv, u_pc, v_pc)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at A grid (contravariant/covariant)
    ! to the C grid (contravariant/covariant) ghost cells only
    real(kind=8), allocatable, intent(inout) :: u_pu(:,:,:), v_pu(:,:,:)
    real(kind=8), allocatable, intent(inout) :: u_pv(:,:,:), v_pv(:,:,:)
    real(kind=8), allocatable, intent(inout) :: u_pc(:,:,:), v_pc(:,:,:)
    integer(i4):: i, j, p, h
    real(kind=8) :: c1, c2
    
    ! cubic interpolation coeffs
    c1 =  9.d0/16.d0
    c2 = -1.d0/16.d0

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
    !$OMP SHARED(c1, c2) &
    !$OMP SHARED(u_pu, u_pv, v_pu, v_pv, u_pc, v_pc)
    ! Let us interpolate the ghost cell edges - cubic interpolation
    ! Panel from south
    u_pu(i0:iend+1,n0:j0-1,:) = &
    c1*(u_pc(i0:iend+1  ,n0:j0-1,:) + u_pc(i0-1:iend  ,n0:j0-1,:)) + &
    c2*(u_pc(i0+1:iend+2,n0:j0-1,:) + u_pc(i0-2:iend-1,n0:j0-1,:))

    v_pu(i0:iend+1,n0:j0-1,:) = &
    c1*(v_pc(i0:iend+1  ,n0:j0-1,:) + v_pc(i0-1:iend  ,n0:j0-1,:)) + &
    c2*(v_pc(i0+1:iend+2,n0:j0-1,:) + v_pc(i0-2:iend-1,n0:j0-1,:))

    ! Panel from north
    u_pu(i0:iend+1, jend+1:nend,:) = &
    c1*(u_pc(i0:iend+1  ,jend+1:nend,:) + u_pc(i0-1:iend  ,jend+1:nend,:)) + &
    c2*(u_pc(i0+1:iend+2,jend+1:nend,:) + u_pc(i0-2:iend-1,jend+1:nend,:))

    v_pu(i0:iend+1,jend+1:nend,:) = &
    c1*(v_pc(i0:iend+1  ,jend+1:nend,:) + v_pc(i0-1:iend  ,jend+1:nend,:)) + &
    c2*(v_pc(i0+1:iend+2,jend+1:nend,:) + v_pc(i0-2:iend-1,jend+1:nend,:))

    ! Panel from west
    u_pv(n0:i0-1,j0:jend+1,:) = &
    c1*(u_pc(n0:i0-1,j0:jend+1,  :) + u_pc(n0:i0-1,j0-1:jend,  :)) + &
    c2*(u_pc(n0:i0-1,j0+1:jend+2,:) + u_pc(n0:i0-1,j0-2:jend-1,:))

    v_pv(n0:i0-1,j0:jend+1,:) = &
    c1*(v_pc(n0:i0-1,j0:jend+1  ,:) + v_pc(n0:i0-1,j0-1:jend  ,:)) + &
    c2*(v_pc(n0:i0-1,j0+1:jend+2,:) + v_pc(n0:i0-1,j0-2:jend-1,:))

    ! Panel from east 
    u_pv(iend+1:nend, j0:jend+1, :) = &
    c1*(u_pc(iend+1:nend,j0:jend+1  ,:) + u_pc(iend+1:nend,j0-1:jend  ,:)) + &
    c2*(u_pc(iend+1:nend,j0+1:jend+2,:) + u_pc(iend+1:nend,j0-2:jend-1,:))

    v_pv(iend+1:nend, j0:jend+1, :) = &
    c1*(v_pc(iend+1:nend,j0:jend+1  ,:) + v_pc(iend+1:nend,j0-1:jend  ,:)) + &
    c2*(v_pc(iend+1:nend,j0+1:jend+2,:) + v_pc(iend+1:nend,j0-2:jend-1,:))


    !-----------------------------------------------------------------------------------
    ! Interpolation needed for RK2 departure point scheme
    ! Panel from west
    u_pu(i0-1,n0:nend,:) = &
    c1*(u_pc(i0-2,n0:nend,:) + u_pc(i0-1,n0:nend,:)) + &
    c2*(u_pc(i0  ,n0:nend,:) + u_pc(i0-3,n0:nend,:))

    v_pu(i0-1,n0:nend,:) = &
    c1*(v_pc(i0-2,n0:nend,:) + v_pc(i0-1,n0:nend,:)) + &
    c2*(v_pc(i0  ,n0:nend,:) + v_pc(i0-3,n0:nend,:))

    ! Panel from east
    u_pu(iend+2,n0:nend,:) = &
    c1*(u_pc(iend+1,n0:nend,:) + u_pc(iend+2,n0:nend,:)) + &
    c2*(u_pc(iend  ,n0:nend,:) + u_pc(iend+3,n0:nend,:))

    v_pu(iend+2,:,:) = &
    c1*(v_pc(iend+1,n0:nend,:) + v_pc(iend+2,n0:nend,:)) + &
    c2*(v_pc(iend  ,n0:nend,:) + v_pc(iend+3,n0:nend,:))
 
    ! Panel from south
    u_pv(n0:nend,j0-1,:) = &
    c1*(u_pc(n0:nend,j0-2,:) + u_pc(n0:nend,j0-1,:)) + &
    c2*(u_pc(n0:nend,j0  ,:) + u_pc(n0:nend,j0-3,:))

    v_pv(n0:nend,j0-1,:) = &
    c1*(v_pc(n0:nend,j0-2,:) + v_pc(n0:nend,j0-1,:)) + &
    c2*(v_pc(n0:nend,j0  ,:) + v_pc(n0:nend,j0-3,:))

    ! Panel from east
    u_pv(n0:nend,jend+2,:) = &
    c1*(u_pc(n0:nend,jend+1,:) + u_pc(n0:nend,jend+2,:)) + &
    c2*(u_pc(n0:nend,jend  ,:) + u_pc(n0:nend,jend+3,:))

    v_pv(n0:nend,jend+2,:) = &
    c1*(v_pc(n0:nend,jend+1,:) + v_pc(n0:nend,jend+2,:)) + &
    c2*(v_pc(n0:nend,jend  ,:) + v_pc(n0:nend,jend+3,:))
    !$OMP END PARALLEL WORKSHARE

end subroutine interp_A2Cduogrid


subroutine interp_A2Cgrid(u_pu, v_pv, u_pc, v_pc, id)
    !---------------------------------------------------
    ! interpolation of vector field (contra/covar/latlon) given at A grid 
    ! (including ghost cells) to the C grid (contra/covar/latlon) inner cells
    real(kind=8), allocatable, intent(inout) :: u_pu(:,:,:), v_pv(:,:,:)
    real(kind=8), allocatable, intent(inout) :: u_pc(:,:,:), v_pc(:,:,:)
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p, h
    real(kind=8) :: c1, c2
    
    ! cubic interpolation coeffs
    c1 =  9.d0/16.d0
    c2 = -1.d0/16.d0

    ! Now, let us interpolate from the A grid to the C-grid points
    ! that are not ghost cells
    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
        !$OMP SHARED(u_pu, v_pv, u_pc, v_pc)
        u_pu(i0-2:iend+3,n0:nend,:) = &
        (u_pc(i0-3:iend+2,n0:nend,:) + u_pc(i0-2:iend+3,n0:nend,:))*0.5d0
        v_pv(n0:nend,j0-2:jend+3,:) = &
        (v_pc(n0:nend,j0-3:jend+2,:) + v_pc(n0:nend,j0-2:jend+3,:))*0.5d0
        !$OMP END PARALLEL WORKSHARE

    elseif(id == 3) then
        ! cubic interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
        !$OMP SHARED(u_pu, v_pv, u_pc, v_pc, c1, c2)
        u_pu(i0-2:iend+3,n0:nend,:) = &
        c1*(u_pc(i0-3:iend+2,n0:nend,:) + u_pc(i0-2:iend+3 ,n0:nend,:)) + &
        c2*(u_pc(i0-4:iend+1,n0:nend,:) + u_pc(i0-1:iend+4 ,n0:nend,:))

        v_pv(n0:nend,j0-2:jend+3,:) = &
        c1*(v_pc(n0:nend,j0-3:jend+2,:) + v_pc(n0:nend,j0-2:jend+3  ,:)) + &
        c2*(v_pc(n0:nend,j0-4:jend+1,:) + v_pc(n0:nend,j0-1:jend+4,:))
        !$OMP END PARALLEL WORKSHARE

    else
        print*, 'ERROR in interp_A2Cgrid: invalid id, ', id
        stop
    end if
end subroutine interp_A2Cgrid




subroutine interp_D2Aduogrid(U_pu, U_pv, U_pc, L, mesh)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at D grid (covariant)
    ! to the A grid (covariant) ghost cell values at centers
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(velocity_field), intent(inout) :: U_pu, U_pv, U_pc
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4):: i, j, p, h
    real(kind=8) :: a1, a2, a3, a4
    real(kind=8) :: b1, b2, b3, b4
    real(kind=8) :: c1, c2

    h = hs-1
    ! cubic interpolation coeffs
    a1 =  5.d0/16.d0
    a2 = 15.d0/16.d0
    a3 = -5.d0/16.d0
    a4 =  1.d0/16.d0

    b1 = -1.d0/16.d0
    b2 =  9.d0/16.d0
    b3 =  9.d0/16.d0
    b4 = -1.d0/16.d0
 
    c1 =  9.d0/16.d0
    c2 = -1.d0/16.d0

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend) &
    !$OMP SHARED(a1, a2, a3, a4) &
    !$OMP SHARED(b1, b2, b3, b4) &
    !$OMP SHARED(h, hs) &
    !$OMP SHARED(U_pc, U_pu, U_pv, mesh)
    ! west boundary
    U_pc%vcovari%f(i0,j0:jend,:) = &
      a1*U_pu%vcovari%f(i0,j0:jend,:) &
    + a2*U_pu%vcovari%f(i0+1,j0:jend,:) &
    + a3*U_pu%vcovari%f(i0+2,j0:jend,:) &
    + a4*U_pu%vcovari%f(i0+3,j0:jend,:)

    U_pc%vcovari%f(i0+1:i0+h,j0:jend,:) = & 
      b1*U_pu%vcovari%f(i0:i0-1+h,j0:jend,:) &
    + b2*U_pu%vcovari%f(i0+1:i0+h,j0:jend,:) &
    + b3*U_pu%vcovari%f(i0+2:i0+1+h,j0:jend,:) &
    + b4*U_pu%vcovari%f(i0+3:i0+2+h,j0:jend,:)

    U_pc%ucovari%f(i0:i0+h,j0+hs:jend-h,:) = &
      b1*U_pv%ucovari%f(i0:j0+h,j0+hs-1:jend-h-1,:) &
    + b2*U_pv%ucovari%f(i0:i0+h,j0+hs:jend-h,:) &
    + b3*U_pv%ucovari%f(i0:i0+h,j0+hs+1:jend-h+1,:) &
    + b4*U_pv%ucovari%f(i0:i0+h,j0+hs+2:jend-h+2,:)

    ! east boundary
    U_pc%vcovari%f(iend,j0:jend,:) = &
      a4*U_pu%vcovari%f(iend-2,j0:jend,:) &
    + a3*U_pu%vcovari%f(iend-1,j0:jend,:) &
    + a2*U_pu%vcovari%f(iend,j0:jend,:) &
    + a1*U_pu%vcovari%f(iend+1,j0:jend,:)

    U_pc%vcovari%f(iend-hs:iend-1,j0:jend,:) = &
      b4*U_pu%vcovari%f(iend-hs-1:iend-2,j0:jend,:) &
    + b3*U_pu%vcovari%f(iend-hs+0:iend-1,j0:jend,:) &
    + b2*U_pu%vcovari%f(iend-hs+1:iend+0,j0:jend,:) &
    + b1*U_pu%vcovari%f(iend-hs+2:iend+1,j0:jend,:)

    U_pc%ucovari%f(iend-hs:iend,j0+hs:jend-h,:) = &
      b1*U_pv%ucovari%f(iend-hs:iend,j0+hs-1:jend-h-1,:) &
    + b2*U_pv%ucovari%f(iend-hs:iend,j0+hs:jend-h,:) &
    + b3*U_pv%ucovari%f(iend-hs:iend,j0+hs+1:jend-h+1,:) &
    + b4*U_pv%ucovari%f(iend-hs:iend,j0+hs+2:jend-h+2,:)

    ! south boundary
    U_pc%ucovari%f(i0:iend,j0,:) = &
      a1*U_pv%ucovari%f(i0:iend,j0,:) &
    + a2*U_pv%ucovari%f(i0:iend,j0+1,:) &
    + a3*U_pv%ucovari%f(i0:iend,j0+2,:) &
    + a4*U_pv%ucovari%f(i0:iend,j0+3,:)

    U_pc%ucovari%f(i0:iend,j0+1:j0+h,:) = &
      b1*U_pv%ucovari%f(i0:iend,j0:j0-1+h,:) &
    + b2*U_pv%ucovari%f(i0:iend,j0+1:j0+h,:) &
    + b3*U_pv%ucovari%f(i0:iend,j0+2:j0+1+h,:) &
    + b4*U_pv%ucovari%f(i0:iend,j0+3:j0+2+h,:)

    U_pc%vcovari%f(i0+hs:iend-h,j0:j0+h,:) = &
      b1*U_pu%vcovari%f(i0+hs-1:iend-h-1,j0:j0+h,:) &
    + b2*U_pu%vcovari%f(i0+hs:iend-h    ,j0:j0+h,:) &
    + b3*U_pu%vcovari%f(i0+hs+1:iend-h+1,j0:j0+h,:) &
    + b4*U_pu%vcovari%f(i0+hs+2:iend-h+2,j0:j0+h,:) 

    ! north boundary
    U_pc%ucovari%f(i0:iend,jend,:) = &
    a4*U_pv%ucovari%f(i0:iend,jend-2,:) + &
    a3*U_pv%ucovari%f(i0:iend,jend-1,:) + &
    a2*U_pv%ucovari%f(i0:iend,jend-0,:) + &
    a1*U_pv%ucovari%f(i0:iend,jend+1,:)

    U_pc%ucovari%f(i0:iend,jend-hs:jend-1,:) = &
      b4*U_pv%ucovari%f(i0:iend,jend-hs-1:jend-2,:) &
    + b3*U_pv%ucovari%f(i0:iend,jend-hs+0:jend-1,:) &
    + b2*U_pv%ucovari%f(i0:iend,jend-hs+1:jend,:) &
    + b1*U_pv%ucovari%f(i0:iend,jend-hs+2:jend+1,:)

    U_pc%vcovari%f(i0+hs:iend-h,jend-hs:jend,:) = & 
      b1*U_pu%vcovari%f(i0+hs-1:iend-h-1,jend-hs:jend,:) &
    + b2*U_pu%vcovari%f(i0+hs:iend-h    ,jend-hs:jend,:) &
    + b3*U_pu%vcovari%f(i0+hs+1:iend-h+1,jend-hs:jend,:) &
    + b4*U_pu%vcovari%f(i0+hs+2:iend-h+2,jend-hs:jend,:)


    ! Convert from contravariant to latlon
    U_pc%u%f(i0:i0+h,j0:jend,:) = &
    U_pc%ucovari%f(i0:i0+h,j0:jend,:)*mesh%covari2ll_pc(i0:i0+h,j0:jend,:)%M(1,1)  + &
    U_pc%vcovari%f(i0:i0+h,j0:jend,:)*mesh%covari2ll_pc(i0:i0+h,j0:jend,:)%M(1,2) 
    U_pc%v%f(i0:i0+h,j0:jend,:) = &
    U_pc%ucovari%f(i0:i0+h,j0:jend,:)*mesh%covari2ll_pc(i0:i0+h,j0:jend,:)%M(2,1) + &
    U_pc%vcovari%f(i0:i0+h,j0:jend,:)*mesh%covari2ll_pc(i0:i0+h,j0:jend,:)%M(2,2) 

    U_pc%u%f(i0:iend,j0:j0+h,:) = &
    U_pc%ucovari%f(i0:iend,j0:j0+h,:)*mesh%covari2ll_pc(i0:iend,j0:j0+h,:)%M(1,1)  + &
    U_pc%vcovari%f(i0:iend,j0:j0+h,:)*mesh%covari2ll_pc(i0:iend,j0:j0+h,:)%M(1,2) 
    U_pc%v%f(i0:iend,j0:j0+h,:) = &
    U_pc%ucovari%f(i0:iend,j0:j0+h,:)*mesh%covari2ll_pc(i0:iend,j0:j0+h,:)%M(2,1) + &
    U_pc%vcovari%f(i0:iend,j0:j0+h,:)*mesh%covari2ll_pc(i0:iend,j0:j0+h,:)%M(2,2) 

    U_pc%u%f(iend-h:iend,j0:jend,:) = &
    U_pc%ucovari%f(iend-h:iend,j0:jend,:)*mesh%covari2ll_pc(iend-h:iend,j0:jend,:)%M(1,1)+&
    U_pc%vcovari%f(iend-h:iend,j0:jend,:)*mesh%covari2ll_pc(iend-h:iend,j0:jend,:)%M(1,2) 
    U_pc%v%f(iend-h:iend,j0:jend,:) = &
    U_pc%ucovari%f(iend-h:iend,j0:jend,:)*mesh%covari2ll_pc(iend-h:iend,j0:jend,:)%M(2,1)+&
    U_pc%vcovari%f(iend-h:iend,j0:jend,:)*mesh%covari2ll_pc(iend-h:iend,j0:jend,:)%M(2,2) 

    U_pc%u%f(i0:iend,jend-h:jend,:) = &
    U_pc%ucovari%f(i0:iend,jend-h:jend,:)*mesh%covari2ll_pc(i0:iend,jend-h:jend,:)%M(1,1)+&
    U_pc%vcovari%f(i0:iend,jend-h:jend,:)*mesh%covari2ll_pc(i0:iend,jend-h:jend,:)%M(1,2) 
    U_pc%v%f(i0:iend,jend-h:jend,:) = &
    U_pc%ucovari%f(i0:iend,jend-h:jend,:)*mesh%covari2ll_pc(i0:iend,jend-h:jend,:)%M(2,1)+&
    U_pc%vcovari%f(i0:iend,jend-h:jend,:)*mesh%covari2ll_pc(i0:iend,jend-h:jend,:)%M(2,2) 
    !$OMP END PARALLEL WORKSHARE

    ! Interpolate latlon to the ghost cell centers
    call dg_interp(U_pc%u%f, L)
    call dg_interp(U_pc%v%f, L)

    ! Convert from latlon to contravariant
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
    !$OMP SHARED(h, hs) &
    !$OMP SHARED(U_pc, U_pu, U_pv, mesh)
    U_pc%ucovari%f(i0-hs:i0-1,n0:nend,:) = &
    U_pc%u%f(i0-hs:i0-1,n0:nend,:)*mesh%ll2covari_pc(i0-hs:i0-1,n0:nend,:)%M(1,1) + &
    U_pc%v%f(i0-hs:i0-1,n0:nend,:)*mesh%ll2covari_pc(i0-hs:i0-1,n0:nend,:)%M(1,2) 
    U_pc%vcovari%f(i0-hs:i0-1,n0:nend,:) =&
    U_pc%u%f(i0-hs:i0-1,n0:nend,:)*mesh%ll2covari_pc(i0-hs:i0-1,n0:nend,:)%M(2,1) + &
    U_pc%v%f(i0-hs:i0-1,n0:nend,:)*mesh%ll2covari_pc(i0-hs:i0-1,n0:nend,:)%M(2,2) 



    U_pc%ucovari%f(iend+1:iend+hs,n0:nend,:) = &
    U_pc%u%f(iend+1:iend+hs,n0:nend,:)*mesh%ll2covari_pc(iend+1:iend+hs,n0:nend,:)%M(1,1) + &
    U_pc%v%f(iend+1:iend+hs,n0:nend,:)*mesh%ll2covari_pc(iend+1:iend+hs,n0:nend,:)%M(1,2) 
    U_pc%vcovari%f(iend+1:iend+hs,n0:nend,:) =&
    U_pc%u%f(iend+1:iend+hs,n0:nend,:)*mesh%ll2covari_pc(iend+1:iend+hs,n0:nend,:)%M(2,1) + &
    U_pc%v%f(iend+1:iend+hs,n0:nend,:)*mesh%ll2covari_pc(iend+1:iend+hs,n0:nend,:)%M(2,2) 



    U_pc%ucovari%f(i0:iend,i0-hs:i0-1,:) = &
    U_pc%u%f(i0:iend,i0-hs:i0-1,:)*mesh%ll2covari_pc(i0:iend,i0-hs:i0-1,:)%M(1,1) + &
    U_pc%v%f(i0:iend,i0-hs:i0-1,:)*mesh%ll2covari_pc(i0:iend,i0-hs:i0-1,:)%M(1,2) 
    U_pc%vcovari%f(i0:iend,i0-hs:i0-1,:) = &
    U_pc%u%f(i0:iend,i0-hs:i0-1,:)*mesh%ll2covari_pc(i0:iend,i0-hs:i0-1,:)%M(2,1) + &
    U_pc%v%f(i0:iend,i0-hs:i0-1,:)*mesh%ll2covari_pc(i0:iend,i0-hs:i0-1,:)%M(2,2) 



    U_pc%ucovari%f(i0:iend,iend+1:iend+hs,:) = &
    U_pc%u%f(i0:iend,iend+1:iend+hs,:)*mesh%ll2covari_pc(i0:iend,iend+1:iend+hs,:)%M(1,1) + &
    U_pc%v%f(i0:iend,iend+1:iend+hs,:)*mesh%ll2covari_pc(i0:iend,iend+1:iend+hs,:)%M(1,2) 
    U_pc%vcovari%f(i0:iend,iend+1:iend+hs,:) = &
    U_pc%u%f(i0:iend,iend+1:iend+hs,:)*mesh%ll2covari_pc(i0:iend,iend+1:iend+hs,:)%M(2,1) + &
    U_pc%v%f(i0:iend,iend+1:iend+hs,:)*mesh%ll2covari_pc(i0:iend,iend+1:iend+hs,:)%M(2,2) 
    !$OMP END PARALLEL WORKSHARE


end subroutine interp_D2Aduogrid



subroutine interp_C2Bgrid(Q_po, Q_pu, Q_pv, id)
    !---------------------------------------------------
    ! C2B grid
    !---------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: Q_po(:,:,:), Q_pu(:,:,:), Q_pv(:,:,:)
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2, a3, a4

    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(Q_po, Q_pu, Q_pv)
        Q_po(i0:iend+1,j0:jend+1,:) = &
        (Q_pv(i0-1:iend,j0:jend+1,:) + Q_pv(i0:iend+1,j0:jend+1,:))*0.25d0 + &
        (Q_pu(i0:iend+1,i0-1:iend,:) + Q_pu(i0:iend+1,j0:jend+1,:))*0.25d0! + &
        !$OMP END PARALLEL WORKSHARE

    elseif(id == 3) then
        ! cubic interpolation coeffs
        a1 = -1.d0/16.d0
        a2 =  9.d0/16.d0
        a3 =  9.d0/16.d0
        a4 = -1.d0/16.d0
     
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(Q_po, Q_pu, Q_pv) &
        !$OMP SHARED(a1, a2, a3, a4)
        Q_po(i0:iend+1,j0:jend+1,:) =(a1*Q_pv(i0-2:iend-1,j0:jend+1,:) &
                                    + a2*Q_pv(i0-1:iend  ,j0:jend+1,:) &
                                    + a3*Q_pv(i0  :iend+1,j0:jend+1,:) &
                                    + a4*Q_pv(i0+1:iend+2,j0:jend+1,:))*0.5d0 &
                                    +(a1*Q_pu(i0:iend+1,j0-2:jend-1,:) &
                                    + a2*Q_pu(i0:iend+1,j0-1:jend  ,:) &
                                    + a3*Q_pu(i0:iend+1,j0  :jend+1,:) &
                                    + a4*Q_pu(i0:iend+1,j0+1:jend+2,:))*0.5d0
 
        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in interp_C2Bgrid: invalid id, ', id
        stop
    end if
end subroutine interp_C2Bgrid

subroutine interp_windC2Bgrid(U_po, V_po, U_pu, V_pv, id)
    !---------------------------------------------------
    ! Wind C2B grid
    !---------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: U_po(:,:,:), V_po(:,:,:)
    real(kind=8), allocatable, intent(inout) :: U_pu(:,:,:), V_pv(:,:,:)
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2, a3, a4

    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_po, V_po, U_pu, V_pv)
        U_po(i0-1:iend+2,j0:jend+1,:) = (U_pu(i0-1:iend+2,j0-1:jend,:) + U_pu(i0-1:iend+2,j0:jend+1,:))*0.5d0
        V_po(i0:iend+1,j0-1:jend+2,:) = (V_pv(i0-1:iend,j0-1:jend+2,:) + V_pv(i0:iend+1,j0-1:jend+2,:))*0.5d0
        !$OMP END PARALLEL WORKSHARE

    elseif(id == 3) then
        ! cubic interpolation coeffs
        a1 = -1.d0/16.d0
        a2 =  9.d0/16.d0
        a3 =  9.d0/16.d0
        a4 = -1.d0/16.d0

        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_po, V_po, U_pu, V_pv) &
        !$OMP SHARED(a1, a2, a3, a4)
        U_po(i0-1:iend+2,j0:jend+1,:) = a1*U_pu(i0-1:iend+2,j0-2:jend-1,:) &
                                      + a2*U_pu(i0-1:iend+2,j0-1:jend  ,:) &
                                      + a3*U_pu(i0-1:iend+2,j0  :jend+1,:) &
                                      + a4*U_pu(i0-1:iend+2,j0+1:jend+2,:)

        V_po(i0:iend+1,j0-1:jend+2,:) = a1*V_pv(i0-2:iend-1,j0-1:jend+2,:) &
                                      + a2*V_pv(i0-1:iend  ,j0-1:jend+2,:) &
                                      + a3*V_pv(i0  :iend+1,j0-1:jend+2,:) &
                                      + a4*V_pv(i0+1:iend+2,j0-1:jend+2,:)
        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in interp_windC2Bgrid: invalid id, ', id
        stop
    end if
end subroutine interp_windC2Bgrid



subroutine compute_lagrange_cs(L, mesh)
    !---------------------------------------------------
    ! compute the lagrange polynomials at ghost cell
    ! points of the cubed-sphere
    !--------------------------------------------------
    type(cubedsphere), intent(inout):: mesh
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4) ::  i, j, p, g, d, k, h, jnearest

    if(mesh%kind .ne. 'equiangular') then
        print*, 'ERROR in compute_lagrange_cs: invalid mesh kind, ', mesh%kind
        stop
    end if

    L%offset = ceiling(L%order*0.5d0)

    if(L%pos==1) then
        ! Compute the nodes
        ! Init local coordinates grid
        L%y_support(n0) = -pio4 - (mesh%halosize-0.5d0)*mesh%dx
        do j = n0+1, nend
            L%y_support(j) = L%y_support(j-1) + mesh%dx
        end do

        ! Invert node points using the cs mapping from panel 2
        p = 2
        do g = 1, nghost
            do j = n0, nend
                ! invert it
                call inverse_equiangular_gnomonic_map2(L%x_nodes(j,g), L%y_nodes(j,g), p, &
                mesh%pc(iend+g,j,1)%p, mesh)
            end do
        end do


        ! Compute stencils
        do g = 1, nghost
           do j = j0-g, jend+g
              do i = n0, nend-1
                 if( (L%y_support(i) <= L%y_nodes(j,g)) .and. (L%y_nodes(j,g) <= L%y_support(i+1) )) then
                    exit
                 end if
              enddo

              ! not corner points
              if(j>j0-g .and. j<jend+g) then
                 L%kend(j,g)   = i + L%offset
                 L%k0(j,g) = L%kend(j,g) - L%order + 1

                 if(L%k0(j,g)<j0)then 
                    L%k0(j,g) = j0
                    L%kend(j,g) = L%k0(j,g) + L%order - 1
                 else if(L%kend(j,g)>jend)then 
                    L%kend(j,g) = jend
                    L%k0(j,g) = L%kend(j,g) - L%order + 1
                 end if

              ! corner points
              else if (j==jend+g) then
                 L%kend(j,g) = min(nend, jend + L%offset)
                 L%k0(j,g) = L%kend(j,g) - L%order + 1
              else
                 L%k0(j,g) = max(n0, j0 - L%offset)
                 L%kend(j,g) = L%k0(j,g) + L%order - 1
              end if
           enddo
        enddo

        ! Debug
        do g = 1, nghost
            do i = i0-g, iend+g
                if(L%kend(i,g)-L%k0(i,g) .ne. L%degree)then
                    print*, 'ERROR in compute_lagrange_polynomials: stencil size is not correct.'
                    stop
                end if

                if(L%kend(i,g)<i0 .or. L%k0(i,g)>iend)then
                    print*, 'ERROR in compute_lagrange_polynomials: stencil bounds are not correct.'
                    stop
                end if

                if(L%y_support(L%k0(i,g))>L%y_nodes(i,g))then
                    print*, 'ERROR in compute_lagrange_polynomials: error in neighboring points. >'
                    stop
                end if

                if(L%y_support(L%kend(i,g))<L%y_nodes(i,g))then
                    print*, 'ERROR in compute_lagrange_polynomials: error in neighboring points. <'
                    stop
                end if
 
            end do
            !read(*,*)
        end do  

        ! Store in f the nearest support points used in Lagrange interpolation
        do g = 1, nghost
            do j = j0-g, jend+g
                L%f_nodes(j,g,:) = L%y_support(L%k0(j,g):L%kend(j,g))
            end do
        end do

        ! Compute the Lagrange nodes at halo region
        do g = 1, nghost
            do j = j0-g, jend+g
                do d = 1, L%order
                    call lagrange_basis(L%y_nodes(j,g), L%f_nodes(j,g,:), L%degree, d, L%p_nodes(j,g,d))
                end do
            end do
        end do

    end if
end subroutine compute_lagrange_cs



subroutine lagrange_basis(x, x_support, N, j, Lj)
    !---------------------------------------------------
    ! Compute the jth Lagrange polynomial of degree N
    !--------------------------------------------------
    integer(i4), intent(in) :: N, j
    real(kind=8), intent(in) :: x, x_support(1:N+1)
    real(kind=8), intent(inout) :: Lj
    integer(i4) :: i

    Lj = 1.d0
    do i = 1, N+1
        if (i .ne. j) then
            Lj = Lj*(x-x_support(i))/(x_support(j)-x_support(i))
        end if
    end do
end subroutine lagrange_basis

subroutine gethalodata_PL07(Qx, Qy)
    !---------------------------------------------------
    ! Fill the center ghost cell values of Qx and Qy that are needed
    ! for reconstruction using the PL07 approach
    !--------------------------------------------------
    real(kind=8), allocatable, intent(inout) :: Qx(:,:,:), Qy(:,:,:)
    integer(i4) :: p, east, north, south, west


    ! --------------------- Panel 1 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 1
    north = 5
    south = 6
    east  = 2
    west  = 4
 
    ! Data of panel 1 from east
    Qx(iend+1:nend,n0:nend,p) = Qx(i0:i0+nghost-1, n0:nend, east) ! Panel 2

    ! Data of panel 1 from  west
    Qx(n0:i0-1,n0:nend,p) = Qx(iend-nghost+1:iend, n0:nend, west) ! Panel 4

    ! Data of panel 1 from north
    Qy(n0:nend,jend+1:nend,p) = Qy(n0:nend, j0:j0+nghost-1, north) ! Panel 5

    ! Data of panel 1 from south
    Qy(n0:nend,n0:j0-1,p) = Qy(n0:nend, jend-nghost+1:jend, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE

    ! --------------------- Panel 2 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 2
    north = 5
    south = 6
    east  = 3
    west  = 1

    ! Data of panel 2 from east
    Qx(iend+1:nend,n0:nend,p) = Qx(i0:i0+nghost-1, n0:nend, east) ! Panel 3

    ! Data of panel 2 from west
    Qx(n0:i0-1,n0:nend,p) = Qx(iend-nghost+1:iend, n0:nend, west) ! Panel 1

    ! Data of panel 2 from north
    Qy(n0:nend,jend+1:nend,p) = transpose(Qx(jend:jend-nghost+1:-1, n0:nend, north)) ! Panel 5

    ! Data of panel 2 from south
    Qy(n0:nend,n0:j0-1,p) = transpose(Qx(jend-nghost+1:jend, nend:n0:-1,south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 3 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 3
    north = 5
    south = 6
    east  = 4
    west  = 2

    ! Data of panel 3 from east
    Qx(iend+1:nend,n0:nend,p) = Qx(i0:i0+nghost-1, n0:nend, east) ! Panel 4

    ! Data of panel 3 from west
    Qx(n0:i0-1,n0:nend,p) = Qx(iend-nghost+1:iend, n0:nend, west) ! Panel 2

    ! Data of panel 3 from north
    Qy(n0:nend,jend+1:nend,p) = Qy(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 5

    ! Data of panel 3 from south
    Qy(n0:nend,n0:j0-1,p) = Qy(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 4 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 4
    north = 5
    south = 6
    east  = 1
    west  = 3

    ! Data of panel 4 from east
    Qx(iend+1:nend,n0:nend,p) = Qx(i0:i0+nghost-1, n0:nend, east) ! Panel 1

    ! Data of panel 4 from west
    Qx(n0:i0-1,n0:nend,p) = Qx(iend-nghost+1:iend, n0:nend, west) ! Panel 3

    ! Data of panel 4 from north
    Qy(n0:nend,jend+1:nend,p) = transpose(Qx(i0:i0+nghost-1, nend:n0:-1, north)) ! Panel 5

    ! Data of panel 4 from south
    Qy(n0:nend,n0:j0-1,p) = transpose(Qx(i0+nghost-1:i0:-1, n0:nend, south)) ! Panel 6
    !$OMP END PARALLEL WORKSHARE



    ! --------------------- Panel 5 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 5
    north = 3
    south = 1
    east  = 2
    west  = 4

    ! Data of panel 5 from east
    Qx(iend+1:nend,n0:nend,p) = transpose(Qy(n0:nend, jend:jend-nghost+1:-1, east)) ! Panel 2

    ! Data of panel 5 from west
    Qx(n0:i0-1,n0:nend,p) = transpose(Qy(nend:n0:-1,jend-nghost+1:jend, west)) ! Panel 4

    ! Data of panel 5 from north
    Qy(n0:nend,jend+1:nend,p) = Qy(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 3

    ! Data of panel 5 from south
    Qy(n0:nend,n0:j0-1,p) = Qy(n0:nend, jend-nghost+1:jend, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE

    ! --------------------- Panel 6 ----------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Qx, Qy, i0, iend, j0, jend, n0, nend, nghost) &
    !$OMP SHARED(north, south, east, west, p)
    p = 6
    north = 1
    south = 3
    east  = 2
    west  = 4

    ! Data of panel 6 from east
    Qx(iend+1:nend,n0:nend,p) = transpose(Qy(nend:n0:-1, j0:j0+nghost-1, east)) ! Panel 2

    ! Data of panel 6 from west
    Qx(n0:i0-1,n0:nend,p) = transpose(Qy(n0:nend,j0+nghost-1:j0:-1, west)) ! Panel 4

    ! Data of panel 6 from north
    Qy(n0:nend,jend+1:nend,p) = Qy(n0:nend, j0:j0+nghost-1, north) ! Panel 3

    ! Data of panel 6 from south
    Qy(n0:nend,n0:j0-1,p) = Qy(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE


end subroutine gethalodata_PL07



end module duogrid_interpolation
