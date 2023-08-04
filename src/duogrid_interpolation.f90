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
    vector_field, &
    simulation, &
    lagrange_poly_cs

! CS mapping
use sphgeo, only: &
    inverse_equiangular_gnomonic_map2, &
    ll2contra, contra2ll, &
    contra2covari, &
    covari2contra

implicit none

contains 

subroutine gethalodata(Q, L)
    !---------------------------------------------------
    ! Fill the center ghost cell values of Q that are needed
    ! for duogrid interpolation
    !--------------------------------------------------
    type(scalar_field), intent(inout) :: Q
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
    L%halodata_east(1:nghost,n0:nend,p) = Q%f(i0:i0+nghost-1, n0:nend, east) ! Panel 2

    ! Data of panel 1 from  west
    L%halodata_west(1:nghost,n0:nend,p) = Q%f(iend-nghost+1:iend, n0:nend, west) ! Panel 4

    ! Data of panel 1 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q%f(n0:nend, j0:j0+nghost-1, north) ! Panel 5

    ! Data of panel 1 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q%f(n0:nend, jend-nghost+1:jend, south) ! Panel 6
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
    L%halodata_east(1:nghost,n0:nend,p) = Q%f(i0:i0+nghost-1, n0:nend, east) ! Panel 3

    ! Data of panel 2 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q%f(iend-nghost+1:iend, n0:nend, west) ! Panel 1

    ! Data of panel 2 from north
    L%halodata_north(n0:nend,1:nghost,p) = transpose(Q%f(jend:jend-nghost+1:-1, n0:nend, north)) ! Panel 5

    ! Data of panel 2 from south
    L%halodata_south(n0:nend,1:nghost,p) = transpose(Q%f(jend-nghost+1:jend, nend:n0:-1,south)) ! Panel 6
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
    L%halodata_east(1:nghost,n0:nend,p) = Q%f(i0:i0+nghost-1, n0:nend, east) ! Panel 4

    ! Data of panel 3 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q%f(iend-nghost+1:iend, n0:nend, west) ! Panel 2

    ! Data of panel 3 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q%f(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 5

    ! Data of panel 3 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q%f(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 6
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
    L%halodata_east(1:nghost,n0:nend,p) = Q%f(i0:i0+nghost-1, n0:nend, east) ! Panel 1

    ! Data of panel 4 from west
    L%halodata_west(1:nghost,n0:nend,p) = Q%f(iend-nghost+1:iend, n0:nend, west) ! Panel 3

    ! Data of panel 4 from north
    L%halodata_north(n0:nend,1:nghost,p) = transpose(Q%f(i0:i0+nghost-1, nend:n0:-1, north)) ! Panel 5

    ! Data of panel 4 from south
    L%halodata_south(n0:nend,1:nghost,p) = transpose(Q%f(i0+nghost-1:i0:-1, n0:nend, south)) ! Panel 6
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
    L%halodata_east(1:nghost,n0:nend,p) = transpose(Q%f(n0:nend, jend:jend-nghost+1:-1, east)) ! Panel 2

    ! Data of panel 5 from west
    L%halodata_west(1:nghost,n0:nend,p) = transpose(Q%f(nend:n0:-1,jend-nghost+1:jend, west)) ! Panel 4

    ! Data of panel 5 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q%f(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 3

    ! Data of panel 5 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q%f(n0:nend, jend-nghost+1:jend, south) ! Panel 1
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
    L%halodata_east(1:nghost,n0:nend,p) = transpose(Q%f(nend:n0:-1, j0:j0+nghost-1, east)) ! Panel 2

    ! Data of panel 6 from west
    L%halodata_west(1:nghost,n0:nend,p) = transpose(Q%f(n0:nend,j0+nghost-1:j0:-1, west)) ! Panel 4

    ! Data of panel 6 from north
    L%halodata_north(n0:nend,1:nghost,p) = Q%f(n0:nend, j0:j0+nghost-1, north) ! Panel 3

    ! Data of panel 6 from south
    L%halodata_south(n0:nend,1:nghost,p) = Q%f(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE


end subroutine gethalodata



subroutine dg_interp(Q, L)
    !---------------------------------------------------
    !   duogrid interpolation of scalar field Q
    ! (ghost cells are defined at cell centers)
    !--------------------------------------------------
    type(scalar_field), intent(inout) :: Q
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
                Q%f(iend+g,j,p) = 0.d0
                do d = 1, L%order
                    Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
                Q%f(i0-g2,j,p) = 0.d0
                do d = 1, L%order
                    Q%f(i0-g2,j,p) = Q%f(i0-g2,j,p) + L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
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
                Q%f(i,jend+g,p) = 0.d0
                do d = 1, L%order
                    Q%f(i,jend+g,p) = Q%f(i,jend+g,p) + L%f_nodes(i,g,d)*L%p_nodes(i,g,d)
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
                Q%f(i,j0-g2,p) = 0.d0
                do d = 1, L%order
                    Q%f(i,j0-g2,p) = Q%f(i,j0-g2,p) + L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
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
            Q%f(iend+g,j,p) = 0.d0
            do d = 1, L%order
                Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
                Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
            Q%f(iend+g,j,p) = 0.d0
            do d = 1, L%order
                Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
            end do
        end do

        ! North
        do g = 1, nghost
            g2 = nghost-g+1
            i = iend+g2
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g2,:)= L%halodata_south(L%k0(i,g2):L%kend(i,g2), g, p)
            do d = 1, L%order
                Q%f(i,j0-g2,p) = Q%f(i,j0-g2,p) + 0.5d0*L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
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
            Q%f(i0-g2,j,p) = 0.d0
            do d = 1, L%order
                Q%f(i0-g2,j,p) = Q%f(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
            end do
        end do

        ! North
        do g = 1, nghost
            i = i0-g
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g,:) = L%halodata_north(L%k0(i,g):L%kend(i,g), g, p)
            ! Does the interpolation
            do d = 1, L%order
                Q%f(i,jend+g,p) = Q%f(i,jend+g,p) + 0.5d0*L%f_nodes(i,g,d)*L%p_nodes(i,g,d)
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
            Q%f(i0-g2,j,p) = 0.d0
            do d = 1, L%order
                Q%f(i0-g2,j,p) = Q%f(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
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
                Q%f(i,j0-g2,p) = Q%f(i,j0-g2,p) + 0.5d0*L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
            end do
        end do
        !--------------------------------------------------------------------------
    end do

end subroutine dg_interp


subroutine dg_interp_C2Agrid(U_pu, U_pv, U_pc, L, mesh, id)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at C grid (contravariant)
    ! to the A grid (contravariant), including its ghost cell values at centers
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(vector_field), intent(inout) :: U_pu, U_pv, U_pc
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2, a3, a4
    real(kind=8) :: b1, b2, b3, b4
    
    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_pu, U_pv, U_pc)
        U_pc%ucontra%f(i0:iend,j0:jend,:) = (U_pu%ucontra%f(i0:iend,j0:jend,:)+U_pu%ucontra%f(i0+1:iend+1,j0:jend,:))*0.5d0
        U_pc%vcontra%f(i0:iend,j0:jend,:) = (U_pv%vcontra%f(i0:iend,j0:jend,:)+U_pv%vcontra%f(i0:iend,j0+1:jend+1,:))*0.5d0
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
        !$OMP SHARED(U_pu, U_pv, U_pc)

        ! west boundary
        U_pc%ucontra%f(i0,j0:jend,:) = &
          a1*U_pu%ucontra%f(i0,j0:jend,:) &
        + a2*U_pu%ucontra%f(i0+1,j0:jend,:) &
        + a3*U_pu%ucontra%f(i0+2,j0:jend,:) &
        + a4*U_pu%ucontra%f(i0+3,j0:jend,:)

        ! east boundary
        U_pc%ucontra%f(iend,j0:jend,:) = &
          a4*U_pu%ucontra%f(iend-2,j0:jend,:) &
        + a3*U_pu%ucontra%f(iend-1,j0:jend,:) &
        + a2*U_pu%ucontra%f(iend,j0:jend,:) &
        + a1*U_pu%ucontra%f(iend+1,j0:jend,:)

        ! south boundary
        U_pc%vcontra%f(i0:iend,j0,:) = &
          a1*U_pv%vcontra%f(i0:iend,j0,:) &
        + a2*U_pv%vcontra%f(i0:iend,j0+1,:) &
        + a3*U_pv%vcontra%f(i0:iend,j0+2,:) &
        + a4*U_pv%vcontra%f(i0:iend,j0+3,:)

        ! north boundary
        U_pc%vcontra%f(i0:iend,jend,:) = &
        a4*U_pv%vcontra%f(i0:iend,jend-2,:) + &
        a3*U_pv%vcontra%f(i0:iend,jend-1,:) + &
        a2*U_pv%vcontra%f(i0:iend,jend-0,:) + &
        a1*U_pv%vcontra%f(i0:iend,jend+1,:)

        ! remaing cells
        U_pc%ucontra%f(i0+1:iend-1,j0:jend,:) = & 
          b1*U_pu%ucontra%f(i0:iend-2,j0:jend,:) &
        + b2*U_pu%ucontra%f(i0+1:iend-1,j0:jend,:) &
        + b3*U_pu%ucontra%f(i0+2:iend  ,j0:jend,:) &
        + b4*U_pu%ucontra%f(i0+3:iend+1,j0:jend,:)

        U_pc%vcontra%f(i0:iend,j0+1:jend-1,:) = & 
          b1*U_pv%vcontra%f(i0:iend,j0:jend-2,:) &
        + b2*U_pv%vcontra%f(i0:iend,j0+1:jend-1,:) &
        + b3*U_pv%vcontra%f(i0:iend,j0+2:jend  ,:) &
        + b4*U_pv%vcontra%f(i0:iend,j0+3:jend+1,:)

        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in dg_interp_C2Agrid: invalid id, ', id
        stop
    end if

    ! Convert from contravariant to latlon
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(U_pc, mesh) & 
    !$OMP SHARED(n0, nend, i0, iend, j0, jend, nbfaces) &
    !$OMP PRIVATE(i, j, p) &
    !$OMP SCHEDULE(static)
    do i = i0, iend
        do j = j0, jend
            do p = 1, nbfaces
                call contra2ll(U_pc%u%f(i,j,p), U_pc%v%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%contra2ll_pc(i,j,p)%M)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! Interpolate latlon to the ghost cell centers
    call dg_interp(U_pc%u, L)
    call dg_interp(U_pc%v, L)

    ! Convert from latlon to contravariant
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(U_pc, mesh) & 
    !$OMP SHARED(i0, iend, j0, jend, n0, nend, nbfaces) &
    !$OMP PRIVATE(i, p) &
    !$OMP SCHEDULE(static)
    do i = n0, nend
        do j = n0, nend
            do p = 1, nbfaces
                call ll2contra(U_pc%u%f(i,j,p), U_pc%v%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%ll2contra_pc(i,j,p)%M)
            end do
        end do
    end do
    !$OMP END PARALLEL DO


end subroutine dg_interp_C2Agrid


subroutine dg_interp_D2Agrid(U_pu, U_pv, U_pc, L, mesh, id)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at D grid (covariant)
    ! to the A grid (contravariant), including its ghost cell values at centers
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(vector_field), intent(inout) :: U_pu, U_pv, U_pc
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4), intent(in) :: id
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2, a3, a4
    real(kind=8) :: b1, b2, b3, b4

    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_pu, U_pv, U_pc)
        U_pc%ucovari%f(i0:iend,j0:jend,:) = (U_pu%ucovari%f(i0:iend,j0:jend,:)+U_pu%ucovari%f(i0+1:iend+1,j0:jend,:))*0.5d0
        U_pc%vcovari%f(i0:iend,j0:jend,:) = (U_pv%vcovari%f(i0:iend,j0:jend,:)+U_pv%vcovari%f(i0:iend,j0+1:jend+1,:))*0.5d0
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
        !$OMP SHARED(U_pu, U_pv, U_pc)

        ! west boundary
        U_pc%vcovari%f(i0,j0:jend,:) = &
          a1*U_pu%vcovari%f(i0,j0:jend,:) &
        + a2*U_pu%vcovari%f(i0+1,j0:jend,:) &
        + a3*U_pu%vcovari%f(i0+2,j0:jend,:) &
        + a4*U_pu%vcovari%f(i0+3,j0:jend,:)

        ! east boundary
        U_pc%vcovari%f(iend,j0:jend,:) = &
          a4*U_pu%vcovari%f(iend-2,j0:jend,:) &
        + a3*U_pu%vcovari%f(iend-1,j0:jend,:) &
        + a2*U_pu%vcovari%f(iend,j0:jend,:) &
        + a1*U_pu%vcovari%f(iend+1,j0:jend,:)

        ! south boundary
        U_pc%ucovari%f(i0:iend,j0,:) = &
          a1*U_pv%ucovari%f(i0:iend,j0,:) &
        + a2*U_pv%ucovari%f(i0:iend,j0+1,:) &
        + a3*U_pv%ucovari%f(i0:iend,j0+2,:) &
        + a4*U_pv%ucovari%f(i0:iend,j0+3,:)

        ! north boundary
        U_pc%ucovari%f(i0:iend,jend,:) = &
        a4*U_pv%ucovari%f(i0:iend,jend-2,:) + &
        a3*U_pv%ucovari%f(i0:iend,jend-1,:) + &
        a2*U_pv%ucovari%f(i0:iend,jend-0,:) + &
        a1*U_pv%ucovari%f(i0:iend,jend+1,:)

        ! remaing cells
        U_pc%vcovari%f(i0+1:iend-1,j0:jend,:) = & 
          b1*U_pu%vcovari%f(i0:iend-2,j0:jend,:) &
        + b2*U_pu%vcovari%f(i0+1:iend-1,j0:jend,:) &
        + b3*U_pu%vcovari%f(i0+2:iend  ,j0:jend,:) &
        + b4*U_pu%vcovari%f(i0+3:iend+1,j0:jend,:)

        U_pc%ucovari%f(i0:iend,j0+1:jend-1,:) = & 
          b1*U_pv%ucovari%f(i0:iend,j0:jend-2,:) &
        + b2*U_pv%ucovari%f(i0:iend,j0+1:jend-1,:) &
        + b3*U_pv%ucovari%f(i0:iend,j0+2:jend  ,:) &
        + b4*U_pv%ucovari%f(i0:iend,j0+3:jend+1,:)

        !$OMP END PARALLEL WORKSHARE
    else
        print*, 'ERROR in dg_interp_D2Agrid: invalid id, ', id
        stop
    end if

    
    ! Convert from covariant to latlon
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(U_pc, mesh) & 
    !$OMP SHARED(n0, nend, i0, iend, j0, jend, nbfaces) &
    !$OMP PRIVATE(i, j, p) &
    !$OMP SCHEDULE(static)
    do i = i0, iend
        do j = j0, jend
            do p = 1, nbfaces
                call covari2contra(U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p),&
                U_pc%ucovari%f(i,j,p), U_pc%vcovari%f(i,j,p), mesh%covari2contra_pc(i,j,p)%M)
 
                call contra2ll(U_pc%u%f(i,j,p), U_pc%v%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%contra2ll_pc(i,j,p)%M)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! Interpolate latlon to the ghost cell centers
    call dg_interp(U_pc%u, L)
    call dg_interp(U_pc%v, L)

    ! Convert from latlon to contravariant and then covariant
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(U_pc, mesh) & 
    !$OMP SHARED(i0, iend, j0, jend, n0, nend, nbfaces) &
    !$OMP PRIVATE(i, p) &
    !$OMP SCHEDULE(static)
    do i = n0, nend
        do j = n0, nend
            do p = 1, nbfaces
                call ll2contra(U_pc%u%f(i,j,p), U_pc%v%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%ll2contra_pc(i,j,p)%M)

                call contra2covari(U_pc%ucovari%f(i,j,p), U_pc%vcovari%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%contra2covari_pc(i,j,p)%M)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

 
end subroutine dg_interp_D2Agrid



subroutine dg_interp_A2Cghostgrid(U_pu, U_pv, U_pc, L, mesh, id)
    !---------------------------------------------------
    ! duogrid interpolation of vector field given at A grid (contravariant)
    ! to the C grid (contravariant) (including ghost cells)
    type(cubedsphere), intent(inout) :: mesh
    type(vector_field), intent(inout) :: U_pu, U_pv, U_pc
    integer(i4), intent(in) :: id
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4):: i, j, p, h
    real(kind=8) :: c1, c2
    
    ! cubic interpolation coeffs
    c1 =  9.d0/16.d0
    c2 = -1.d0/16.d0

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(i0, iend, j0, jend, n0, nend) &
    !$OMP SHARED(c1, c2) &
    !$OMP SHARED(U_pc, U_pu, U_pv)
    ! Let us interpolate the ghost cell edges - cubic interpolation
    ! Panel from south
    U_pu%ucontra%f(i0:iend+1,n0:j0-1,:) = &
    c1*(U_pc%ucontra%f(i0:iend+1  ,n0:j0-1,:) + U_pc%ucontra%f(i0-1:iend  ,n0:j0-1,:)) + &
    c2*(U_pc%ucontra%f(i0+1:iend+2,n0:j0-1,:) + U_pc%ucontra%f(i0-2:iend-1,n0:j0-1,:))

    U_pu%vcontra%f(i0:iend+1,n0:j0-1,:) = &
    c1*(U_pc%vcontra%f(i0:iend+1  ,n0:j0-1,:) + U_pc%vcontra%f(i0-1:iend  ,n0:j0-1,:)) + &
    c2*(U_pc%vcontra%f(i0+1:iend+2,n0:j0-1,:) + U_pc%vcontra%f(i0-2:iend-1,n0:j0-1,:))

    ! Panel from north
    U_pu%ucontra%f(i0:iend+1, jend+1:nend,:) = &
    c1*(U_pc%ucontra%f(i0:iend+1  ,jend+1:nend,:) + U_pc%ucontra%f(i0-1:iend  ,jend+1:nend,:)) + &
    c2*(U_pc%ucontra%f(i0+1:iend+2,jend+1:nend,:) + U_pc%ucontra%f(i0-2:iend-1,jend+1:nend,:))

    U_pu%vcontra%f(i0:iend+1,jend+1:nend,:) = &
    c1*(U_pc%vcontra%f(i0:iend+1  ,jend+1:nend,:) + U_pc%vcontra%f(i0-1:iend  ,jend+1:nend,:)) + &
    c2*(U_pc%vcontra%f(i0+1:iend+2,jend+1:nend,:) + U_pc%vcontra%f(i0-2:iend-1,jend+1:nend,:))

    ! Panel from west
    U_pv%ucontra%f(n0:i0-1,j0:jend+1,:) = &
    c1*(U_pc%ucontra%f(n0:i0-1,j0:jend+1,  :) + U_pc%ucontra%f(n0:i0-1,j0-1:jend,  :)) + &
    c2*(U_pc%ucontra%f(n0:i0-1,j0+1:jend+2,:) + U_pc%ucontra%f(n0:i0-1,j0-2:jend-1,:))

    U_pv%vcontra%f(n0:i0-1,j0:jend+1,:) = &
    c1*(U_pc%vcontra%f(n0:i0-1,j0:jend+1  ,:) + U_pc%vcontra%f(n0:i0-1,j0-1:jend  ,:)) + &
    c2*(U_pc%vcontra%f(n0:i0-1,j0+1:jend+2,:) + U_pc%vcontra%f(n0:i0-1,j0-2:jend-1,:))

    ! Panel from east 
    U_pv%ucontra%f(iend+1:nend, j0:jend+1, :) = &
    c1*(U_pc%ucontra%f(iend+1:nend,j0:jend+1  ,:) + U_pc%ucontra%f(iend+1:nend,j0-1:jend  ,:)) + &
    c2*(U_pc%ucontra%f(iend+1:nend,j0+1:jend+2,:) + U_pc%ucontra%f(iend+1:nend,j0-2:jend-1,:))

    U_pv%vcontra%f(iend+1:nend, j0:jend+1, :) = &
    c1*(U_pc%vcontra%f(iend+1:nend,j0:jend+1  ,:) + U_pc%vcontra%f(iend+1:nend,j0-1:jend  ,:)) + &
    c2*(U_pc%vcontra%f(iend+1:nend,j0+1:jend+2,:) + U_pc%vcontra%f(iend+1:nend,j0-2:jend-1,:))


    !-----------------------------------------------------------------------------------
    ! Interpolation needed for RK2 departure point scheme
    ! Panel from west
    U_pu%ucontra%f(i0-1,n0:nend,:) = &
    c1*(U_pc%ucontra%f(i0-2,n0:nend,:) + U_pc%ucontra%f(i0-1,n0:nend,:)) + &
    c2*(U_pc%ucontra%f(i0,n0:nend,:) + U_pc%ucontra%f(i0-3,n0:nend,:))

    U_pu%vcontra%f(i0-1,n0:nend,:) = &
    c1*(U_pc%vcontra%f(i0-2,n0:nend,:) + U_pc%vcontra%f(i0-1,n0:nend,:)) + &
    c2*(U_pc%vcontra%f(i0,n0:nend,:) + U_pc%vcontra%f(i0-3,n0:nend,:))


    ! Panel from east
    U_pu%ucontra%f(iend+2,n0:nend,:) = &
    c1*(U_pc%ucontra%f(iend+1,:,:) + U_pc%ucontra%f(iend+2,:,:)) + &
    c2*(U_pc%ucontra%f(iend,:,:) + U_pc%ucontra%f(iend+3,:,:))

    U_pu%vcontra%f(iend+2,:,:) = &
    c1*(U_pc%vcontra%f(iend+1,:,:) + U_pc%vcontra%f(iend+2,:,:)) + &
    c2*(U_pc%vcontra%f(iend,:,:) + U_pc%vcontra%f(iend+3,:,:))
 
    ! Panel from south
    U_pv%ucontra%f(n0:nend,j0-1,:) = &
    c1*(U_pc%ucontra%f(n0:nend,j0-2,:) + U_pc%ucontra%f(n0:nend,j0-1,:)) + &
    c2*(U_pc%ucontra%f(n0:nend,j0,:) + U_pc%ucontra%f(n0:nend,j0-3,:))

    U_pv%vcontra%f(n0:nend,j0-1,:) = &
    c1*(U_pc%vcontra%f(n0:nend,j0-2,:) + U_pc%vcontra%f(n0:nend,j0-1,:)) + &
    c2*(U_pc%vcontra%f(n0:nend,j0,:) + U_pc%vcontra%f(n0:nend,j0-3,:))

    ! Panel from east
    U_pv%ucontra%f(:,jend+2,:) = &
    c1*(U_pc%ucontra%f(:,jend+1,:) + U_pc%ucontra%f(:,jend+2,:)) + &
    c2*(U_pc%ucontra%f(:,jend,:) + U_pc%ucontra%f(:,jend+3,:))

    U_pv%vcontra%f(:,jend+2,:) = &
    c1*(U_pc%vcontra%f(:,jend+1,:) + U_pc%vcontra%f(:,jend+2,:)) + &
    c2*(U_pc%vcontra%f(:,jend,:) + U_pc%vcontra%f(:,jend+3,:))
    !$OMP END PARALLEL WORKSHARE



    ! Now, let us interpolate from the A grid to the C-grid points
    ! that are not ghost cells
    if(id == 1) then
        ! linear interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_pu, U_pv, U_pc)
        U_pu%ucontra%f(i0:iend+1,j0:jend,:) = (U_pc%ucontra%f(i0-1:iend,j0:jend,:)+U_pc%ucontra%f(i0:iend+1,j0:jend,:))*0.5d0
        U_pv%vcontra%f(i0:iend,j0:jend+1,:) = (U_pc%vcontra%f(i0:iend,j0-1:jend,:)+U_pc%vcontra%f(i0:iend,j0:jend+1,:))*0.5d0
        !$OMP END PARALLEL WORKSHARE

    elseif(id == 3) then
        ! cubic interpolation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(i0, iend, j0, jend) &
        !$OMP SHARED(U_pu, U_pv, U_pc, c1, c2)
        U_pu%ucontra%f(i0:iend+1,j0:jend,:) = &
        c1*(U_pc%ucontra%f(i0-1:iend,j0:jend,:)+U_pc%ucontra%f(i0:iend+1,j0:jend,:)) + &
        c2*(U_pc%ucontra%f(i0-2:iend-1,j0:jend,:)+U_pc%ucontra%f(i0+1:iend+2,j0:jend,:))

        U_pv%vcontra%f(i0:iend,j0:jend+1,:) = &
        c1*(U_pc%vcontra%f(i0:iend,j0-1:jend,:)+U_pc%vcontra%f(i0:iend,j0:jend+1,:)) + &
        c2*(U_pc%vcontra%f(i0:iend,j0-2:jend-1,:)+U_pc%vcontra%f(i0:iend,j0+1:jend+2,:))
        !$OMP END PARALLEL WORKSHARE

    else
        print*, 'ERROR in dg_interp_C2Agrid: invalid id, ', id
        stop
    end if
end subroutine dg_interp_A2Cghostgrid


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

        L%k0(:,:) = 0
        L%kend(:,:) =0

        ! Compute stencils
        do g = 1, nghost
            do j = j0-g, jend+g
                ! Point that are not at the corners
                if (j>j0-g .and. j<jend+g)then
                    jnearest =  minloc(abs(L%y_support(j0:jend)-L%y_nodes(j,g)),DIM=1)
                    if(L%y_nodes(j,g) > L%y_support(jnearest))then
                        L%kend(j,g) = jnearest  + ceiling(L%order*0.5d0)
                        L%k0(j,g) = L%kend(j,g) - L%order + 1 
                    else 
                        L%k0(j,g) = jnearest - ceiling(L%order*0.5d0)
                        L%kend(j,g) = L%k0(j,g) + L%order - 1 
                    end if
                
                    if(L%k0(j,g)<j0)then
                        L%k0(j,g) = j0
                        L%kend(j,g) = L%k0(j,g) + L%order - 1 
                    end if

                    if(L%kend(j,g)>jend)then
                        L%kend(j,g) = jend
                        L%k0(j,g) = L%kend(j,g) - L%order + 1 
                    end if

                ! Corner points
                else if (j==jend+g) then
                    L%kend(j,g) = min(nend, jend + ceiling(L%order*0.5d0))
                    L%k0(j,g)   = L%kend(j,g) - L%order + 1
                else 
                    L%k0(j,g)   = max(n0, j0-ceiling(L%order*0.5d0))
                    L%kend(j,g) = L%k0(j,g) + L%order - 1
                end if
            end do
            !read(*,*)
            !stop
        end do

        do g = 1, nghost
            do i = i0-g, iend+g
                if(L%kend(i,g)-L%k0(i,g) .ne. L%degree)then
                    print*, 'error1'
                    stop
                end if

                !if(L%kend(i,g)<i0 .or. L%k0(i,g)>iend)then
                !    print*, 'error2'
                !    stop
                !end if

                if(L%y_support(L%k0(i,g))>L%y_nodes(i,g))then
                    print*, 'error3'
                    print*, i,g,L%k0(i,g), L%kend(i,g)
                    read(*,*)
                    stop
                end if

                if(L%y_support(L%kend(i,g))<L%y_nodes(i,g))then
                    print*, 'error4'
                    print*, i,g,L%k0(i,g), L%kend(i,g)
                    read(*,*)
                    stop
                end if
 
                !print*, i, g, L%y_support(L%k0(i,g)), L%y_nodes(i,g), L%y_support(L%kend(i,g))
                !print*, i, g,L%k0(i,g), L%kend(i,g)
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
    type(scalar_field), intent(inout) :: Qx, Qy
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
    Qx%f(iend+1:nend,n0:nend,p) = Qx%f(i0:i0+nghost-1, n0:nend, east) ! Panel 2

    ! Data of panel 1 from  west
    Qx%f(n0:i0-1,n0:nend,p) = Qx%f(iend-nghost+1:iend, n0:nend, west) ! Panel 4

    ! Data of panel 1 from north
    Qy%f(n0:nend,jend+1:nend,p) = Qy%f(n0:nend, j0:j0+nghost-1, north) ! Panel 5

    ! Data of panel 1 from south
    Qy%f(n0:nend,n0:j0-1,p) = Qy%f(n0:nend, jend-nghost+1:jend, south) ! Panel 6
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
    Qx%f(iend+1:nend,n0:nend,p) = Qx%f(i0:i0+nghost-1, n0:nend, east) ! Panel 3

    ! Data of panel 2 from west
    Qx%f(n0:i0-1,n0:nend,p) = Qx%f(iend-nghost+1:iend, n0:nend, west) ! Panel 1

    ! Data of panel 2 from north
    Qy%f(n0:nend,jend+1:nend,p) = transpose(Qx%f(jend:jend-nghost+1:-1, n0:nend, north)) ! Panel 5

    ! Data of panel 2 from south
    Qy%f(n0:nend,n0:j0-1,p) = transpose(Qx%f(jend-nghost+1:jend, nend:n0:-1,south)) ! Panel 6
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
    Qx%f(iend+1:nend,n0:nend,p) = Qx%f(i0:i0+nghost-1, n0:nend, east) ! Panel 4

    ! Data of panel 3 from west
    Qx%f(n0:i0-1,n0:nend,p) = Qx%f(iend-nghost+1:iend, n0:nend, west) ! Panel 2

    ! Data of panel 3 from north
    Qy%f(n0:nend,jend+1:nend,p) = Qy%f(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 5

    ! Data of panel 3 from south
    Qy%f(n0:nend,n0:j0-1,p) = Qy%f(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 6
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
    Qx%f(iend+1:nend,n0:nend,p) = Qx%f(i0:i0+nghost-1, n0:nend, east) ! Panel 1

    ! Data of panel 4 from west
    Qx%f(n0:i0-1,n0:nend,p) = Qx%f(iend-nghost+1:iend, n0:nend, west) ! Panel 3

    ! Data of panel 4 from north
    Qy%f(n0:nend,jend+1:nend,p) = transpose(Qx%f(i0:i0+nghost-1, nend:n0:-1, north)) ! Panel 5

    ! Data of panel 4 from south
    Qy%f(n0:nend,n0:j0-1,p) = transpose(Qx%f(i0+nghost-1:i0:-1, n0:nend, south)) ! Panel 6
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
    Qx%f(iend+1:nend,n0:nend,p) = transpose(Qy%f(n0:nend, jend:jend-nghost+1:-1, east)) ! Panel 2

    ! Data of panel 5 from west
    Qx%f(n0:i0-1,n0:nend,p) = transpose(Qy%f(nend:n0:-1,jend-nghost+1:jend, west)) ! Panel 4

    ! Data of panel 5 from north
    Qy%f(n0:nend,jend+1:nend,p) = Qy%f(nend:n0:-1, jend:jend-nghost+1:-1, north) ! Panel 3

    ! Data of panel 5 from south
    Qy%f(n0:nend,n0:j0-1,p) = Qy%f(n0:nend, jend-nghost+1:jend, south) ! Panel 1
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
    Qx%f(iend+1:nend,n0:nend,p) = transpose(Qy%f(nend:n0:-1, j0:j0+nghost-1, east)) ! Panel 2

    ! Data of panel 6 from west
    Qx%f(n0:i0-1,n0:nend,p) = transpose(Qy%f(n0:nend,j0+nghost-1:j0:-1, west)) ! Panel 4

    ! Data of panel 6 from north
    Qy%f(n0:nend,jend+1:nend,p) = Qy%f(n0:nend, j0:j0+nghost-1, north) ! Panel 3

    ! Data of panel 6 from south
    Qy%f(n0:nend,n0:j0-1,p) = Qy%f(nend:n0:-1, j0+nghost-1:j0:-1, south) ! Panel 1
    !$OMP END PARALLEL WORKSHARE


end subroutine gethalodata_PL07



end module duogrid_interpolation
