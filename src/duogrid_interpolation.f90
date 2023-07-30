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
    ll2contra, contra2ll

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
            !h = g-1
            do j = n0, nend
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
            !Q%f(iend+g,j,p) = 0.d0
            do d = 1, L%order
            !    Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
            !    Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
            !Q%f(iend+g,j,p) =0.d0
            do d = 1, L%order
            !    Q%f(iend+g,j,p) = Q%f(iend+g,j,p) + 0.5d0*L%f_nodes(j,g,d)*L%p_nodes(j,g,d)
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
            !Q%f(i0-g2,j,p) = 0.d0
            do d = 1, L%order
            !    Q%f(i0-g2,j,p) = Q%f(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
            end do
        end do

        ! North
        do g = 1, nghost
            i = i0-g
            ! Store in f the support points used in Lagrange interpolation
            L%f_nodes(i,g,:) = L%halodata_north(L%k0(i,g):L%kend(i,g), g, p)
            ! Does the interpolation
            do d = 1, L%order
            !    Q%f(i,jend+g,p) = Q%f(i,jend+g,p) + 0.5d0*L%f_nodes(i,g,d)*L%p_nodes(i,g,d)
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
            !Q%f(i0-g2,j,p) = 0.d0
            do d = 1, L%order
            !    Q%f(i0-g2,j,p) = Q%f(i0-g2,j,p) + 0.5d0*L%f_nodes(j,g2,d)*L%p_nodes(j,g2,d)
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
            !   Q%f(i,j0-g2,p) = Q%f(i,j0-g2,p) + 0.5d0*L%f_nodes(i,g2,d)*L%p_nodes(i,g2,d)
            end do
        end do
        !--------------------------------------------------------------------------
    end do

end subroutine dg_interp

subroutine dg_vf_interp(U_pu, U_pv, U_pc, L, mesh)
    !---------------------------------------------------
    !   duogrid interpolation of vector field 
    ! (ghost cells are defined at cell centers)
    !--------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(vector_field), intent(inout) :: U_pu, U_pv, U_pc
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4):: i, j, p
    real(kind=8) :: a1, a2

    ! cubic interpolation coeffs
    a1 =  9.d0/16.d0
    a2 = -1.d0/16.d0

    ! Let us populate the center values needed for the ghost cell interpolation
    ! using a centered averaged (2nd order accurate)
    U_pc%ucontra%f(i0:iend,j0:jend,:) = (U_pu%ucontra%f(i0:iend,j0:jend,:) + U_pu%ucontra%f(i0+1:iend+1,j0:jend,:))*0.5d0
    U_pc%vcontra%f(i0:iend,j0:jend,:) = (U_pv%vcontra%f(i0:iend,j0:jend,:) + U_pv%vcontra%f(i0:iend,j0+1:jend+1,:))*0.5d0

    ! Convert from contravariant to latlon
    do i = i0, iend
        do j = j0, jend
            do p = 1, nbfaces
                call contra2ll(U_pc%u%f(i,j,p), U_pc%v%f(i,j,p), &
                U_pc%ucontra%f(i,j,p), U_pc%vcontra%f(i,j,p), mesh%contra2ll_pc(i,j,p)%M)
            end do
        end do
    end do

    ! Interpolate latlon to the ghost cell centers
    call dg_interp(U_pc%u, L)
    call dg_interp(U_pc%v, L)

    ! Now, let us interpolate the ghost cell edges - cubic interpolation
    ! Panel from south
    U_pu%u%f(i0:iend+1,n0:j0-1,:) = &
    a1*(U_pc%u%f(i0:iend+1  ,n0:j0-1,:) + U_pc%u%f(i0-1:iend  ,n0:j0-1,:)) + &
    a2*(U_pc%u%f(i0+1:iend+2,n0:j0-1,:) + U_pc%u%f(i0-2:iend-1,n0:j0-1,:))

    U_pu%v%f(i0:iend+1,n0:j0-1,:) = &
    a1*(U_pc%v%f(i0:iend+1  ,n0:j0-1,:) + U_pc%v%f(i0-1:iend  ,n0:j0-1,:)) + &
    a2*(U_pc%v%f(i0+1:iend+2,n0:j0-1,:) + U_pc%v%f(i0-2:iend-1,n0:j0-1,:))

    ! Convert from latlon to contravariant
    do i = i0, iend+1
        do p = 1, nbfaces
            do j = n0, j0-1
                call ll2contra(U_pu%u%f(i,j,p), U_pu%v%f(i,j,p), &
                U_pu%ucontra%f(i,j,p), U_pu%vcontra%f(i,j,p), mesh%ll2contra_pu(i,j,p)%M)
            end do
        end do
    end do

    ! Panel from north
    U_pu%u%f(i0:iend+1, jend+1:nend,:) = &
    a1*(U_pc%u%f(i0:iend+1  ,jend+1:nend,:) + U_pc%u%f(i0-1:iend  ,jend+1:nend,:)) + &
    a2*(U_pc%u%f(i0+1:iend+2,jend+1:nend,:) + U_pc%u%f(i0-2:iend-1,jend+1:nend,:))

    U_pu%v%f(i0:iend+1,jend+1:nend,:) = &
    a1*(U_pc%v%f(i0:iend+1  ,jend+1:nend,:) + U_pc%v%f(i0-1:iend  ,jend+1:nend,:)) + &
    a2*(U_pc%v%f(i0+1:iend+2,jend+1:nend,:) + U_pc%v%f(i0-2:iend-1,jend+1:nend,:))

    ! Convert from latlon to contravariant
    do i = i0, iend+1
        do p = 1, nbfaces
            do j = jend+1, nend
                call ll2contra(U_pu%u%f(i,j,p), U_pu%v%f(i,j,p), &
                U_pu%ucontra%f(i,j,p), U_pu%vcontra%f(i,j,p), mesh%ll2contra_pu(i,j,p)%M)
            end do
        end do
    end do
 
    ! Now, let us interpolate the ghost cell edges - cubic interpolation
    ! Panel from south
    U_pv%u%f(n0:i0-1,j0:jend+1,:) = &
    a1*(U_pc%u%f(n0:i0-1,j0:jend+1,  :) + U_pc%u%f(n0:i0-1,j0-1:jend,  :)) + &
    a2*(U_pc%u%f(n0:i0-1,j0+1:jend+2,:) + U_pc%u%f(n0:i0-1,j0-2:jend-1,:))

    U_pv%v%f(n0:i0-1,j0:jend+1,:) = &
    a1*(U_pc%v%f(n0:i0-1,j0:jend+1  ,:) + U_pc%v%f(n0:i0-1,j0-1:jend  ,:)) + &
    a2*(U_pc%v%f(n0:i0-1,j0+1:jend+2,:) + U_pc%v%f(n0:i0-1,j0-2:jend-1,:))

    ! Convert from latlon to contravariant
    do j = j0, jend+1
        do p = 1, nbfaces
            do i = n0, i0-1
                call ll2contra(U_pv%u%f(i,j,p), U_pv%v%f(i,j,p), &
                U_pv%ucontra%f(i,j,p), U_pv%vcontra%f(i,j,p), mesh%ll2contra_pv(i,j,p)%M)
            end do
        end do
    end do

    ! Panel from north
    U_pv%u%f(iend+1:nend, j0:jend+1, :) = &
    a1*(U_pc%u%f(iend+1:nend,j0:jend+1  ,:) + U_pc%u%f(iend+1:nend,j0-1:jend  ,:)) + &
    a2*(U_pc%u%f(iend+1:nend,j0+1:jend+2,:) + U_pc%u%f(iend+1:nend,j0-2:jend-1,:))

    U_pv%v%f(iend+1:nend, j0:jend+1, :) = &
    a1*(U_pc%v%f(iend+1:nend,j0:jend+1  ,:) + U_pc%v%f(iend+1:nend,j0-1:jend  ,:)) + &
    a2*(U_pc%v%f(iend+1:nend,j0+1:jend+2,:) + U_pc%v%f(iend+1:nend,j0-2:jend-1,:))

    ! Convert from latlon to contravariant
    do j = j0, jend+1
        do p = 1, nbfaces
            do i = iend+1, nend
                call ll2contra(U_pv%u%f(i,j,p), U_pv%v%f(i,j,p), &
                U_pv%ucontra%f(i,j,p), U_pv%vcontra%f(i,j,p), mesh%ll2contra_pv(i,j,p)%M)
            end do
        end do
    end do
 

    !-----------------------------------------------------------------------------------
    ! Interpolation needed for RK2 departure point scheme
    ! Panel from west
    U_pu%u%f(i0-1,n0:nend,:) = a1*(U_pc%u%f(i0-2,n0:nend,:) + U_pc%u%f(i0-1,n0:nend,:))&
    + a2*(U_pc%u%f(i0,n0:nend,:) + U_pc%u%f(i0-3,n0:nend,:))
    U_pu%v%f(i0-1,n0:nend,:) = a1*(U_pc%v%f(i0-2,n0:nend,:) + U_pc%v%f(i0-1,n0:nend,:))&
    + a2*(U_pc%v%f(i0,n0:nend,:) + U_pc%v%f(i0-3,n0:nend,:))

    ! Convert from latlon to contravariant
    do j = n0, nend
        do p = 1, nbfaces
            call ll2contra(U_pu%u%f(i0-1,j,p), U_pu%v%f(i0-1,j,p), &
        U_pu%ucontra%f(i0-1,j,p), U_pu%vcontra%f(i0-1,j,p), mesh%ll2contra_pu(i0-1,j,p)%M)
        end do
    end do

    ! Panel from east
    U_pu%u%f(iend+2,n0:nend,:) = a1*(U_pc%u%f(iend+1,:,:) + U_pc%u%f(iend+2,:,:))&
    + a2*(U_pc%u%f(iend,:,:) + U_pc%u%f(iend+3,:,:))
    U_pu%v%f(iend+2,:,:) = a1*(U_pc%v%f(iend+1,:,:) + U_pc%v%f(iend+2,:,:))&
    + a2*(U_pc%v%f(iend,:,:) + U_pc%v%f(iend+3,:,:))

    ! Convert from latlon to contravariant
    do j = n0, nend
        do p = 1, nbfaces
            call ll2contra(U_pu%u%f(iend+2,j,p), U_pu%v%f(iend+2,j,p), &
            U_pu%ucontra%f(iend+2,j,p), U_pu%vcontra%f(iend+2,j,p), mesh%ll2contra_pu(iend+2,j,p)%M)
        end do
    end do
 
    !-----------------------------------------------------------------------------------
    ! Panel from south
    U_pv%u%f(n0:nend,j0-1,:) = a1*(U_pc%u%f(n0:nend,j0-2,:) + U_pc%u%f(n0:nend,j0-1,:))&
    + a2*(U_pc%u%f(n0:nend,j0,:) + U_pc%u%f(n0:nend,j0-3,:))
    U_pv%v%f(n0:nend,j0-1,:) = a1*(U_pc%v%f(n0:nend,j0-2,:) + U_pc%v%f(n0:nend,j0-1,:))&
    + a2*(U_pc%v%f(n0:nend,j0,:) + U_pc%v%f(n0:nend,j0-3,:))

    ! Convert from latlon to contravariant
    do i = n0, nend
        do p = 1, nbfaces
            call ll2contra(U_pv%u%f(i,j0-1,p), U_pv%v%f(i,j0-1,p), &
            U_pv%ucontra%f(i,j0-1,p), U_pv%vcontra%f(i,j0-1,p), mesh%ll2contra_pv(i,i0-1,p)%M)
        end do
    end do

    ! Panel from east
    U_pv%u%f(:,jend+2,:) = a1*(U_pc%u%f(:,jend+1,:) + U_pc%u%f(:,jend+2,:))&
    + a2*(U_pc%u%f(:,jend,:) + U_pc%u%f(:,jend+3,:))
    U_pv%v%f(:,jend+2,:) = a1*(U_pc%v%f(:,jend+1,:) + U_pc%v%f(:,jend+2,:))&
    + a2*(U_pc%v%f(:,jend,:) + U_pc%v%f(:,jend+3,:))

    ! Convert from latlon to contravariant
    do i = n0, nend
        do p = 1, nbfaces
            call ll2contra(U_pv%u%f(i,jend+2,p), U_pv%v%f(i,jend+2,p), &
            U_pv%ucontra%f(i,jend+2,p), U_pv%vcontra%f(i,jend+2,p), mesh%ll2contra_pv(i,jend+2,p)%M)
        end do
    end do
 
end subroutine dg_vf_interp



subroutine compute_lagrange_cs(L, mesh)
    !---------------------------------------------------
    ! compute the lagrange polynomials at ghost cell
    ! points of the cubed-sphere
    !--------------------------------------------------
    type(cubedsphere), intent(inout):: mesh
    type(lagrange_poly_cs), intent(inout):: L
    integer(i4) ::  i, j, p, g, d, k, h, jnearest
    real(kind=8):: dist, kk

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
        print*
        do g = 1, nghost
            do j = n0, nend
            !    print*, L%y_nodes(j,g),  L%x_nodes(j,g), mesh%dx
            end do
            !read(*,*)
        end do
        !stop
        ! Compute stencils for all ghost cells expect at corners
        print*, mesh%n
        do g = 1, nghost
            !h = g-1
            do j = n0, nend
                kk = ((L%y_nodes(j,g)-L%y_support(n0))/mesh%dx)
                L%kend(j,g) = kk+ceiling(L%order*0.5d0)
                L%k0(j,g) = L%kend(j,g) - L%order + 1
                L%kend(j,g) = L%kend(j,g)+n0
                L%k0(j,g)   = L%k0(j,g)+n0
 
                if (j>=j0 .and. j<=jend)then
                    if(L%kend(j,g)>=jend+1)then
                        L%kend(j,g) = jend
                        L%k0(j,g) = L%kend(j,g) - L%order + 1
                    else if(L%k0(j,g)<j0)then
                        L%k0(j,g) = j0
                        L%kend(j,g) = L%k0(j,g) + L%order - 1
                    end if
                end if

                if (j>jend .and. L%kend(j,g)>=nend+1)then
                    L%kend(j,g) = nend
                    L%k0(j,g) = L%kend(j,g) - L%order + 1
                end if

                if (j<j0 .and. L%k0(j,g)<n0)then
                    L%k0(j,g) = n0
                    L%kend(j,g) = L%k0(j,g) + L%order - 1
                end if
     
                !print*,L%k0(j,g)+3,L%kend(j,g)+3
            end do
            !read(*,*)
        end do


        !stop
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
                    print*, i,g, 'error3'
                    print*, L%k0(i,g), L%kend(i,g)
                    !read(*,*)
                    stop
                end if

                if(L%y_support(L%kend(i,g))<L%y_nodes(i,g))then
                    print*, i,g, 'error4'
                    !print*, L%k0(i,g), L%kend(i,g)
                    !read(*,*)
                    stop
                end if
 
                !print*, i, g, L%y_support(L%k0(i,g)), L%y_nodes(i,g), L%y_support(L%kend(i,g))
                !print*, i, g,L%k0(i,g), L%kend(i,g)
            end do
            !read(*,*)
        end do  

        ! Store in f the nearest support points used in Lagrange interpolation
        do g = 1, nghost
            do j = n0, nend
                L%f_nodes(j,g,:) = L%y_support(L%k0(j,g):L%kend(j,g))
            end do
        end do
        print*
        do g = 1, nghost
            do j = n0, nend
            !    print*, L%f_nodes(j,g,:)
            end do
            !read(*,*)
        end do

        !stop
        ! Compute the Lagrange nodes at halo region
        do g = 1, nghost
            do j = n0, nend
                do d = 1, L%order
                    call lagrange_basis(L%y_nodes(j,g), L%f_nodes(j,g,:), L%degree, d, L%p_nodes(j,g,d))
                end do
            end do
        end do

        do g = 1, nghost
            do j = n0, nend
            !    print*, L%p_nodes(j,g,:)
            end do
            !read(*,*)
        end do
        !stop
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

end module duogrid_interpolation
