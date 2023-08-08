module output
!===============================================================================================
!   This module contains all the output routines
!   Routines based on iModel (https://github.com/pedrospeixoto/iModel)
!
!   Luan Santos 2022
!===============================================================================================

!Global constants
use constants, only: &
   datadir, &
   deg2rad, &
   eps, &
   eps2, &
   griddir, &
   i4, &
   pardir, &
   pi, &
   pi2, &
   pio2, &
   rad2deg, &
   showonscreen, &
   nbfaces, &
   erad, &
   n0, nend, &
   i0, iend, &
   j0, jend

!Data structures
use datastruct, only: &
   cubedsphere, &
   scalar_field, &
   simulation

use miscellaneous, only: &
   getunit

use allocation, only: &
  r8_2darray_allocation 

! Advection model vars
use advection_vars

! Diagnostics 
use diagnostics, only: &
    adv_diagnostics, &
    swm_diagnostics

! Linear algebra
use linear_algebra, only: &
    error_norm_max_rel, &
    error_norm_1_rel, &
    error_norm_2_rel, &
    norm_max, &
    norm_1, &
    norm_2

! Spherical geometry
use sphgeo, only: &
  sph2cart
 
implicit none

contains 

!---------------------------------------
!Simple header printing routine
!---------------------------------------
subroutine printheader()
    print*,"-----------------------------------------------------------"
    print*,"  Numerical Analysis on a Cubed Sphere grid                "
    print*,"                                                           "
    print*
    print*,"  Luan Santos                                              "
    print*,"  December 2022                                            "
    print*,"-----------------------------------------------------------"
    print*
end subroutine printheader
    
!---------------------------------------
!Simple header printing routine
!------------------------------------
subroutine printending()
    print*
    print*,"-----------------------------------------------------------"
    print*,"  End of program  "
    print*,"-----------------------------------------------------------"
    print*  
end subroutine printending

subroutine printmesh(mesh)
    !-------------------------------------------
    !PRINTMESH
    ! Print main mesh caracteristics on screen
    ! Do not make confusion with meshprint routine, that
    !  writes all mesh structure onto txt files
    !-------------------------------------------
    type(cubedsphere):: mesh


    print*
    select case(trim(mesh%kind))
        case("equiangular")
            print'(a)', " Mesh                        : Equiangular cubed-sphere"
        case("read")
            print'(a)', " Mesh                        : Read from file "//trim(mesh%filename)
            return
        case default
            print'(a)', " PRINTMESH ERROR: Invalid mesh kind : ", mesh%kind
            stop
    end select

    print '(a,i8)',        " N                           : ", mesh%n
    print '(a,i8)',        " Number of cells (6*N*N)     : ", mesh%nbcells
    print '(a,l8)',        " Loadable    ?                 " , mesh%loadable
    print '(a33, 3e16.8)', " Min Max Mean distances (km) : ",  &
    mesh%mindist*erad/1e3, mesh%maxdist*erad/1e3, mesh%meandist*erad/1e3
    print '(a33, 3e16.8)', " Min Max Mean areas (km^2)     : ",  &
    mesh%minarea*erad**2/1e6, mesh%maxarea*erad**2/1e6, mesh%meanarea*erad**2/1e6
    print*
    print*
return
end subroutine printmesh

subroutine meshstore(mesh, header)
    !--------------------------------------------------------------
    ! MESHSTORE
    !  Subroutine that writes on files the mesh structure
    !  It is written on binary files .dat
    !--------------------------------------------------------------
    type(cubedsphere), intent(in) :: mesh
    character (len=128)::  header

    !Names for files and numbers
    character (len=128):: filename

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: p

    !Save mesh data for future uses
    print*
    print*, "Saving generated grid in binary format (.dat)"
    print*, "Directory: ", trim(griddir)
    print*

    !---------------------------------------
    !Write vertices
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_po.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                write(iunit,*) mesh%po(i,j,p)%lat, mesh%po(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write center
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%pc(i,j,p)%lat, mesh%pc(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write mipoints pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%pu(i,j,p)%lat, mesh%pu(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write mipoints pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%pv(i,j,p)%lat, mesh%pv(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgx at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%tgx_pu(i,j,p)%v(1), mesh%tgx_pu(i,j,p)%v(2), mesh%tgx_pu(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgy at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%tgy_pu(i,j,p)%v(1), mesh%tgy_pu(i,j,p)%v(2), mesh%tgy_pu(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgx at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%tgx_pv(i,j,p)%v(1), mesh%tgx_pv(i,j,p)%v(2), mesh%tgx_pv(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgy at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%tgy_pv(i,j,p)%v(1), mesh%tgy_pv(i,j,p)%v(2), mesh%tgy_pv(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgx at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%tgx_pc(i,j,p)%v(1), mesh%tgx_pc(i,j,p)%v(2), mesh%tgx_pc(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write tgy at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%tgy_pc(i,j,p)%v(1), mesh%tgy_pc(i,j,p)%v(2), mesh%tgy_pc(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)


    !---------------------------------------
    ! Write latlon grid indexes
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ll_grid.dat"
    open(iunit, file=filename, status='replace')
    do i = 0, mesh%nlon
        do j = 0, mesh%nlat
            write(iunit,*)  mesh%ix_ll(i,j),mesh%jy_ll(i,j), mesh%panels_ll(i,j)
        end do
    end do
    close(iunit)
end subroutine meshstore

subroutine plot_scalarfield(var, mesh)
    !----------------------------------------------
    ! PLOT_SCALARFIELD
    !   Writes in 'var%name file' a uniform mesh with
    !   nlat latitudes and nlon longitudes
    !   using nearest neighbour interpolation
    !---------------------------------------------

    !Mesh structure
    type(cubedsphere), intent(in) :: mesh

    !Variable to be plotted
    type(scalar_field), intent(inout) :: var

    !Variable for uniform grid
    real (kind=8), allocatable :: varll(:, :)
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ix, jy, p ! cubedsphere coordinates
    integer (i4):: iunit

    !File name for output
    character (len=256):: filename

    ! Variable on latlon grid
    call r8_2darray_allocation(varll, 0, mesh%nlon, 0, mesh%nlat )

    filename=trim(datadir)//trim(var%name)//"_"//trim(mesh%name)//".dat"
    print*, " Plotting scalarfield  ", trim(filename)

    select case(var%pos)
    case(0) !Values on centers
        ! no interpolation is needed
        case(1) !Values on vertices
            print*, 'ERROR in plot_scalarfield: Not implemented yet'
            stop
        case(2) !Values on midpoint u
            print*, 'ERROR in plot_scalarfield: Not implemented yet'
            stop
    case(3) !Values on midpoint v
        print*, 'ERROR in plot_scalarfield: Not implemented yet'
        stop
    case default
        print*, 'ERROR in plot_scalarfield: invalid var position ', var%pos
        stop
    end select

    ! Nearest neighbour
    do i = 0, mesh%nlon
        do j = 0, mesh%nlat
            ix = mesh%ix_ll(i,j)
            jy = mesh%jy_ll(i,j)
            p  = mesh%panels_ll(i,j)
            varll(i,j) = var%f(ix,jy,p)
        end do
    end do

    !Write values on file
    call getunit(iunit)
    !print*
    open(iunit,file=filename, status='replace', access='stream', form='unformatted')

    !Write whole block to file (much faster)
    write(iunit) varll

    close(iunit)

    deallocate(varll)
    return
end subroutine plot_scalarfield

subroutine write_final_errors_adv(advsimul, mesh, filename) 
    !----------------------------------------------------------
    !  write the final errors of advection simulation in a file
    !----------------------------------------------------------
    type(simulation), intent(in) :: advsimul
    type(cubedsphere), intent(in) :: mesh

    !File name for output
    character (len=256), intent(inout) :: filename
    character (len=256):: name

    !File units
    integer (i4):: iunit
    logical::  iopen

    !File for errors
    filename=trim(datadir)//trim(filename)//".txt"
    call getunit(iunit)

    open(iunit,file=filename, status='replace')
    write(iunit, *) advsimul%linf_error
    write(iunit, *) advsimul%l1_error
    write(iunit, *) advsimul%l2_error
    write(iunit, *) advsimul%cfl
    write(iunit, *) advsimul%mass_variation
    close(iunit)
end subroutine write_final_errors_adv

subroutine write_final_errors_swm(swm_simul, mesh, filename) 
    !----------------------------------------------------------
    !  write the final errors of shallow water model simulation in a file
    !----------------------------------------------------------
    type(simulation), intent(in) :: swm_simul
    type(cubedsphere), intent(in) :: mesh

    !File name for output
    character (len=256), intent(inout) :: filename
    character (len=256):: name

    !File units
    integer (i4):: iunit
    logical::  iopen

    !File for errors
    filename=trim(datadir)//trim(filename)//".txt"
    call getunit(iunit)

    open(iunit,file=filename, status='replace')
    write(iunit, *) swm_simul%linf_error
    write(iunit, *) swm_simul%l1_error
    write(iunit, *) swm_simul%l2_error
    write(iunit, *) swm_simul%cfl
    write(iunit, *) swm_simul%mass_variation
    close(iunit)
end subroutine write_final_errors_swm



subroutine write_final_errors_interp(filename, error_q, error_u, error_v)
    !----------------------------------------------------------
    !  write the final errors of advection simulation in a file
    !----------------------------------------------------------
    real(kind=8), intent(in) :: error_q, error_u, error_v

    !File name for output
    character (len=256), intent(inout) :: filename

    !File units
    integer (i4):: iunit
    logical::  iopen

    !File for errors
    filename=trim(datadir)//trim(filename)//".txt"
    call getunit(iunit)

    open(iunit,file=filename, status='replace')
    write(iunit, *) error_q
    write(iunit, *) error_u
    write(iunit, *) error_v
    close(iunit)
end subroutine write_final_errors_interp

subroutine output_adv(mesh)
    use advection_vars
    !----------------------------------------------------------
    !  output of the advection model
    !----------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    real(kind=8) :: lat, lon
    integer(i4):: i, j, p
    character (len=60):: an

    if(advsimul%n>0 .and. (showonscreen .or. advsimul%n==advsimul%nsteps)) then
        ! Compute diagnostics
        call adv_diagnostics(advsimul, mesh, Q)

        if(advsimul%n==advsimul%nsteps)then
            advsimul%exactsolution = .true.
        end if

        ! Screen output
        print*
        print*, "Step = ", advsimul%n, " of ", advsimul%nsteps
 
        ! Compute exact solution and errors
        if (advsimul%exactsolution .or. advsimul%n==advsimul%nsteps) then
            ! this is the only case where we need to update the exact solution
            if(advsimul%ic==2)then
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(NONE) & 
                !$OMP SHARED(Q_exact, Q, Q_error, mesh, advsimul) & 
                !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
                !$OMP PRIVATE(i, j, p, lat, lon) &
                !$OMP SCHEDULE(static) 
                do p = 1, nbfaces
                    do i = i0, iend
                        do j = j0, jend
                            lat  = mesh%pc(i,j,p)%lat
                            lon  = mesh%pc(i,j,p)%lon
                            Q_exact%f(i,j,p) = qexact_adv(lat, lon, advsimul%ic, advsimul%t)
                            Q_error%f(i,j,p) = Q_exact%f(i,j,p) - Q%f(i,j,p)
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

            call compute_errors_field(Q, Q_exact, Q_error, &
            advsimul%linf_error, advsimul%l1_error, advsimul%l2_error, mesh)

            print '(a22, 3e16.8)','linf, l1, l2 errors:', &
            advsimul%linf_error, advsimul%l1_error, advsimul%l2_error

        else if(.not. advsimul%exactsolution)then
            call compute_norms_field(Q,&
            advsimul%linf_error, advsimul%l1_error, advsimul%l2_error, mesh)

            print '(a22, 3e16.8)','linf, l1, l2 Q norms:', &
            advsimul%linf_error, advsimul%l1_error, advsimul%l2_error

        end if

        print '(a22, 1e16.8)','mass change:', advsimul%mass_variation
    end if


    ! Plot the solution
    if (mod(advsimul%n,advsimul%nplot)==0 .or. advsimul%n==advsimul%nsteps .or. &
        advsimul%n==0) then
        ! Plot scalar fields
        write(an,'(i8)') advsimul%plotcounter

        Q%name = "adv_"//trim(advsimul%name)//"_Q_t"//trim(adjustl(an))
        call plot_scalarfield(Q, mesh)

        if(advsimul%exactsolution .and. advsimul%n>0)then
            if(advsimul%ic==2)then
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(NONE) & 
                !$OMP SHARED(Q_exact, Q, Q_error, mesh, advsimul) & 
                !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
                !$OMP PRIVATE(i, j, p, lat, lon) &
                !$OMP SCHEDULE(static) 
                do p = 1, nbfaces
                    do i = i0, iend
                        do j = j0, jend
                            lat  = mesh%pc(i,j,p)%lat
                            lon  = mesh%pc(i,j,p)%lon
                            Q_exact%f(i,j,p) = qexact_adv(lat, lon, advsimul%ic, advsimul%t)
                            Q_error%f(i,j,p) = Q_exact%f(i,j,p) - Q%f(i,j,p)
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end if
            Q_error%name = "adv_"//trim(advsimul%name)//"_Q_error_t"//trim(adjustl(an))
            call plot_scalarfield(Q_error, mesh)
        end if
        advsimul%plotcounter = advsimul%plotcounter + 1
 
    end if
end subroutine output_adv


subroutine output_swm(mesh)
    use swm_vars
    !----------------------------------------------------------
    !  output of the shallow water model
    !----------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    real(kind=8) :: lat, lon
    integer(i4):: i, j, p
    character (len=60):: an

    if(swm_simul%n>0 .and. (showonscreen .or. swm_simul%n==swm_simul%nsteps)) then
        ! Compute diagnostics
        call swm_diagnostics(swm_simul, mesh, H)

        if(swm_simul%n==swm_simul%nsteps)then
            swm_simul%exactsolution = .true.
        end if

        ! Screen output
        print*
        print*, "Step = ", swm_simul%n, " of ", swm_simul%nsteps
 
        ! Compute exact solution and errors
        if (swm_simul%exactsolution .or. swm_simul%n==swm_simul%nsteps) then
            ! this is the only case where we need to update the exact solution
            if(swm_simul%ic==2)then
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(NONE) & 
                !$OMP SHARED(H_exact, H, H_error, mesh, swm_simul) & 
                !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
                !$OMP PRIVATE(i, j, p, lat, lon) &
                !$OMP SCHEDULE(static) 
                do p = 1, nbfaces
                    do i = i0, iend
                        do j = j0, jend
                            lat  = mesh%pc(i,j,p)%lat
                            lon  = mesh%pc(i,j,p)%lon
                            H_exact%f(i,j,p) = qexact_adv(lat, lon, swm_simul%ic, swm_simul%t)
                            H_error%f(i,j,p) = H_exact%f(i,j,p) - H%f(i,j,p)
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

            call compute_errors_field(H, H_exact, H_error, &
            swm_simul%linf_error, swm_simul%l1_error, swm_simul%l2_error, mesh)

            print '(a22, 3e16.8)','linf, l1, l2 errors:', &
            swm_simul%linf_error, swm_simul%l1_error, swm_simul%l2_error

        else if(.not. swm_simul%exactsolution)then
            call compute_norms_field(H,&
            swm_simul%linf_error, swm_simul%l1_error, swm_simul%l2_error, mesh)

            print '(a22, 3e16.8)','linf, l1, l2 Q norms:', &
            swm_simul%linf_error, swm_simul%l1_error, swm_simul%l2_error

        end if

        print '(a22, 1e16.8)','mass change:', swm_simul%mass_variation
    end if


    ! Plot the solution
    if (mod(swm_simul%n,swm_simul%nplot)==0 .or. swm_simul%n==swm_simul%nsteps .or. &
        swm_simul%n==0) then
        ! Plot scalar fields
        write(an,'(i8)') swm_simul%plotcounter

        H%name = "adv_"//trim(swm_simul%name)//"_H_t"//trim(adjustl(an))
        call plot_scalarfield(H, mesh)

        if(swm_simul%exactsolution .and. swm_simul%n>0)then
            if(swm_simul%ic==2)then
                !$OMP PARALLEL DO &
                !$OMP DEFAULT(NONE) & 
                !$OMP SHARED(Q_exact, Q, Q_error, mesh, swm_simul) & 
                !$OMP SHARED(i0, iend, j0, jend, nbfaces) &
                !$OMP PRIVATE(i, j, p, lat, lon) &
                !$OMP SCHEDULE(static) 
                do p = 1, nbfaces
                    do i = i0, iend
                        do j = j0, jend
                            lat  = mesh%pc(i,j,p)%lat
                            lon  = mesh%pc(i,j,p)%lon
                            !Q_exact%f(i,j,p) = qexact_adv(lat, lon, swm_simul%ic, swm_simul%t)
                            !Q_error%f(i,j,p) = Q_exact%f(i,j,p) - Q%f(i,j,p)
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
            end if
            !Q_error%name = "swm_"//trim(swm_simul%name)//"_Q_error_t"//trim(adjustl(an))
            call plot_scalarfield(H_error, mesh)
        end if
        swm_simul%plotcounter = swm_simul%plotcounter + 1
 
    end if
end subroutine output_swm


subroutine compute_errors_field(Q, Q_ref, Q_error, linf, l1, l2, mesh)
    !---------------------------------------------------
    ! COMPUTE_ERRORS_FIELD
    !
    ! Given the scalar field reference field Q_ref,
    ! this routine compute the l1, l2, linf errors
    ! of Q
    !--------------------------------------------------
    type(cubedsphere),intent(in) :: mesh
    type(scalar_field), intent(inout) :: Q      ! numerical approximation
    type(scalar_field), intent(inout) :: Q_ref   ! reference solution
    type(scalar_field), intent(inout) :: Q_error ! error
    real(kind=8), intent(out) :: linf, l1, l2
    ! aux vars
    integer (i4) :: x0, xend
    integer (i4) :: y0, yend

    select case(Q%pos)
        case(0) ! pc
            x0 = i0
            y0 = j0
            xend = iend
            yend = jend

        case(1) ! po
            x0 = i0+1
            y0 = j0+1
            xend = iend+1
            yend = jend+1

        case(2) ! pu
            x0 = i0+1
            y0 = j0
            xend = iend+1
            yend = jend

        case(3) ! pv
            x0 = i0
            y0 = j0+1
            xend = iend
            yend = jend+1

        case default
            print*, 'ERROR on compute_errors_field: invalid position', Q%pos

    end select

    ! Compute the errors
    Q_error%f = Q_ref%f - Q%f
    linf = error_norm_max_rel(Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 
    l1   = error_norm_1_rel  (Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 
    l2   = error_norm_2_rel  (Q%f(x0:xend,y0:yend,:), Q_ref%f(x0:xend,y0:yend,:)) 


end subroutine compute_errors_field

subroutine compute_norms_field(Q, linf, l1, l2, mesh)
    !---------------------------------------------------
    ! COMPUTE_NORMS_FIELD
    !
    ! Given the scalar field reference field Q_ref,
    ! this routine compute the l1, l2, linf norms
    ! of Q
    !--------------------------------------------------
    type(cubedsphere),intent(in) :: mesh
    type(scalar_field), intent(inout) :: Q
    real(kind=8), intent(out) :: linf, l1, l2
    ! aux vars
    integer (i4) :: x0, xend
    integer (i4) :: y0, yend

    select case(Q%pos)
        case(0) ! pc
            x0 = i0
            y0 = j0
            xend = iend
            yend = jend

        case(1) ! po
            x0 = i0+1
            y0 = j0+1
            xend = iend+1
            yend = jend+1

        case(2) ! pu
            x0 = i0+1
            y0 = j0
            xend = iend+1
            yend = jend

        case(3) ! pv
            x0 = i0
            y0 = j0+1
            xend = iend
            yend = jend+1

        case default
            print*, 'ERROR on compute_norms_field: invalid position', Q%pos

    end select

    ! Compute the errors
    linf = norm_max(Q%f(x0:xend,y0:yend,:)) 
    l1   = norm_1(Q%f(x0:xend,y0:yend,:)) 
    l2   = norm_2(Q%f(x0:xend,y0:yend,:)) 

end subroutine compute_norms_field


function qexact_adv(lat, lon, ic, t)
    !--------------------------------------------------
    ! Compute the exact solution of the advection
    ! problem on the sphere
    ! 
    ! Possible initial conditions (ic)
    ! 1 - constant field
    ! 2 - one gaussian hill
    ! 3 - two gaussian hills
    !--------------------------------------------------
    integer(i4), intent(in) :: ic
    real(kind=8), intent(in) :: lat, lon, t
    real(kind=8) :: qexact_adv

    ! aux vars
    real(kind=8) :: alpha ! rotation angle
    real(kind=8) :: x, y, z ! r3 coordinates
    real(kind=8) :: x0, y0, z0 ! r3 coordinates of center point
    real(kind=8) :: x1, y1, z1 ! r3 coordinates of center point
    real(kind=8) :: lat0, lon0 ! latlon coordinates of center point
    real(kind=8) :: lat1, lon1 ! latlon coordinates of center point
    real(kind=8) :: b0 ! Gaussian width
    real(kind=8) :: u0, f, ws, wt
    real(kind=8) :: cosa, cos2a, sina, sin2a, coswt, sinwt
    real(kind=8) :: rotX, rotY, rotZ

    integer(i4):: m, n

    select case(ic)
        case(1) ! constant scalar field
            qexact_adv = 1.d0

        case(2) ! one Gaussian hill
            u0 =  2.d0*pi/5.d0 ! Wind speed
            ws = -u0
            wt = ws*t
            alpha = -45.d0*deg2rad ! Rotation angle

            !Rotation parameters
            cosa  = dcos(alpha)
            cos2a = cosa*cosa
            sina  = dsin(alpha)
            sin2a = sina*sina
            coswt = dcos(wt)
            sinwt = dsin(wt)

            call sph2cart(lon, lat, x, y, z)

            rotX = (coswt*cos2a+sin2a)*x -sinwt*cosa*y + (coswt*cosa*sina-cosa*sina)*z
            rotY =  sinwt*cosa*x + coswt*y + sina*sinwt*z
            rotZ = (coswt*sina*cosa-sina*cosa)*x -sinwt*sina*y + (coswt*sin2a+cos2a)*z

            ! Gaussian center
            lon0 = pi*0.25d0
            lat0 = pi/6.d0
 
            call sph2cart(lon0, lat0, x0, y0, z0)
            b0 = 10.d0
            qexact_adv = dexp(-b0*((rotX-x0)**2+ (rotY-y0)**2 + (rotZ-z0)**2))
   
        case(3) ! two Gaussian hills
            call sph2cart(lon, lat, x, y, z)
            ! Gaussian hill centers
            lon0 = -pi/6.d0
            lat0 = 0.d0
            lon1 = pi/6.d0
            lat1 = 0.d0
            call sph2cart(lon0, lat0, x0, y0, z0)
            call sph2cart(lon1, lat1, x1, y1, z1)
            b0 = 5.d0
            qexact_adv = dexp(-b0*((x-x1)**2+ (y-y1)**2 + (z-z1)**2)) + &
                     dexp(-b0*((x-x0)**2+ (y-y0)**2 + (z-z0)**2))

        case(4) ! steady state from will92
            alpha = -45.d0*deg2rad ! Rotation angle
            u0 = 2.d0*pi/5.d0     ! Wind speed
            f = (-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))
            qexact_adv = 1.0 - f*f

        case default
            print*, "ERROR on qexact_adv: invalid initial condition."
            stop
    end select

    return
end function qexact_adv

subroutine print_advparameters(advsimul)
    !-------------------------------------------
    ! PRINT_ADVPARAMETERS
    ! Prints advection parameters from file named "advection.par"
    !--------------------------------------------------
    type(simulation), intent(inout):: advsimul

    print*
    print '(a,i8)',      " Initial condition            : ", advsimul%ic
    print '(a,i8)',      " Velocity field               : ", advsimul%vf
    print '(a, 3e16.8)', " dt                           : ", advsimul%dt
    print '(a, 3e16.8)', " CFL                          : ", advsimul%cfl
    print '(a, a21)',    " 1D reconstruction            : ", advsimul%recon1d
    print '(a, a21)',    " Operator spltting            : ", advsimul%opsplit
    print '(a, a21)',    " Metric tensor                : ", advsimul%mt
    print '(a, a21)',    " Departure point              : ", advsimul%dp
    print '(a, a21)',    " Mass fixer                   : ", advsimul%dp
    print '(a, a21)',    " Edge treatment               : ", advsimul%et
    print '(a, i8)',     " Duo grid interpolation degree: ", advsimul%id
    print*
return
end subroutine print_advparameters

subroutine print_swmparameters(swm_simul)
    !-------------------------------------------
    ! PRINT_SWMPARAMETERS
    ! Prints sw parameters from file named "swm.par"
    !--------------------------------------------------
    type(simulation), intent(inout):: swm_simul

    print*
    print '(a,i8)',      " Test case                    : ", swm_simul%ic
    print '(a, 3e16.8)', " dt                           : ", swm_simul%dt
    print '(a, 3e16.8)', " CFL                          : ", swm_simul%cfl
    print '(a, a21)',    " 1D reconstruction            : ", swm_simul%recon1d
    print '(a, a21)',    " Operator spltting            : ", swm_simul%opsplit
    print '(a, a21)',    " Metric tensor                : ", swm_simul%mt
    print '(a, a21)',    " Departure point              : ", swm_simul%dp
    print '(a, a21)',    " Mass fixer                   : ", swm_simul%dp
    print '(a, a21)',    " Edge treatment               : ", swm_simul%et
    print '(a, i8)',     " Duo grid interpolation degree: ", swm_simul%id
    print*
return
end subroutine print_swmparameters


end module output
