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
   r8, &
   r4, &
   r16, &
   rad2deg, &
   showonscreen, &
   nbfaces, &
   erad, &
   n0, nend

!Data structures
use datastruct, only: &
   cubedsphere, &
   scalar_field, &
   simulation

use miscellaneous, only: &
   getunit

use allocation, only: &
  r8_2darray_allocation 

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
    mesh%mindist*erad/1e3_r8, mesh%maxdist*erad/1e3_r8, mesh%meandist*erad/1e3_r8
    print '(a33, 3e16.8)', " Min Max Mean areas (km^2)     : ",  &
    mesh%minarea*erad**2/1e6_r8, mesh%maxarea*erad**2/1e6_r8, mesh%meanarea*erad**2/1e6_r8
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
    !Write metric tensor at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_mt_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%mt_pc(i,j,p)
            end do 
        end do
    end do

    !---------------------------------------
    !Write metric tensor at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_mt_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%mt_pu(i,j,p)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write metric tensor at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_mt_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%mt_pv(i,j,p)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write ll2contra at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ll2contra_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%ll2contra_pu(i,j,p)%M(1,1), mesh%ll2contra_pu(i,j,p)%M(1,2), &
                             mesh%ll2contra_pu(i,j,p)%M(2,1), mesh%ll2contra_pu(i,j,p)%M(2,2)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write contra2ll at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_contra2ll_pu.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                write(iunit,*) mesh%contra2ll_pu(i,j,p)%M(1,1), mesh%contra2ll_pu(i,j,p)%M(1,2), &
                             mesh%contra2ll_pu(i,j,p)%M(2,1), mesh%contra2ll_pu(i,j,p)%M(2,2)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write ll2contra at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ll2contra_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%ll2contra_pv(i,j,p)%M(1,1), mesh%ll2contra_pv(i,j,p)%M(1,2), &
                             mesh%ll2contra_pv(i,j,p)%M(2,1), mesh%ll2contra_pv(i,j,p)%M(2,2)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write contra2ll at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_contra2ll_pv.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                write(iunit,*) mesh%contra2ll_pv(i,j,p)%M(1,1), mesh%contra2ll_pv(i,j,p)%M(1,2), &
                             mesh%contra2ll_pv(i,j,p)%M(2,1), mesh%contra2ll_pv(i,j,p)%M(2,2)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write ll2contra at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ll2contra_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%ll2contra_pc(i,j,p)%M(1,1), mesh%ll2contra_pc(i,j,p)%M(1,2), &
                             mesh%ll2contra_pc(i,j,p)%M(2,1), mesh%ll2contra_pc(i,j,p)%M(2,2)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    !Write contra2ll at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_contra2ll_pc.dat"
    open(iunit, file=filename, status='replace')
    !Write coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                write(iunit,*) mesh%contra2ll_pc(i,j,p)%M(1,1), mesh%contra2ll_pc(i,j,p)%M(1,2), &
                             mesh%contra2ll_pc(i,j,p)%M(2,1), mesh%contra2ll_pc(i,j,p)%M(2,2)
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
    real (r8), allocatable :: varll(:, :)
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

subroutine write_final_errors(advsimul, mesh, filename) 
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
end subroutine write_final_errors


end module output
