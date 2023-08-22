module input
!========================================================================
!
! Module for input routines
!
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
  pio2, &
  n0, nend, nbfaces

!Data structures
use datastruct, only: &
  cubedsphere, &
  simulation

! Spherical geometry
use sphgeo, only: &
  sph2cart, &
  cart2sph

! Miscellaneous
use miscellaneous, only: &
  getunit

implicit none

contains 

subroutine getparameters(mesh)
    !---------------------------------------------------
    ! GETPARAMETERS
    !    Reads mesh parameters from file named "mesh.par"
    !    Saves parameters on mesh structure
    !--------------------------------------------------
    type(cubedsphere), intent(inout):: mesh
    character (len=60):: filename
    character (len=300):: buffer
    integer (i4):: fileunit
    integer:: i
    integer:: n

    !Variables for line argument reading
    character(len=60):: argv
    integer:: iargc
    integer:: nargs

    !Standard definition of the mesh caracteristics    
    mesh%loadable=0

    !Standard mesh parameters file
    filename=trim(pardir)//"mesh.par"

    n=0
    nargs=iargc()
    select case (nargs)
    case(0) !If no arguments are given, read file "mesh.par"
    case(1) !If a file is given, use it
        call getarg(1, argv)
        write(filename,'(a)') argv
    case(2) !If a file is given, and a number of mesh points
        call getarg(1, argv)
        write(filename,'(a)') argv
        !This number of grid points will overlap the file given one
        call getarg(2, argv)
        read (argv, *) n
    end select
    print*,"Mesh parameters: ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%n
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%kind
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%loadable
    read(fileunit,*)  buffer
    read(fileunit,*)  i  !showonscreen
    read(fileunit,*)  buffer
    read(fileunit,*)  simulcase
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%name

    close(fileunit)

    if(i==1) showonscreen=.true.
    if(n>0) mesh%n=n
    return
end subroutine getparameters

subroutine meshload(mesh, header)
    !--------------------------------------------------------------
    !meshload
    !  Subroutine that reads the mesh structure from grid files
    !--------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    character (len=128), intent(in)::  header

    !Names for files and numbers
    character (len=128):: filename
    character (len=16):: buffer

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: p

    print*
    print*, "Loading mesh..." 
    print*

    ! Local coordinates grid size
    if (mesh%kind=='equiangular') then
        mesh%dx = pio2/mesh%n
        mesh%dy = pio2/mesh%n
    else
        mesh%dx = 2.d0/mesh%n
        mesh%dy = 2.d0/mesh%n
    end if

    !---------------------------------------
    ! Read vertices
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_po.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                read(iunit,*) mesh%po(i,j,p)%lat, mesh%po(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read center
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pc.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                read(iunit,*) mesh%pc(i,j,p)%lat, mesh%pc(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read mipoints pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pu.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                read(iunit,*) mesh%pu(i,j,p)%lat, mesh%pu(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read mipoints pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_pv.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                read(iunit,*) mesh%pv(i,j,p)%lat, mesh%pv(i,j,p)%lon
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgx at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pu.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                read(iunit,*) mesh%tgx_pu(i,j,p)%v(1), mesh%tgx_pu(i,j,p)%v(2), mesh%tgx_pu(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgy at pu
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pu.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                read(iunit,*) mesh%tgy_pu(i,j,p)%v(1), mesh%tgy_pu(i,j,p)%v(2), mesh%tgy_pu(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgx at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pv.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                read(iunit,*) mesh%tgx_pv(i,j,p)%v(1), mesh%tgx_pv(i,j,p)%v(2), mesh%tgx_pv(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgy at pv
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pv.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                read(iunit,*) mesh%tgy_pv(i,j,p)%v(1), mesh%tgy_pv(i,j,p)%v(2), mesh%tgy_pv(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgx at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_pc.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                read(iunit,*) mesh%tgx_pc(i,j,p)%v(1), mesh%tgx_pc(i,j,p)%v(2), mesh%tgx_pc(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgy at pc
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_pc.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                read(iunit,*) mesh%tgy_pc(i,j,p)%v(1), mesh%tgy_pc(i,j,p)%v(2), mesh%tgy_pc(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgx at po
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgx_po.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                read(iunit,*) mesh%tgx_po(i,j,p)%v(1), mesh%tgx_po(i,j,p)%v(2), mesh%tgx_po(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read tgy at po
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tgy_po.dat"
    open(iunit, file=filename, status='old')
    ! Read coordinates
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                read(iunit,*) mesh%tgy_po(i,j,p)%v(1), mesh%tgy_po(i,j,p)%v(2), mesh%tgy_po(i,j,p)%v(3)
            end do 
        end do
    end do
    close(iunit)

    !---------------------------------------
    ! Read latlon grid indexes
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ll_grid.dat"
    open(iunit, file=filename, status='old')
    do i = 0, mesh%nlon
        do j = 0, mesh%nlat
            read(iunit,*)  mesh%ix_ll(i,j),mesh%jy_ll(i,j), mesh%panels_ll(i,j)
        end do
    end do
    close(iunit)

    print*, "Converting latlon to xyz..."


    !---------------------------------------
    ! Convert po to xyz
    !---------------------------------------
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend+1
                call sph2cart(mesh%po(i,j,p)%lon, mesh%po(i,j,p)%lat, &
                mesh%po(i,j,p)%p(1), mesh%po(i,j,p)%p(2), mesh%po(i,j,p)%p(3))
            end do
        end do
    end do

    !---------------------------------------
    ! Convert pc to xyz
    !---------------------------------------
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend
                call sph2cart(mesh%pc(i,j,p)%lon , mesh%pc(i,j,p)%lat, &
                mesh%pc(i,j,p)%p(1), mesh%pc(i,j,p)%p(2), mesh%pc(i,j,p)%p(3))
            end do
        end do
    end do

    !---------------------------------------
    ! Convert pu to xyz
    !---------------------------------------
    do p = 1, nbfaces
        do i = n0, nend+1
            do j = n0, nend
                call sph2cart(mesh%pu(i,j,p)%lon, mesh%pu(i,j,p)%lat, &
                mesh%pu(i,j,p)%p(1), mesh%pu(i,j,p)%p(2), mesh%pu(i,j,p)%p(3))
            end do
        end do
    end do
 
    !---------------------------------------
    ! Convert pv to xyz
    !---------------------------------------
    do p = 1, nbfaces
        do i = n0, nend
            do j = n0, nend+1
                call sph2cart(mesh%pv(i,j,p)%lon, mesh%pv(i,j,p)%lat, &
                mesh%pv(i,j,p)%p(1), mesh%pv(i,j,p)%p(2), mesh%pv(i,j,p)%p(3))
            end do
        end do
    end do


end subroutine

subroutine meshread(mesh)
    !--------------------------------------------------------------
    !meshread
    !  Subroutine that reads the mesh nodes from a file
    !  The file must be in grid/
    !--------------------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh

    !Names for files and numbers
    character (len=128):: filename
    character (len=128):: tmp
    character (len=16):: ext
    !character (len=256):: buffer

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: n
    integer (i4):: m
    integer (i4):: ios
    real (kind=8):: lon
    real (kind=8):: lat
    logical:: ifile
    integer (i4), dimension(1:6):: nbtmp, nbtmp2

    print*
    print*, "Reading mesh nodes from file..."
    print*, "Not implemented yet."
    print*
    stop
end subroutine

subroutine getadvparameters(advsimul)
    !---------------------------------------------------
    ! GETADVPARAMETERS
    !    Reads advection parameters from file named "advection.par"
    !--------------------------------------------------
    type(simulation), intent(inout):: advsimul
    character (len=60):: filename
    character (len=300):: buffer
    integer (i4):: fileunit
    integer:: i
    integer:: n

    !Standard advection parameters file
    filename=trim(pardir)//"advection.par"

    print*,"Advection parameters: ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%ic
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%vf
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%dt
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%nplot
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%recon1d
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%opsplit
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%mt
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%dp
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%mf
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%et
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%id

    close(fileunit)
    return
end subroutine getadvparameters


subroutine getswmparameters(swm_simul)
    !---------------------------------------------------
    ! GETSWMPARAMETERS
    ! Reads swm parameters from file named "swm.par"
    !--------------------------------------------------
    type(simulation), intent(inout):: swm_simul
    character (len=60):: filename
    character (len=300):: buffer
    integer (i4):: fileunit
    integer:: i
    integer:: n

    !Standard swmection parameters file
    filename=trim(pardir)//"swm.par"

    print*,"Shallow water parameters: ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%ic
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%tf
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%dt
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%nplot
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%recon1d
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%opsplit
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%mt
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%dp
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%mf
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%et
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%id
    read(fileunit,*)  buffer
    read(fileunit,*)  swm_simul%avd


    close(fileunit)
    return
end subroutine getswmparameters


subroutine getinterp_parameters(advsimul)
    !---------------------------------------------------
    ! GETINTERP_PARAMETERS
    !    Reads advection parameters from file named "advection.par"
    !--------------------------------------------------
    type(simulation), intent(inout):: advsimul
    character (len=60):: filename
    character (len=300):: buffer
    integer (i4):: fileunit
    integer:: i
    integer:: n

    !Standard advection parameters file
    filename=trim(pardir)//"interpolation.par"

    print*,"Interpolation parameters: ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%ic
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%vf
    read(fileunit,*)  buffer
    read(fileunit,*)  advsimul%id
    read(fileunit,*)  buffer
 
    close(fileunit)
    return
end subroutine getinterp_parameters

end module input 

