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
      simulcase

  !Data structures
  use datastruct, only: &
      cubedsphere

  ! Spherical geometry
  use sphgeo, only: &
      sph2cart

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
    type(cubedsphere):: mesh
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
    read(fileunit,*)  mesh%midpos
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

    mesh%resolution = "unif"
    !mesh%resolution = "var"

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
    integer (i4):: k
    integer (i4):: l
    !integer (i4):: m
    integer (i4):: n
    integer (i4):: nnb

    print*
    print*, "Loading mesh..." 
    print*, "Not implemented yet."
    print*
    stop
  end subroutine

  subroutine meshread(mesh)
    !--------------------------------------------------------------
    !meshread
    !  Subroutine that reads the mesh nodes from a file
    !  The file must be in grid/
    !  The file name is given implicitily via mesh%readpath
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
    real (r8):: lon
    real (r8):: lat
    logical:: ifile
    integer (i4), dimension(1:6):: nbtmp, nbtmp2

    print*
    print*, "Reading mesh nodes from file..."
    print*, "Not implemented yet."
    print*
    stop
  end subroutine

end module input 

