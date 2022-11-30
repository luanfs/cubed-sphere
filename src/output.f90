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
       nbfaces

  !Data structures
  use datastruct, only: &
       cubedsphere

  use miscellaneous, only: &
       getunit

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

    !print*
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
    integer (i4):: i0, iend, j0, jend

    !Save mesh data for future uses
    print*
    print*, "Saving generated grid in binary format (.dat)"
    print*, "Directory: ", trim(griddir)
    print*


    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    !---------------------------------------
    !Write vertices
    !---------------------------------------
    call getunit(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_vert.dat"
    open(iunit, file=filename, status='replace', access='stream', form='unformatted')

    !Write coordinates
    do p = 1, nbfaces
      do i = i0-1, iend
        do j = j0-1, jend
          write(iunit) mesh%po(i,j,p)%lat, mesh%po(i,j,p)%lon
        end do 
      end do
    end do
    close(iunit)

    !---------------------------------------
    !Write center
    !---------------------------------------
    call getunit(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_center.dat"
    open(iunit, file=filename, status='replace', access='stream', form='unformatted')

    !Write coordinates
    do p = 1, nbfaces
      do i = i0, iend
        do j = j0, jend
          write(iunit) mesh%pc(i,j,p)%lat, mesh%pc(i,j,p)%lon
        end do 
      end do
    end do

    close(iunit)

    !---------------------------------------
    !Write mipoints pu
    !---------------------------------------
    call getunit(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_midu.dat"
    open(iunit, file=filename, status='replace', access='stream', form='unformatted')

    !Write coordinates
    do p = 1, nbfaces
      do i = i0-1, iend
        do j = j0, jend
          write(iunit) mesh%pu(i,j,p)%lat, mesh%pu(i,j,p)%lon
        end do 
      end do
    end do

    close(iunit)

    !---------------------------------------
    !Write mipoints pv
    !---------------------------------------
    call getunit(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_midv.dat"
    open(iunit, file=filename, status='replace', access='stream', form='unformatted')

    !Write coordinates
    do p = 1, nbfaces
      do i = i0, iend
        do j = j0-1, jend
          write(iunit) mesh%pv(i,j,p)%lat, mesh%pv(i,j,p)%lon
        end do 
      end do
    end do

    close(iunit)



  end subroutine meshstore

  subroutine meshwrite(mesh, header)
   type(cubedsphere), intent(in) :: mesh
    character (len=128), intent(in)::  header


  end subroutine meshwrite

end module output
