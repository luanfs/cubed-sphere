program main
  !--------------------------------------------------
  ! Cubed sphere main program
  !
  ! based on iModel (https://github.com/pedrospeixoto/iModel)
  !
  ! Luan Santos 2022
  !--------------------------------------------------

  !Global constants
  use constants, only: &
       simulcase

  !Data structures
  use datastruct, only: &
       cubedsphere

  !Main mesh build
  use cubed_sphere, only: &
       meshbuild

  !Output routines
  use output, only: &
       printending, &
       printheader

  !Input routines
  use input, only: &
       getparameters

  !Variable declaration
  implicit none

  type(cubedsphere) :: mesh

  !Print a header on screen
  call printheader()

  !Create/Load mesh

  !Read mesh user defined parameters and simulation case
  call getparameters(mesh)

  !Call mesh generation algorithm
  call meshbuild(mesh)

  !Do a simulation/test with the mesh loaded
  select case(simulcase)
  case(1) !Grid generation and storage
    print*
  case default
    print*, "Please select a proper simulation case ...:", simulcase
  end select
  !Print finishing line
  call printending()

end program main

