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

!Simulation routines
use simulpack, only: &
   grid_quality, &
   div_test

!Deallocation routines
use deallocation, only: &
   meshdeallocation

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
    case(1) ! Grid generation, storage and quality
        call grid_quality(mesh)

    case(2) ! Divergence test
        call div_test(mesh)

    case default
        print*, "Please select a proper simulation case ...:", simulcase
end select

! Deallocation
call meshdeallocation(mesh)

!Print finishing line
call printending()

end program main

