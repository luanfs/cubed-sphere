module diagnostics
  !===============================================================================================
  !   This module contains all the routines to compute the diagnostics
  !   for the advection and shallow-water models
  !===============================================================================================

  !Global constants
  use constants, only: &
        i4, &
        r8, &
        nbfaces

  ! Data structures
  use datastruct, only: &
      cubedsphere, &
      scalar_field, &
      simulation

 implicit none

contains 

  function mass_computation(Q, mesh)
    type(cubedsphere)  :: mesh
    type(scalar_field) :: Q
    real(r8) :: mass_computation
    integer :: i0, iend, j0, jend

    ! Interior of panel grid
    i0 = mesh%i0
    j0 = mesh%j0
    iend = mesh%iend
    jend = mesh%jend

    !---------------------------------------------------------
    ! Computes the mass of the scalar field Q
    !---------------------------------------------------------
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(Q,mesh) &
    !$OMP SHARED(i0,iend,j0,jend)
    mass_computation = sum(Q%f(i0:iend,j0:jend,:)*mesh%area(i0:iend,j0:jend,:))
    !$OMP END PARALLEL WORKSHARE


  end function mass_computation

  subroutine adv_diagnostics(advsimul, mesh, Q)
    !---------------------------------------------------------
    ! Computes the diagnostics for the advection model
    !---------------------------------------------------------
    type(cubedsphere),  intent(in) :: mesh
    type(simulation),   intent(inout) :: advsimul
    type(scalar_field), intent(in) :: Q

    ! compute mass
    advsimul%mass = mass_computation(Q,mesh)

    ! mass change
    advsimul%mass_variation = abs(advsimul%mass-advsimul%mass0)/abs(advsimul%mass0)
  end subroutine adv_diagnostics


end module diagnostics
