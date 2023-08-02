module mass_fixer
!===============================================================================================
!
!
! Luan Santos  - 2023
! Module for routines dedicated to ensure mass conservation on the duogrid
!
!===============================================================================================
! Constants
use constants, only: &
  i4, &
  nbfaces, &
  i0, iend, &
  j0, jend, &
  n0, nend

!Data structures
use datastruct, only: &
  ppm_parabola, &
  scalar_field, &
  simulation, &
  cubedsphere

! Diagnostics
use diagnostics, only: &
  mass_computation
 
implicit none
contains 

subroutine average_flux_at_cube_intefaces(px, py, dx, dy, dt)
    !---------------------------------------------------------------------------------
    ! AVERAGE_FLUX_AT_CUBE_INTERFACES
    !
    ! Mass fixer that average the flux values at cube interfaces
    !--------------------------------------------------------------------------------
    type(ppm_parabola), intent(inout) :: px ! parabola
    type(ppm_parabola), intent(inout) :: py ! parabola
    real(kind=8), intent(in) :: dx, dy, dt ! parabola
    real(kind=8) :: a, b

    a = 0.5d0
    b = 0.5d0
    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(px, py, i0, iend, j0, jend, dx, dy, dt, a, b)
    ! Average panels 1-2,2-3,3-4,4-1
    px%f_upw(iend+1,j0:jend,1:3) = a*px%f_upw(iend+1,j0:jend,1:3) + b*px%f_upw(i0,j0:jend,2:4)
    px%f_upw(i0,j0:jend,2:4)     = px%f_upw(iend+1,j0:jend,1:3)

    px%f_upw(iend+1,j0:jend,4) = a*px%f_upw(iend+1,j0:jend,4) + b*px%f_upw(i0,j0:jend,1)
    px%f_upw(i0,j0:jend,1)  = px%f_upw(iend+1,j0:jend,4)

    ! Average panels 1-5
    py%f_upw(i0:iend,j0,5)   = a*py%f_upw(i0:iend,j0,5) + b*py%f_upw(i0:iend,jend+1,1)
    py%f_upw(i0:iend,jend+1,1) = py%f_upw(i0:iend,j0,5)

    ! Average panels 2-5
    px%f_upw(iend+1,j0:jend,5) = a*px%f_upw(iend+1,j0:jend,5) - b*py%f_upw(i0:iend,jend+1,2)
    py%f_upw(i0:iend,jend+1,2) = -px%f_upw(iend+1,j0:jend,5)

    ! Average panels 3-5
    py%f_upw(i0:iend,jend+1,5) = a*py%f_upw(i0:iend,jend+1,5) - b*py%f_upw(iend:i0:-1,jend+1,3)
    py%f_upw(i0:iend,jend+1,3) = -py%f_upw(iend:i0:-1,jend+1,5)

    ! Average panels 4-5
    px%f_upw(i0,j0:jend,5)    = a*px%f_upw(i0,j0:jend,5) + b*(py%f_upw(iend:i0:-1,jend+1,4))
    py%f_upw(i0:iend,jend+1,4) = (px%f_upw(i0,jend:j0:-1,5))

    ! Average panels 1-6
    py%f_upw(i0:iend,jend+1,6) = a*py%f_upw(i0:iend,jend+1,6) + b*py%f_upw(i0:iend,j0,1)
    py%f_upw(i0:iend,j0,1)   = py%f_upw(i0:iend,jend+1,6)

    ! Average panels 2-6
    py%f_upw(i0:iend,j0,2)     = a*py%f_upw(i0:iend,j0,2) + b*px%f_upw(iend+1,jend:j0:-1,6)
    px%f_upw(iend+1,j0:jend,6) = py%f_upw(iend:i0:-1,j0,2)

    ! Average panels 3-6
    py%f_upw(i0:iend,j0,3) = a*py%f_upw(i0:iend,j0,3) - b*py%f_upw(iend:i0:-1,j0,6)
    py%f_upw(i0:iend,j0,6) = -py%f_upw(iend:i0:-1,j0,3)

    ! Average panels 4-6
    py%f_upw(i0:iend,j0,4) = -a*px%f_upw(i0,j0:jend,6) + b*py%f_upw(i0:iend,j0,4)
    px%f_upw(i0,j0:jend,6) = -py%f_upw(i0:iend,j0,4)

    ! Update divergence
    px%df(i0,:,:)   = -(dt/dx)*(px%f_upw(i0+1,:,:)   - px%f_upw(i0,:,:))
    px%df(iend,:,:) = -(dt/dx)*(px%f_upw(iend+1,:,:) - px%f_upw(iend,:,:))
    py%df(:,j0,:)   = -(dt/dy)*(py%f_upw(:,j0+1,:)   - py%f_upw(:,j0,:))
    py%df(:,jend,:) = -(dt/dy)*(py%f_upw(:,jend+1,:) - py%f_upw(:,jend,:))

    !$OMP END PARALLEL WORKSHARE
end subroutine average_flux_at_cube_intefaces


subroutine divergence_projection(div_ugq,  advsimul, mesh)
    !---------------------------------------------------
    ! Uses the L2 divergence projection on the zero average grid
    ! function space to ensure mass conservation
    !---------------------------------------------------
    type(cubedsphere), intent(inout) :: mesh
    type(simulation), intent(inout) :: advsimul
    type(scalar_field), intent(inout) :: div_ugq
    real(kind=8) :: l

    advsimul%mass = mass_computation(div_ugq, mesh)

    ! sum of metric tensor^2 of cells at boundary
    l = advsimul%mass/advsimul%a2

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(div_ugq, l, i0, iend, j0, jend, advsimul, mesh)
    !div_ugq%f(i0  ,j0:jend,:) = div_ugq%f(i0  ,j0:jend,:) - l*mesh%mt_pc(i0  ,j0:jend,:)
    !div_ugq%f(iend,j0:jend,:) = div_ugq%f(iend,j0:jend,:) - l*mesh%mt_pc(iend,j0:jend,:)
    !div_ugq%f(i0+1:iend-1,j0  ,:) = div_ugq%f(i0+1:iend-1,j0  ,:) - l*mesh%mt_pc(i0+1:iend-1,j0  ,:)
    !div_ugq%f(i0+1:iend-1,jend,:) = div_ugq%f(i0+1:iend-1,jend,:) - l*mesh%mt_pc(i0+1:iend-1,jend,:)
    div_ugq%f(i0:iend,j0:jend,:) = div_ugq%f(i0:iend,j0:jend,:) - &
    mesh%mt_pc(i0:iend,j0:jend,:)*l
    !$OMP END PARALLEL WORKSHARE

end subroutine divergence_projection


end module mass_fixer
