module advection_vars
!===============================================================================================
! Module the variables in the advection model on the sphere
!
! Luan da Fonseca Santos - 2022
! (luan.santos@usp.br)
!===============================================================================================

! Data structures
use datastruct, only: &
  scalar_field, &
  velocity_field, &
  ppm_parabola, &
  lagrange_poly_cs, &
  simulation

! Advection class
type(simulation) :: advsimul

! Lagrange polynomials at ghost cell centers
type(lagrange_poly_cs) :: L_pc 

! Scalar fields
type(scalar_field) :: Q       ! average values of the advected scalar field
type(scalar_field) :: gQ      ! metrictensor*Q
type(scalar_field) :: Q_exact ! exact solution
type(scalar_field) :: Q_error ! error

! Vector field
type(velocity_field) :: wind_pu, wind_pv, wind_pc, wind_po

! CFL
type(scalar_field) :: cx_pu ! cfl x direction at pu
type(scalar_field) :: cy_pv ! cfl y direction at pv

! Divergence
type(scalar_field) :: div_ugq
type(scalar_field) :: div_ugq_exact
type(scalar_field) :: div_ugq_error

! Dimension splitting vars
type(scalar_field) :: Qx
type(scalar_field) :: Qy

! PPM fields
type(ppm_parabola) :: px ! ppm in x direction
type(ppm_parabola) :: py ! ppm in y direction

contains 

end module advection_vars
