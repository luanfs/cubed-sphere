module swm_vars
!===============================================================================================
! Module the variables in the shallow water model on the sphere
!
! Luan da Fonseca Santos - 2023
! (luan.santos@usp.br)
!===============================================================================================

! Data structures
use datastruct, only: &
  scalar_field, &
  vector_field, &
  ppm_parabola, &
  lagrange_poly_cs, &
  simulation

! Simulation class
type(simulation) :: swm_simul

! Lagrange polynomials at ghost cell centers
type(lagrange_poly_cs) :: L_pc 

! Scalar fields
type(scalar_field) :: H       ! average values of the fluid depth
type(scalar_field) :: H_exact ! exact fluid depth
type(scalar_field) :: H_error ! fluid depth error

! Vector field
type(vector_field) :: wind_pu, wind_pv, wind_pc

! CFL
type(scalar_field) :: cx_pu ! cfl x direction at pu
type(scalar_field) :: cy_pv ! cfl y direction at pv

! Divergence
type(scalar_field) :: div_ugH
type(scalar_field) :: div_ugH_exact
type(scalar_field) :: div_ugH_error

! Dimension splitting vars
type(scalar_field) :: Qx
type(scalar_field) :: Qy

! PPM fields
type(ppm_parabola) :: px ! ppm in x direction
type(ppm_parabola) :: py ! ppm in y direction

contains 

end module swm_vars
