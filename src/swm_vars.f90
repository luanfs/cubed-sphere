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

! Scalar fields for depth
type(scalar_field) :: H       ! average values of the fluid depth
type(scalar_field) :: H_exact ! exact fluid depth
type(scalar_field) :: H_error ! fluid depth error
type(scalar_field) :: H_po
type(scalar_field) :: H_pu
type(scalar_field) :: H_pv
type(scalar_field) :: H_pu_exact
type(scalar_field) :: H_pv_exact
type(scalar_field) :: H_po_exact

! Coriolis force at pc
type(scalar_field) :: fcoriolis_pc

! Relative vorticity
type(scalar_field) :: rel_vort
type(scalar_field) :: rel_vort_exact
type(scalar_field) :: rel_vort_error

! Absolute vorticity
type(scalar_field) :: abs_vort
type(scalar_field) :: abs_vort_exact
type(scalar_field) :: abs_vort_error
type(scalar_field) :: div_abs_vort


! Absolute vorticity fluxes
type(scalar_field) :: abs_vort_flux_pu
type(scalar_field) :: abs_vort_flux_pv
type(scalar_field) :: abs_vort_flux_exact_pu
type(scalar_field) :: abs_vort_flux_exact_pv
type(scalar_field) :: abs_vort_flux_error_pu
type(scalar_field) :: abs_vort_flux_error_pv


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

! logical var for duogrid interpolation
logical :: dginterp

contains 

end module swm_vars
