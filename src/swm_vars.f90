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
  velocity_field, &
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

! Gradient of h
type(scalar_field) :: dy_H_pu
type(scalar_field) :: dx_H_pv

! Gradient of divuh
type(scalar_field) :: dy_div_ugh_pu
type(scalar_field) :: dx_div_ugh_pv

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

! Kinectic energy
type(scalar_field) :: Ku_po ! u part at po
type(scalar_field) :: Kv_po ! v part at po
type(scalar_field) :: dy_K_pu ! derivative in x direction of (Ku^2+Kv^2)/2
type(scalar_field) :: dx_K_pv ! derivative in y direction of (Ku^2+Kv^2)/2
type(scalar_field) :: Ku_po_exact ! u part at po
type(scalar_field) :: Kv_po_exact ! v part at po


! Vector field
type(velocity_field) :: wind_pu, wind_pv, wind_pc, wind_po

! CFL
type(scalar_field) :: cx_pu ! cfl x direction at pu
type(scalar_field) :: cy_pv ! cfl y direction at pv

! Divergence
type(scalar_field) :: div_ugH
type(scalar_field) :: div_ugH_exact
type(scalar_field) :: div_ugH_error
type(scalar_field) :: div_ugH_pu
type(scalar_field) :: div_ugH_pv
type(scalar_field) :: div_ugH_po

! Dimension splitting vars
type(scalar_field) :: Qx
type(scalar_field) :: Qy

! PPM fields
type(ppm_parabola) :: px ! ppm in x direction
type(ppm_parabola) :: py ! ppm in y direction
type(ppm_parabola) :: Kv_px ! ppm in x direction for kinetic energy
type(ppm_parabola) :: Ku_py ! ppm in y direction for kinetic energy



contains 

end module swm_vars
