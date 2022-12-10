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
      vector_field, &
      ppm_parabola, &
      simulation

  ! Advection class
  type(simulation):: advsimul

  ! Scalar fields
  type(scalar_field) :: Q      ! average values of the advected scalar field
  type(scalar_field) :: Q_exact ! exact solution

  ! Vector field
  type(vector_field) :: wind_pu, wind_pv

  ! CFL
  type(scalar_field) :: cx_pu ! cfl x direction at pu
  type(scalar_field) :: cy_pv ! cfl y direction at pv

  ! Divergence
  type(scalar_field) :: div_ugq
  type(scalar_field) :: div_ugq_exact
  type(scalar_field) :: div_ugq_error

  ! Metric tensor x Q
  type(scalar_field) :: gQ

  ! Dimension splitting vars
  type(scalar_field) :: Qx
  type(scalar_field) :: Qy
  type(scalar_field) :: F_gQ
  type(scalar_field) :: G_gQ
  type(scalar_field) :: FG_gQ
  type(scalar_field) :: GF_gQ

  ! PPM fields
  type(ppm_parabola) :: px ! ppm in x direction
  type(ppm_parabola) :: py ! ppm in y direction
contains 

end module advection_vars
