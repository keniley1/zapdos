# This input file sets up a pipe of length L = 1 m and width W = 0.25 m
# Two immiscible fluids flow due to a pressure gradient in the pipe, going from
# p = 1e-3 at the left and p = 0 at the right.
# The boundary between the fluids is located at W/2.
# 
# A more viscous and denser fluid flows in the bottom half of the region
# (block 0). 

[GlobalParams]
  offset = 20
  #supg = true
  #pspg = true
  potential_units = kV
  use_moles = true
  convective_term = true
  #transient_term = true
  integrate_p_by_parts = false
  #supg = true
  #integrate_p_by_parts = false
  gravity = '0 0 0'
[]

[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'two_fluid_mesh.msh'
  []

  [interface1]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'interface0'
    input = geo
  []
  [interface2]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'interface1'
    input = interface1
  []
[]

[Problem]
  type = FEProblem
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  #type = Transient
  #end_time = 1e6
  #automatic_scaling = true
  #compute_scaling_once = false
  #line_search = 'basic'
  #petsc_options = '-snes_converged_reason'
  #solve_type = newton
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  #petsc_options_value = 'lu NONZERO 1.e-10 0'
  #nl_rel_tol = 1e-6
  #nl_abs_tol = 1e-10
  #dtmin = 1e-17
  #l_max_its = 100
  #nl_max_its = 20
  #steady_state_detection = true
  #steady_state_tolerance = 1e-6
  #[TimeStepper]
  #  type = IterationAdaptiveDT
  #  cutback_factor = 0.4
  #  dt = 1e-4
  #  growth_factor = 1.4
  #  optimal_iterations = 10
  #[]
  type = Steady
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
[]

[Outputs]
  [out_01]
    type = Exodus
    interval = 1
  []
[]

[Variables]
  [vel_x]
    # Velocity in radial (r) direction
    family = LAGRANGE
    order = SECOND
    block = '0 1'
    initial_condition = 0
  []
  [vel_y]
    # Velocity in axial (z) direction
    family = LAGRANGE
    order = SECOND
    block = '0 1'
    initial_condition = 0
  []
  [p]
    family = LAGRANGE
    order = FIRST
    block = '0 1'
  []
[]

[Kernels]
  # Navier-Stokes problem test
  #[x_momentum_time]
  #  type = INSMomentumTimeDerivative
  #  variable = vel_x
  #  block = '0 1'
  #[]
  #[y_momentum_time]
  #  type = INSMomentumTimeDerivative
  #  variable = vel_y
  #  block = '0 1'
  #[]
  [mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    p = p
    block = '0 1'
  []
  [x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
    block = '0 1'
  []
  [y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
    block = '0 1'
  []
[]

[BCs]
  # Left BCs -- setting pressure difference
  [pressure_boundary_top_left]
    type = DirichletBC
    variable = p
    boundary = 'upper_left'
    value = 0.001
  []
  [pressure_boundary_bottom_left]
    type = DirichletBC
    variable = p
    boundary = 'lower_left'
    value = 0.001
  []
  [pressure_boundary_top_right]
    type = DirichletBC
    variable = p
    boundary = 'upper_right'
    value = 0
  []
  [pressure_boundary_bottom_right]
    type = DirichletBC
    variable = p
    boundary = 'lower_right'
    value = 0
  []

  [u_no_slip]
    type = DirichletBC
    boundary = 'top'
    variable = vel_x
    value = 0
  []
  [v_no_slip]
    type = DirichletBC
    boundary = 'upper_left top'
    variable = vel_y
    value = 0
  []

  #############################
  # BOTTOM REGION
  #############################
  [w_vel_y_noslip]
    type = DirichletBC
    boundary = 'lower_left bottom'
    variable = vel_y
    value = 0
  []

  # Axis boundary conditions - no radial velocity (due to symmetry)
  [water_u_axis]
    type = DirichletBC
    boundary = 'bottom'
    variable = vel_x 
    value = 0
  []
[]

[Functions]
  [pressure_function]
    type = ParsedFunction
    vars = 'p0'
    vals = '0.001'
    value = 'p0*(1-x)'
  []
[]

[ICs]
  [p_ic]
    type = FunctionIC
    variable = p
    function = pressure_function
  []
[]

[Materials]
  # INS stuff
  [const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '2 3'
  []

  [const_w]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'rho mu'
    prop_values = '1 1'
  []
[]
