dom0Scale=1.0
dom1Scale=1.0
dom2Scale=1.0

[GlobalParams]
  offset = 20
  #supg = true
  #pspg = true
  potential_units = kV
  use_moles = true
  convective_term = true
  integrate_p_by_parts = true
  gravity = '0 0 0'
[]

[Mesh]
  [./geo]
    type = FileMeshGenerator
    file = 'sankaran_pin.msh'
  [../]

  [./interface1]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'gas_bottom'
    input = geo
  [../]
  [./interface2]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'water_top'
    input = interface1
  [../]

  # The next two definitions create boundary conditions named
  # 'left' and 'right', where 'left' is at x = 0 and 'right' is at x = 1.1 mm.
  #[./left]
  #  type = SideSetsFromNormalsGenerator
  #  normals = '-1 0 0'
  #  new_boundary = 'left'
  #  input = interface2
  #[../]
  #[./right]
  #  type = SideSetsFromNormalsGenerator
  #  normals = '1 0 0'
  #  new_boundary = 'right'
  #  input = left
  #[../]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    #ksp_norm = none
  [../]
[]

[Executioner]
  #type = Steady
  #automatic_scaling = true
  type = Transient
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  dtmin = 1e-17
  l_max_its = 100
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-6
    growth_factor = 1.2
    optimal_iterations = 10
  [../]
[]

[Outputs]
  #[out_no_water_convection]
  [out_04]
    type = Exodus
    interval = 1
  []
[]

[Debug]
  #show_var_residual_norms = true
[]

[Variables]
  [./vel_x]
    # Velocity in radial (r) direction
    family = LAGRANGE
    order = SECOND
    block = 0
    initial_condition = 0
  [../]
  [./vel_y]
    # Velocity in axial (z) direction
    family = LAGRANGE
    order = SECOND
    block = 0
    initial_condition = 0
  [../]
  [./p]
    family = LAGRANGE
    order = FIRST
    block = 0
    #initial_condition = 101325
    initial_condition = 0
  [../]
  [./w_vel_x]
    # Velocity in radial (r) direction
    family = LAGRANGE
    order = SECOND
    block = 1
    initial_condition = 0
  [../]
  [./w_vel_y]
    # Velocity in axial (z) direction
    family = LAGRANGE
    order = SECOND
    block = 1
    initial_condition = 0
  [../]
  [./w_p]
    family = LAGRANGE
    order = FIRST
    block = 1
    #initial_condition = 101325
    initial_condition = 0
  [../]
  [T]
    initial_condition = 300
  []
[]

[Kernels]
  # Navier-Stokes problem test
  [x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  []
  [y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  []
  [mass]
    type = INSMassRZ
    variable = p
    u = vel_x
    v = vel_y
    p = p
  []
  [x_momentum_space]
    type = INSMomentumLaplaceFormRZ
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  []
  [y_momentum_space]
    type = INSMomentumLaplaceFormRZ
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  []

  [x_momentum_time_w]
    type = INSMomentumTimeDerivative
    variable = w_vel_x
    block = 1 
  []
  [y_momentum_time_w]
    type = INSMomentumTimeDerivative
    variable = w_vel_y
    block = 1 
  []
  [mass_w]
    type = INSMassRZ
    variable = w_p
    u = w_vel_x
    v = w_vel_y
    p = w_p
    block = 1 
  []
  [x_momentum_space_w]
    type = INSMomentumLaplaceFormRZ
    variable = w_vel_x
    u = w_vel_x
    v = w_vel_y
    p = w_p
    component = 0
    block = 1 
  []
  [y_momentum_space_w]
    type = INSMomentumLaplaceFormRZ
    variable = w_vel_y
    u = w_vel_x
    v = w_vel_y
    p = w_p
    component = 1
    block = 1 
  []

  # Block 0: Gas
  # Block 1: Water
  # Block 2: Dielectric tubing
  # Block 4: electrode
  #[HeatDiff]
  #  type = ADHeatConduction
  #  variable = T
  #  block = 0
  #[]
  #[HeatTdot]
  #  type = ADHeatConductionTimeDerivative
  #  variable = T
  #  block = 0
  #[]
   # temperature
 [temperature_time]
   type = INSTemperatureTimeDerivative
   variable = T
   block = 0
 []
 [temperature_space]
   type = INSTemperature
   variable = T
   u = vel_x
   v = vel_y
   block = 0
 []

  [temperature_time_w]
   type = INSTemperatureTimeDerivative
   variable = T
   block = 1
  []
  [temperature_space_w]
   type = INSTemperature
   variable = T
   u = w_vel_x
   v = w_vel_y
   block = 1
  []
  

  #[HeatDiff_water]
  #  type = ADHeatConduction
  #  variable = T
  #  block = 1
  #[]
  #[HeatTdot_water]
  #  type = ADHeatConductionTimeDerivative
  #  variable = T
  #  block = 1
  #[]

  [HeatDiff_dielectric]
    type = ADHeatConduction
    variable = T
    block = 2
  []
  [HeatTdot_dielectric]
    type = ADHeatConductionTimeDerivative
    variable = T
    block = 2
  []

  [HeatDiff_cathode]
    type = ADHeatConduction
    variable = T
    block = 3
  []
  [HeatTdot_cathode]
    type = ADHeatConductionTimeDerivative
    variable = T
    block = 3
  []

[]

[AuxVariables]
  [heat_x]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
  [heat_y]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
[]

[AuxKernels]
  [heat_x_calc]
    type = HeatFlux
    variable = heat_x
    u = vel_x 
    v = vel_y
    temperature = T
    component = 0
    position_units = ${dom0Scale}
    block = 0
  []
  [heat_y_calc]
    type = HeatFlux
    variable = heat_y
    u = vel_x 
    v = vel_y
    temperature = T
    component = 1
    position_units = ${dom0Scale}
    block = 0
  []
[]

[BCs]
  [pressure_boundary]
    type = DirichletBC
    variable = p
    boundary = 'top_right'
    value = 0
  []
  [u_in]
    type = DirichletBC
    boundary = inlet
    variable = vel_x
    value = 0
  []
  [v_in]
    #type = FunctionDirichletBC
    #function = 'inlet_func'
    type = DirichletBC
    boundary = inlet
    variable = vel_y
    value = -0.088
    #value = -1
  []
  [u_axis_and_walls]
    type = DirichletBC
    #boundary = 'dielectric_left dielectric_right dielectric_tip electrode_side electrode_needle electrode_tip water_surface'
    boundary = 'dielectric_left dielectric_right dielectric_tip electrode_side electrode_needle electrode_tip'
    variable = vel_x
    value = 0
  []
  [v_no_slip]
    type = DirichletBC
    boundary = 'dielectric_left dielectric_right dielectric_tip electrode_side electrode_needle electrode_tip water_surface'
    variable = vel_y
    value = 0
  []

  [tip_temperature]
    type = DirichletBC
    variable = T
    value = 2500
    #value = 400
    boundary = 'electrode_tip'
  []
  [needle_top]
    type = DirichletBC
    variable = T
    value = 300
    boundary = 'needle_top'
  []

  #############################
  # WATER REGION
  #############################
  [temperature_water_bulk]
    type = DirichletBC
    variable = T 
    value = 300
    boundary = 'water_bulk'
  []

  [u_w_noslip]
    type = DirichletBC
    boundary = 'water_surface'
    variable = w_vel_y
    value = 0
  []

  [y_w_eq_y_bc] 
    type = MatchedValueBC
    boundary = 'water_surface'
    variable = w_vel_x
    v = vel_x
  []

  # pressure continuity
  [pressure_continuity]
    type = MatchedValueBC
    variable = w_p
    v = p
    boundary = 'water_surface'
  []
  
  ## Water no slip boundaries
  [w_vel_y_no_slip]
    type = DirichletBC
    #boundary = 'water_bulk ground'
    boundary = 'water_bulk'
    variable = w_vel_y
    value = 0
  []
  [w_vel_x_no_slip]
    type = DirichletBC
    #boundary = 'water_bulk ground'
    boundary = 'water_bulk'
    variable = w_vel_x
    value = 0
  []
[]

[Postprocessors]
  [flux_down]
    type = SideFluxIntegral
    variable = T
    boundary = 'gas_bottom'
    diffusivity = k
  []
[]

#[InterfaceKernels]
#  #[OH_diff]
#  #  type = InterfaceDiffusionLinear
#  #  variable = w_p
#  #  neighbor_var = p
#  #  h = 6.48e3
#  #  position_units = ${dom1Scale}
#  #  neighbor_position_units = ${dom0Scale}
#  #  boundary = 'water_surface'
#  #[]
#  [vel_x_continuity]
#    type = StressContinuity
#    variable = w_vel_x
#    neighbor_var = vel_x
#    u = vel_x
#    u_neighbor = w_vel_x
#    h = 6.48e3
#    position_units = ${dom1Scale}
#    neighbor_position_units = ${dom0Scale}
#    component = 0
#    boundary = 'water_top'
#  []
#  [vel_y_continuity]
#    type = StressContinuity
#    variable = w_vel_y
#    neighbor_var = vel_y
#    u = vel_x
#    u_neighbor = w_vel_x
#    h = 6.48e3
#    position_units = ${dom1Scale}
#    neighbor_position_units = ${dom0Scale}
#    component = 1
#    boundary = 'water_top'
#  []
#[]

#[ICs]
#  [./potential_ic]
#    type = FunctionIC
#    variable = potential
#    function = potential_ic_func
#    #block = 0
#  [../]
#[]

[Functions]
  [water_ic_func]
    type = ParsedFunction
    value = 'log(8.6949e23/6.022e23)'
    #value = '-2.15'
  []
  [./potential_bc_func]
    type = ParsedFunction
    value = 3
    #value = 10
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  [./em_aq_ic_func]
    type = ParsedFunction
    #value = 'log(exp(-22)*exp(-x*1e-5))'
    #value = 'log(exp(-16)*exp(-(x-1e-3)*7e5))'
    value = '1778 - 1.8e6*x'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-0.01 * (1.001e-3 - x)'
  [../]
  [./charged_gas_ic]
    type = ParsedFunction
    vars = 'sigma mu'
    vals = '25e-6 10e-5'
    value = 'max(log(1/(sigma*sqrt(2*3.14159)) * exp(-0.5*((x - mu)/sigma)^2) / (500000*sqrt(2/3.14159)) * 1e8/6.022e23), -48)'
  [../]
  [./mean_en_ic]
    type = ParsedFunction
    vars = 'sigma mu'
    vals = '25e-6 10e-5'
    value = '(log(1/(sigma*sqrt(2*3.14159)) * exp(-0.5*((x - mu)/sigma)^2)/16000 + 0.0258)) + max(log(1/(sigma*sqrt(2*3.14159)) * exp(-0.5*((x - mu)/sigma)^2) / (500000*sqrt(2/3.14159)) * 1e8/6.022e23), -48)'
  [../]
[]

[Materials]
  [test_0]
    type = ADGenericConstantMaterial
    block = 0
    prop_names = 'diffp' 
    prop_values = '1'
  []
  [test_1]
    type = ADGenericConstantMaterial
    block = 1
    prop_names = 'diffw_p' 
    prop_values = '1'
  []
  # INS stuff
  [const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu cp k'
    prop_values = '1.784 3.77e-5 312.2 0.016'
  []

  [const_w]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'rho mu cp k'
    prop_values = '977 8.9e-4 4182 0.6'
  []

  # Gas phase
  [density_argon]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.784' # kg m^-3
    block = 0
  []
  [thermal_conductivity_argon]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.016' # W m^-1 K^-1
    block = 0
  []
  [specific_heat_argon]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '312.2' # J kg^-1 K^-1
    block = 0
  []


  # Water phase
  [density_water]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '997' # kg m^-3
    block = 1
  []
  [thermal_conductivity_water]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.6' # W m^-1 K^-1
    block = 1
  []
  [specific_heat_water]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '4182' # J kg^-1 K^-1
    block = 1
  []

  # Dielectric
  [density_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '2500' # kg m^-3
    block = 2
  []
  [thermal_conductivity_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '1' # W m^-1 K^-1
    block = 2
  []
  [specific_heat_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '800' # J kg^-1 K^-1
    block = 2
  []

  # Cathode
  [density_tungsten]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '19300' # kg m^-3
    block = 3
  []
  [thermal_conductivity_cathode]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '173' # W m^-1 K^-1
    block = 3
  []
  [specific_heat_cathode]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '134' # J kg^-1 K^-1
    block = 3
  []

[]
