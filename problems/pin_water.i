dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  offset = 20
  potential_units = kV
  use_moles = true
[]

[Mesh]
  type = FileMesh
  #file = 'pin_water.msh'
  file = 'pin_water02.msh'
[]

[MeshModifiers]
  [./subdomain0]
    type = ParsedSubdomainMeshModifier
    combinatorial_geometry = 'x < 1e-3'
    #combinatorial_geometry = 'x <= 1'
    block_id = 0 
  [../]
  [./subdomain1]
    type = ParsedSubdomainMeshModifier
    combinatorial_geometry = 'x >= 1e-3'
    #combinatorial_geometry = 'x > 1'
    block_id = 1
    depends_on = subdomain0
  [../]

  [./interface]
    type = SideSetsBetweenSubdomains
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    depends_on = 'subdomain0 subdomain1' 
  [../]
  [./interface_again]
    type = SideSetsBetweenSubdomains
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
    depends_on = 'subdomain0 subdomain1'
  [../]

  [./left]
    type = SideSetsFromNormals
    normals = '-1 0 0'
    new_boundary = 'left'
  [../]
  [./right]
    type = SideSetsFromNormals
    normals = '1 0 0'
    new_boundary = 'right'
  [../]
[]

[Problem]
  type = FEProblem
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  #end_time = 5e-9
  #num_steps = 10
  end_time = 1e-1
  #line_search = 'none'
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = newton 
  petsc_options_iname = '-pc_type -snes_linesearch_minlambda -ksp_type'
  petsc_options_value = 'lu 1e-3 fgmres'
  #petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_minlambda -ksp_type'
  #petsc_options_value = 'bjacobi lu 1e-3 gmres'
  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -snes_linesearch_minlambda'
  #petsc_options_value = 'asm 2 ilu NONZERO 1e-10 1e-3'
  nl_rel_tol = 1e-4
  nl_abs_tol = 7.6e-5
  dtmin = 1e-14
  l_max_its = 20
  #dt = 1e-11
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    # dt = 1.1
    growth_factor = 1.2
   optimal_iterations = 15
  [../]
[]

[Outputs]
  perf_graph = true
  # print_linear_residuals = false
  #exodus = false
  exodus = true
  #interval = 50
[]

[Debug]
  show_var_residual_norms = true
[]

[UserObjects]
  [./data_provider]
    type = ProvideMobility
    electrode_area = 5.02e-7 # Formerly 3.14e-6
    ballast_resist = 1e6
    e = 1.6e-19
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[Kernels]
  [./em_time_deriv]
    type = ElectronTimeDerivative
    variable = em
    block = 0
  [../]
  [./em_advection]
    type = EFieldAdvectionElectrons
    variable = em
    potential = potential
    mean_en = mean_en
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = CoeffDiffusionElectrons
    variable = em
    mean_en = mean_en
    block = 0
    position_units = ${dom0Scale}
  [../]
  #[./em_ionization]
  #  type = ElectronsFromIonization
  #  variable = em
  #  potential = potential
  #  mean_en = mean_en
  #  em = em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    block = 0
  [../]
  # [./em_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = em
  #   potential = potential
  #   block = 0
  # [../]

  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]

  [./Arp_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = Arp
    block = 0
  [../]
  [./em_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = em
    block = 0
  [../]

  [./Arp_time_deriv]
    type = ElectronTimeDerivative
    variable = Arp
    block = 0
  [../]
  [./Arp_advection]
    type = EFieldAdvection
    variable = Arp
    potential = potential
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Arp_diffusion]
    type = CoeffDiffusion
    variable = Arp
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_log_stabilization]
    type = LogStabilizationMoles
    variable = Arp
    block = 0
  [../]

  [./mean_en_time_deriv]
    type = ElectronTimeDerivative
    variable = mean_en
    block = 0
  [../]
  [./mean_en_advection]
    type = EFieldAdvectionEnergy
    variable = mean_en
    potential = potential
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_diffusion]
    type = CoeffDiffusionEnergy
    variable = mean_en
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_joule_heating]
    type = JouleHeating
    variable = mean_en
    potential = potential
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_log_stabilization]
    type = LogStabilizationMoles
    variable = mean_en
    block = 0
    offset = 15
  [../]
[]

[Variables]
  [./potential]
    block = 0
    #scaling = 1e6
  [../]
  #[./potential_liq]
  #  block = 1
  #  #scaling = 1e6
  #[../]
  [./em]
    block = 0
    #scaling = 1e3
  [../]
  #[./emliq]
  #  block = 1
  #  #scaling = 1e5
  #[../]

  [./Arp]
    block = 0
    #scaling = 1e4
  [../]

  [./mean_en]
    block = 0
    #scaling = 1e3
  [../]

  #[./OHm]
  #  block = 1
  #  #initial_condition = -31
  #  #scaling = 2
  #[../]
[]

[AuxVariables]
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.701920755421197
    #initial_condition = 3.704261
  [../]

  #[./H2O]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_condition = 10.92252
  #  block = 1
  #[../]
  #[./O2]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_condition = -0.609203
  #  block = 1
  #[../]
  [./e_temp]
    block = 0
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./x]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./x_node]
  #[../]
  #[./rho]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./rholiq]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  [./em_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  #[./emliq_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  [./Arp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  #[./OHm_lin]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./Efield]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./Current_em]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./Current_emliq]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  #[./Current_Arp]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./Current_OHm]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./tot_gas_current]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./tot_liq_current]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./tot_flux_OHm]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./EFieldAdvAux_em]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./DiffusiveFlux_em]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./EFieldAdvAux_emliq]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  #[./DiffusiveFlux_emliq]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  #[./PowerDep_em]
  # order = CONSTANT
  # family = MONOMIAL
  # block = 0
  #[../]
  #[./PowerDep_Arp]
  # order = CONSTANT
  # family = MONOMIAL
  # block = 0
  #[../]
[]

[AuxKernels]
  [./e_temp]
    type = ElectronTemperature
    variable = e_temp
    electron_density = em
    mean_en = mean_en
    block = 0
  [../]
  #[./x_g]
  #  type = Position
  #  variable = x
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]
  #[./x_ng]
  #  type = Position
  #  variable = x_node
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]
  #[./rho]
  #  type = ParsedAux
  #  variable = rho
  #  args = 'em_lin Arp_lin'
  #  function = 'Arp_lin - em_lin'
  #  execute_on = 'timestep_end'
  #  block = 0
  #[../]
  #[./tot_gas_current]
  #  type = ParsedAux
  #  variable = tot_gas_current
  #  args = 'Current_em Current_Arp'
  #  function = 'Current_em + Current_Arp'
  #  execute_on = 'timestep_end'
  #  block = 0
  #[../]
  [./em_lin]
    type = DensityMoles
    convert_moles = true
    variable = em_lin
    density_log = em
    block = 0
  [../]
  #[./emliq_lin]
  #  type = DensityMoles
  #  convert_moles = true
  #  variable = emliq_lin
  #  density_log = emliq
  #  block = 1
  #[../]
  [./Arp_lin]
    type = DensityMoles
    convert_moles = true
    variable = Arp_lin
    density_log = Arp
    block = 0
  [../]
  #[./OHm_lin]
  #  type = DensityMoles
  #  convert_moles = true
  #  variable = OHm_lin
  #  density_log = OHm
  #  block = 1
  #[../]
  #[./Efield_g]
  #  type = Efield
  #  component = 0
  #  potential = potential
  #  variable = Efield
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]
  #[./Efield_l]
  #  type = Efield
  #  component = 0
  #  potential = potential_liq
  #  variable = Efield
  #  position_units = ${dom1Scale}
  #  block = 1
  #[../]
  #[./Current_em]
  #  type = Current
  #  potential = potential
  #  density_log = em
  #  variable = Current_em
  #  art_diff = false
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Current_emliq]
  #  type = Current
  #  potential = potential_liq
  #  density_log = emliq
  #  variable = Current_emliq
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom1Scale}
  #[../]
  #[./Current_Arp]
  #  type = Current
  #  potential = potential
  #  density_log = Arp
  #  variable = Current_Arp
  #  art_diff = false
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Current_OHm]
  #  block = 1
  #  type = Current
  #  potential = potential_liq
  #  density_log = OHm
  #  variable = Current_OHm
  #  art_diff = false
  #  position_units = ${dom1Scale}
  #[../]
  #[./tot_flux_OHm]
  #  block = 1
  #  type = TotalFlux
  #  potential = potential_liq
  #  density_log = OHm
  #  variable = tot_flux_OHm
  #[../]
  #[./EFieldAdvAux_em]
  #  type = EFieldAdvAux
  #  potential = potential
  #  density_log = em
  #  variable = EFieldAdvAux_em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./DiffusiveFlux_em]
  #  type = DiffusiveFlux
  #  density_log = em
  #  variable = DiffusiveFlux_em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./EFieldAdvAux_emliq]
  #  type = EFieldAdvAux
  #  potential = potential_liq
  #  density_log = emliq
  #  variable = EFieldAdvAux_emliq
  #  block = 1
  #  position_units = ${dom1Scale}
  #[../]
  #[./DiffusiveFlux_emliq]
  #  type = DiffusiveFlux
  #  density_log = emliq
  #  variable = DiffusiveFlux_emliq
  #  block = 1
  #  position_units = ${dom1Scale}
  #[../]
[]

[BCs]
  [./potential_right]
    type = DirichletBC
    variable = potential
    boundary = bottom
    value = 0
  [../]
  #[./potential_left]
  #  type = NeumannCircuitVoltageMoles_KV
  #  variable = potential
  #  boundary = left
  #  function = potential_bc_func
  #  ip = Arp
  #  data_provider = data_provider
  #  em = em
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./potential_electrode]
    type = DirichletBC
    variable = potential
    boundary = electrode
    value = -1.0
  [../]
  [./em_physical_bottom]
    type = HagelaarElectronBC
    variable = em
    boundary = bottom
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0.99
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_bottom_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = bottom 
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_bottom_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = bottom
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_bottom]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = bottom
    potential = potential
    em = em
    ip = Arp
    r = 0.99
    position_units = ${dom0Scale}
  [../]
  [./em_physical_electrode]
    type = HagelaarElectronBC
    variable = em
    boundary = electrode
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = SecondaryElectronBC
    variable = em
    boundary = electrode
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_inlet_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = electrode
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_inlet_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = electrode
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_inlet]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = electrode
    potential = potential
    em = em
    ip = Arp
    r = 0
    position_units = ${dom0Scale}
  [../]
[]

[ICs]
  [./em_ic]
    type = ConstantIC
    variable = em
    value = -21
    block = 0
  [../]
  [./Arp_ic]
    type = ConstantIC
    variable = Arp
    value = -21
    block = 0
  [../]
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -20
    block = 0
  [../]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
    block = 0
  [../]
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    #value = '1.25*tanh(1e6*t)'
    value = 1.25
    #value = '1.25*tanh(1e9*t)'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-1.25 * (1.0001e-3 - x)'
  [../]
[]

[Materials]
 [./GasBasics]
   type = GasBase
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   # user_p_gas = 1.01325e5
   user_p_gas = 1.01e5
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   property_tables_file = cpc_test/e_vals_test.txt
   #property_tables_file = electron_moments_test.txt
   # property_tables_file = dc_stankov/electron_moments.txt
   position_units = ${dom0Scale}
   block = 0
 [../]
 
  [./gas_species_0]
    type = HeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    block = 0
  [../]
  [./Ar_species]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  [../]
[]

[Reactions]
  [./gas_phase_reactions]
    species = 'em Arp'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    electron_energy = 'mean_en'
    electron_density = 'em'
    file_location = 'cpc_test'
    potential = 'potential'
    position_units = ${dom0Scale}
    use_log = true
    block = 0

    reactions = 'em + Ar -> em + Ar : EEDF [elastic]
                 em + Ar -> em + Ar*  : EEDF [-11.5]
                 em + Ar -> em + em + Arp  : EEDF [-1.576e1]'
  [../]
[]
