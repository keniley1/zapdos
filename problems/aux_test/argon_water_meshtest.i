dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  #offset = 20
  #offset = 24
  #offset = 48
  offset = 20
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  [./geo]
    type = FileMeshGenerator
    #file = 'gopalakrishnan_100um.msh'
    file = 'gopalakrishnan_1um.msh'
  [../]

  [./interface1]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    input = geo
  [../]
  [./interface2]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
    input = interface1
  [../]

  # The next two definitions create boundary conditions named
  # 'left' and 'right', where 'left' is at x = 0 and 'right' is at x = 1.1 mm.
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = interface2
  [../]
  [./right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  [../]
[]

[Problem]
  type = FEProblem
  # kernel_coverage_check = false
  #type = ReferenceResidualProblem
  #extra_tag_vectors = 'ref'
  #reference_vector = 'ref'
[]

[Preconditioning]
  active = smp

  [./smp]
    type = SMP
    full = true
    #ksp_norm = none
  [../]

  [./fsp]
    type = FSP
    topsplit = 'pc'
    #full = 'true'
    [./pc]
      splitting = 'p c'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type -pc_hypre_type'
      petsc_options_value = 'hypre boomeramg'
    [../]
    [./c]
      vars = 'em Arp mean_en emliq OH- O- O2- O3- HO2- H+ O O2_1 O3 H H2 HO2 HO3 OH H2O2'
      petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      petsc_options_value = 'asm 2 ilu NONZERO 1e-10'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  #end_time = 1e-1
  #end_time = 1e10
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  #resid_vs_jac_scaling_param = 1
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  #solve_type = pjfnk
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmin = 1e-17
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-14
    growth_factor = 1.4
    optimal_iterations = 10
  [../]
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  [out_test]
    type = Exodus
  []
[]

[Debug]
  #show_var_residual_norms = true
[]

[UserObjects]
  [./data_provider]
    type = ProvideMobility
    electrode_area = 5.02e-7 # Formerly 3.14e-6
    ballast_resist = 1e6
    #ballast_resist = 2.5e5
    e = 1.6e-19
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[Variables]
  [./potential]
  [../]
  [./em]
    block = 0
    initial_condition = -24
  [../]
  [./emliq]
    block = 1
    initial_condition = -16
  [../]

  [./Arp]
    block = 0
    initial_condition = -24
  [../]

  [./mean_en]
    block = 0
    initial_condition = -24
    # scaling = 1e-1
  [../]
[]

[Kernels]
  # Potential
  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential
    block = 1
    position_units = ${dom1Scale}
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
  [./emliq_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = emliq
    block = 1
  [../]
  
  # Gas Phase
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
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
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

  # Liquid Phase
  [./emliq_time_deriv]
    type = ElectronTimeDerivative
    variable = emliq
    block = 1
  [../]
  [./emliq_advection]
    type = EFieldAdvection
    variable = emliq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./emliq_diffusion]
    type = CoeffDiffusion
    variable = emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./emliq_log_stabilization]
    type = LogStabilizationMoles
    variable = emliq
    block = 1
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

[AuxVariables]
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.701920755421197
    #initial_condition = 3.704261
  [../]
  [./em_density]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Arp_density]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]

  [./H2O]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 10.92252
    block = 1
  [../]
  [./H2O_density]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./emliq_density]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]

  [./e_temp]
    block = 0
    order  = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./x0]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./x1]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
[]

[AuxKernels]
  [./em_density]
    type = DensityMoles
    variable = em_density
    density_log = em
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./Arp_density]
    type = DensityMoles
    variable = Arp_density
    density_log = Arp
    execute_on = 'initial timestep_end'
    block = 0
  [../]

  [./H2O_density]
    type = DensityMoles
    variable = H2O_density
    density_log = H2O
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./emliq_density]
    type = DensityMoles
    variable = emliq_density
    density_log = emliq
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./e_temp]
    type = ElectronTemperature
    variable = e_temp
    electron_density = em
    execute_on = 'initial timestep_end'
    mean_en = mean_en
    block = 0
  [../]
  [./x_g]
    type = Position
    variable = x0
    position_units = ${dom0Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_l]
    type = Position
    variable = x1
    position_units = ${dom1Scale}
    execute_on = 'initial timestep_end'
    block = 1
  [../]
[]

[InterfaceKernels]
  [./em_advection]
    type = InterfaceAdvection
    mean_en_neighbor = mean_en
    potential_neighbor = potential
    neighbor_var = em
    variable = emliq
    boundary = master1_interface
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = InterfaceLogDiffusionElectrons
    mean_en_neighbor = mean_en
    neighbor_var = em
    variable = emliq
    boundary = master1_interface
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
[]

[BCs]
  [./potential_left]
    type = NeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    ip = Arp
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_dirichlet_right]
    type = DirichletBC
    variable = potential
    boundary = right
    value = 0
  [../]
  [./em_physical_right]
    type = HagelaarElectronBC
    variable = em
    boundary = 'master0_interface'
    potential = potential
    #ip = Arp
    mean_en = mean_en
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = 'master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = 'master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_right]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master0_interface'
    potential = potential
    em = em
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./em_physical_left]
    type = HagelaarElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = SecondaryElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = SecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    ip = Arp
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = 'left'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./emliq_right]
    type = DCIonBC
    variable = emliq
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
[]

[ICs]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
    #block = 0
  [../]
  #[./em_ic]
  #  type = FunctionIC
  #  variable = em
  #  function = charged_gas_ic
  #[../]
  #[./Arp_ic]
  #  type = FunctionIC
  #  variable = Arp
  #  function = charged_gas_ic
  #[../]
  #[./mean_en_ic]
  #  type = FunctionIC
  #  variable = mean_en
  #  #function = charged_gas_ic
  #  function = mean_en_ic
  #[../]
  #[./emliq_ic]
  #  type = FunctionIC
  #  variable = emliq
    #function = emliq_ic_func
  #[../]
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    #value = 1.25
    #value = 5
    #value = 4
    value = 3
    #value = 2
    #value = 1
    #value = 0.8
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  [./emliq_ic_func]
    type = ParsedFunction
    #value = 'log(exp(-22)*exp(-x*1e-5))'
    #value = 'log(exp(-16)*exp(-(x-1e-3)*7e5))'
    value = '1778 - 1.8e6*x'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-1.25 * (1.0001e-3 - x)'
    #value = '-5 * (2e-3 - x)'
    #value = '-4 * (2e-3 - x)'
    #value = '-3 * (2e-3 - x)'
    #value = '-2 * (2e-3 - x)'
    #value = '-1 * (2e-3 - x)'
    #value = '-0.8 * (2e-3 - x)'
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
  [./se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'left'
  [../]
 [./GasBasics]
   type = GasElectronMoments
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   user_p_gas = 1.01e5
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   property_tables_file = '/home/shane/projects/zapdos/problems/argon_water_prelim_files/electron_mobility_diffusion.txt'
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
  [./gas_species_1]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar2p
    heavy_species_mass = 13.28e-26
    heavy_species_charge = 1.0
    mobility = 0.0
    block = 0
  [../]

  [./electron_data]
    type = GenericConstantMaterial
    prop_names = 'diffemliq muemliq sgnemliq'
    prop_values = '4.5e-9 0.000173913 -1'
    block = 1
  [../]

  [./N_A_mat]
    type = GenericConstantMaterial
    prop_names = 'N_A e diffpotential diffpotential_liq T_gas p_gas'
    prop_values = '6.022e23 1.602e-19 7.0832e-10 7.0832e-10 300 1.013e5'
    block = 1
  [../]
[]
