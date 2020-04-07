dom0Scale=1.0

[GlobalParams]
  offset = 20
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  [./file]
    type = FileMeshGenerator
    file = 'argon_csv_bc.msh'
  [../]
  [./dielectric_left]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master01_interface'
    input = file
  [../]
  [./plasma_left]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master10_interface'
    input = dielectric_left
  [../]
  [./plasma_right]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '2'
    new_boundary = 'master12_interface'
    input = plasma_left
  [../]
  [./dielectric_right]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '2'
    paired_block = '1'
    new_boundary = 'master21_interface'
    input = plasma_right
  [../]
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = dielectric_right
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
  compute_scaling_once = false
  end_time = 1
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = NEWTON
  line_search = 'basic'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
  dtmin = 1e-18
  dtmax = 1e-6
  l_max_its = 20
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    growth_factor = 1.2
   optimal_iterations = 30
  [../]
[]

[Outputs]
  perf_graph = true
  # print_linear_residuals = false
  [./out]
    type = Exodus
    #execute_on = 'final'
  [../]
  #[./dof_map]
  #  type = DOFMap
  #[../]
[]

[Debug]
  #show_var_residual_norms = true
[]

[UserObjects]
  [./data_provider]
    type = ProvideMobility
    electrode_area = 5.02e-7 # Formerly 3.14e-6
    ballast_resist = 1e6
    e = 1.6e-19
  [../]
[]

[DriftDiffusion]
  #[./Gas]
  #  species =      'Arp  Ar2p  Ar*'
  #  mass =         '1    1     1'
  #  charge =       '1    1     0'
  #  offset_vals =  '20   20    20'
  #  no_log_stabilization = 'Ar*'
  #  boundaries = 'left right'
  #  #offset = 20
  #  potential = potential
  #  block = 1
  #[../]
  [./Gas]
    species = 'Arp  Ar2p'
    mass =    '1    1'
    charge =  '1    1'
    offset_vals =  '20   20'
    #boundaries = 'left right'
    boundaries = 'master10_interface master12_interface'
    #offset = 20
    potential = potential
    block = 1
  [../]
  [./Excited]
    species = 'Ar*'
    mass = '1'
    charge = '0'
    #offset = 40
    #add_bcs = false
    #boundaries = 'left right'
    boundaries = 'master10_interface master12_interface'
    skip_advection = true
    potential = potential
    block = 1
  [../]
[]

[Kernels]
  [./em_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = em
    block = 1
  [../]
  [./em_advection]
    type = ADEFieldAdvection
    variable = em
    potential = potential
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = ADCoeffDiffusion
    variable = em
    mean_en = mean_en
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    #offset = 40
    block = 1
  [../]

  [./potential_diffusion_dom0]
    type = ADCoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom1]
    type = ADCoeffDiffusionLin
    variable = potential
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = ADCoeffDiffusionLin
    variable = potential
    block = 2
    position_units = ${dom0Scale}
  [../]
  #[./Arp_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential
  #  charged = Arp
  #  block = 1
  #[../]
  #[./Ar2p_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential
  #  charged = Ar2p
  #  block = 1
  #[../]
  [./em_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = em
    block = 1
  [../]

  # Heavy Species Log Stabilization
  #[./Arp_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Arp
  #  #offset = 40
  #  block = 1
  #[../]
  #[./Ar2p_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Ar2p
  #  #offset = 40
  #  block = 1
  #[../]
  #[./Arex_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Ar*
  #  block = 1
  #  offset = 40
  #  #offset = 25
  #[../]

  [./mean_en_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = mean_en
    block = 1
  [../]
  [./mean_en_advection]
    type = EFieldAdvectionEnergy
    variable = mean_en
    potential = potential
    em = em
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./mean_en_diffusion]
    type = ADCoeffDiffusion
    variable = mean_en
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./mean_en_joule_heating]
    #type = JouleHeating
    type = ADJouleHeating
    variable = mean_en
    potential = potential
    em = em
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./mean_en_log_stabilization]
    type = LogStabilizationMoles
    variable = mean_en
    block = 1
    offset = 15
    #offset = 40
  [../]
[]

[Variables]
  [./potential]
  [../]

  [./em]
    initial_condition = -21
    block = 1
  [../]

  [./Arp]
    initial_condition = -21
    block = 1
  [../]

  [./Ar*]
    initial_condition = -21
    block = 1
  [../]

  [./Ar2p]
    initial_condition = -21
    block = 1
  [../]

  [./mean_en]
    block = 1
    initial_condition = -20
  [../]
[]

[AuxVariables]
  [./Arex_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ar2p_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ar]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.70109
  [../]
  [./e_temp]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x_node]
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./em_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Arp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Efield]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Current_em]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Current_Arp]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Current_Ar2p]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./tot_gas_current]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./EFieldAdvAux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./DiffusiveFlux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./PowerDep_em]
   order = CONSTANT
   family = MONOMIAL
   block = 1
  [../]
  [./PowerDep_Arp]
   order = CONSTANT
   family = MONOMIAL
   block = 1
  [../]
  #[./ProcRate_el]
  # order = CONSTANT
  # family = MONOMIAL
  # block = 1
  #[../]
  #[./ProcRate_ex]
  # order = CONSTANT
  # family = MONOMIAL
  # block = 1
  #[../]
  #[./ProcRate_iz]
  # order = CONSTANT
  # family = MONOMIAL
  # block = 1
  #[../]
[]

[AuxKernels]
  [./Arex_lin]
    type = DensityMoles
    variable = Arex_lin
    density_log = Ar*
    block = 1
  [../]
  [./Ar2p_lin]
    type = DensityMoles
    variable = Ar2p_lin
    density_log = Ar2p
    block = 1
  [../]
  [./PowerDep_em]
    type = PowerDep
    density_log = em
    potential = potential
    art_diff = false
    potential_units = kV
    variable = PowerDep_em
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./PowerDep_Arp]
    type = PowerDep
    density_log = Arp
    potential = potential
    art_diff = false
    potential_units = kV
    variable = PowerDep_Arp
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./e_temp]
    type = ElectronTemperature
    variable = e_temp
    electron_density = em
    mean_en = mean_en
    block = 1
  [../]
  [./x_g]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./rho]
    type = ParsedAux
    variable = rho
    args = 'em_lin Arp_lin Ar2p_lin'
    function = 'Arp_lin + Ar2p_lin - em_lin'
    execute_on = 'timestep_end'
    block = 1
  [../]
  [./tot_gas_current]
    type = ParsedAux
    variable = tot_gas_current
    args = 'Current_em Current_Arp Current_Ar2p'
    function = 'Current_em + Current_Arp + Current_Ar2p'
    execute_on = 'timestep_end'
    block = 1
  [../]
  [./em_lin]
    type = DensityMoles
    variable = em_lin
    density_log = em
    block = 1
  [../]
  [./Arp_lin]
    type = DensityMoles
    variable = Arp_lin
    density_log = Arp
    block = 1
  [../]
  [./Efield_g]
    type = Efield
    component = 0
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./Current_em]
    type = Current
    potential = potential
    density_log = em
    variable = Current_em
    art_diff = false
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./Current_Arp]
    type = Current
    potential = potential
    density_log = Arp
    variable = Current_Arp
    art_diff = false
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./Current_Ar2p]
    type = Current
    potential = potential
    density_log = Ar2p
    variable = Current_Ar2p
    art_diff = false
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./EFieldAdvAux_em]
    type = EFieldAdvAux
    potential = potential
    density_log = em
    variable = EFieldAdvAux_em
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./DiffusiveFlux_em]
    type = DiffusiveFlux
    density_log = em
    variable = DiffusiveFlux_em
    block = 1
    position_units = ${dom0Scale}
  [../]
[]

[BCs]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    variable = mean_en
    #boundary = 'right'
    boundary = 'master12_interface'
    potential = potential
    em = em
    ip = 'Arp Ar2p'
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = ADHagelaarEnergyBC
    variable = mean_en
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential
    em = em
    ip = 'Arp Ar2p'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./secondary_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential
    em = em
    ip = 'Arp Ar2p'
    r = 0
    position_units = ${dom0Scale}
  [../]

  #[./potential_left]
  #  type = ADNeumannCircuitVoltageMoles_KV
  #  variable = potential
  #  boundary = left
  #  function = potential_bc_func
  #  ip = 'Arp Ar2p'
  #  data_provider = data_provider
  #  em = em
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./potential_left]
    type = ADCSVPotentialBC
    variable = potential
    file_name = 'voltage_data_01.txt'
    boundary = 'left'
  [../]
  [./potential_dirichlet_right]
    type = ADDirichletBC
    variable = potential
    boundary = right
    value = 0
  [../]
  [./em_physical_right]
    type = ADHagelaarElectronBC
    variable = em
    #boundary = 'right'
    boundary = 'master12_interface'
    potential = potential
    mean_en = mean_en
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]


  [./em_physical_left]
    type = ADHagelaarElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential
    ip = 'Arp Ar2p'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
[]

[ICs]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
  [../]
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    # value = '1.25*tanh(1e6*t)'
    value = 0.8
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-0.1 * (1.0001e-3 - x)'
  [../]
[]

[Materials]
  [./electron_moments]
    type = ADGasElectronMoments
    block = 1
    em = em
    mean_en = mean_en
    #property_tables_file = 'argon_chemistry_rates/electron_moments.txt'
    property_tables_file = 'argon_cm_test/electron_mobility_diffusion.txt'
  [../]

  [./dielectric_left_side]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'diffpotential'
    prop_values = '8.85e-11'
  [../]
  [./dielectric_right_side]
    type = GenericConstantMaterial
    block = 2
    prop_names = 'diffpotential'
    prop_values = '8.85e-11'
  [../]
  [./gas_constants]
    type = GenericConstantMaterial
    block = 1
    prop_names = ' e         N_A      diffpotential k_boltz eps  se_coeff se_energy T_gas massem   p_gas  n_gas'
    prop_values = '1.6e-19 6.022e23 8.85e-12      1.38e-23 8.854e-12 0.05     3.        300   9.11e-31 1.01e5 40.4915'
  [../]

  [./gas_species_0]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    #mobility = 0.144409938
    #diffusivity = 6.428571e-3
    block = 1
  [../]
  [./gas_species_1]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar2p
    heavy_species_mass = 1.328e-25
    heavy_species_charge = 1.0
    block = 1
  [../]

  [./gas_species_2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]

  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]
[]

#break libmesh_handleFPE
#  run ...
#  bt

[Reactions]
  # This argon reaction network based on a ZDPlasKin example:
  # zdplaskin.laplace.univ-tlse.fr
  # Example: "Micro-cathode sustained discharged in Ar"
  [./Argon]
    species = 'em Arp Ar2p Ar*'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    #file_location = 'argon_chemistry_rates'
    file_location = 'argon_cm_test'
    equation_constants = 'Tgas e_temp'
    equation_values = '300 34800'
    potential = 'potential'
    use_log = true
    position_units = ${dom0Scale}
    use_ad = true
    #convert_to_moles = true
    #convert_to_meters = 1e-2
    block = 1

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Effective_momentum)
                 em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.54_eV)
                 em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
                 em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ar*_De-excitation_11.54_eV)
                 em + Ar* -> em + em + Arp        : EEDF [-4.3] (C5_Ar*_Ionization_4.20_eV)
                 Ar2p + em -> Ar* + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
                 Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tgas * exp(-15130/Tgas)}
                 Ar* + Ar* -> Ar2p + em           : 3.6132e8
                 Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
                 Ar* + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
                 Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}'
    #reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Effective_momentum)
    #             em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.54_eV)
    #             em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
    #             em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ar*_De-excitation_11.54_eV)
    #             em + Ar* -> em + em + Arp        : EEDF [-4.3] (C5_Ar*_Ionization_4.20_eV)
    #             Ar2p + em -> Ar* + Ar            : {5.1187e17 * (e_temp/300)^(-0.67)}
    #             Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e18 / Tgas * exp(-15130/Tgas)}
    #             Ar* + Ar* -> Ar2p + em           : 3.6132e14
    #             Arp + em + em -> Ar + em         : {3.17314235e21 * (e_temp/11600)^(-4.5)}
    #             Ar* + Ar + Ar -> Ar + Ar + Ar    : 5.07702776e15
    #             Arp + Ar + Ar -> Ar2p + Ar       : {8.1595089e16 * (Tgas/300)^(-0.4)}'
  [../]
[]
