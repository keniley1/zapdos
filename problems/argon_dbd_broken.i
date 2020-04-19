dom0Scale=1
dom1Scale=1
dom2Scale=1

[GlobalParams]
  #offset = 20
  offset = 40
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  #final_generator = dielectric_left_mesh
  [./file]
    type = FileMeshGenerator
    #file = 'argon_csv_bc.msh'
    file = 'argon_dbd_mesh.msh'
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

  #[./dielectric_left_mesh]
  #  type = MeshSideSetGenerator
  #  input = right
  #  boundaries = 'master10_interface'
  #  block_id = 101
  #  block_name = dielectric_left_surface
  #  depends_on = 'plasma_left'
  #[../]
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
  #compute_scaling_once = false
  #end_time = 5e-3
  end_time = 8e-5
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  #solve_type = NEWTON
  solve_type = PJFNK
  line_search = 'basic'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol -ksp_gmres_restart'
  petsc_options_value = 'lu NONZERO 1.e-10 0 40'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  dtmin = 1e-18
  dtmax = 1e-7
  l_max_its = 40
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    growth_factor = 1.2
   optimal_iterations = 30
  [../]
  [./TimeIntegrator]
    type = BDF2
  [../]
[]

[Outputs]
  perf_graph = true
  [./out]
    type = Exodus
    #execute_on = 'final'
  [../]
  #[./sigma_output]
  #  type = CSV
  #  show = 'sigma_left'
  #[../]
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

[Kernels]
  [./em_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = em
    block = 1
  [../]
  [./em_advection]
    type = ADEFieldAdvection
    #type = EFieldAdvectionElectrons
    variable = em
    mean_en = mean_en
    potential = potential
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = ADCoeffDiffusion
    #type = CoeffDiffusionElectrons
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
    type = CoeffDiffusionLin
    variable = potential
    #variable = potential_dom0
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential
    #variable = potential_dom2
    block = 2
    position_units = ${dom0Scale}
  [../]
  [./Arp_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = Arp
    block = 1
  [../]
  [./em_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = em
    block = 1
  [../]

  [./Arp_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = Arp
    block = 1
  [../]
  [./Arp_advection]
    type = ADEFieldAdvection
    #type = EFieldAdvection
    variable = Arp
    potential = potential
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion]
    #type = CoeffDiffusion
    type = ADCoeffDiffusion
    variable = Arp
    position_units = ${dom0Scale}
  [../]
  [./Arp_log_stabilization]
    type = LogStabilizationMoles
    variable = Arp
    block = 1
  [../]

  [./Arex_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = Ar*
    block = 1
  [../]
  [./Arex_diffusion]
    #type = CoeffDiffusion
    type = ADCoeffDiffusion
    variable = Ar*
    position_units = ${dom0Scale}
  [../]
  [./Arex_log_stabilization]
    type = LogStabilizationMoles
    variable = Ar*
    block = 1
  [../]

  [./Ar_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = Ar
    block = 1
  [../]
  #[./Ar_diffusion]
  #  #type = CoeffDiffusion
  #  type = ADCoeffDiffusion
  #  variable = Ar
  #  position_units = ${dom0Scale}
  #[../]
  [./Ar_log_stabilization]
    type = LogStabilizationMoles
    variable = Ar
    block = 1
  [../]

  [./mean_en_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = mean_en
    block = 1
  [../]
  [./mean_en_advection]
    #type = EFieldAdvectionEnergy
    type = ADEFieldAdvection
    variable = mean_en
    potential = potential
    em = em
    block = 1
    position_units = ${dom0Scale}
  [../]
  [./mean_en_diffusion]
    type = ADCoeffDiffusion
    #type = CoeffDiffusion
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
    #offset = 15
    #offset = 40
  [../]
[]

[Variables]
  [./Ar]
    initial_condition = 3.70109
    block = 1
  [../]
  [./potential]
  [../]

  [./em]
    initial_condition = -24
    block = 1
  [../]

  [./Arp]
    initial_condition = -24
    block = 1
  [../]

  [./Ar*]
    initial_condition = -21
    block = 1
  [../]

  [./mean_en]
    block = 1
    initial_condition = -24
  [../]
[]

[AuxVariables]
  [./Ar_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Arex_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./Ar2p_lin]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./Ar]
  #  block = 1
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_condition = 3.70109
  #[../]
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
  #[./Current_Ar2p]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
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
  [./Ar_lin]
    type = DensityMoles
    variable = Ar_lin
    density_log = Ar
    block = 1
  [../]
  [./Arex_lin]
    type = DensityMoles
    variable = Arex_lin
    density_log = Ar*
    block = 1
  [../]
  #[./Ar2p_lin]
  #  type = DensityMoles
  #  variable = Ar2p_lin
  #  density_log = Ar2p
  #  block = 1
  #[../]
  #[./PowerDep_em]
  #  type = PowerDep
  #  density_log = em
  #  potential = potential
  #  art_diff = false
  #  potential_units = kV
  #  variable = PowerDep_em
  #  position_units = ${dom0Scale}
  #  block = 1
  #[../]
  #[./PowerDep_Arp]
  #  type = PowerDep
  #  density_log = Arp
  #  potential = potential
  #  art_diff = false
  #  potential_units = kV
  #  variable = PowerDep_Arp
  #  position_units = ${dom0Scale}
  #  block = 1
  #[../]
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
    #args = 'em_lin Arp_lin Ar2p_lin'
    #function = 'Arp_lin + Ar2p_lin - em_lin'
    args = 'em_lin Arp_lin'
    function = 'Arp_lin - em_lin'
    execute_on = 'timestep_end'
    block = 1
  [../]
  #[./tot_gas_current]
  #  type = ParsedAux
  #  variable = tot_gas_current
  #  #args = 'Current_em Current_Arp Current_Ar2p'
  #  #function = 'Current_em + Current_Arp + Current_Ar2p'
  #  args = 'Current_em Current_Arp'
  #  function = 'Current_em + Current_Arp'
  #  execute_on = 'timestep_end'
  #  block = 1
  #[../]
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
  [./Efield_d0]
    type = Efield
    component = 0
    #potential = potential_dom0
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Efield_d1]
    type = Efield
    component = 0
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./Efield_d2]
    type = Efield
    component = 0
    #potential = potential_dom2
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 2
  [../]
  #[./Current_em]
  #  type = Current
  #  potential = potential
  #  density_log = em
  #  variable = Current_em
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./Current_Arp]
  #  type = Current
  #  potential = potential
  #  density_log = Arp
  #  variable = Current_Arp
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./Current_Ar2p]
  #  type = Current
  #  potential = potential
  #  density_log = Ar2p
  #  variable = Current_Ar2p
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./EFieldAdvAux_em]
  #  type = EFieldAdvAux
  #  potential = potential
  #  density_log = em
  #  variable = EFieldAdvAux_em
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./DiffusiveFlux_em]
  #  type = DiffusiveFlux
  #  density_log = em
  #  variable = DiffusiveFlux_em
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
[]

[BCs]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    #type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master12_interface'
    potential = potential
    em = em
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = ADHagelaarEnergyBC
    #type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master10_interface'
    potential = potential
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./secondary_energy_left]
  #  type = ADSecondaryElectronEnergyBC
  #  variable = mean_en
  #  #boundary = 'left'
  #  boundary = 'master10_interface'
  #  potential = potential
  #  em = em
  #  #ip = 'Arp Ar2p'
  #  ip = 'Arp'
  #  se_coeff = 0.01
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./secondary_energy_right]
  #  type = ADSecondaryElectronEnergyBC
  #  variable = mean_en
  #  #boundary = 'left'
  #  boundary = 'master12_interface'
  #  potential = potential
  #  em = em
  #  ip = 'Arp'
  #  se_coeff = 0.01
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]

  [./em_physical_right]
    type = ADHagelaarElectronBC
    #type = HagelaarElectronBC
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
    #type = HagelaarElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./sec_electrons_left]
  #  type = ADSecondaryElectronBC
  #  variable = em
  #  #boundary = 'left'
  #  boundary = 'master10_interface'
  #  potential = potential
  #  #ip = 'Arp Ar2p'
  #  ip = 'Arp'
  #  se_coeff = 0.01
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./sec_electrons_right]
  #  type = ADSecondaryElectronBC
  #  variable = em
  #  #boundary = 'left'
  #  boundary = 'master12_interface'
  #  potential = potential
  #  #ip = 'Arp Ar2p'
  #  ip = 'Arp'
  #  se_coeff = 0.01
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]

  [./potential_left]
    type = FunctionDirichletBC
    #variable = potential
    #variable = potential_dom0
    variable = potential
    function = potential_input
    boundary = 'left'
  [../]

  [./potential_dirichlet_right]
    type = DirichletBC
    #variable = potential_dom2
    variable = potential
    boundary = right
    value = 0
  [../]

  [./Arp_advection_left]
    type = ADHagelaarIonAdvectionBC
    #type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = master10_interface
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_advection_right]
    type = ADHagelaarIonAdvectionBC
    #type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = master12_interface
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion_left]
    type = ADHagelaarIonDiffusionBC
    #type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = master10_interface
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion_right]
    type = ADHagelaarIonDiffusionBC
    #type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = master12_interface
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./Arex_diffusion_left]
    type = ADHagelaarIonDiffusionBC
    #type = HagelaarIonDiffusionBC
    variable = Ar*
    boundary = master10_interface
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex_diffusion_right]
    type = ADHagelaarIonDiffusionBC
    #type = HagelaarIonDiffusionBC
    variable = Ar*
    boundary = master12_interface 
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
  #[./potential_input]
  #  type = PiecewiseLinear
  #  #data_file = 'voltage_data_01.txt'
  #  data_file = 'voltage_data_temp01.txt'
  #  format = columns
  #[../]
  [./potential_input]
    type = ParsedFunction
    vars = 'f0'
    vals = '50e3'
    value = '-0.75*sin(2*3.1415926*f0*t)'
  [../]
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
    property_tables_file = 'argon_cm_test/electron_mobility_diffusion.txt'
  [../]

  [./dielectric_left_side]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'diffpotential'
    #prop_names = 'diffpotential_dom0'
    prop_values = '8.85e-11'
  [../]
  [./dielectric_right_side]
    type = GenericConstantMaterial
    block = 2
    prop_names = 'diffpotential'
    #prop_names = 'diffpotential_dom2'
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
    #type = HeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    block = 1
  [../]

  [./gas_species_2]
    type = ADHeavySpeciesMaterial
    #type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]

  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    #type = HeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]
[]

[Reactions]
  # This argon reaction network based on a ZDPlasKin example:
  # zdplaskin.laplace.univ-tlse.fr
  # Example: "Micro-cathode sustained discharged in Ar"
  [./Argon]
    #species = 'em Arp Ar2p Ar*'
    species = 'em Arp Ar* Ar'
    #aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    #file_location = 'argon_chemistry_rates'
    #file_location = 'argon_cm_test/townsend/'
    file_location = 'argon_cm_test/'
    equation_constants = 'Tgas e_temp'
    equation_values = '300 34800'
    potential = 'potential'
    use_log = true
    position_units = ${dom1Scale}
    use_ad = true
    #convert_to_moles = true
    #convert_to_meters = 1e-2
    block = 1

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
                 em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
                 em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ars_De-excitation_11.50_eV)
                 em + Ar* -> em + em + Arp        : EEDF [-4.43] (C5_Ars_Ionization_4.43_eV)
                 Ar* + Ar -> Ar + Ar              : 1807
                 Ar* + Ar* -> em + Ar + Arp       : 3.3734e8'
  [../]
[]
