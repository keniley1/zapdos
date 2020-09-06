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
    dt = 1e-9
    growth_factor = 1.4
    optimal_iterations = 10
  [../]
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  #[./m5_kV_1um_02]
  #[./m4_kV_1um_02]
  #[./m3_kV_1um_02]
  #[./m2_kV_1um_02]
  #[./m1_kV_1um_02]
  [m1_kV_1um_03]
    type = Exodus
  [../]
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

[DriftDiffusionAction]
  [./Plasma]
    electrons = em
    charged_particle = 'Arp'
    mean_energy = mean_en
    potential = potential
    is_potential_unique = false
    using_offset = true
    offset = 20
    position_units = ${dom0Scale}
    block = 0
  [../]

  [./Water]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    #charged_particle = 'emliq OH- H2O+ O- H3O+ HO2- O2- O3- H+'
    charged_particle = 'emliq OH- H2O+ O- H3O+ HO2- O2- O3- H+ Na+ Cl-'
    Neutrals = 'H OH H2 H2O2 O2 O HO2 O3'
    is_potential_unique = false
    potential = potential
    using_offset = true
    offset = 30
    position_units = ${dom1Scale}
    block = 1
  [../]
[]

[Variables]
  [./Na+]
    block = 1
    initial_condition = 4.60517
  [../]
  [./Cl-]
    block = 1
    initial_condition = 4.60517
  [../]
  [./potential]
    #block = 0
  [../]
  #[./potential_liq]
  #  block = 1
  #[../]
  [./em]
    block = 0
    initial_condition = -24
  [../]
  [./emliq]
    block = 1
    #initial_condition = -24
    #initial_condition = -21
    # scaling = 1e-5
    initial_condition = -22
    #initial_condition = -24
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

  [./O2]
    block = 1
    #initial_condition = -0.609203
    initial_condition = -24
  [../]

  [./OH-]
    block = 1
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    initial_condition = -9.210340
  [../]

  [./O-]
    block = 1
    initial_condition = -31
  [../]

  [./O2-]
    block = 1
    initial_condition = -26
  [../]

  [./O3-]
    block = 1
    initial_condition = -31
  [../]

  [./HO2-]
    block = 1
    initial_condition = -31
  [../]

  [./H+]
    block = 1
    initial_condition = -21
    #initial_condition = -21
    #initial_condition = -21.872637287474074
    #initial_condition = -9.210340
  [../]

  [./H3O+]
    block = 1
    initial_condition = -9.210340
  [../]

  [./H2O+]
    block = 1
    initial_condition = -31
  [../]

  [./O]
    block = 1
    initial_condition = -31
  [../]

  [./O3]
    block = 1
    initial_condition = -31
  [../]

  [./H]
    block = 1
    initial_condition = -25
  [../]

  [./H2]
    block = 1
    initial_condition = -31
  [../]

  [./HO2]
    block = 1
    initial_condition = -31
  [../]

  [./OH]
    block = 1
    initial_condition = -25
  [../]

  [./H2O2]
    block = 1
    initial_condition = -31
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

  [./H2O]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 10.92252
    block = 1
  [../]
  #[./O2]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_condition = -0.609203
  #  block = 1
  #[../]
  [./H2O_density]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]

  [./e_temp]
    block = 0
    order  = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./x]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./x_node]
    initial_condition = 0
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./rholiq]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./Efield]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./Current_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./Current_Arp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_OH-]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./tot_gas_current]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./tot_liq_current]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./tot_flux_OH-]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./EFieldAdvAux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./DiffusiveFlux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./EFieldAdvAux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./DiffusiveFlux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./PowerDep_em]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
  [./PowerDep_Arp]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
  [./ProcRate_el]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
  [./ProcRate_ex]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
  [./ProcRate_iz]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
[]

[AuxKernels]
  [./H2O_density]
    type = DensityMoles
    variable = H2O_density
    density_log = H2O
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./PowerDep_em]
    type = PowerDep
    density_log = em
    potential = potential
    art_diff = false
    execute_on = 'initial timestep_end'
    potential_units = kV
    variable = PowerDep_em
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./PowerDep_Arp]
    type = PowerDep
    density_log = Arp
    potential = potential
    art_diff = false
    potential_units = kV
    variable = PowerDep_Arp
    execute_on = 'initial timestep_end'
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ProcRate_el]
    type = ProcRate
    em = em
    potential = potential
    proc = el
    variable = ProcRate_el
    execute_on = 'initial timestep_end'
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ProcRate_ex]
    type = ProcRate
    em = em
    potential = potential
    execute_on = 'initial timestep_end'
    proc = ex
    variable = ProcRate_ex
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ProcRate_iz]
    type = ProcRate
    em = em
    potential = potential
    execute_on = 'initial timestep_end'
    proc = iz
    variable = ProcRate_iz
    position_units = ${dom0Scale}
    block = 0
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
    variable = x
    position_units = ${dom0Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_l]
    type = Position
    variable = x
    position_units = ${dom1Scale}
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_nl]
    type = Position
    variable = x_node
    position_units = ${dom1Scale}
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./rho]
    type = ParsedAux
    variable = rho
    args = 'em_density Arp_density'
    function = 'Arp_density - em_density'
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./Efield_g]
    type = Efield
    component = 0
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Efield_l]
    type = Efield
    component = 0
    #potential = potential_liq
    potential = potential
    variable = Efield
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./Current_em]
    type = Current
    potential = potential
    density_log = em
    variable = Current_em
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Current_emliq]
    type = Current
    #potential = potential_liq
    potential = potential
    density_log = emliq
    variable = Current_emliq
    art_diff = false
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Current_Arp]
    type = Current
    potential = potential
    density_log = Arp
    variable = Current_Arp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Current_OH-]
    block = 1
    type = Current
    #potential = potential_liq
    potential = potential
    density_log = OH-
    variable = Current_OH-
    art_diff = false
    position_units = ${dom1Scale}
  [../]
  [./tot_flux_OH-]
    block = 1
    type = TotalFlux
    #potential = potential_liq
    potential = potential
    density_log = OH-
    variable = tot_flux_OH-
  [../]
  [./EFieldAdvAux_em]
    type = EFieldAdvAux
    potential = potential
    density_log = em
    variable = EFieldAdvAux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./DiffusiveFlux_em]
    type = DiffusiveFlux
    density_log = em
    variable = DiffusiveFlux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./EFieldAdvAux_emliq]
    type = EFieldAdvAux
    #potential = potential_liq
    potential = potential
    density_log = emliq
    variable = EFieldAdvAux_emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./DiffusiveFlux_emliq]
    type = DiffusiveFlux
    density_log = emliq
    variable = DiffusiveFlux_emliq
    block = 1
    position_units = ${dom1Scale}
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

  #[./water_interface]
  #  type = InterfaceFluxConservation
  #  neighbor_var = potential
  #  variable = potential_liq
  #  boundary = 'master1_interface'
  #  region_name = 'diffpotential_liq'
  #  neighbor_region_name = 'diffpotential'
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #[../]
[]

[BCs]
  #[./emliq_physical_left]
  #  type = HagelaarElectronBC
  #  variable = emliq
  #  boundary = 'master1_interface'
  #  potential = potential
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Nap_right]
  #  type = DCIonBC
  #  variable = Na+
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./Clm_right]
  #  type = DCIonBC
  #  variable = Cl-
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  [Nap_right]
    type = DirichletBC
    variable = Na+
    value = 4.60517
    boundary = right
  []
  [Clm_right]
    type = DirichletBC
    variable = Cl-
    value = 4.60517
    boundary = right
  []
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
  #[./potential_left]
  #  type = DirichletBC
  #  variable = potential
  #  value = -2.5
  #  boundary = left
  #[../]
  #[./potential_left]
  #  type = FunctionDirichletBC
  #  variable = potential
  #  boundary = left
  #  function = test_bc
  #[../]
  [./potential_dirichlet_right]
    type = DirichletBC
    variable = potential
    boundary = right
    value = 0
  [../]
  #[./em_physical_right]
  #  type = HagelaarElectronBC
  #  variable = em
  #  boundary = 'master0_interface'
  #  potential = potential
  #  #ip = Arp
  #  mean_en = mean_en
  #  #r = 0.99
  #  r = 0.0
  #  position_units = ${dom0Scale}
  #[../]
  [em_physical_right]
    type = ADDriftDiffusionOpenBC
    variable = em
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'master0_interface'
  []
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
  #[./emliq_right]
  #  type = DirichletBC
  #  variable = emliq
  #  boundary = 'right'
  #  value = -20
  #[../]
  [./OH-_physical]
    type = DCIonBC
    variable = OH-
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]

  [./O-_physical]
    type = DCIonBC
    variable = O-
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./O-_right]
  #  type = DirichletBC
  #  variable = O- 
  #  boundary = 'right'
  #  value = -20
  #[../]
  [./O2-_physical]
    type = DCIonBC
    variable = O2-
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./O2-_right]
  #  type = DirichletBC
  #  variable = O2- 
  #  boundary = 'right'
  #  value = -20
  #[../]
  [./O3-_physical]
    type = DCIonBC
    variable = O3-
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./O3-_right]
  #  type = DirichletBC
  #  variable = O3- 
  #  boundary = 'right'
  #  value = -20
  #[../]
  [./HO2-_physical]
    type = DCIonBC
    variable = HO2-
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./HO2-_right]
  #  type = DirichletBC
  #  variable = HO2- 
  #  boundary = 'right'
  #  value = -20
  #[../]
  [./Hp_physical]
    type = DCIonBC
    variable = H+
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./Hp_right]
  #  type = DirichletBC
  #  variable = H+ 
  #  boundary = 'right'
  #  value = -20
  #[../]
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
#  [./gas_block]
#    type = Gas
#    interp_trans_coeffs = true
#    interp_elastic_coeff = true
#    ramp_trans_coeffs = false
#    em = em
#    potential = potential
#    ip = Arp
#    mean_en = mean_en
#    user_se_coeff = .05
#    block = 0
#    property_tables_file = td_argon_mean_en.txt
#  [../]
 # [./water_block]
 #   type = Water
 #   block = 1
 #   potential = potential
 # [../]
 [./GasBasics]
   #type = GasBase
   type = GasElectronMoments
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   # user_p_gas = 1.01325e5
   user_p_gas = 1.01e5
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   #property_tables_file = cpc_test/e_vals_test.txt
   property_tables_file = '/home/shane/projects/zapdos/problems/argon_water_prelim_files/electron_mobility_diffusion.txt'
   block = 0
 [../]

  [./Nap_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = Na+
    heavy_species_mass = 3.816e-26
    heavy_species_charge = 1.0
    diffusivity = 2e-9
    block = 1
  [../]
  [./Clm_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = Cl-
    heavy_species_mass = 5.887e-26
    heavy_species_charge = -1.0
    diffusivity = 2e-9
    block = 1
  [../]
  [./NO2-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = NO2-
    heavy_species_mass = 7.6e-26 
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]
  [./NO2_2m_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = NO2_2- 
    heavy_species_mass = 7.6e-26 
    heavy_species_charge = -2
    diffusivity = 5e-9
    block = 1
  [../]
  [./NO3-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = NO3-
    heavy_species_mass = 1e-25
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]
  [./NO3_2m_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = NO3_2-
    heavy_species_mass = 1e-25 
    heavy_species_charge = -2
    diffusivity = 5e-9 
    block = 1
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

 # [./water_block]
 #   type = Water
 #   block = 1
 #   potential = potential
 # [../]


  [./OH-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = OH-
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]

  [./O_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./O3_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O3
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = 0
    diffusivity = 5e-9
    block = 1
  [../]

  [./OH_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = OH
    heavy_species_mass = 2.82431e-26
    heavy_species_charge = 0
    diffusivity = 4.5e-9
    block = 1
  [../]

  [./HO2_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = HO2
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./H2O2_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H2O2
    heavy_species_mass = 5.64840e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./H2Op_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H2O+
    heavy_species_mass = 2.992e-26 
    heavy_species_charge = 1
  [../]

  [./H2_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H2
    heavy_species_mass = 3.34752e-26
    heavy_species_charge = 0
    diffusivity = 4.50e-9
    block = 1
  [../]

  [./H_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 0
    diffusivity = 5.0e-9
    block = 1
  [../]

  [./H+_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H+
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 1
  [../]

  [./H3O+_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = H3O+
    heavy_species_mass = 3e-26
    # ^ just estimated...whatever
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 1
  [../]

  [./HO2-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = HO2-
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./O-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O-
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./O2-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O2-
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = -1
    diffusivity = 5e-9 
    block = 1
  [../]

  [./O3-_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O3-
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./O2_mat]
    type = HeavySpeciesMaterial
    heavy_species_name = O2
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 2e-9
    block = 1
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

[Reactions]
  [./gas_phase_reactions]
    species = 'em Arp'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    #reaction_coefficient_format = 'rate'
    electron_energy = 'mean_en'
    electron_density = 'em'
    #file_location = 'cpc_test'
    file_location = '/home/shane/projects/zapdos/problems/argon_water_prelim_files/townsend'
    potential = 'potential'
    position_units = ${dom0Scale}
    use_log = true
    track_rates = false
    block = 0

    reactions = 'em + Ar -> em + Ar : EEDF [elastic] (C1_Ar_Effective_momentum)
                 em + Ar -> em + Ar*  : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
                 em + Ar -> em + em + Arp  : EEDF [-1.576e1] (C3_Ar_Ionization_15.80_eV)'
  [../]


  [./liquid_phase_reactions]
    # Missing NO2-, NO2_2-, NO3-, NO3_2-
    species = 'emliq OH- H2O+ O- H3O+ HO2- O2- O3- H OH H2 H2O2 O2 O HO2 O3 H+'
    aux_species = 'H2O'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'emliq + H2O -> H + OH-               : 1.9e-2
                 emliq + H2O+ -> H + OH               : 6e8
                 #emliq + emliq -> H2 + OH- + OH-      : 5.5e9
                 emliq + emliq -> H2 + OH- + OH-      : 3.0703e8
                 emliq + H + H2O -> H2 + OH-          : 2.5e4
                 emliq + OH -> OH-                    : 3e7
                 emliq + O- + H2O -> OH + OH          : 2.2e4
                 emliq + H3O+ -> H + H2O              : 2.3e7
                 emliq + H2O2 -> OH + OH-             : 1.1e7
                 emliq + HO2- + H2O -> OH + OH- + OH- : 3.5e3
                 emliq + O2 -> O2-                    : 1.9e7
                 emliq + O -> O-                      : 1.9e7
                 H + H2O -> H2 + OH                   : 1e-2
                 H + H -> H2                          : 7.5e6
                 H + OH -> H2O                        : 7e6
                 H + OH- -> H2O + emliq               : 2.2e4
                 H + H2O2 -> H2O + OH                 : 9e4
                 H2 + H2O2 -> H2O + OH + H            : 6e3
                 H + O2 -> HO2                        : 2.1e7
                 H + HO2 -> H2O2                      : 1e7
                 O + H2O -> OH + OH                   : 1.3e1
                 O + O2 -> O3                         : 3e6
                 OH + OH -> H2O2                      : 5.5e6
                 OH + O- -> HO2-                      : 2e8
                 OH + H2 -> H + H2O                   : 4.2e4
                 OH + OH- -> O- + H2O                 : 1.3e7
                 OH + HO2 -> O2 + H2O                 : 6e6
                 OH + O2- -> O2 + OH-                 : 8e6
                 O- + H2O -> OH- + OH                 : 1.8e3
                 O- + H2 -> OH- + H                   : 8e4
                 O- + H2O2 -> O2- + H2O               : 5e5
                 O- + HO2- -> O2- + OH-               : 4e5
                 O- + O2 -> O3-                       : 3.6e6
                 O- + O2- + H2O -> OH- + OH- + O2     : 6e2
                 OH + H2O2 -> H2O + HO2               : 2.7e4
                 OH + HO2- -> OH- + HO2               : 7.5e6
                 H2O+ + H2O -> H3O+ + OH              : 6
                 H3O+ + OH- -> H2O + H + OH           : 6e7
                 HO2 + H2O -> H3O+ + O2-              : 2
                 H3O+ + O2- -> H2O + HO2              : 6e-2
                 #emliq + NO2- -> NO2_2-               : 5.2e6
                 #emliq + NO3- -> NO3_2-               : 7e6
                 H2O -> H+ + OH-                      : 1.3963e-3
                 H+ + OH- -> H2O                      : 1.3973e8'
  [../]
[]
