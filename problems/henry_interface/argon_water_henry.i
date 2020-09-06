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
    file = 'mesh.msh'
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
  petsc_options_iname = '-pc_type -pc_factor_shift-type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu NONZERO 1.e-10 superlu_dist'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_package -snes_stol'
  #petsc_options_value = 'lu NONZERO 1.e-10 mumps 0'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  #petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  #nl_abs_tol = 1e-8
  dtmin = 1e-18
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    #dt = 2.56e-13
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
  #[m1_kV_1um_02]
  [out_01]
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

[DriftDiffusionActionAD]
  [./Plasma]
    electrons = em
    charged_particle = 'Arp Ar2p H2Op'
    Neutrals = 'Ar* Ar** Ar*** Ar2* H2O OH'
    mean_energy = mean_en
    potential = potential
    #Is_potential_unique = false
    using_offset = true
    offset = 30
    use_ad = true
    position_units = ${dom0Scale}
    block = 0
  [../]

  [./Water]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    charged_particle = 'emliq OH-'
    Neutrals = 'OH_aq'
    #Is_potential_unique = false
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    position_units = ${dom1Scale}
    block = 1
  [../]
[]

[Variables]
  [H2O]
    block = 0
    #initial_condition = 0.0367321
  []
  [H2Op]
    block = 0
    initial_condition = -25
  []
  [./potential]
    #block = 0
  [../]
  #[./potential_liq]
  #  block = 1
  #[../]
  [./em]
    block = 0
    initial_condition = -20
  [../]
  [./emliq]
    block = 1
    #initial_condition = -24
    #initial_condition = -21
    # scaling = 1e-5
    initial_condition = -14
    #initial_condition = -24
  [../]

  [./Arp]
    block = 0
    initial_condition = -20.693147
  [../]
  [./Ar2p]
    block = 0
    initial_condition = -20.693147
  [../]
  [OH]
    block = 0
    initial_condition = -20
  []
  [./OH_aq]
    block = 1
    initial_condition = -20
  [../]

  [./Ar*]
    block = 0
    initial_condition = -25
  [../]
  [./Ar**]
    block = 0
    initial_condition = -25
  [../]
  [./Ar***]
    block = 0
    initial_condition = -25
  [../]
  [./Ar2*]
    block = 0
    initial_condition = -25
  [../]
  #[./Om]
  #  block = 0
  #  initial_condition = -25
  #[../]
  #[./OHm]
  #  block = 0
  #  initial_condition = -25
  #[../]

  [./mean_en]
    block = 0
    initial_condition = -20
    # scaling = 1e-1
  [../]

  [./OH-]
    block = 1
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    #initial_condition = -9.210340
    initial_condition = -14
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

  [./H2O_aq]
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
  [./H2O_aq_density]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]

  [./Te]
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
  [./ADCurrent_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./ADCurrent_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./ADCurrent_Arp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./ADCurrent_OH-]
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
  [./ADEFieldAdvAux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./ADDiffusiveFlux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./ADEFieldAdvAux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./ADDiffusiveFlux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./ADPowerDep_em]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
  [./ADPowerDep_Arp]
   order = CONSTANT
   family = MONOMIAL
   block = 0
    initial_condition = 0
  [../]
[]

[AuxKernels]
  [./H2O_aq_density]
    type = DensityMoles
    variable = H2O_aq_density
    density_log = H2O_aq
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./ADPowerDep_em]
    type = ADPowerDep
    density_log = em
    potential = potential
    art_diff = false
    execute_on = 'initial timestep_end'
    potential_units = kV
    variable = ADPowerDep_em
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ADPowerDep_Arp]
    type = ADPowerDep
    density_log = Arp
    potential = potential
    art_diff = false
    potential_units = kV
    variable = ADPowerDep_Arp
    execute_on = 'initial timestep_end'
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./e_temp]
    type = ElectronTemperature
    variable = Te
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
  [./ADCurrent_em]
    type = ADCurrent
    potential = potential
    density_log = em
    variable = ADCurrent_em
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_emliq]
    type = ADCurrent
    #potential = potential_liq
    potential = potential
    density_log = emliq
    variable = ADCurrent_emliq
    art_diff = false
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./ADCurrent_Arp]
    type = ADCurrent
    potential = potential
    density_log = Arp
    variable = ADCurrent_Arp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_OH-]
    block = 1
    type = ADCurrent
    #potential = potential_liq
    potential = potential
    density_log = OH-
    variable = ADCurrent_OH-
    art_diff = false
    position_units = ${dom1Scale}
  [../]
  #[./tot_flux_OH-]
  #  block = 1
  #  type = TotalFlux
  #  #potential = potential_liq
  #  potential = potential
  #  density_log = OH-
  #  variable = tot_flux_OH-
  #[../]
  [./ADEFieldAdvAux_em]
    type = ADEFieldAdvAux
    potential = potential
    density_log = em
    variable = ADEFieldAdvAux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADDiffusiveFlux_em]
    type = ADDiffusiveFlux
    density_log = em
    variable = ADDiffusiveFlux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADEFieldAdvAux_emliq]
    type = ADEFieldAdvAux
    #potential = potential_liq
    potential = potential
    density_log = emliq
    variable = ADEFieldAdvAux_emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./ADDiffusiveFlux_emliq]
    type = ADDiffusiveFlux
    density_log = emliq
    variable = ADDiffusiveFlux_emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
[]

[InterfaceKernels]
  [./em_advection]
    type = ADInterfaceAdvection
    mean_en_neighbor = mean_en
    potential_neighbor = potential
    neighbor_var = em
    variable = emliq
    boundary = master1_interface
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    #type = InterfaceLogDiffusionElectrons
    type = ADInterfaceLogDiffusion
    mean_en_neighbor = mean_en
    neighbor_var = em
    variable = emliq
    boundary = master1_interface
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [henry_test]
    type = ADHenryInterface
    variable = OH_aq
    neighbor_var = OH
    h = 6.92e2
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'master1_interface'
  []
[]

[BCs]
  # H2O+ boundary conditions
  [./H2Op_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2Op 
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./H2Op_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = H2Op
    boundary = 'left master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  # O- and OH- boundary conditions
  #[./Om_physical_diffusion]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = Om 
  #  boundary = 'left master0_interface'
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Om_physical_advection]
  #  type = ADHagelaarIonAdvectionBC
  #  variable = Om
  #  boundary = 'left master0_interface'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./OHm_physical_diffusion]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = OHm 
  #  boundary = 'left master0_interface'
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./OHm_physical_advection]
  #  type = ADHagelaarIonAdvectionBC
  #  variable = OHm
  #  boundary = 'left master0_interface'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]

  # H2O evaporation boundary condition
  [H2O_interface]
    type = DirichletBC
    variable = H2O
    value = 0.367321
    boundary = 'master0_interface'
  []
  #[OH_left]
  #  type = DirichletBC
  #  variable = OH
  #  value = -10
  #  boundary = 'left'
  #[]
  [OH_right]
    type = ADHagelaarIonDiffusionBC
    variable = OH
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  []
  #[H2O_physical_left]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = H2O
  #  boundary = 'left'
  #  r = 0
  #  position_units = ${dom0Scale}
  #[]

  [./Arex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar*
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex2_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar**
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex3_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar***
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2ex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2*
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./Arex_henry]
  #  type = MatchedValueLogBC
  #  variable = Ar_aq
  #  v = Ar*
  #  H = 0.1
  #  boundary = 'master1_interface'
  #[../]
  [./Ar2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2p
    boundary = 'left master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Ar2p
    boundary = 'left master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_left]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    ip = 'Arp Ar2p H2Op'
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
    type = ADHagelaarElectronBC
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
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'master0_interface'
    potential = potential
    em = em
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./em_physical_left]
    type = ADHagelaarElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    ip = 'Arp Ar2p H2Op'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    ip = 'Arp Ar2p H2Op'
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'left'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./emliq_right]
    type = ADDCIonBC
    variable = emliq
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./OH-_physical]
    type = ADDCIonBC
    variable = OH-
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
  [H2O_ic]
    type = FunctionIC
    variable = H2O
    function = water_ic_func
  []
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
  [./water_ic_func]
    type = ParsedFunction
    value = 'log(8.6949e23/6.022e23)'
    # close to the boundary condition, essentially
  [../]
  [./potential_bc_func]
    type = ParsedFunction
    value = 0.8
    #value = 1.0
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  [./emliq_ic_func]
    type = ParsedFunction
    value = '1778 - 1.8e6*x'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-0.8 * (1.001e-3 - x)'
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
    boundary = 'left master0_interface'
  [../]
 [./GasBasics]
   #type = GasBase
   type = ADGasElectronMoments
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   # user_p_gas = 1.01325e5
   user_p_gas = 1.01e5
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   #property_tables_file = '/home/shane/projects/zapdos/problems/argon_cm_test/electron_mobility_diffusion.txt'
   property_tables_file = 'lumped_argon_rates/electron_mobility_diffusion.txt'
   block = 0
 [../]
 [gas_constants]
   type = GenericConstantMaterial
   prop_names = 'e    N_A    k_boltz    eps     se_energy    T_gas    massem    p_gas diffpotential'
   prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12 1 400 9.11e-31 1.01e5 8.854e-12'
   block = 0
 []

  [H2O_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2O
    heavy_species_mass = 2.9907e-26
    heavy_species_charge = 0
    diffusivity = 2.3e-5
  []
  [H2Op_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2Op
    heavy_species_mass = 2.9907e-26
    heavy_species_charge = 1
    diffusivity = 2.3e-5
  []

  [./gas_species_0]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    block = 0
  [../]
  [./Ar_species]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  [../]
  [./gas_species_1]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar2p
    heavy_species_mass = 13.28e-26
    heavy_species_charge = 1.0
    block = 0
  [../]
  [./gas_species_2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar**
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_4]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar***
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_5]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar2*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]

  [./OH-_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH-
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]

  [OH_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = 0 
    diffusivity = 4e-5
    block = 0
  []
  [OHm_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OHm
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 4e-5
    block = 0
  []
  [Om_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Om
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1 
    diffusivity = 5e-5
    block = 0
  []

  [./OH_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH_aq
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = 0 
    diffusivity = 5.29e-9
    block = 1
  [../]

  [./O2_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 2e-9
    block = 1
  [../]

  [./electron_data]
    type = ADGenericConstantMaterial
    prop_names = 'diffemliq muemliq Temliq'
    prop_values = '4.5e-9 0.000173913 300'
    block = 1
  [../]
  [./electron_sign]
    # I increased the electron mass by a factor of 10 here
    # Not sure what it's really supposed to be
    type = GenericConstantMaterial
    prop_names = 'sgnemliq massemliq'
    prop_values = '-1 9.11e-30'
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
  [./Argon]
    species = 'em Arp Ar2p Ar* Ar** Ar*** OH H2O'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    #file_location = 'argon_chemistry_rates'
    file_location = 'lumped_argon_rates'
    equation_constants = 'Tgas Tn'
    equation_values = '300 1'
    equation_variables = 'Te'
    potential = 'potential'
    use_log = true
    position_units = ${dom0Scale}
    use_ad = true
    block = 0

    # Ar*    - Ar(1s1), Ar(1s3) (metastable states of Ar(1s) manifold)
    # Ar**   - Ar(1s2), Ar(1s4) (radiative states of Ar(1s) manifold)
    # Ar***  - Ar(4p) and higher states
    #reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (reaction1)
    #             em + Ar -> em + Ar*              : EEDF [-11.5]   (reaction2)
    #             em + Ar -> em + em + Arp         : EEDF [-15.76]  (reaction3)
    #             em + Ar* -> em + Ar              : EEDF [11.5]    (reaction4)
    #             em + Ar* -> em + em + Arp        : EEDF [-4.3]    (reaction5)
    #             Ar2p + em -> Ar* + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
    #             Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tgas * exp(-15130/Tgas)}
    #             Ar* + Ar* -> Ar2p + em           : {3.6132e8}
    #             Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
    #             Ar* + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
    #             Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}
    reactions = 'em + Ar -> em + Ar                       : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ar*                      : EEDF [-11.549] (C2_Ar_Excitation_11.55_eV)
                 em + Ar -> em + Ar*                      : EEDF [-11.723] (C4_Ar_Excitation_11.72_eV)
                 em + Ar -> em + Ar**                     : EEDF [-11.624] (C3_Ar_Excitation_11.62_eV)
                 em + Ar -> em + Ar**                     : EEDF [-11.828] (C5_Ar_Excitation_11.83_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.074] (C6_Ar_Excitation_13.07_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.204] (C7_Ar_Excitation_13.20_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.227] (C8_Ar_Excitation_13.23_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.276] (C9_Ar_Excitation_13.28_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.299] (C10_Ar_Excitation_13.30_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.396] (C11_Ar_Excitation_13.40_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.397] (C12_Ar_Excitation_13.40_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.418] (C13_Ar_Excitation_13.42_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.442] (C14_Ar_Excitation_13.44_eV)
                 em + Ar -> em + Ar***                    : EEDF [-13.594] (C15_Ar_Excitation_13.59_eV)
                 em + Ar* -> em + Ar**                    : EEDF [-0.075] (C27_Ar1s5_Excitation_0.075_eV)
                 em + Ar* -> em + Ar**                    : EEDF [-0.279] (C28_Ar1s5_Excitation_0.28_eV)
                 #em + Ar** -> em + Ar*                    : EEDF [] (C_Ar_Excitation__eV)
                 em + Ar* -> em + Ar***                   : EEDF [-1.525] (C29_Ar1s5_Excitation_1.52_eV)
                 em + Ar* -> em + Ar***                   : EEDF [-1.655] (C30_Ar1s5_Excitation_1.66_eV)
                 em + Ar* -> em + Ar***                   : EEDF [-1.678] (C31_Ar1s5_Excitation_1.68_eV)
                 em + Ar* -> em + Ar***                   : EEDF [-1.727] (C32_Ar1s5_Excitation_1.73_eV)
                 em + Ar* -> em + Ar***                   : EEDF [-1.75] (C33_Ar1s5_Excitation_1.75_eV)
                 #em + Ar*** -> em + Ar*                   : EEDF [] (C_Ar_Excitation__eV)
                 em + Ar** -> em + Ar***                  : EEDF [-1.45] (C17_Ar1s4_Excitation_1.45_eV)
                 em + Ar** -> em + Ar***                  : EEDF [-1.58] (C18_Ar1s4_Excitation_1.58_eV)
                 em + Ar** -> em + Ar***                  : EEDF [-1.60] (C19_Ar1s4_Excitation_1.60_eV)
                 em + Ar** -> em + Ar***                  : EEDF [-1.65] (C20_Ar1s4_Excitation_1.65_eV)
                 em + Ar** -> em + Ar***                  : EEDF [-1.68] (C21_Ar1s4_Excitation_1.68_eV)
                 #em + Ar*** -> em + Ar**                  : EEDF [] (C_Ar_Excitation__eV)
                 em + Ar -> em + em + Arp                 : EEDF [-15.76] (C16_Ar_Ionization_15.76_eV)
                 em + Ar* -> em + em + Arp                : EEDF [-4.21] (C38_Ar1s5_Ionization_4.21_eV)
                 em + Ar** -> em + em + Arp               : EEDF [-4.14] (C26_Ar1s4_Ionization_4.14_eV)
                 #em + Ar*** -> em + em + Arp              : EEDF [] (C_Ar_Ionization__eV)
                 em + Arp -> Ar***                        : {4e-13*Te^(-0.5)}
                 em + em + Arp -> Ar*** + em              : {5e-27 * Te^(-4.5)}
                 em + Ar2* -> Ar2p + em + em              : {9e-8*Te^0.7}
                 em + Ar2* -> Ar + Ar + em                : 1e-7
                 em + Ar2p -> Ar*** + Ar                  : {5.38e-8*Te^(-0.66)}
                 Ar* + Ar* -> Ar + Arp + em               : {5e-10*Tn^0.5}
                 Ar** + Ar** -> Ar + Arp + em             : {5e-10*Tn^0.5}
                 Ar*** + Ar*** -> Ar + Arp + em           : {5e-10*Tn^0.5}
                 #Arp + Ar -> Ar + Arp                     : {5.66e-10*Tn^0.5}
                 Arp + Ar + Ar -> Ar + Ar2p               : {1.41e-31*Tn^(-0.5)}
                 Ar* + Ar + Ar -> Ar + Ar2*               : {1.14e-32}
                 Ar** + Ar + Ar -> Ar + Ar2*              : {1.14e-32}
                 Ar*** + Ar + Ar -> Ar + Ar2*             : {1.14e-32}
                 #########################
                 # Ar-H2O reactions
                 #########################
                 Arp + H2O -> Ar + H2Op                  : 1.5e-10
                 Ar2p + H2O -> Ar + Ar + H2Op            : 1.5e-10
                 Ar* + H2O -> Ar + OH + H                : 4.8e-10
                 Ar** + H2O -> Ar + OH + H              : 4.8e-10
                 Ar*** + H2O -> Ar + OH* + H            : 4.8e-10
                 Ar2* + H2O  -> Ar + OH + H             : 4.8e-10
                 Ar2* + H2O  -> Ar + OH* + H            : 4.8e-10
                 #Arp + Om + Ar -> Ar + Ar + Ar          : {2e-25 * (Tgas/300)^(-2.5)}
                 #Arp + Om + H2O -> Ar + Ar + H2O        : {2e-25 * (Tgas/300)^(-2.5)}
                 #Ar2p + Om + Ar -> Ar + Ar + Ar         : {2e-25 * (Tgas/300)^(-2.5)} 
                 #Ar2p + Om + H2O -> Ar + Ar + H2O       : {2e-25 * (Tgas/300)^(-2.5)} 
                 #Arp + OHm + Ar -> Ar + Ar + OH         : {2e-25 * (Tgas/300)^(-2.5)}
                 #Arp + OHm + H2O -> Ar + H2O + OH       : {2e-25 * (Tgas/300)^(-2.5)}
                 #Ar2p + OHm + Ar -> Ar + Ar + OH        : {2e-25 * (Tgas/300)^(-2.5)}
                 #Ar2p + OHm + H2O -> Ar + Ar + H2O + OH : {2e-25 * (Tgas/300)^(-2.5)}
                 #########################
                 # H2O, OH, and H reactions
                 # (only a select few for now)
                 #########################
                 #OHm + H2Op + Ar -> OH + H2O + Ar       : {2e-25*Tn^(-2.5)}
                 #OHm + H2Op + H2O -> OH + H2O + H2O     : {2e-25*Tn^(-2.5)}
                 #Om + H2Op + Ar -> O + H2O + Ar         : {2e-25*Tn^(-2.5)}
                 #Om + H2Op + H2O -> O + H2O + H2O       : {2e-25*Tn^(-2.5)}
                 #########################
                 # Radiative Transitions
                 #########################
                 Ar*** -> Ar*                           : 3.3e7
                 Ar*** -> Ar                            : 3.1e5
                 Ar** -> Ar                             : 5.3e5
                 Ar2* -> Ar + Ar                        : 6e7'
  [../]


  [./water_reactions]
    species = 'emliq OH- OH_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'emliq + H2O_aq -> H + OH-               : 1.9e-2
                 #emliq + H2O+ -> H + OH               : 6e8
                 #emliq + emliq -> H2 + OH- + OH-      : 5.5e9
                 emliq + emliq -> H2 + OH- + OH-      : 3.0703e8'
                 #emliq + H + H2O_aq -> H2 + OH-          : 2.5e4
                 #emliq + OH -> OH-                    : 3e7
                 #emliq + O- + H2O_aq -> OH + OH          : 2.2e4
                 #emliq + H3O+ -> H + H2O_aq              : 2.3e7
                 #emliq + H2O2 -> OH + OH-             : 1.1e7
                 #emliq + HO2- + H2O_aq -> OH + OH- + OH- : 3.5e3
                 #emliq + O2 -> O2-                    : 1.9e7
                 #emliq + O -> O-                      : 1.9e7
                 #H + H2O_aq -> H2 + OH                   : 1e-2
                 #H + H -> H2                          : 7.5e6
                 #H + OH -> H2O_aq                        : 7e6
                 #H + OH- -> H2O_aq + emliq               : 2.2e4
                 #H + H2O2 -> H2O_aq + OH                 : 9e4
                 #H2 + H2O2 -> H2O_aq + OH + H            : 6e3
                 #H + O2 -> HO2                        : 2.1e7
                 #H + HO2 -> H2O2                      : 1e7
                 #O + H2O_aq -> OH + OH                   : 1.3e1
                 #O + O2 -> O3                         : 3e6
                 #OH + OH -> H2O2                      : 5.5e6
                 #OH + O- -> HO2-                      : 2e8
                 #OH + H2 -> H + H2O_aq                   : 4.2e4
                 #OH + OH- -> O- + H2O_aq                 : 1.3e7
                 #OH + HO2 -> O2 + H2O_aq                 : 6e6
                 #OH + O2- -> O2 + OH-                 : 8e6
                 #O- + H2O_aq -> OH- + OH                 : 1.8e3
                 #O- + H2 -> OH- + H                   : 8e4
                 #O- + H2O2 -> O2- + H2O_aq               : 5e5
                 #O- + HO2- -> O2- + OH-               : 4e5
                 #O- + O2 -> O3-                       : 3.6e6
                 #O- + O2- + H2O_aq -> OH- + OH- + O2     : 6e2
                 #OH + H2O2 -> H2O_aq + HO2               : 2.7e4
                 #OH + HO2- -> OH- + HO2               : 7.5e6
                 #H2O+ + H2O_aq -> H3O+ + OH              : 6
                 #H3O+ + OH- -> H2O_aq + H + OH           : 6e7
                 #HO2 + H2O_aq -> H3O+ + O2-              : 2
                 #H3O+ + O2- -> H2O_aq + HO2              : 6e-2
                 ##emliq + NO2- -> NO2_2-               : 5.2e6
                 ##emliq + NO3- -> NO3_2-               : 7e6
                 #H2O_aq -> H+ + OH-                      : 1.3963e-3
                 #H+ + OH- -> H2O_aq                      : 1.3973e8'
  [../]
[]
