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
    file = 'mesh01_gas.msh'
    #file = 'mesh02_gas.msh'
  [../]
  # The next two definitions create boundary conditions named
  # 'left' and 'right', where 'left' is at x = 0 and 'right' is at x = 1.1 mm.
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = geo
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
[]

[Executioner]
  type = Transient
  end_time = 1e6
  automatic_scaling = true
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  #solve_type = pjfnk
  petsc_options_iname = '-pc_type -pc_factor_shift-type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-9
  dtmin = 1e-18
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 2e-15
    growth_factor = 1.4
    optimal_iterations = 10
  [../]
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  #[out_2kV_330kOhm]
  [out]
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
    #ballast_resist = 1e6
    #ballast_resist = 9.0e5
    #ballast_resist = 8.5e5
    #ballast_resist = 7.0e5
    #ballast_resist = 5.5e5
    #ballast_resist = 4.0e5
    #ballast_resist = 3.5e5
    #ballast_resist = 2.5e5
    #ballast_resist = 2.0e5
    #ballast_resist = 1.5e5
    #ballast_resist = 1.0e5
    ballast_resist = 2.2e5
    e = 1.6e-19
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[DriftDiffusionActionAD]
  [./Plasma]
    electrons = em
    #charged_particle = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    #Neutrals = 'Ars Arss Arsss Ar2s H2O OH O2 O2s O Os H2 H Hs O3 HO2 H2O2 OHs'
    charged_particle = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
    Neutrals = 'Ars Arss Arsss Ar2s H2O OH O2 O2s O Os H2 H Hs O3 HO2 H2O2 OHs'
    mean_energy = mean_en
    potential = potential
    #Is_potential_unique = false
    using_offset = true
    offset = 20
    use_ad = true
    order = FIRST 
    position_units = ${dom0Scale}
    block = 0
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
  [Hp]
    block = 0
    initial_condition = -30
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

  [./Arp]
    block = 0
    initial_condition = -20.693147
  [../]
  [./Ar2p]
    block = 0
    initial_condition = -20.693147
  [../]
  [O2]
    block = 0
    initial_condition = -30
  []
  [OHp]
    block = 0
    initial_condition = -30
  []
  [OHm]
    block = 0
    initial_condition = -30
  []

  [./Ars]
    block = 0
    initial_condition = -25
  [../]
  [./Arss]
    block = 0
    initial_condition = -25
  [../]
  [./Arsss]
    block = 0
    initial_condition = -25
  [../]
  [./Ar2s]
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

  # Humid Species (H2O, H2, O2, and their derivatives)
  [O2s]
    block = 0
    initial_condition = -30
  []
  [O]
    block = 0
    initial_condition = -30
  []
  [Om]
    block = 0
    initial_condition = -30
  []
  [O2p]
    block = 0
    initial_condition = -30
  []
  [O2m]
    block = 0
    initial_condition = -30
  []
  [Os]
    block = 0
    initial_condition = -30
  []
  [OH]
    block = 0
    initial_condition = -30
  []
  [OHs]
    block = 0
    initial_condition = -30
  []
  [H2]
    block = 0
    initial_condition = -30
  []
  [H]
    block = 0
    initial_condition = -30
  []
  [Hs]
    block = 0
    initial_condition = -30
  []
  [O3]
    block = 0
    initial_condition = -30
  []
  [HO2]
    block = 0
    initial_condition = -30
  []
  [H2O2]
    block = 0
    initial_condition = -30
  []


[]

[AuxVariables]
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.701920755421197
    #initial_condition = 3.704261
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
  [./Current_Arp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_Ar2p]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_H2Op]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_OHp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_OHm]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_Om]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_O2m]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_O2p]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]
  [./Current_Hp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 0
  [../]

  [./tot_gas_current]
    order = CONSTANT
    family = MONOMIAL
    block = 0
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
[]

[AuxKernels]
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
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [rho_calc]
    type = ChargeDensity
    variable = rho
    charged = 'em Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  []
  [./Efield_g]
    type = Efield
    component = 0
    potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ADCurrent_em]
    type = ADCurrent
    potential = potential
    density_log = em
    variable = Current_em
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_Arp]
    type = ADCurrent
    potential = potential
    density_log = Arp
    variable = Current_Arp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_Ar2p]
    type = ADCurrent
    potential = potential
    density_log = Ar2p
    variable = Current_Ar2p
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_H2Op]
    type = ADCurrent
    potential = potential
    density_log = H2Op
    variable = Current_H2Op
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_OHp]
    type = ADCurrent
    potential = potential
    density_log = OHp
    variable = Current_OHp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_OHm]
    type = ADCurrent
    potential = potential
    density_log = OHm
    variable = Current_OHm
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_Om]
    type = ADCurrent
    potential = potential
    density_log = Om
    variable = Current_Om
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_O2m]
    type = ADCurrent
    potential = potential
    density_log = O2m
    variable = Current_O2m
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_O2p]
    type = ADCurrent
    potential = potential
    density_log = O2p
    variable = Current_O2p
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADCurrent_Hp]
    type = ADCurrent
    potential = potential
    density_log = Hp
    variable = Current_Hp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [tot_gas_current]
    type = ParsedAux
    variable = tot_gas_current
    args = 'Current_em Current_Arp Current_Ar2p Current_H2Op Current_OHp Current_OHm Current_Om Current_O2m Current_O2p Current_Hp'
    function = 'Current_em + Current_Arp + Current_Ar2p + Current_H2Op + Current_OHp + Current_OHm + Current_Om + Current_O2m + Current_O2p + Current_Hp'
    execute_on = 'timestep_end'
    block = 0
  []

  [./ADEFieldAdvAux_em]
    type = ADEFieldAdvAux
    potential = potential
    density_log = em
    variable = EFieldAdvAux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./ADDiffusiveFlux_em]
    type = ADDiffusiveFlux
    density_log = em
    variable = DiffusiveFlux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
[]

[BCs]
  [H2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2
    boundary = 'left right'
    #boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  # H2O+ boundary conditions
  [./H2Op_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2Op 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./H2Op_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = H2Op
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Om_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Om 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Om_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Om
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHm_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHm 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHm_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = OHm
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHp_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHp 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHp_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = OHp
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2m_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2m 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2m_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = O2m
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2p 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = O2p
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./Hp_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Hp 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Hp_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Hp
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]

  # H2O evaporation boundary condition
  [H2O_interface]
    type = DirichletBC
    variable = H2O
    value = 0.367321
    boundary = 'right'
  []
  #[OH_left]
  #  type = DirichletBC
  #  variable = OH
  #  value = -10
  #  boundary = 'left'
  #[]
  [OH_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OH
    boundary = 'left right'
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
  [O2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [O2s_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2s
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [O_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Os_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Os
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [H_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Hs_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Hs
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [O3_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O3
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [HO2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = HO2
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [H2O2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2O2
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [OHs_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHs 
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  []

  [./Arex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ars
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex2_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arss
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex3_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arsss
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2ex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2s
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2p
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Ar2p
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_left]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
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
    boundary = 'right'
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
    boundary = 'right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'right'
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
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
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
  #[./em_aq_ic]
  #  type = FunctionIC
  #  variable = em_aq
    #function = em_aq_ic_func
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
    #value = -1.25
    #value = 1.0
    #value = -2.5
    value = -0.8
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  [./em_aq_ic_func]
    type = ParsedFunction
    value = '1778 - 1.8e6*x'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-1.25 * (1.001e-3 - x)'
    #value = '-2.5 * (1.001e-3 - x)'
    value = '-0.8 * (1.010e-3 - x)'
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
    boundary = 'left right'
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
   property_tables_file = 'data_files/electron_mobility_diffusion.txt'
   block = 0
 [../]
 [gas_constants]
   type = GenericConstantMaterial
   prop_names = 'e    N_A    k_boltz    eps     se_energy    T_gas    massem    p_gas diffpotential'
   prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12 1 400 9.11e-31 1.01e5 8.854e-12'
   block = 0
 []

  [./gas_species_O2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
    block = 0
  [../]

  [./gas_species_O2s]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2s
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
    block = 0
  [../]

  [./gas_species_Os]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Os
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 0
    diffusivity = 6e-5
    block = 0
  [../]

  [./gas_species_O]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 0
    diffusivity = 6e-5
    block = 0
  [../]

  [./gas_species_HO2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO2
    heavy_species_mass = 5.481069e-26
    heavy_species_charge = 0
    diffusivity = 2e-5
    block = 0
  [../]

  [./gas_species_H2O2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2O2
    heavy_species_mass = 5.64829e-26
    heavy_species_charge = 0
    diffusivity = 2e-5
    block = 0
  [../]

  [./gas_species_H]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H
    heavy_species_mass = 1.673597e-27
    heavy_species_charge = 0
    diffusivity = 8.8e-5
    block = 0
  [../]

  [./gas_species_Hs]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Hs
    heavy_species_mass = 1.673597e-27
    heavy_species_charge = 0
    diffusivity = 8.8e-5
    block = 0
  [../]

  [./gas_species_Hp]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Hp
    heavy_species_mass = 1.673597e-27
    heavy_species_charge = 1
    diffusivity = 8.8e-5
    block = 0
  [../]

  [./gas_species_H2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2
    heavy_species_mass = 3.347526e-27
    heavy_species_charge = 0
    diffusivity = 7.8e-5
    block = 0
  [../]

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
    heavy_species_name = Ars
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Arss
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_4]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Arsss
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./gas_species_5]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar2s
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
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

  [./gas_species_OHp]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OHp
    heavy_species_mass = 2.824311e-26
    heavy_species_charge = 1
    diffusivity = 4e-5
    block = 0
  [../]

  [./gas_species_OHs]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OHs
    heavy_species_mass = 2.824311e-26
    heavy_species_charge = 0
    diffusivity = 4e-5
    block = 0
  [../]

  [./gas_species_O3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O3
    heavy_species_mass = 7.970475e-26
    heavy_species_charge = 0
    diffusivity = 2e-5 
    block = 0
  [../]

  [./gas_species_Op]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Op
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 1
    diffusivity = 5.8e-5
    block = 0
  [../]

  [./gas_species_O2p]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2p
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 1
    diffusivity = 5.6e-5
    block = 0
  [../]

  [./gas_species_Om]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Om
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = -1
    diffusivity = 7.0e-5
    block = 0
  [../]

  [./gas_species_O2m]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2m
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = -1
    diffusivity = 5.6e-5
    block = 0
  [../]
[]

[Reactions]
  [./Argon]
    #species = 'em Arp Ar2p Ars Arss Arsss OH H2O O2 O2p O2m O2s Om O Os H2O H2Op H2 H Hs O3 HO2 H2O2 OHm OHp OHs'
    species = 'em Arp Ar2p Ars Arss Arsss Ar2s OH H2O O2 O2p O2m O2s Om O Os H2Op H2 H Hs O3 HO2 H2O2 OHm OHp OHs Hp'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'data_files'
    equation_constants = 'Tgas Tn'
    equation_values = '300 1'
    equation_variables = 'Te'
    potential = 'potential'
    use_log = true
    convert_to_moles = true
    convert_to_meters = 1e-2
    position_units = ${dom0Scale}
    use_ad = true
    block = 0

    # Note that the monatomic hydrogen reactions were taken from 
    # Van Gaens et al. "Kinetic modeling for an atmospheric pressure argon 
    # plasma jet in humid air." J. Phys. D: Appl. Phys. 46 275201 (2013)
    #
    # Ars    - Ar(1s1), Ar(1s3) (metastable states of Ar(1s) manifold)
    # Arss   - Ar(1s2), Ar(1s4) (radiative states of Ar(1s) manifold)
    # Arsss  - Ar(4p) and higher states
    #reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (reaction1)
    #             em + Ar -> em + Ars              : EEDF [-11.5]   (reaction2)
    #             em + Ar -> em + em + Arp         : EEDF [-15.76]  (reaction3)
    #             em + Ars -> em + Ar              : EEDF [11.5]    (reaction4)
    #             em + Ars -> em + em + Arp        : EEDF [-4.3]    (reaction5)
    #             Ar2p + em -> Ars + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
    #             Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tgas * exp(-15130/Tgas)}
    #             Ars + Ars -> Ar2p + em           : {3.6132e8}
    #             Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
    #             Ars + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
    #             Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}
    reactions = 'em + Ar -> em + Ar                       : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ars                      : EEDF [-11.549] (C2_Ar_Excitation_11.55_eV)
                 em + Ar -> em + Ars                      : EEDF [-11.723] (C4_Ar_Excitation_11.72_eV)
                 em + Ar -> em + Arss                     : EEDF [-11.624] (C3_Ar_Excitation_11.62_eV)
                 em + Ar -> em + Arss                     : EEDF [-11.828] (C5_Ar_Excitation_11.83_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.074] (C6_Ar_Excitation_13.07_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.204] (C7_Ar_Excitation_13.20_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.227] (C8_Ar_Excitation_13.23_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.276] (C9_Ar_Excitation_13.28_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.299] (C10_Ar_Excitation_13.30_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.396] (C11_Ar_Excitation_13.40_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.397] (C12_Ar_Excitation_13.40_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.418] (C13_Ar_Excitation_13.42_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.442] (C14_Ar_Excitation_13.44_eV)
                 em + Ar -> em + Arsss                    : EEDF [-13.594] (C15_Ar_Excitation_13.59_eV)
                 em + Ars -> em + Arss                    : EEDF [-0.075] (C27_Ar1s5_Excitation_0.075_eV)
                 em + Ars -> em + Arss                    : EEDF [-0.279] (C28_Ar1s5_Excitation_0.28_eV)
                 #em + Arss -> em + Ars                    : EEDF [] (C_Ar_Excitation__eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.525] (C29_Ar1s5_Excitation_1.52_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.655] (C30_Ar1s5_Excitation_1.66_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.678] (C31_Ar1s5_Excitation_1.68_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.727] (C32_Ar1s5_Excitation_1.73_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.75] (C33_Ar1s5_Excitation_1.75_eV)
                 #em + Arsss -> em + Ars                   : EEDF [] (C_Ar_Excitation__eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.45] (C17_Ar1s4_Excitation_1.45_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.58] (C18_Ar1s4_Excitation_1.58_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.60] (C19_Ar1s4_Excitation_1.60_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.65] (C20_Ar1s4_Excitation_1.65_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.68] (C21_Ar1s4_Excitation_1.68_eV)
                 #em + Arsss -> em + Arss                  : EEDF [] (C_Ar_Excitation__eV)
                 em + Ar -> em + em + Arp                 : EEDF [-15.76] (C16_Ar_Ionization_15.76_eV)
                 em + Ars -> em + em + Arp                : EEDF [-4.21] (C38_Ar1s5_Ionization_4.21_eV)
                 em + Arss -> em + em + Arp               : EEDF [-4.14] (C26_Ar1s4_Ionization_4.14_eV)
                 #em + Arsss -> em + em + Arp              : EEDF [] (C_Ar_Ionization__eV)
                 em + Arp -> Arsss                        : {4e-13*Te^(-0.5)}
                 em + em + Arp -> Arsss + em              : {5e-27*Te^(-4.5)}
                 em + Ar2s -> Ar2p + em + em              : {9e-8*Te^0.7}
                 em + Ar2s -> Ar + Ar + em                : 1e-7
                 em + Ar2p -> Arsss + Ar                  : {5.38e-8*Te^(-0.66)}
                 Ars + Ars -> Ar + Arp + em               : {5e-10*Tn^0.5}
                 Arss + Arss -> Ar + Arp + em             : {5e-10*Tn^0.5}
                 Arsss + Arsss -> Ar + Arp + em           : {5e-10*Tn^0.5}
                 #Arp + Ar -> Ar + Arp                     : {5.66e-10*Tn^0.5}
                 Arp + Ar + Ar -> Ar + Ar2p               : {1.41e-31*Tn^(-0.5)}
                 Ars + Ar + Ar -> Ar + Ar2s               : {1.14e-32}
                 Arss + Ar + Ar -> Ar + Ar2s              : {1.14e-32}
                 Arsss + Ar + Ar -> Ar + Ar2s             : {1.14e-32}
                 #########################
                 # Ar-H2O reactions
                 #########################
                 Arp + H2O -> Ar + H2Op                  : 1.5e-10
                 Ar2p + H2O -> Ar + Ar + H2Op            : 1.5e-10
                 Ars + H2O -> Ar + OH + H                : 4.8e-10
                 Arss + H2O -> Ar + OH + H              : 4.8e-10
                 Arsss + H2O -> Ar + OHs + H            : 4.8e-10
                 Ar2s + H2O  -> Ar + OH + H             : 4.8e-10
                 Ar2s + H2O  -> Ar + OHs + H            : 4.8e-10
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
                 #########################
                 em + H2O -> em + H2O             : EEDF [elastic] (C56_H2O_Elastic)
                 #em + H2O -> H2Ov + em            : EEDF [-0.198] (C61_H2O_Excitation_0.20_eV)
                 em + H2O -> OHm + H              : EEDF (C54_H2O_Attachment)
                 # It is unclear what Hs* and Hs** refer to in the thesis
                 #em + H2O -> em + OH + Hs         : EEDF [] (C_H2O_Excitation__eV)
                 #em + H2O -> em + OH + Hs*        : EEDF [] (C_H2O_Excitation__eV)
                 #em + H2O -> em + OH + Hs**       : EEDF [] (C_H2O_Excitation__eV)
                 em + H2O -> em + OHs + H         : EEDF [-9.00] (C62_H2O_Excitation_9.00_eV) 
                 em + H -> em + Hs                : EEDF [-10.2] (C65_H_Excitation_10.20_eV)
                 em + H -> em + Hs                : EEDF [-10.2] (C66_H_Excitation_10.20_eV)
                 em + Hs -> em + H                : EEDF [10.2] (C68_H2p10.2eV_De-excitation_10.20_eV)
                 em + Hs -> em + H                : EEDF [10.2] (C69_H2s10.2eV_De-excitation_10.20_eV)
                 em + O2 -> em + O2               : EEDF [elastic] (C40_O2_Elastic)
                 em + O2 -> O2s + em              : EEDF [-0.98] (C42_O2_Excitation_0.98_eV)
                 em + O2 -> O2s + em              : EEDF [-1.63] (C43_O2_Excitation_1.63_eV)
                 em + O2s -> em + O2              : EEDF [0.98] (C46_O20.98_De-excitation_0.98_eV)
                 em + O2s -> em + O2              : EEDF [1.63] (C50_O21.63_De-excitation_1.63_eV)
                 em + O2 -> em + em + O2p         : EEDF [-12.1] (C44_O2_Ionization_12.10_eV)
                 em + O2s -> em + em + O2p        : EEDF [-11.10] (C48_O20.98_Ionization_11.10_eV)
                 em + O2s -> em + em + O2p        : EEDF [-10.45] (C52_O21.63_Ionization_10.45_eV)
                 #em + O2s -> em + em + Op + O     : EEDF [-18.52] (C49_O20.98_Ionization_18.52_eV)
                 #em + O2s -> em + em + Op + O     : EEDF [-17.87] (C53_O21.63_Ionization_17.87_eV)
                 em + O2 -> Om + O                : EEDF (C39_O2_Attachment)
                 em + O2s -> Om + O               : EEDF (C47_O20.98_Attachment)
                 em + O2s -> Om + O               : EEDF (C51_O21.63_Attachment)
                 em + O2p -> O + O                : {1.2e-8*Te^(-0.7)} 
                 #em + H -> em + Hs*               : EEDF [] (C_H_Excitation__eV)
                 #em + Hs* -> em + H               : EEDF [] (C_H_Excitation__eV)
                 #em + H -> em + Hs**              : EEDF [] (C_H_Excitation__eV)
                 #em + Hs** -> em + H              : EEDF [] (C_H_Excitation__eV)
                 #em + Hs -> em + Hs*              : EEDF [] (C_H_Excitation__eV)
                 #em + Hs* -> em + Hs              : EEDF [] (C_H_Excitation__eV)
                 #em + Hs -> em + Hs**             : EEDF [] (C_H_Excitation__eV)
                 #em + Hs** -> em + Hs             : EEDF [] (C_H_Excitation__eV)
                 #em + Hs* -> em + Hs**            : EEDF [] (C_H_Excitation__eV)
                 #em + Hs** -> em + Hs*            : EEDF [] (C_H_Excitation__eV)
                 # Cannot find cross section of OH
                 #em + OH -> em + OHs              : EEDF [] (C_OH_Excitation__eV)
                 em + OHs -> em + O + H           : {2.7e-10*Te^(0.5)} 
                 Hs + H2O -> H + H2O              : 9.1e-9
                 #Hs* + H2O -> H + H2O             : 9.1e-9 
                 #Hs** + H2O -> H + H2O            : 9.1e-9
                 #Hs* + H2O -> Hs + H2O            : 9.1e-9
                 #Hs** + H2O -> Hs + H2O           : 9.1e-9
                 #Hs** + H2O -> Hs* + H2O          : 9.1e-9
                 OHs + H2O -> OH + H2O            : 9.1e-9
                 OHm + H -> H2O + em              : 1.8e-9
                 OH + H -> H2O                    : {6.87e-31*Tn^(-2)}
                 OHs + H -> H2O                   : {6.87e-31*Tn^(-2)}
                 OHs + OHs + O2 -> H2O2 + O2      : {6.9e-31*Tn^(-0.8)} 
                 OHs + OHs + H2O -> H2O2 + H2O    : {6.9e-31*Tn^(-0.8)} 
                 H2 + HO2 -> H2O2 + H             : {5e-11*exp(-Tgas/11310)}
                 HO2 + HO2 -> H2O2 + O2           : {8.05e-11*Tn^(-1)}
                 HO2 + HO2 + O2 -> H2O2 + O2 + O2 : {1.9e-33*exp(980/Tgas)}
                 HO2 + HO2 + H2O -> H2O2 + O2 + H2O : {1.9e-33*exp(980/Tgas)}
                 HO2 + H2O -> H2O2 + OH           : {4.65e-11*exp(-11647/Tgas)} 
                 H + H2O2 -> HO2 + H2             : {8e-11*exp(-4000/Tgas)}
                 H + H2O2 -> OH + H2O             : {4e-11*exp(-2000/Tgas)}
                 O2 + H2O2 -> HO2 + HO2           : {9e-11*exp(-19965/Tgas)}
                 O + H2O2 -> HO2 + OH             : {1.4e-12*exp(-2000/Tgas)}
                 Os + H2O2 -> O2 + H2O            : 5.2e-10
                 OH + H2O2 -> HO2 + H2O           : {2.9e-12*exp(-160/Tgas)}
                 H2O2 -> OH + OH                  : {1.96e-9*Tn^(-4.86)*exp(-26800/Tgas)}
                 H + H2O2 -> OH + H2O             : {4e-11*exp(-2000/Tgas)}
                 Hs + H2O2 -> OH + H2O            : {4e-11*exp(-2000/Tgas)}
                 #Hs* + H2O2 -> OH + H2O           : {4e-11*exp(-2000/Tgas)}
                 #Hs** + H2O2 -> OH + H2O          : {4e-11*exp(-2000/Tgas)}
                 OHm + OHp + O2 -> H2O2 + O2      : {2e-25*Tn^(-2.5)}
                 OHm + OHp + H2O -> H2O2 + H2O    : {2e-25*Tn^(-2.5)}
                 OHm + H2Op + O2 -> OH + H2O + O2 : {2e-25*Tn^(-2.5)}
                 OHm + H2Op + H2O -> OH + H2O + H2O : {2e-25*Tn^(-2.5)}
                 O2m + O2p -> O2 + O2             : 2e-6
                 O2m + H2Op -> O2 + H2O           : 2e-6
                 Om + O2p -> O + O2               : 3e-6
                 Om + OHp + O2 -> HO2 + O2        : {2e-25*Tn^(-2.5)}
                 Om + OHp + H2O -> HO2 + H2O      : {2e-25*Tn^(-2.5)}
                 Om + H2Op + O2 -> O + H2O + O2   : {2e-25*Tn^(-2.5)}
                 Om + H2Op + H2O -> O + H2O + H2O : {2e-25*Tn^(-2.5)}
                 Om + O2 -> O2m + O               : 1.5e-20
                 Om + O3 -> O2 + O2m              : 1e-11
                 O2m + O -> O3 + em               : 1.5e-10
                 O2m + O -> Om + O2               : 1.5e-10
                 O2m + O2s -> O2 + O2 + em        : 2e-10
                 O2s + O2 -> O2 + O2              : 2.2e-18
                 O2s + H2O -> O2 + H2O            : 2.2e-18
                 O2s + O2 -> O + O3               : 2.9e-21
                 O2s + O3 -> O2 + O2 + O          : 9.9e-11
                 ##########################################
                 # Radiative Transitions Ar
                 ##########################################
                 Arsss -> Ars                           : 3.3e7
                 Arsss -> Ar                            : 3.1e5
                 Arss -> Ar                             : 5.3e5
                 Ar2s -> Ar + Ar                        : 6e7
                 ##########################################
                 # Radiative Transitions H2O, OH, etc.
                 ##########################################
                 OHs -> OH                              : 1.3e6
                 Hs -> H                                : 4.7e8
                 ##########################################
                 # Additional reactions (from Van Gaens et al)
                 # Includes Ar and H reactions
                 ##########################################
                 Ars + OH -> Ar + OHs           : 6.6e-11 
                 Arp + H -> Ar + Hp             : 1e-10
                 Arp + H2 -> Ar + H2Op          : 1.1e-9
                 Ar2p + H -> Ar + Ar + Hp       : 5e-11
                 Ar2p + OHm -> Ar + Ar + OH     : 1e-7
                 Ar2p + OHm -> Ar + Ar + O + H  : 1e-7
                 H2Op + Ar -> Arp + H2          : 2.2e-10
                 # Next one is commented because Op is not tracked
                 #H + Op -> Hp + O               : {5.66e-10*(Tgas/300)^0.36*exp(8.6/Tgas)}
                 H + Om -> OH + em              : 5e-10
                 H + O2 + Ar -> HO2 + Ar        : {6.09e-32*(Tgas/300)^(-0.8)}
                 H + H + Ar -> H2 + Ar          : {2e-32*(Tgas/300)^(-1)}
                 Hp + OH -> H + OHp             : 2.1e-9
                 Hp + H2O -> H2Op + H           : 6.9e-9
                 H2Op + O2 -> H2 + O2p          : 8e-10
                 OH + O -> H + O2               : {2.08e-11*(Tgas/300)^(-0.186)*exp(-154/Tgas)}
                 OH + O2m -> OHm + O2           : 1e-10
                 OH + O3 -> HO2 + O2            : 1.69e-12
                 #################################
                 # H2O+ Rxns
                 # (Van Gaens, page 40)
                 #################################
                 H2Op + Om + Ar -> H2O2 + Ar    : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + H2O -> H2O2 + H2O  : {2e-25*(Tgas/300)^(-2.5)}
                 em + H2Op -> O + H2            : {6.27e-9*(Tgas/300)^(-0.5)}
                 em + H2Op -> O + H + H         : {4.9e-8*(Tgas/300)^(-0.5)}
                 ##################################
                 # Electron-impact of Hydrogen
                 ##################################
                 em + H -> em + em + Hp         : EEDF [-13.60] (C67_H_Ionization_13.60_eV)
                 ###################################
                 # OH, OHp reactions
                 # (Van Gaens)
                 ###################################
                 em + OHp -> O + H              : {6.03e-9*(Tgas/300)^(-0.5)}
                 em + OH -> em + O + H          : {2.08e-7*(Tgas/300)^(-0.76)*exp(-6.9/Tgas)}'
  [../]
[]
