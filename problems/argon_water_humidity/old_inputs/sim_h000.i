dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  offset = 20
  potential_units = kV
  use_moles = true
[]

[Mesh]
  [./geo]
    type = FileMeshGenerator
    file = 'mesh01.msh'
  [../]

  [./interface1]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'gas_right'
    input = geo
  [../]
  [./interface2]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'water_left'
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
[]

[Executioner]
  type = Transient
  end_time = 3600 # 1 hour
  automatic_scaling = true
  compute_scaling_once = false
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift-type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-9
  dtmin = 1e-18
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 2e-16
    growth_factor = 1.4
    optimal_iterations = 10
  []
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  [out_02]
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
    ballast_resist = 4.0e5
    #ballast_resist = 3.5e5
    #ballast_resist = 2.5e5
    #ballast_resist = 2.0e5
    #ballast_resist = 1.5e5
    #ballast_resist = 1.0e5
    e = 1.6e-19
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[DriftDiffusionAction]
  [./Plasma]
    electrons = em
    charged_particle = 'Arp Ar2p'
    Neutrals = 'Ars Arss Arsss Ar2s H2 O3'
    mean_energy = mean_en
    potential = potential
    #Is_potential_unique = false
    using_offset = true
    offset = 30
    use_ad = true
    order = FIRST 
    position_units = ${dom0Scale}
    block = 0
  [../]

  [./Water]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    #charged_particle = 'em_aq OHm_aq'
    #Neutrals = 'OH_aq'
    charged_particle = 'em_aq H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
    Neutrals = 'H_aq H2O2_aq OH_aq O2_aq O_aq H2_aq HO2_aq O3_aq HO3_aq'
    #Is_potential_unique = false
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    order = FIRST 
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./Salt]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    charged_particle = 'Nap_aq Clm_aq'
    First_DriftDiffusionAction_in_block = false
    potential = potential
    using_offset = true
    offset = -2.3026
    use_ad = true
    order = FIRST 
    position_units = ${dom1Scale}
    block = 1
  [../]
[]

[Variables]
  #[H2O]
  #  block = 0 
  #  initial_condition = -0.900556 # 1 % humidity
  #  #initial_condition = -0.2074 # 2 % humidity
  #  #initial_condition = 0.70881 # 5 % humidity
  #  #initial_condition = 1.402028 # 10 % humidity
  #[]

  [Nap_aq]
    block = 1
    # 100 mM
    #initial_condition = 4.60517
    # 10 mM
    initial_condition = 2.3026
  [] 
  [Clm_aq]
    block = 1
    # 100 mM
    #initial_condition = 4.60517
    # 10 mM
    initial_condition = 2.3026
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
  [H2]
    block = 0
    initial_condition = -30
  []
  [O3]
    block = 0
    initial_condition = -30
  []

  # Water Species
  [./em_aq]
    block = 1
    #initial_condition = -24
    #initial_condition = -21
    # scaling = 1e-5
    #initial_condition = -14
    #initial_condition = -24
    #initial_condition = -20
    initial_condition = -16
  [../]

  [./OHm_aq]
    block = 1
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    #initial_condition = -9.210340
    initial_condition = -14
  [../]

  [./OH_aq]
    block = 1
    initial_condition = -20
  [../]
  [./H3Op_aq]
    block = 1
    #initial_condition = -20
    initial_condition = -14
  [../]
  [./O2m_aq]
    block = 1
    initial_condition = -20
  [../]
  [./HO2m_aq]
    block = 1
    initial_condition = -20
  [../]
  [./H2Op_aq]
    block = 1
    initial_condition = -20
  [../]
  [./O3m_aq]
    block = 1
    initial_condition = -20
  [../]
  [./H_aq]
    block = 1
    initial_condition = -20
  [../]
  [./H2O2_aq]
    block = 1
    initial_condition = -20
  [../]
  [./O2_aq]
    block = 1
    initial_condition = -20
  [../]
  [./O_aq]
    block = 1
    initial_condition = -20
  [../]
  [./H2_aq]
    block = 1
    initial_condition = -20
  [../]
  [./HO2_aq]
    block = 1
    initial_condition = -20
  [../]
  [./HO3_aq]
    block = 1
    initial_condition = -20
  [../]
  [./O3_aq]
    block = 1
    initial_condition = -20
  [../]
  [./Om_aq]
    block = 1
    initial_condition = -20
  [../]
[]

[AuxVariables]
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL 
    initial_condition = 3.69456 # 1 % humidity 
    #initial_condition = 3.68441 # 2 % humidity 
    #initial_condition = 3.65332 # 5 % humidity 
    #initial_condition = 3.59925 # 10 % humidity 
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
  [./ADCurrent_em_aq]
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
  [./ADCurrent_OHm_aq]
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
  [./tot_flux_OHm_aq]
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
  [./ADEFieldAdvAux_em_aq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
    initial_condition = 0
  [../]
  [./ADDiffusiveFlux_em_aq]
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
  #[./rho]
  #  type = ParsedAux
  #  variable = rho
  #  args = 'em_density Arp_density'
  #  function = 'Arp_density - em_density'
  #  execute_on = 'timestep_end'
  #  block = 0
  #[../]
  [rho_calc]
    type = ChargeDensity
    variable = rho
    #charged = 'em Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
    charged = 'em Arp Ar2p'
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  []
  [rholiq_calc]
    type = ChargeDensity
    variable = rholiq
    charged = 'em_aq H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
    execute_on = 'INITIAL TIMESTEP_END'
    block = 1
  []
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
  [./ADCurrent_em_aq]
    type = ADCurrent
    #potential = potential_liq
    potential = potential
    density_log = em_aq
    variable = ADCurrent_em_aq
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
  [./ADCurrent_OHm_aq]
    block = 1
    type = ADCurrent
    #potential = potential_liq
    potential = potential
    density_log = OHm_aq
    variable = ADCurrent_OHm_aq
    art_diff = false
    position_units = ${dom1Scale}
  [../]
  #[./tot_flux_OHm_aq]
  #  block = 1
  #  type = TotalFlux
  #  #potential = potential_liq
  #  potential = potential
  #  density_log = OHm_aq
  #  variable = tot_flux_OHm_aq
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
  [./ADEFieldAdvAux_em_aq]
    type = ADEFieldAdvAux
    #potential = potential_liq
    potential = potential
    density_log = em_aq
    variable = ADEFieldAdvAux_em_aq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./ADDiffusiveFlux_em_aq]
    type = ADDiffusiveFlux
    density_log = em_aq
    variable = ADDiffusiveFlux_em_aq
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
    variable = em_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    #type = InterfaceLogDiffusionElectrons
    type = ADInterfaceLogDiffusion
    mean_en_neighbor = mean_en
    neighbor_var = em
    variable = em_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]

  # Henry coefficients (From Lietz, 2016)
  # H2O2  -  1.92e6
  # HO2   -  1.32e5
  # OH    -  6.2e2
  # H     -  6.48e3
  # H2    -  1.8e-2
  # O     -  1
  # O2    -  3.24e-2
  # O3    -  3e-1
  # N2    -  1.6e-2
  # N2O3  -  6.0e2
  # N2O4  -  3.69e1
  # N2O5  -  4.85e1
  # N2O   -  5.99e-1
  # HO2NO2 - 2.99e5
  # NO    -  4.4e-2
  # NO2   -  2.8e-1
  # NO3   -  4.15e1
  # HNO2, HNO  - 1.15e3
  # HNO3, ONOOH - 4.8e6
  # CO, CO(v)  - 2.42e-2
  # CO2, CO2(v)  - 8.23e-1
  # NH        - 1.47e3
#########################
#########################
  [H2_diff]
    type = InterfaceDiffusionTest
    variable = H2_aq
    neighbor_var = H2
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [H2_henry]
    type = InterfaceReactionTest
    variable = H2_aq
    neighbor_var = H2
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 1.8e-2 
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []

  [O3_diff]
    type = InterfaceDiffusionTest
    variable = O3_aq
    neighbor_var = O3
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [O3_henry]
    type = InterfaceReactionTest
    variable = O3_aq
    neighbor_var = O3
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 3e-1 
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
[]

[BCs]
  [./Nap_aq_physical]
    type = ADDCIonBC
    variable = Nap_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Clm_aq_physical]
    type = ADDCIonBC
    variable = Clm_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[H2_henry_test]
  #  type = MatchedValueLogBC
  #  variable = H2
  #  v = H2_aq
  #  H = 55
  #  boundary = 'gas_right'
  #[]
  #[H2_aq_open]
  #  type = ADDriftDiffusionOpenBC
  #  variable = H2_aq
  #  potential = potential
  #  position_units = ${dom1Scale}
  #  boundary = 'right'
  #[]
  [./Arex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ars
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex2_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arss
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex3_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arsss
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2ex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2s
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./Arex_henry]
  #  type = MatchedValueLogBC
  #  variable = Ar_aq
  #  v = Ars
  #  H = 0.1
  #  boundary = 'water_left'
  #[../]
  [./Ar2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2p
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Ar2p
    boundary = 'left gas_right'
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
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp'
    ip = 'Arp Ar2p'
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
    boundary = 'gas_right'
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
    boundary = 'gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'gas_right'
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
    ip = 'Arp Ar2p'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    ip = 'Arp Ar2p'
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
  [./em_aq_right]
    type = ADDCIonBC
    variable = em_aq
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./OHm_aq_physical]
    type = ADDCIonBC
    variable = OHm_aq
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H3Op_aq_physical]
    type = ADDCIonBC
    variable = H3Op_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O2m_aq_physical]
    type = ADDCIonBC
    variable = O2m_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Om_aq_physical]
    type = ADDCIonBC
    variable = Om_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./HO2m_aq_physical]
    type = ADDCIonBC
    variable = HO2m_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2Op_aq_physical]
    type = ADDCIonBC
    variable = H2Op_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O3m_aq_physical]
    type = ADDCIonBC
    variable = O3m_aq
    boundary = 'right'
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
[]

[Functions]
  [./water_ic_func]
    type = ParsedFunction
    value = 'log(8.6949e23/6.022e23)'
    # close to the boundary condition, essentially
  [../]
  [./potential_bc_func]
    type = ParsedFunction
    value = -0.8
    #value = 1.0
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
  # Salts!
  [./Nap_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Nap_aq
    heavy_species_mass = 3.816e-26
    heavy_species_charge = 1.0
    diffusivity = 2e-9
    block = 1
  [../]
  [./Clm_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Clm_aq
    heavy_species_mass = 5.887e-26
    heavy_species_charge = -1.0
    diffusivity = 2e-9
    block = 1
  [../]
  [./se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'left gas_right'
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
   property_tables_file = 'bolsig_files_01/electron_mobility_diffusion.txt'
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

  [./electron_data]
    type = ADGenericConstantMaterial
    prop_names = 'diffem_aq muem_aq Tem_aq'
    prop_values = '4.5e-9 0.000173913 300'
    block = 1
  [../]
  [./electron_sign]
    # I increased the electron mass by a factor of 10 here
    # Not sure what it's really supposed to be
    type = GenericConstantMaterial
    prop_names = 'sgnem_aq massem_aq'
    prop_values = '-1 9.11e-30'
    block = 1
  [../]

  [./N_A_mat]
    type = GenericConstantMaterial
    prop_names = 'N_A e diffpotential diffpotential_liq T_gas p_gas'
    prop_values = '6.022e23 1.602e-19 7.0832e-10 7.0832e-10 300 1.013e5'
    block = 1
  [../]

  ###################################
  # WATER SPECIES
  ###################################
  [./OHm_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OHm_aq
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]

  [./O_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O_aq
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./O3_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O3_aq
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = 0
    diffusivity = 5e-9
    block = 1
  [../]

  [./O3m_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O3m_aq
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./OH_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH_aq
    heavy_species_mass = 2.82431e-26
    heavy_species_charge = 0
    diffusivity = 4.5e-9
    block = 1
  [../]

  [./HO2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO2_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./HO3_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO3_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./H2O2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2O2_aq
    heavy_species_mass = 5.64840e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./H2Op_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2Op_aq
    heavy_species_mass = 2.992e-26 
    heavy_species_charge = 1
  [../]

  [./H2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2_aq
    heavy_species_mass = 3.34752e-26
    heavy_species_charge = 0
    diffusivity = 4.50e-9
    block = 1
  [../]

  [./H_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H_aq
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 0
    diffusivity = 5.0e-9
    block = 1
  [../]

  [./H+_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Hp_aq
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 1
  [../]

  [./H3O+_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H3Op_aq
    heavy_species_mass = 3e-26
    # ^ just estimated...whatever
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 1
  [../]

  [./HO2-_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO2m_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./Om_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Om_aq
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./O2m_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2m_aq
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = -1
    diffusivity = 5e-9 
    block = 1
  [../]

  [./O2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2_aq
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 2e-9
    block = 1
  [../]
[]

[Reactions]
  [./Argon]
    #species = 'em Arp Ar2p Ars Arss Arsss OH H2O O2 O2p O2m O2s Om O Os H2O H2Op H2 H Hs O3 HO2 H2O2 OHm OHp OHs Hp'
    species = 'em Arp Ar2p Ars Arss Arsss H2'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'bolsig_files_01'
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
                 em + Ars -> em + Arss                    : EEDF [-0.075] (C31_Ar1s5_Excitation_0.075_eV)
                 em + Ars -> em + Arss                    : EEDF [-0.279] (C32_Ar1s5_Excitation_0.28_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.525] (C33_Ar1s5_Excitation_1.52_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.655] (C34_Ar1s5_Excitation_1.66_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.678] (C35_Ar1s5_Excitation_1.68_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.727] (C36_Ar1s5_Excitation_1.73_eV)
                 em + Ars -> em + Arsss                   : EEDF [-1.75] (C37_Ar1s5_Excitation_1.75_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.45] (C20_Ar1s4_Excitation_1.45_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.58] (C21_Ar1s4_Excitation_1.58_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.60] (C22_Ar1s4_Excitation_1.60_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.65] (C23_Ar1s4_Excitation_1.65_eV)
                 em + Arss -> em + Arsss                  : EEDF [-1.68] (C24_Ar1s4_Excitation_1.68_eV)
                 #em + Arsss -> em + Arss                  : EEDF [] (C_Ar_Excitation__eV)
                 em + Ar -> em + em + Arp                 : EEDF [-15.76] (C16_Ar_Ionization_15.76_eV)
                 em + Ars -> em + em + Arp                : EEDF [-4.21] (C42_Ar1s5_Ionization_4.21_eV)
                 em + Arss -> em + em + Arp               : EEDF [-4.14] (C29_Ar1s4_Ionization_4.14_eV)
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
                 ##########################################
                 # Radiative Transitions Ar
                 ##########################################
                 Arsss -> Ars                           : 3.3e7
                 Arsss -> Ar                            : 3.1e5
                 Arss -> Ar                             : 5.3e5
                 Ar2s -> Ar + Ar                        : 6e7
                 ##########################################
                 # Additional reactions (from Van Gaens et al)
                 # Includes Ar and H reactions
                 ##########################################
                 #Arp + H2 -> Ar + H2Op          : 1.1e-9
                 # ^ Changed to the below reaction to conserve charge and still degrade H2
                 Arp + H2 -> Arp                : 1.1e-9'
  [../]


  # Rate coefficients are in m^3 s^-1
  # Taken from Wei Tian's thesis
  # Note the difference in values (1 L = 1000 m^-3)
  [water2]
    species = 'em_aq H3Op_aq H_aq H2O2_aq OH_aq OHm_aq O2_aq O2m_aq O_aq Om_aq H2_aq HO2_aq HO2m_aq H2Op_aq O3m_aq O3_aq HO3_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    reaction_coefficient_format = 'rate'
    block = 1

    reactions = 'em_aq + H2O_aq -> H_aq + OHm_aq            : 1.9e-2
                 #em_aq + H2Op_aq -> H_aq + OH_aq            : 6e-8 
                 em_aq + H2Op_aq -> H_aq + OH_aq            : 6e8 
                 # The units on this next one make no sense and have no consistency 
                 # across literature. This is the Buxton value (halved)
                 #em_aq + H_aq + H2O_aq -> H2_aq + OHm_aq    : 2.5e7
                 em_aq + OH_aq -> OHm_aq                    : 3e7
                 # where does this one come from???
                 #em_aq + em_aq + H2O_aq -> OHm_aq + OHm_aq  : 2.2e4
                 em_aq + H3Op_aq -> H_aq + H2O_aq                 : 2.3e7
                 em_aq + H2O2_aq -> OH_aq + OHm_aq                : 1.1e7 
                 #em_aq + HO2m_aq + H2O_aq -> OH_aq + OHm_aq + OHm_aq     : 3.5e3
                 em_aq + HO2m_aq -> OH_aq + OHm_aq + OHm_aq   : 3.5e6
                 em_aq + O2_aq -> O2m_aq                        : 1.9e7
                 # Next one is approximated by analogy (probably wrong...)
                 em_aq + O_aq -> Om_aq                          : 1.9e7
                 # This one is listed as 1e10 in Chens work. Completely different.
                 # I am going with this value because I have seen it in multiple references.
                 H_aq + OH_aq -> H2O_aq                        : 7e6
                 H_aq + OHm_aq -> em_aq + H2O_aq                  : 2.2e4
                 H_aq + H2O2_aq -> OH_aq + H2O_aq                 : 9e4
                 H_aq + O2_aq -> HO2_aq                        : 2.1e7
                 H_aq + HO2_aq -> H2O2_aq                      : 1e7
                 O_aq + H2O_aq -> OH_aq + OH_aq                   : 1.3e1
                 O_aq + O2_aq -> O3_aq                         : 3e6
                 OH_aq + OH_aq -> H2O2_aq    : 5.5e6
                 OH_aq + Om_aq -> HO2m_aq    : 2e7
                 OH_aq + OHm_aq -> Om_aq + H2O_aq   : 1.3e7
                 OH_aq + HO2_aq -> H2O_aq + O2_aq   : 6e6
                 OH_aq + O2m_aq -> OHm_aq + O2_aq     : 8e6
                 Om_aq + H2O_aq -> OHm_aq + OH_aq     : 1.8e3
                 Om_aq + H2O2_aq -> O2m_aq + H2O_aq   : 5e5
                 Om_aq + HO2m_aq -> O2m_aq + OHm_aq   : 4e5
                 Om_aq + O2_aq -> O3m_aq           : 3.6e6
                 #Om_aq + O2m_aq + H2O_aq -> OHm_aq + OHm_aq + O2_aq   : 6e2
                 Om_aq + O2m_aq -> OHm_aq + OHm_aq + O2_aq   : 6e5
                 OH_aq + H2O2_aq -> H2O_aq + HO2_aq     : 2.7e4
                 OH_aq + HO2m_aq -> OHm_aq + HO2_aq     : 7.5e6
                 # What the fuck
                 #H2Op_aq + H2O_aq -> OHm_aq + HO2_aq    : 6
                 H2Op_aq + H2O_aq -> H3Op_aq + OH_aq    : 6
                 H3Op_aq + OHm_aq -> H_aq + OH_aq + H2O_aq     : 6e7
                 HO2_aq + H2O_aq -> H3Op_aq + O2m_aq        : 2
                 H3Op_aq + O2m_aq -> HO2_aq + H2O_aq        : 6e-2
                 ##################################
                 # Additional reactions from Chen
                 #
                 # Some of these are from: 
                 # Elliot, A. John and McCracken, David R. "Computer modelling 
                 # of the radiolysis in an aqueous lithium salt blanket: 
                 # Suppression of radiolysis by addition of hydrogen." 
                 # Fusion Engineering and Design 13 (1990) 21-27
                 # doi: 10.1016/0920-3796(90)90028-5 
                 # 
                 # Note the reactions with H2O are often given with wrong
                 # (or confusing) rate coefficients. 
                 # e.g. k = 1.3e10 / [H2O] - 
                 # this means that the rate coefficient is essentially
                 # for a two body reaction since H2O is already included
                 #################################
                 O_aq + O_aq -> O2_aq          : 2.8e7
                 #em_aq + O2m_aq + H2O_aq -> HO2m_aq + OHm_aq   : 1.3e4
                 em_aq + O2m_aq -> HO2m_aq + OHm_aq   : 1.3e7
                 em_aq + HO2_aq -> HO2m_aq     : 2e7
                 #em_aq + Om_aq + H2O_aq -> OHm_aq + OHm_aq     : 2.2e4
                 # This one is listed with conflicting units in literature. 
                 # (Three body reaction with a two body reaction rate coefficient.)
                 em_aq + Om_aq -> OHm_aq + OHm_aq       : 2.2e7
                 #em_aq + O3m_aq + H2O_aq -> O2_aq + OHm_aq + OHm_aq   : 1.6e4 
                 em_aq + O3m_aq -> O2_aq + OHm_aq + OHm_aq : 1.6e7
                 em_aq + O3_aq -> O3m_aq     : 3.6e7
                 H_aq + Om_aq -> OHm_aq      : 1.1e7
                 H_aq + HO2m_aq -> OHm_aq + OH_aq   : 9e7 
                 H_aq + O3m_aq -> OHm_aq + O2_aq    : 1e7
                 H_aq + O2m_aq -> HO2m_aq        : 1.8e7
                 # Include HO3_aq or no?
                 H_aq + O3_aq -> HO3_aq          : 3.8e7 
                 OH_aq + O3m_aq ->  O3_aq + OHm_aq    : 2.6e6
                 #OH_aq + O3m_aq -> O2_aq + O2_aq + Hp     : 6e6
                 OH_aq + O3_aq -> HO2_aq + O2_aq          : 1.1e5
                 HO2_aq + O2m_aq -> HO2m_aq + O2_aq       : 8e4
                 HO2_aq + HO2_aq -> H2O2_aq + O2_aq       : 7e2
                 HO2_aq + Om_aq -> O2_aq + OHm_aq         : 6e6
                 HO2_aq + H2O2_aq -> OH_aq + O2_aq + H2O_aq    : 5e-4
                 HO2_aq + HO2m_aq -> OHm_aq + OH_aq + O2_aq    : 5e-4
                 O2m_aq + O2m_aq -> H2O2_aq + O2_aq + OHm_aq + OHm_aq : 1e-1
                 Om_aq + O2m_aq -> OHm_aq + OHm_aq + O2_aq       : 6e5
                 O2m_aq + H2O2_aq -> OH_aq + O2_aq + OHm_aq            : 1.3e-4
                 O2m_aq + HO2m_aq -> Om_aq + O2_aq + OHm_aq            : 1.3e-4
                 #O2m_aq + O3m_aq + H2O_aq -> OHm_aq + OHm_aq + O2_aq + O2_aq : 1e1 
                 O2m_aq + O3m_aq -> OHm_aq + OHm_aq + O2_aq + O2_aq   : 1e1
                 O2m_aq + O3_aq -> O3m_aq + O2_aq      : 1.5e6
                 Om_aq + Om_aq -> OHm_aq + HO2m_aq    : 1e6 
                 O3m_aq + H_aq -> O2_aq + OH_aq       : 9e7
                 O_aq + OHm_aq -> HO2m_aq          : 1.1e2
                 O_aq + H2O2_aq -> OH_aq + HO2_aq     : 1.6e5
                 O_aq + HO2m_aq -> OH_aq + O2m_aq     : 5.3e6
                 O3_aq + H2O2_aq -> OH_aq + HO2_aq + O2_aq   : 3e9
                 #HO3_aq -> O2_aq + OH_aq           : 1e2
                 O3m_aq + H3Op_aq -> O2_aq + OH_aq + H2O_aq      : 9e7
                 HO3_aq -> O2_aq + OH_aq          : 1.1e5
                 H2O2_aq -> OH_aq + OH_aq               : 4.4e-9
                 HO2m_aq -> Om_aq + OH_aq         : 1e-5'
  []

  [water3]
    name = 'test'
    species = 'em_aq H2_aq OHm_aq H_aq Om_aq OH_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = true
    reaction_coefficient_format = 'rate'
    block = 1

    reactions = 'em_aq + em_aq -> H2_aq + OHm_aq + OHm_aq   : 5.5e6 
                 em_aq + H_aq -> H2_aq + OHm_aq             : 2.5e7
                 H_aq + H2O_aq -> H2_aq + OH_aq                   : 1e-2
                 H_aq + H_aq -> H2_aq                          : 7.5e6
                 H2_aq + H2O2_aq -> H_aq + OH_aq + H2O_aq            : 6e3
                 OH_aq + H2_aq -> H_aq + H2O_aq     : 4.2e4
                 Om_aq + H2_aq -> OHm_aq + H_aq       : 8e7'
  []

  # More comprehensive reaction network is available in "Air plasma treatment
  # of liquid covered tissue: long timescale chemistry" 
  # Lietz et al. J. Phys. D: Appl. Phys. 49 425204 (2016)
  #[water_nitrogen]
  #  species = 'OH NO2m OHm NO2 H2O2 OONm H 
  #  aux_species = 'H2O'
  #  reactions = 'OH + NO2m -> OHm + NO2     : 1e7
  #               H2O2 = NO2m -> OONOm + H2O   : 4.5e5
  #               H + NO2m -> OHm + NO         : 1.2e6
  #               Om + NO2m + H2O -> OHm + OHm + NO2     : 3.6e5
  #               NO + NO + O2 -> NO2 + NO2              : 2.3e3
  #               NO + NO2 + H2O -> HNO2 + HNO2          : 2e5
  #               H + HNO2 -> NO + H2O         : 4.5e5
  #               NO + OH -> HNO2        : 2e7
  #               NO + H -> HNO2       : 1e7
  #               HNO3 + OH -> NO3 + H2O   : 1.2e5
  #               NO + HO2 -> HNO3     : 8e6
  #               NO2 + OH -> HNO3     : 3e7
  #               O2m + NO -> NO3m     : 1.6e7
  #               NO + HO2 -> HOONO      : 3.2e6
  #               NO2 + OH -> HOONO      : 1.2e7
  #               O2m + NO -> OONOm      : 6.6e6
  #               HNO2 + H2O -> H3Op + NO2m    : 1.8e-2
  #               H3Op + NO2m -> HNO2 + H2O    : 2
  #               H3Op + NO3m -> HNO3 + H2O    : 2e-1
  #               # The next three have strange units. Why M^-2 (m^6) when there
  #               # are only two reactants?
  #               # I only reduced by a factor of 1e-3 instead of 1e-6 to match
  #               # the reactants. 
  #               #N2O3 + H2O -> HNO2 + HNO2      : 1.1e-1
  #               #N2O4 + H2O -> HNO2 + HNO3      : 8e-1
  #               #N2O5 + H2O -> HNO3 + HNO3      : 1.2e-3
  #               NO2 + NO2 + H2O -> HNO2 + H3Op + NO3m    : 1.5e2
  #               NO2 + NO2 + H2O -> H3Op + NO2m + H3Op + NO3m   : 5e1'
  #[]
[]
