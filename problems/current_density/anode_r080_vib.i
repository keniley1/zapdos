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
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
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
  [out_01]
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
    #ballast_resist = 651e3   # This is the experimental value
    #ballast_resist = 9.0e5
    ballast_resist = 8e5
    e = 1.6e-19
  [../]
[]

[DriftDiffusionAction]
  [./Plasma]
    electrons = em
    charged_particle = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
    Neutrals = 'Ars Arss Arsss Ar2s H2O OH O2 O2s O Os H2 H Hs O3 HO2 H2O2 OHs H2Ov'
    mean_energy = mean_en
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    order = FIRST 
    position_units = ${dom0Scale}
    block = 0
  [../]

  [./Water]
    charged_particle = 'H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
    Neutrals = 'H_aq H2O2_aq OH_aq O2_aq O_aq H2_aq HO2_aq O3_aq HO3_aq'
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    order = FIRST 
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./Salt]
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
  [H2Ov]
    block = 0
    initial_condition = -30
  []
  [H2O]
    block = 0 
    #initial_condition = -1.26279318 # 20 % humidity 
    #initial_condition = -0.70316062 # 35 % humidity 
    #initial_condition = -0.34647896 # 50 % humidity 
    #initial_condition = -0.08411109 # 65 % humidity 
    initial_condition = 0.12351586 # 80 % humidity 
  []

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
  [H2Op]
    block = 0
    initial_condition = -25
  []
  [Hp]
    block = 0
    initial_condition = -30
  []
  [Op]
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


  # Water Species
  #[./em_aq]
  #  block = 1
  #  #initial_condition = -24
  #  #initial_condition = -21
  #  # scaling = 1e-5
  #  #initial_condition = -14
  #  #initial_condition = -24
  #  #initial_condition = -20
  #  initial_condition = -16
  #[../]

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
  [Tgas]
    #order = CONSTANT
    #family = MONOMIAL
    #order = LINEAR
    family = LAGRANGE
    block = 0
  []
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL 
    #initial_condition = 3.69781886 # 20 % humidity 
    #initial_condition = 3.69254831 # 35 % humidity 
    #initial_condition = 3.68724983 # 50 % humidity 
    #initial_condition = 3.68192314 # 65 % humidity 
    initial_condition = 3.67656834 # 80 % humidity 
  [../]
  [Ar_density]
    block = 0
    order = CONSTANT
    family = MONOMIAL
  []

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
  #[./ADCurrent_em_aq]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #  initial_condition = 0
  #[../]
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
  #[./ADEFieldAdvAux_em_aq]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #  initial_condition = 0
  #[../]
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
  [Tgas_interpolate]
    type = InterpolateAlongAxis
    variable = Tgas
    axis = 0
    property_file = 'T_gas.csv'
    block = 0
  []
  [./H2O_aq_density]
    type = DensityMoles
    variable = H2O_aq_density
    density_log = H2O_aq
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [Ar_density_calc]
    type = DensityMoles
    variable = Ar_density
    density_log = Ar
    execute_on = 'initial'
    block = 0
  []
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
    charged = 'em Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  []
  [rholiq_calc]
    type = ChargeDensity
    variable = rholiq
    charged = 'H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
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
  #[./ADCurrent_em_aq]
  #  type = ADCurrent
  #  #potential = potential_liq
  #  potential = potential
  #  density_log = em_aq
  #  variable = ADCurrent_em_aq
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom1Scale}
  #[../]
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
[]

[InterfaceKernels]
  #[H2_g_to_l]
  #  type = ADHenryDiffusion
  #  variable = H2_aq
  #  neighbor_var = H2 
  #  h = 1.8e-2
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
  #[H2_l_to_g]
  #  type = ADHenryDiffusion
  #  variable = H2
  #  neighbor_var = H2_aq 
  #  h = 55.5556
  #  position_units = ${dom0Scale}
  #  neighbor_position_units = ${dom1Scale}
  #  boundary = 'gas_right'
  #[]
  #[H2_henry]
  #  type = ADHenryInterface
  #  variable = H2_aq
  #  neighbor_var = H2 
  #  h = 1.8e-2
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
  #[H2_henry_liquid_to_gas]
  #  type = ADHenryInterface
  #  variable = H2
  #  neighbor_var = H2_aq 
  #  h = 55
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
  #[./em_advection]
  #  type = ADInterfaceAdvection
  #  mean_en_neighbor = mean_en
  #  potential_neighbor = potential
  #  neighbor_var = em
  #  variable = em_aq
  #  boundary = water_left
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #[../]
  #[./em_diffusion]
  #  #type = InterfaceLogDiffusionElectrons
  #  type = ADInterfaceLogDiffusion
  #  mean_en_neighbor = mean_en
  #  neighbor_var = em
  #  variable = em_aq
  #  boundary = water_left
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #[../]
  [./H2Op_advection]
    type = ADInterfaceAdvection
    potential_neighbor = potential
    neighbor_var = H2Op
    variable = H2Op_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./H2Op_diffusion]
    type = ADInterfaceLogDiffusion
    neighbor_var = H2Op
    variable = H2Op_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./Arp_advection]
    type = ADInterfaceAdvection
    potential_neighbor = potential
    neighbor_var = Arp
    variable = H2Op_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion]
    type = ADInterfaceLogDiffusion
    neighbor_var = Arp
    variable = H2Op_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./Ar2p_advection]
    type = ADInterfaceAdvection
    potential_neighbor = potential
    neighbor_var = Ar2p
    variable = H2Op_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./Ar2p_diffusion]
    type = ADInterfaceLogDiffusion
    neighbor_var = Ar2p
    variable = H2Op_aq
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

  [OH_diff]
    type = InterfaceDiffusionTest
    variable = OH_aq
    neighbor_var = OH
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [OH_henry]
    type = InterfaceReactionTest
    variable = OH_aq
    neighbor_var = OH
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 6.2e2 
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []

  [H2O2_diff]
    type = InterfaceDiffusionTest
    variable = H2O2_aq
    neighbor_var = H2O2
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [H2O2_henry]
    type = InterfaceReactionTest
    variable = H2O2_aq
    neighbor_var = H2O2
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 1.92e6 
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []

  [HO2_diff]
    type = InterfaceDiffusionTest
    variable = HO2_aq
    neighbor_var = HO2
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [HO2_henry]
    type = InterfaceReactionTest
    variable = HO2_aq
    neighbor_var = HO2
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 1.32e5 
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []


  [H_diff]
    type = InterfaceDiffusionTest
    variable = H_aq
    neighbor_var = H
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [H_henry]
    type = InterfaceReactionTest
    variable = H_aq
    neighbor_var = H
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []

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

  [O_diff]
    type = InterfaceDiffusionTest
    variable = O_aq
    neighbor_var = O
    h = 6.48e3
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []
  [O_henry]
    type = InterfaceReactionTest
    variable = O_aq
    neighbor_var = O
    #kf = 6.48e3
    #kb = 1
    kf = 1
    kb = 1 
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
#########################
#########################
[]

[BCs]
  # H2O evaporation boundary condition
  #[H2O_interface]
  #  type = DirichletBC
  #  variable = H2O
  #  #value = -1.26279318 # 20 % humidity 
  #  #value = -0.70316062 # 35 % humidity 
  #  #value = -0.34647896 # 50 % humidity 
  #  #value = -0.08411109 # 65 % humidity 
  #  value = 0.12351586 # 80 % humidity 
  #  boundary = 'gas_right'
  #[]
  [H2Ov_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2Ov
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []
  [vapor_pressure]
    type = VaporValueBC
    variable = H2O
    gas_temperature = Tgas
    boundary = 'gas_right'
  []

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
  [H2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [OH_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OH
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [H2O2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2O2
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [HO2_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = HO2
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [H_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [O_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  [O3_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O3
    #boundary = 'left gas_right'
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  []
  # H2O+ boundary conditions
  [./H2Op_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = H2Op 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./H2Op_physical_right_advection]
    type = ADHagelaarIonAdvectionBC
    variable = H2Op
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  # O- and OHm_aq boundary conditions
  [./Om_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Om 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Om_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Om
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHm_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHm 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHm_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = OHm
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHp_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHp 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./OHp_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = OHp
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2m_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2m 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2m_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = O2m
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2p 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./O2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = O2p
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./Hp_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Hp 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Hp_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Hp
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Op_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Op 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Op_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Op
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]

  #[OH_left]
  #  type = DirichletBC
  #  variable = OH
  #  value = -10
  #  boundary = 'left'
  #[]
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
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []
  [O2s_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = O2s
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Os_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Os
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Hs_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Hs
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []
  [OHs_bc_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = OHs 
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  []

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
  #[./potential_left]
  #  type = ADNeumannCircuitVoltageMoles_KV
  #  variable = potential
  #  boundary = left
  #  function = potential_bc_func
  #  #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
  #  ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
  #  data_provider = data_provider
  #  em = em
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [potential_left]
    type = CircuitBC
    variable = potential
    boundary = left
    function = potential_bc_func
    data_provider = data_provider
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
    em = em
    mean_en = mean_en
    r_ion = 0
    r_electron = 0
    position_units = ${dom0Scale}
  []
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
    boundary = 'gas_right'
    potential = potential
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'gas_right'
    potential = potential
    #ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p'
    ip = 'Arp Ar2p H2Op OHp OHm Om O2m O2p Hp Op'
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

  # Water BCs
  #  charged_particle = 'H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
  [./H3Op_aq_right]
    type = ADDCIonBC
    variable = H3Op_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2Op_aq_right]
    type = ADDCIonBC
    variable = H2Op_aq
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]

  #[./em_aq_right]
  #  type = ADDCIonBC
  #  variable = em_aq
  #  boundary = 'right'
  #  #potential = potential_liq
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./OHm_aq_physical]
  #  type = ADDCIonBC
  #  variable = OHm_aq
  #  boundary = 'right'
  #  #potential = potential_liq
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./H3Op_aq_physical]
  #  type = ADDCIonBC
  #  variable = H3Op_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./O2m_aq_physical]
  #  type = ADDCIonBC
  #  variable = O2m_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./Om_aq_physical]
  #  type = ADDCIonBC
  #  variable = Om_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./HO2m_aq_physical]
  #  type = ADDCIonBC
  #  variable = HO2m_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./H2Op_aq_physical]
  #  type = ADDCIonBC
  #  variable = H2Op_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./O3m_aq_physical]
  #  type = ADDCIonBC
  #  variable = O3m_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
[]

[ICs]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
    #block = 0
  [../]
  #[H2O_ic]
  #  type = FunctionIC
  #  variable = H2O
  #  function = water_ic_func
  #[]
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
    value = 1.65
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  #[./em_aq_ic_func]
  #  type = ParsedFunction
  #  value = '1778 - 1.8e6*x'
  #[../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '1.65 * (1.001e-3 - x)'
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
  [H2Ov_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2Ov
    heavy_species_mass = 2.9907e-26
    heavy_species_charge = 0
    diffusivity = 2.3e-5
  []
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
  [./se_coefficient_left]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0'
    boundary = 'left'
  [../]
  [./se_coefficient_right]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'gas_right'
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
   property_tables_file = 'bolsig_files_080/electron_mobility_diffusion.txt'
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

  #[./electron_data]
  #  type = ADGenericConstantMaterial
  #  prop_names = 'diffem_aq muem_aq Tem_aq'
  #  prop_values = '4.5e-9 0.000173913 300'
  #  block = 1
  #[../]
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
    #species = 'em Arp Ar2p Ars Arss Arsss OH H2O O2 O2p O2m O2s Om O Os H2O H2Op H2 H Hs O3 HO2 H2O2 OHm OHp OHs'
    species = 'em Arp Ar2p Ars Arss Arsss OH H2O O2 O2p O2m O2s Om O Os H2O H2Op H2 H Hs O3 HO2 H2O2 OHm OHp OHs Hp Op H2Ov'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'bolsig_files_080'
    #equation_constants = 'Tgas Tn'
    #equation_values = '300 1'
    #equation_variables = 'Te'
    equation_constants = 'Tn'
    equation_values = '1'
    equation_variables = 'Tgas Te'
    potential = 'potential'
    use_log = true
    convert_to_moles = true
    convert_to_meters = 1e-2
    position_units = ${dom0Scale}
    use_ad = true
    block = 0

    # NOTE: Most reactions are taken from Van Gaens, specifically from the B reaction list.
    # Some additional ones are included, specifically those including H2 and H2O (since H2 
    # is produced in more significant quantities in this work due to electrolysis, and 
    # H2O is introduced in significant quantities through humidity)

    # Note that the monatomic hydrogen reactions were taken from 
    # Van Gaens et al. "Kinetic modeling for an atmospheric pressure argon 
    # plasma jet in humid air." J. Phys. D: Appl. Phys. 46 275201 (2013)
    #
    # Ars    - Ar(1s1), Ar(1s3) (metastable states of Ar(1s) manifold)
    # Arss   - Ar(1s2), Ar(1s4) (radiative states of Ar(1s) manifold)
    # Arsss  - Ar(4p) and higher states
    # Os: 
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
                 # Note that some rate coefficients are multiplied by a factor of 2
                 # since Ars and Arss account for multiple excited states
                 Ars + Ar2s -> Arp + Ar + Ar + em         : 5e-10
                 Ars + Ars -> Ar + Arp + em               : 1.28e-9
                 Arss + Arss -> Ar + Arp + em             : 4.2e-9 
                 Arsss + Arsss -> Ar + Arp + em           : 1e-9 
                 #Arp + Ar -> Ar + Arp                     : {5.66e-10*Tn^0.5}
                 # Van gaens reactions 309-320
                 Ars + Ar + Ar -> Ar + Ar2s               : 6.6e-32 
                 Arss + Ar + Ar -> Ar + Ar2s              : 6.6e-32 
                 Arsss + Ar + Ar -> Ar + Ar2s             : 1e-31
                 Arp + Ar + Ar -> Ar + Ar2p               : 2.5e-31 
                 Ar2s + Ar2s -> Ar2p + Ar + Ar + em       : 5e-10
                 Ar2s + Ar -> Ar + Ar + Ar                : 5e-15
                 # Van Gaens reactions 317-320 
                 Ar2s + Ar -> Ars + Ar + Ar               : 1e-14 
                 Ar2s + Ar -> Arss + Ar + Ar              : 1e-14 
                 #########################
                 # H2O electron-impact reactions
                 #########################
                 em + H2O -> em + H2O             : EEDF [elastic] (C60_H2O_Elastic)
                 em + H2O -> H2Ov + em            : EEDF [-0.198] (C65_H2O_Excitation_0.20_eV)
                 H2Ov + Ar -> H2O + Ar            : 1e-13
                 H2Ov + H2O -> H2O + H2O          : 1e-13
                 H2Ov + O2 -> H2O + O2            : 1e-13
                 em + H2O -> OHm + H              : EEDF (C58_H2O_Attachment)
                 em + H2O -> em + OHs + H         : EEDF [-9.00] (C66_H2O_Excitation_9.00_eV) 
                 em + H2O -> em + em + H2Op       : EEDF [-13.5] (C67_H2O_Ionization_13.50_eV)
                 em + H -> em + Hs                : EEDF [-10.2] (C69_H_Excitation_10.20_eV)
                 em + H -> em + Hs                : EEDF [-10.2] (C70_H_Excitation_10.20_eV)
                 em + Hs -> em + H                : EEDF [10.2] (C72_H2p10.2eV_De-excitation_10.20_eV)
                 em + Hs -> em + H                : EEDF [10.2] (C73_H2s10.2eV_De-excitation_10.20_eV)
                 em + O2 -> em + O2               : EEDF [elastic] (C44_O2_Elastic)
                 em + O2 -> O2s + em              : EEDF [-0.98] (C46_O2_Excitation_0.98_eV)
                 em + O2 -> O2s + em              : EEDF [-1.63] (C47_O2_Excitation_1.63_eV)
                 em + O2 -> em + em + O2p         : EEDF [-12.1] (C48_O2_Ionization_12.10_eV)
                 em + O2 -> em + em + O + Op      : EEDF [-19.5] (C49_O2_Ionization_19.50_eV)
                 em + O2s -> em + O2              : EEDF [0.98] (C50_O20.98_De-excitation_0.98_eV)
                 em + O2s -> em + O2              : EEDF [1.63] (C54_O21.63_De-excitation_1.63_eV)
                 em + O2s -> em + em + O2p        : EEDF [-11.10] (C52_O20.98_Ionization_11.10_eV)
                 em + O2s -> em + em + O2p        : EEDF [-10.45] (C56_O21.63_Ionization_10.45_eV)
                 em + O2s -> em + em + Op + O     : EEDF [-18.52] (C53_O20.98_Ionization_18.52_eV)
                 em + O2s -> em + em + Op + O     : EEDF [-17.87] (C57_O21.63_Ionization_17.87_eV)
                 em + O2 -> Om + O                : EEDF (C43_O2_Attachment)
                 em + O2s -> Om + O               : EEDF (C51_O20.98_Attachment)
                 em + O2s -> Om + O               : EEDF (C55_O21.63_Attachment)
                 em + O2p -> O + O                : {1.2e-8*Te^(-0.7)} 
                 em + H -> em + em + Hp           : EEDF [-13.60] (C71_H_Ionization_13.60_eV)
                 # Van Gaens reactions 266-268 -- seem important!
                 em + H2Op -> O + H2              : {6.27e-9*Te^(-0.5)}
                 em + H2Op -> O + H + H           : {4.90e-8*Te^(-0.5)}
                 em + H2Op -> OH + H              : {1.38e-8*Te^(-0.5)}
                 # H2O2 dissociative attachment
                 # Cross sections originally from D. Nandi, Chem. Phys. Lett. 373 (2003), \cite{Nandi2003}
                 # But the form used here is from 
                 # Liu, D. PSST 19 (2010), \cite{Liu2010a}
                 # (reactions 83-84)
                 em + H2O2 -> H2O + Om            : {1.57e-10*Te^(-0.55)}
                 em + H2O2 -> OH + OHm            : {2.7e-10*Te^(-0.5)}  
                 #########################
                 # Argon-heavy particle collisions
                 #########################
                 Ars + O2 -> Ar + O + O                   : 4.2e-10
                 Ars + O2s -> Ar + O + O                  : 2.1e-10
                 Arss + O2 -> Ar + O + O                  : 4.2e-10
                 Arsss + O2 -> Ar + O + O                 : 2.1e-10
                 Arsss + O2s -> Ar + O + O                : 2.1e-10
                 Arp + O -> Ar + Op                       : 6.4e-12
                 Arp + O2 -> Ar + O2p                     : 5e-11
                 Arp + O2m -> Ar + O2                     : {2e-7*(Tgas/300)^(-0.5)}
                 Arp + O2m -> Ar + O + O                  : 1e-7
                 Arp + Om -> Ar + O                       : {2e-7*(Tgas/300)^(-0.5)}
                 # Van Gaens, Reactions 427, 428, 433-437
                 Ar2s + O2 -> Ar + Ar + O + O             : 4.6e-11
                 Ar2s + O3 -> Ar + Ar + O2 + O            : 2.1e-10
                 Ar2p + Om -> Ar + Ar + O                 : 1e-7
                 Ar2p + O2 -> Ar + Ar + O2p               : 1e-10
                 Ar2p + O2m -> Ar + Ar + O + O            : 1e-7
                 Ar2p + O2m -> Ar + Ar + O2               : 1e-7
                 Ar2p + O3 -> Ar + Ar + O2p + O           : 5e-9
                 #########################
                 # Oxygen-heavy particle collisions
                 #########################
                 O + O + Ar -> O2 + Ar                    : {5.21e-35*exp(900/Tgas)}
                 O + O + O2 -> O2 + O2                    : {5.21e-35*exp(900/Tgas)}
                 O + O + H2 -> O2 + H2                    : {1.3e-34*exp(900/Tgas)}
                 O + O + H2O -> O2 + H2O                  : {2.6e-34*exp(900/Tgas)}
                 # Reaction 451-452
                 O + O + Ar -> O2s + Ar                   : {6.93e-35*(Tgas/300)^(-0.63)}
                 # Reaction 454
                 O + Op + Ar -> O2p + Ar                  : 1e-29
                 O + Om -> O2 + em                        : 1.5e-10 
                 O + O2 + Ar -> O3 + Ar                   : {6.4e-35*exp(663/Tgas)} 
                 # Reaction 459-460 - included for H2 and H2O
                 O + O2 + H2 -> O3 + H2                   : {1.3e-34*exp(663/Tgas)}
                 O + O2 + H2O -> O3 + H2O                 : {3.2e-34*exp(663/Tgas)} 
                 # oxygen charge exchange reactions 
                 O + O2m -> O2 + Om                       : {1.5e-10*(Tgas/300)^(0.5)}
                 O + O2m -> O3 + em                       : 1.5e-10
                 O + O3 -> O2 + O2                        : {8e-12*(-2060/Tgas)}
                 # Ignored because O3m not included
                 #O + O3m -> O2m + O2                      : 3.2e-10
                 #O + O3m -> O2 + O2 + em                  : 3e-10
                 Os + O2 -> O2s + O                       : {2.6e-11*(67/Tgas)}
                 Os + O3 -> O2s + O2                      : 1.2e-10
                 Os + Ar -> Ar + O                        : 5e-12
                 Op + O2 -> O + O2p                       : 2.1e-11
                 Om + O2 -> O2m + O                       : 1e-10
                 Om + O2s -> O3 + em                      : 3e-10
                 Om + O2p + Ar -> O + O2 + Ar             : {2e-25*(Tgas/300)^(-2.5)}

                 # Next is commented because O4p not tracked
                 #O2 + O2p + Ar -> O4p + Ar                : {2.4e-30*(Tgas/300)^(-3.2)}
                 # Rxns 505, 508, 512, 515 ignored because vibrational levels of O2 are not tracked
                 O2s + O2s -> O2s + O2                    : {1.8e-19*(Tgas/300)^(3.8)*exp(700/Tgas)}
                 O2s + O2m -> O2 + O2 + em                : 2e-10
                 O2s + O3 -> O2 + O2 + O                  : 1e-14
                 O2s + Ar -> O2 + Ar                      : {3e-18*exp(-200/Tgas)}
                 # Next included for H2 
                 O2s + H2 -> O2 + H2                      : {6e-18*exp(-200/Tgas)}
                 O2s + H2O -> O2 + H2O                    : 2.3e-12
                 O2s + O3 -> O2 + O3                      : 2.2e-11
                 # Reactions 546-550 ignored because O2s levels are lumped together
                 O2p + O2m -> O2 + O + O                  : 1e-7
                 O2p + O2m + Ar -> O2 + O2 + Ar           : {1e-24*(Tgas/300)^(-2.5)}
                 # Ignored because O3m not included
                 #O2p + O3m + Ar -> O3 + O2 + Ar           : {2e-25*(Tgas/300)^(-2.5)}
                 #O2m + O3 -> O2 + O3m                     : 4e-10
                 O2m + Ar -> O2 + em + Ar                 : {2.7e-10*(Tgas/300)^(0.5)*exp(-5590/Tgas)}
                 O3 + Ar -> O + O2 + Ar                   : {7.3e-10*exp(-11400/Tgas)}
                 O3 + H2 -> O + O2 + H                    : {1.4e-9*exp(-11400/Tgas)}
                 O3 + H2O -> O + O2 + H2O                 : {4e-9*exp(-11400/Tgas)}
                 #########################
                 # Argon-humid air reactions
                 # All of the B reactions from Van Gaens, plus all H2 reactions
                 #########################
                 Ars + H2 -> Ar + H + H                   : 1.32e-10
                 Ars + OH -> Ar + O + H                   : 6.6e-11
                 Ars + OH -> Ar + OHs                     : 6.6e-11
                 Ars + H2O -> Ar + OH + H                 : 2.1e-10
                 Arss + H2 -> Ar + H + H                  : 1.32e-10 
                 Arss + OH -> Ar + O + H                  : 6.6e-11
                 Arss + OH -> Ar + OHs                    : 6.6e-11
                 Arss + H2O -> Ar + OH + H                : 2.1e-10
                 Arss + H2 -> Ar + H + H                  : 6.6e-11
                 Arsss + H2 -> Ar + H + H                 : 4e-10
                 Arsss + OH -> Ar + O + H                 : 2.1e-10
                 Arsss + H2O -> Ar + OH + H               : 2.1e-10
                 Arp + H -> Ar + Hp                       : 1e-10
                 #Arp + H2 -> Ar + H2Op                    : 1.1e-9
                 Arp + H2O -> Ar + H2Op                   : 7e-10
                 Arp + OHm -> Ar + OH                     : {2e-7*(Tgas/300)^(-0.5)}
                 # Next ones excluded because ArHp is not tracked
                 # I should look into ArH molecules -- significant or not?
                 #Arp + H2O -> ArHp + OH                   : 3e-10
                 #Ar2p + H2O -> ArHp + Ar + OH             : 4e-10
                 Ar2s + H2 -> Ar + Ar + H + H             : 6.6e-11
                 Ar2s + H2O -> Ar + Ar + OH + H           : 1e-10
                 Ar2p + H -> Ar + Ar + Hp                 : 5e-11
                 Ar2p + H2O -> Ar + Ar + H2Op             : 1.6e-9
                 Ar2p + OHm -> Ar + Ar + OH               : 1e-7
                 Ar2p + OHm -> Ar + Ar + O + H            : 1e-7
                 #########################
                 # Humid air-heavy particle collisions
                 #########################
                 H + O + Ar -> OH + Ar                    : {1.5e-32*(300/Tgas)}
                 H + O + O2 -> OH + O2                    : {1.5e-32*(300/Tgas)}
                 H + O + H2 -> OH + H2                    : {3e-32*(300/Tgas)}
                 H + O + H2O -> OH + H2O                  : {7.5e-32*(300/Tgas)}
                 H + Op -> Hp + O                         : {5.66e-10*(Tgas/300)^(0.36)*exp(8.6/Tgas)}
                 H + O2 + Ar -> HO2 + Ar                  : {6.09e-32*(Tgas/300)^(-0.8)}
                 H + O2 + O2 -> HO2 + O2                  : {6.09e-32*(Tgas/300)^(-0.8)}
                 H + O2 + H2O -> HO2 + H2O                : {3.05e-31*(Tgas/300)^(-0.8)}
                 H + O2 + H2 -> HO2 + H2                  : {1.53e-31*(Tgas/300)^(-0.8)}
                 H + O2s -> OH + O                        : {1.83e-13*exp(-1550/Tgas)}
                 # Next excluded because Hm is not tracked
                 #H + O2m -> Hm + O2                       : 7e-10
                 H + O2m -> HO2 + em                      : 7e-10
                 H + O3 -> OH + O2                        : {2.71e-11*(Tgas/300)^0.75}
                 H + O3 -> HO2 + O                        : 7.51e-13
                 H + H + Ar -> H2 + Ar                    : {2e-32*(Tgas/300)^(-1)}
                 H + H + H2 -> H2 + H2                    : {4e-32*(Tgas/300)^(-1)}
                 H + H + H2O -> H2 + H2O                  : {9.2e-32*(Tgas/300)^(-1)}
                 H + OH + Ar -> H2O + Ar                  : {8e-31*(Tgas/300)^(-2.6)}
                 H + OH + O3 -> H2O + O3                  : {1.8e-30*(Tgas/300)^(-2.6)}
                 H + OHm -> H2O + em                      : 1.8e-9
                 H + HO2 -> H2 + O2                       : {2.06e-11*(Tgas/300)^0.84*exp(-277/Tgas)}
                 H + HO2 -> OH + OH                       : {1.66e-10*exp(-413/Tgas)}
                 H + HO2 -> H2O + O                       : {5e-11*exp(-866/Tgas)}
                 H + HO2 -> H2O + Os                      : {2.32e-12*(Tgas/300)^1.55*exp(80.85/Tgas)}
                 H + H2O -> OH + H + H                    : {5.8e-9*exp(-52900/Tgas)}
                 H + H2O -> OH + H2                       : {6.89e-12*exp(-9720/Tgas)}
                 H + H2O2 -> HO2 + H2                     : {8e-11*exp(-4000/Tgas)}
                 H + H2O2 -> H2O + OH                     : {4e-11*exp(-2000/Tgas)}
                 # Skipping nitrogen, moving to rxn number 1104
                 Hs + O2 -> H + O2                        : 1e-13
                 Hs + O3 -> H + O3                        : 1e-13
                 Hs + O -> H + O                          : 1e-13
                 Hs + H2O -> H + H2O                      : 1e-13
                 Hs + H2 -> H + H2                        : 1e-13
                 Hs + Ar -> H + Ar                        : 1e-13
                 Hp + O -> H + Op                         : {3.04e-10*(Tgas/300)^0.47*exp(-11.5/Tgas)}
                 Hp + Om -> H + O                         : {2e-7*(Tgas/300)^(-0.5)}
                 Hp + Om + Ar -> H + O + Ar               : {2e-25*(Tgas/300)^(-2.5)}
                 Hp + Om + Ar -> OH + Ar                  : {2e-25*(Tgas/300)^(-2.5)}
                 Hp + O2 -> H + O2p                       : 2e-9
                 Hp + O2m -> H + O2                       : {2e-7*(Tgas/300)^(-0.5)}
                 # skipping rxns 1116-1121 because O2m and O3m are not terribly relevant
                 # Also skipping up through 1133 since we have no nitrogen
                 # next three are all rxn 1136
                 # commented out because of erratum - H2Op should be H2p, which I ignore
                 #Hp + H + Ar -> H2Op + Ar                 : 1e-34
                 #Hp + H + H2 -> H2Op + H2                 : 1e-34
                 #Hp + H + H2O -> H2Op + H2O               : 1e-34
                 Hp + OH -> H + OHp                       : 2.1e-9
                 Hp + OHm -> OH + H                       : {2e-7*(Tgas/300)^(-0.5)}
                 Hp + OHm -> H + H + O                    : 1e-7
                 # skipping 1141-1142
                 Hp + H2O -> H + H2Op                     : 6.9e-9
                 # Skipping up through 1193 
                 H2 + O -> OH + H                         : {1.6e-11*exp(-4570/Tgas)}
                 H2 + Op -> OHp + H                       : 1.62e-9
                 H2 + Om -> H2O + em                      : 6.72e-10
                 H2 + Om -> OHm + H                       : 2.8e-11
                 H2 + O2 -> H + HO2                       : {2.4e-10*exp(-28500/Tgas)}
                 H2 + O2 -> O + O + H2                    : {6e-9*exp(-52300/Tgas)}
                 H2 + O2 -> OH + OH                       : {3.16e-10*exp(-21890/Tgas)}
                 H2 + O2m -> OHm + OH                     : 5e-13
                 # T too small
                 #H2 + H2Op -> H2 + Hp + H                 : {1e-8*exp(-84100/Tgas)} 
                 H2 + OH -> H2O + H                       : {9.54e-13*exp(-1490/Tgas)}
                 H2 + OHp -> H2Op + H                     : 1.3e-9
                 # Next is excluded because Hm is not tracked
                 #H2 + OHm -> H2O + Hm                     : 5e-12
                 ##############################
                 ##############################
                 H2 + HO2 -> H2O2 + H                     : {4.38e-12*exp(-10751/Tgas)}
                 ##############################
                 ##############################
                 H2 + H2O -> OH + H2 + H                  : {4.38e-12*exp(-52900/Tgas)}
                 #H2Op + O -> H + OHp                      : 1.5e-9
                 #H2Op + Om -> H2 + O                      : {2e-7*(Tgas/300)^(-0.5)}
                 OH + O -> H + O2                         : {1.81e-11*(Tgas/300)^(-0.31)*exp(177/Tgas)}
                 OH + Os -> H + O2                        : {2.08e-11*(Tgas/300)^(-0.186)*exp(-154/Tgas)}
                 OH + Op -> OHp + O                       : 3.6e-10
                 OH + Op -> O2p + H                       : 3.6e-10
                 OH + O2m -> OHm + O2                     : 1e-10
                 OH + O3 -> HO2 + O2                      : 1.69e-12
                 OH + OH -> H2O + O                       : {5.49e-14*(Tgas/300)^(2.42)*exp(970/Tgas)}
                 ##############################
                 ##############################
                 OH + OH + Ar -> H2O2 + Ar                : {8e-31*(Tgas/300)^(-0.8)}
                 OH + OH + O2 -> H2O2 + O2                : {8e-31*(Tgas/300)^(-0.9)} 
                 OH + OH + H2 -> H2O2 + H2                : {1.6e-30*(Tgas/300)^(-0.9)} 
                 OH + OH + H2O -> H2O2 + H2O              : {4e-30*(Tgas/300)^(-0.9)} 
                 OH + OH + O3 -> H2O2 + O3                : {1.6e-30*(Tgas/300)^(-0.9)}
                 # second order OH recombination to H2O2 is strange. See:
                 # D. Liu. "Global model of low-temperature atmospheric-pressure He + H2O
                 #   plasmas," Plasma Sources Sci. Technol. 19 (2010)
                 # 
                 # But this is apparently the "high pressure limit" -- see
                 # Atkinson, Atmos. Chem. Phys. Discussions 3 (2003)
                 # I am excluding for now. Not entirely sure what to make of this.
                 #OH + OH -> H2O2                          : {1.5e-11*(Tgas/300)^(-0.37)}
                 ##############################
                 ##############################
                 OH + OHs -> H2O + O                      : {1.5e-11*(Tgas/300)^(-0.37)}
                 OH + OHp -> H2Op + O                     : 7e-10
                 OH + HO2 -> H2O + O2                     : {4.8e-11*exp(250/Tgas)}
                 OH + H2O2 -> H2O + HO2                   : {4.53e-12*exp(-288.9/Tgas)}
                 OHs + H2O2 -> HO2 + H2O                  : 1.47e-10
                 OHs + Ar -> OH + Ar                      : {1.44e-15*(Tgas/300)^(-0.5)}
                 OHs + H2O -> OH + H2O                    : {7.61e-10*(Tgas/300)^(-0.5)}
                 OHs + H2 -> OH + H2                      : 7.1e-11
                 OHp + O -> O2p + H2                      : 7.1e-10
                 OHp + Om -> O + OH                       : {2e-7*(Tgas/300)^(-0.5)}
                 OHp + Om -> O + O + H                    : 1e-7
                 OHp + Om + Ar -> OH + O + Ar             : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + Om + H2 -> OH + O + H2             : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + Om + Ar -> HO2 + Ar                : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + Om + H2 -> HO2 + H2                : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + O2 -> OH + O2p                     : 3.8e-10
                 # Next 3 are all rxn number 1399
                 OHp + OHm + Ar -> H2O2 + Ar              : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + OHm + H2O -> H2O2 + H2O            : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + OHm + H2 -> H2O2 + H2              : {2e-25*(Tgas/300)^(-2.5)}
                 OHp + H2O -> H2Op + OH                   : 1.6e-9
                 # Skipped because H3Op not included
                 #OHp + H2O -> H3Op + O                    : 1.3e-9
                 OHm + O -> HO2 + em                      : {4e-10*(Tgas/300)^0.5}
                 OHm + Op + Ar -> HO2 + Ar                : {2e-25*(Tgas/300)^(-2.5)}
                 OHm + O2p -> OH + O2                     : {2e-7*(Tgas/300)^(-0.5)}
                 OHm + O2p + Ar -> OH + O2 + Ar           : {2e-25*(Tgas/300)^(-2.5)}
                 OHm + O3 -> O2m + HO2                    : 1.08e-11
                 ###################################################
                 # OHm and H2O+ recombination
                 # Van gaens reactions 1443-1447
                 OHm + H2Op -> OH + H2O                   : {2e-7*(Tgas/300)^-2.5}
                 OHm + H2Op -> O + H + H2O                : 1e-7
                 OHm + H2Op -> OH + OH + H                : 1e-7
                 OHm + H2Op -> O + H + OH + H             : 1e-7
                 OHm + H2Op + Ar -> OH + H2O + Ar         : {2e-25*(Tgas/300)^(-2.5)}
                 OHm + H2Op + H2O -> OH + H2O + H2O        : {2e-25*(Tgas/300)^(-2.5)}
                 OHm + H2Op + H2 -> OH + H2O + H2         : {2e-25*(Tgas/300)^(-2.5)}
                 ###################################################
                 HO2 + O -> OH + O2                       : {2.71e-11*exp(224/Tgas)}
                 HO2 + O2s -> OH + O2 + O                 : 1.66e-11
                 HO2 + HO2 -> H2O2 + O2                   : {2.2e-13*exp(600/Tgas)}
                 HO2 + HO2 + Ar -> H2O2 + O2 + Ar         : {1.9e-33*exp(980/Tgas)}
                 HO2 + HO2 + H2 -> H2O2 + O2 + H2         : {4e-33*exp(980/Tgas)}
                 HO2 + HO2 + H2O -> H2O2 + O2 + H2O       : {1e-32*exp(980/Tgas)}
                 HO2 + H2O -> H2O2 + OH                   : {6.45e-11*exp(-16500/Tgas)}
                 HO2 + H2O2 -> OH + H2O + O2              : 1e-16
                 H2O + O -> OH + OH                       : {1.67e-11*(Tgas/300)^(1.14)*exp(-8680/Tgas)}
                 H2O + Os -> H2 + O2                      : 2.2e-12
                 H2O + Op -> H2Op + O                     : 2.6e-9
                 H2O + Om -> H2O2 + em                    : 6e-13
                 H2O + O2m -> H2O + O2 + em               : {5e-9*exp(-5000/Tgas)}
                 H2Op + Om -> OH + H + O                  : 1e-7
                 # Next six reactions are Van Gaens, reactions 1521-1522 (split M)
                 H2Op + Om + Ar -> H2O + O + Ar           : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + H2O -> H2O + O + H2O         : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + H2 -> H2O + O + H2           : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + Ar -> H2O2 + Ar              : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + H2O -> H2O2 + H2O            : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + Om + H2 -> H2O2 + H2              : {2e-25*(Tgas/300)^(-2.5)}
                 H2Op + O2 -> H2O + O2p                   : 3.3e-10
                 H2Op + O2m -> OH + H + O2                : 1e-7
                 H2Op + O2m + Ar -> H2O + O2 + Ar         : {2e-25*(Tgas/300)^(-0.5)} 
                 # Reactions 1255 and 1256 are WRONG -- see Erratum
                 #H2Op + O2m + Ar -> H2O2 + Ar             : {2e-25*(Tgas/300)^(-2.5)}
                 #H2Op + O2m + H2O -> H2O2 + H2O           : {2e-25*(Tgas/300)^(-2.5)}
                 #H2Op + O2m + H2 -> H2O2 + H2             : {2e-25*(Tgas/300)^(-2.5)}
                 #H2Op + O2m + Ar -> H2 + O2 + Ar          : {2e-25*(Tgas/300)^(-2.5)}
                 #H2Op + O2m + H2O -> H2 + O2 + H2O        : {2e-25*(Tgas/300)^(-2.5)}
                 #H2Op + O2m + H2 -> H2 + O2 + H2          : {2e-25*(Tgas/300)^(-2.5)}
                 H2O2 + O -> HO2 + OH                     : {1.79e-13*(Tgas/300)^(2.92)*exp(-1394/Tgas)}
                 H2O2 + O -> H2O + O2                     : 1.45e-15
                 H2O2 + Os -> H2O + O2                    : 5.2e-10
                 ##############################
                 # Radiation terms
                 # Note that Arss is a photoemission state
                 # Ars is a metastable, and accordingly has a low decay rate
                 # Arss is a photon emission state (as is Arsss, I presume -- must check)
                 ##############################
                 Arsss -> Ars                             : 2e6
                 Arsss -> Arss                            : 2e6
                 Ars -> Ar                                : 2e-2
                 Arss -> Ar                               : 1e8
                 Ar2s -> Ar + Ar                          : 3.17e5
                 O2s -> O2                                : 2.7e-4
                 OHs -> OH                                : 1.25e6'
                 
  [../]


  # Rate coefficients are in m^3 s^-1
  # Taken from Wei Tian's thesis
  # Note the difference in values (1 L = 1000 m^-3)
  [water2]
    species = 'H3Op_aq H_aq H2O2_aq OH_aq OHm_aq O2_aq O2m_aq O_aq Om_aq H2_aq HO2_aq HO2m_aq H2Op_aq O3m_aq O3_aq HO3_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    reaction_coefficient_format = 'rate'
    block = 1

    #reactions = 'em_aq + H2O_aq -> H_aq + OHm_aq            : 1.9e-2
    #             em_aq + H2Op_aq -> H_aq + OH_aq            : 6e8 
    #             em_aq + OH_aq -> OHm_aq                    : 3e7
    #             em_aq + H3Op_aq -> H_aq + H2O_aq                 : 2.3e7
    #             em_aq + H2O2_aq -> OH_aq + OHm_aq                : 1.1e7 
    #             em_aq + HO2m_aq -> OH_aq + OHm_aq + OHm_aq   : 3.5e6
    #             em_aq + O2_aq -> O2m_aq                        : 1.9e7
    #             # Next one is approximated by analogy (probably wrong...)
    #             em_aq + O_aq -> Om_aq                          : 1.9e7
                 # This one is listed as 1e10 in Chens work. Completely different.
                 # I am going with this value because I have seen it in multiple references.
    reactions = 'H_aq + OH_aq -> H2O_aq                        : 7e6
                 #H_aq + OHm_aq -> em_aq + H2O_aq                  : 2.2e4
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
                 # What is with this reaction?
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
                 ##em_aq + O2m_aq + H2O_aq -> HO2m_aq + OHm_aq   : 1.3e4
                 #em_aq + O2m_aq -> HO2m_aq + OHm_aq   : 1.3e7
                 #em_aq + HO2_aq -> HO2m_aq     : 2e7
                 ##em_aq + Om_aq + H2O_aq -> OHm_aq + OHm_aq     : 2.2e4
                 # This one is listed with conflicting units in literature. 
                 # (Three body reaction with a two body reaction rate coefficient.)
                 #em_aq + Om_aq -> OHm_aq + OHm_aq       : 2.2e7
                 ##em_aq + O3m_aq + H2O_aq -> O2_aq + OHm_aq + OHm_aq   : 1.6e4 
                 #em_aq + O3m_aq -> O2_aq + OHm_aq + OHm_aq : 1.6e7
                 #em_aq + O3_aq -> O3m_aq     : 3.6e7
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
                 HO2m_aq -> Om_aq + OH_aq         : 1e-5
                 H_aq + H2O_aq -> H2_aq + OH_aq                   : 1e-2
                 H_aq + H_aq -> H2_aq                          : 7.5e6
                 H2_aq + H2O2_aq -> H_aq + OH_aq + H2O_aq            : 6e3
                 OH_aq + H2_aq -> H_aq + H2O_aq     : 4.2e4
                 Om_aq + H2_aq -> OHm_aq + H_aq       : 8e7'
  []

  #[water3]
  #  name = 'test'
  #  species = 'H2_aq OHm_aq H_aq Om_aq OH_aq'
  #  aux_species = 'H2O_aq'
  #  use_log = true
  #  position_units = ${dom1Scale}
  #  track_rates = true
  #  reaction_coefficient_format = 'rate'
  #  block = 1
  #
  #  #reactions = 'em_aq + em_aq -> H2_aq + OHm_aq + OHm_aq   : 5.5e6 
  #  #             em_aq + H_aq -> H2_aq + OHm_aq             : 2.5e7
  #[]
[]
