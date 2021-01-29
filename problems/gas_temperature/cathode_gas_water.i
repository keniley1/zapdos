dom0Scale=1.0
dom1Scale=1.0
dom2Scale=1.0

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
    #file = 'mesh_nowater.msh'
    #file = 'gas_water_mesh.msh'
    file = 'cathode_gas_water_mesh.msh'
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

  [./interface3]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '2'
    new_boundary = 'gas_left'
    input = interface2
  [../]
  [./interface4]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '2'
    paired_block = '0'
    new_boundary = 'cathode_right'
    input = interface3
  [../]

  # The next two definitions create boundary conditions named
  # 'left' and 'right', where 'left' is at x = 0 and 'right' is at x = 1.1 mm.
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = interface4
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
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmin = 1e-18
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
  [out_05]
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
    ballast_resist = 651e3
    e = 1.6e-19
  [../]
[]

[DriftDiffusionActionAD]
  [./Plasma]
    electrons = em
    charged_particle = 'Arp Ar2p'
    Neutrals = 'Ar* H2O OH'
    mean_energy = mean_en
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    position_units = ${dom0Scale}
    block = 0
  [../]

  [./Water]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    charged_particle = 'emliq OH-'
    #Neutrals = 'OH_aq'
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
  [Tg]
    block = 0
    initial_condition = 300
  []
  [Tw]
    block = 1
    initial_condition = 300
  []
  [Tc]
    block = 2
    initial_condition = 300
  []
  [H2O]
    block = 0
    #initial_condition = 0.0367321
  []
  [./potential]
    #block = 0
    block = '0 1'
  [../]
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
  #[OH]
  #  block = 0
  #  initial_condition = -20
  #[]
  [./Ar*]
    block = 0
    initial_condition = -25
  [../]

  [./mean_en]
    block = 0
    initial_condition = -20
    # scaling = 1e-1
  [../]

  [./emliq]
    block = 1
    #initial_condition = -24
    #initial_condition = -21
    # scaling = 1e-5
    initial_condition = -14
    #initial_condition = -24
  [../]
  #[./OH_aq]
  #  block = 1
  #  #initial_condition = -20
  #  initial_condition = -5
  #[../]
  [./OH-]
    block = 1
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    #initial_condition = -9.210340
    initial_condition = -14
  [../]
[]

[Kernels]
  [HeatDiff]
    type = ADHeatConduction
    variable = Tg
    block = 0
  []
  [HeatTdot]
    type = ADHeatConductionTimeDerivative
    variable = Tg
    block = 0
  []
  [GasJouleHeating]
    type = JouleHeatingIons
    variable = Tg
    ions = 'Arp Ar2p'
    potential = potential
    position_units = ${dom0Scale}
    block = 0
  []

  [HeatDiff_water]
    type = ADHeatConduction
    variable = Tw
    block = 1
  []
  [HeatTdot_water]
    type = ADHeatConductionTimeDerivative
    variable = Tw
    block = 1
  []

  [HeatDiff_cathode]
    type = ADHeatConduction
    variable = Tc
    block = 2
  []
  [HeatTdot_cathode]
    type = ADHeatConductionTimeDerivative
    variable = Tc
    block = 2
  []

  #[ElasticTest]
  #  type = GasTemperatureElastic
  #  variable = Tg
  #  potential = potential
  #  electrons = em
  #  mean_energy = mean_en
  #  target = Ar
  #  reaction = 'em + Ar -> em + Ar'
  #  position_units = ${dom0Scale}
  #  number = 0
  #[]
[]

[AuxVariables]
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.701920755421197
    #initial_condition = 3.704261
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
  [./ADCurrent_Arp]
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
  [./ADCurrent_em]
    type = ADCurrent
    potential = potential
    density_log = em
    variable = ADCurrent_em
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
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
  [./em_advection]
    type = ADInterfaceAdvection
    potential_neighbor = potential
    neighbor_var = em
    variable = emliq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    #type = InterfaceLogDiffusionElectrons
    type = ADInterfaceLogDiffusion
    neighbor_var = em
    variable = emliq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]

  # Now we test the henry interface condition
  #[OH_diff]
  #  type = InterfaceDiffusionTest
  #  variable = OH_aq
  #  neighbor_var = OH
  #  h = 6.48e3
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
  #[OH_henry]
  #  type = InterfaceReactionTest
  #  variable = OH_aq
  #  neighbor_var = OH
  #  #kf = 6.48e3
  #  #kb = 1
  #  kf = 1
  #  kb = 6.2e2 
  #  #kb = 1
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]

  [temp]
    type = TemperatureTest
    variable = Tw
    neighbor_var = Tg
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'water_left'
  []

  [temp2]
    type = TemperatureTestOneWay
    variable = Tc
    neighbor_var = Tg
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = 'cathode_right'
  []
  #[temp2]
  #  type = IonBombardment
  #  variable = Tc
    
  #[]
[]

[BCs]
  # nodal BC to ensure temperatures are the same
  [tg_tw_bc]
    type = MatchedValueBC
    variable = Tg
    v = Tw
    boundary = gas_right
  []
  [cathode_Tc]
    type = MatchedValueBC
    variable = Tg
    v = Tc
    boundary = gas_left
  []
  [cathode_test]
    type = IonHeating
    variable = Tg
    ions = 'Arp Ar2p'
    potential = potential
    r = 0
    position_units = ${dom2Scale}
    boundary = 'gas_left'
  []
  # Gas temperature boundary conditions
  #[lefttemp]
  #  type = ADDirichletBC
  #  variable = Tg
  #  boundary = gas_left
  #  value = 350
  #[]
  #[lefttemp]
  #  type = ADConvectiveHeatFluxBC
  #  variable = Tg
  #  boundary = left
  #  T_infinity = Tinf
  #  heat_transfer_coefficient = htc
  #[]
  [righttemp]
    type = ADDirichletBC
    variable = Tw
    boundary = right
    value = 300
  []
  # H2O evaporation boundary condition
  #[H2O_interface]
  #  type = DirichletBC
  #  variable = H2O
  #  value = 0.367321
  #  boundary = 'gas_right'
  [H2O_interface]
    type = VaporPressureBC
    variable = H2O
    gas_temp = Tg
    boundary = gas_right
  []
  #[OH_bc]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = OH
  #  boundary = 'gas_left gas_right'
  #  r = 0
  #  position_units = ${dom0Scale}
  #[]
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
    boundary = 'gas_left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar2p
    boundary = 'gas_left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Ar2p_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Ar2p
    boundary = 'gas_left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_left]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = gas_left
    function = potential_bc_func
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
  [./em_physical_bc]
    type = ADHagelaarElectronBC
    variable = em
    boundary = 'gas_left gas_right'
    potential = potential
    mean_en = mean_en
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'gas_left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_advection]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'gas_left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'gas_left gas_right'
    potential = potential
    em = em
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    boundary = 'gas_left'
    potential = potential
    ip = 'Arp Ar2p'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'gas_left'
    potential = potential
    ip = 'Arp Ar2p'
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
    block = '0 1'
  [../]
  [H2O_ic]
    type = FunctionIC
    variable = H2O
    function = water_ic_func
    block = 0
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
    value = -1.5
  [../]
  [./test_bc]
    type = ParsedFunction
    value = '-2.5*tanh(1e9*t)'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-1.5 * (1.001e-3 - x)'
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
  [convective_bc_props]
    type = ADGenericConstantMaterial
    prop_names = 'Tinf htc'
    prop_values = '300 10'
    block = 0
  []
  [density_argon]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.784' # kg m^-3
    block = 0
  []

  # Gas phase
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
  [density_water]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '997' # kg m^-3
    block = 1
  []

  # Water phase
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
  [density_tungsten]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '19300' # kg m^-3
    block = 2
  []
  
  # Cathode
  [thermal_conductivity_cathode]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '173' # W m^-1 K^-1
    block = 2
  []
  [specific_heat_cathode]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '134' # J kg^-1 K^-1
    block = 2
  []
  [work_function_cathode]
    type = GenericConstantMaterial
    prop_names = 'work_function'
    prop_values = '4.5' # W m^-1 K^-1
    block = 2
  []
  
  [./se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'gas_left gas_right'
  [../]
 [./GasBasics]
   type = ADGasElectronMoments
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   user_p_gas = 1.01e5
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   property_tables_file = 'argon_chemistry_rates/electron_moments.txt'
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

  [OH_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = 0
    diffusivity = 4e-5
    block = 0
  []

  [./OH-_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH-
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]

  [./OH_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH_aq
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = 0 
    diffusivity = 5.29e-9
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

  # Why do these need to be declared? Where are they used?
  [cathode_constants]
   type = GenericConstantMaterial
   prop_names = 'p_gas T_gas'
   prop_values = '101325 300'
   block = 2
  []
[]

[Reactions]
  [./Argon]
    #species = 'em Arp Ar2p Ar* OH'
    species = 'em Arp Ar2p Ar*'
    aux_species = 'Ar H2O'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'argon_chemistry_rates'
    #equation_constants = 'Tgas'
    #equation_values = '300'
    equation_variables = 'e_temp Tg'
    potential = 'potential'
    use_log = true
    position_units = ${dom0Scale}
    use_ad = true
    block = 0

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (reaction1)
                 em + Ar -> em + Ar*              : EEDF [-11.5]   (reaction2)
                 em + Ar -> em + em + Arp         : EEDF [-15.76]  (reaction3)
                 em + Ar* -> em + Ar              : EEDF [11.5]    (reaction4)
                 em + Ar* -> em + em + Arp        : EEDF [-4.3]    (reaction5)
                 Ar2p + em -> Ar* + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
                 #Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tgas * exp(-15130/Tgas)}
                 Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tg * exp(-15130/Tg)}
                 Ar* + Ar* -> Ar2p + em           : 3.6132e8
                 Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
                 Ar* + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
                 #Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}
                 Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tg/300)^(-0.4)}
                 #Ar* + H2O -> Ar + OH + H         : 2.89056e8'
                 #Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}'
  [../]

  [./liquid_phase_reactions]
    #species = 'emliq OH- OH_aq'
    species = 'emliq OH-'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    block = 1
    reaction_coefficient_format = 'rate'
    #reactions = 'emliq + H2O -> H + OH-               : 1.9e-2
    #             emliq + emliq -> H2 + OH- + OH-      : 3.0703e8'
    reactions = 'emliq + emliq -> H2 + OH- + OH-       : 3.0703e8'
                 #emliq + OH_aq -> OH-                     : 3e7'
  [../]
[]
