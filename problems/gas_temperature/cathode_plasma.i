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
  [mesh]
    type = FileMeshGenerator
    file = 'cathode_plasma.msh'
  []

  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'plasma_left'
    input = mesh
  []

  [interface_again]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'electrode_right'
    input = interface
  []

  # The next two definitions create boundary conditions named
  # 'left' and 'right', where 'left' is at x = 0 and 'right' is at x = 1.1 mm.
  [left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = interface_again
  []
  [right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  []
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
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  #resid_vs_jac_scaling_param = 1
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  #solve_type = pjfnk
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmin = 1e-18
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-14
    #dt = 2.56e-13
    growth_factor = 1.4
    optimal_iterations = 10
  []
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  [out_01]
    type = Exodus
  [../]
[]

[Debug]
  #show_var_residual_norms = true
[]

[UserObjects]
  [data_provider]
    type = ProvideMobility
    electrode_area = 5.02e-7 # Formerly 3.14e-6
    ballast_resist = 651e3
    e = 1.6e-19
  []
[]

[DriftDiffusionAction]
  [Plasma]
    electrons = em
    charged_particle = 'Arp Ar2p'
    Neutrals = 'Ars'
    mean_energy = mean_en
    potential = potential
    using_offset = true
    offset = 30
    use_ad = true
    position_units = ${dom0Scale}
    block = 0
  []
[]

[Variables]
  # Potential variable
  [potential]
    block = 0
  []

  # Gas temperature variable
  # Note that it exists in both blocks (plasma and cathode regions)
  [Tg]
    block = '0 1'
    initial_condition = 300
  []

  # Plasma species variables
  [em]
    block = 0
    initial_condition = -20
  []
  [Arp]
    block = 0
    initial_condition = -20.693147
  []
  [Ar2p]
    block = 0
    initial_condition = -20.693147
  []
  [Ars]
    block = 0
    initial_condition = -25
  []
  [mean_en]
    block = 0
    initial_condition = -20
  []
[]

[Kernels]
  # Diffusion and time deriviatve happen in both blocks, but joule heating only occurs in the plasma phase
  [HeatDiff]
    type = ADHeatConduction
    variable = Tg
    block = '0 1'
  []
  [HeatTdot]
    type = ADHeatConductionTimeDerivative
    variable = Tg
    block = '0 1'
  []

  [GasJouleHeating]
    type = JouleHeatingIons
    variable = Tg
    ions = 'Arp Ar2p'
    potential = potential
    position_units = ${dom0Scale}
    block = 0
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
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    execute_on = 'initial timestep_end'
    block = 0
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

[BCs]
  # Gas temperature boundary conditions
  [lefttemp]
    type = ADDirichletBC
    variable = Tg
    boundary = left
    value = 300
  []

  [Tgas_left]
    type = DirichletBC
    variable = Tg
    boundary = 'plasma_left'
    value = 350
  []

  [Arex_physical_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Ars
    boundary = 'plasma_left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Ar2p_physical_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Ar2p
    boundary = 'plasma_left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Ar2p_physical_advection]
    type = HagelaarIonAdvectionBC
    variable = Ar2p
    boundary = 'plasma_left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  []
  [potential_left]
    type = NeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = plasma_left
    function = potential_bc_func
    ip = 'Arp Ar2p'
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  []
  [potential_dirichlet_right]
    type = DirichletBC
    variable = potential
    boundary = right
    value = 0
  []
  [em_physical_bc]
    type = HagelaarElectronBC
    variable = em
    boundary = 'plasma_left right'
    potential = potential
    mean_en = mean_en
    r = 0.0
    position_units = ${dom0Scale}
  []
  [Arp_physical_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = 'plasma_left right'
    r = 0
    position_units = ${dom0Scale}
  []
  [Arp_physical_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = 'plasma_left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  []
  [mean_en_physical]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'plasma_left right'
    potential = potential
    em = em
    r = 0.0
    position_units = ${dom0Scale}
  []
  [sec_electrons_left]
    type = SecondaryElectronBC
    variable = em
    boundary = 'plasma_left'
    potential = potential
    ip = 'Arp Ar2p'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  []
  [sec_electrons_energy_left]
    type = SecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'plasma_left'
    potential = potential
    ip = 'Arp Ar2p'
    em = em
    r = 0
    position_units = ${dom0Scale}
  []
[]

[ICs]
  [potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
  []
[]

[Functions]
  [potential_bc_func]
    type = ParsedFunction
    value = -1.0
  []
  [potential_ic_func]
    type = ParsedFunction
    value = '-1.0 * (1.1e-2 - x)'
  []
[]

[Materials]
  [convective_bc_props]
    type = ADGenericConstantMaterial
    prop_names = 'Tinf htc'
    prop_values = '300 10' 
  []
  
  # Heat conduction properties -- argon (block 0)
  [density_argon]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.784' # kg m^-3
    block = 0
  []
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

  # Heat conduction properties -- electrode (block 1)
  [density_ss420]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '7850' # kg m^-3
    block = 1
  []
  [thermal_conductivity_ss420]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '16.26' # W m^-1 K^-1
    block = 1
  []
  [specific_heat_ss420]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '502.1' # J kg^-1 K^-1
    block = 1
  []


  [se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'plasma_left right'
  []
  [GasBasics]
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
  []
  [gas_constants]
    type = GenericConstantMaterial
    prop_names = 'e    N_A    k_boltz    eps     se_energy    T_gas    massem    p_gas'
    prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12 1 400 9.11e-31 1.01e5'
    block = 0
  []

  [potential_permittivity]
    type = ADGenericConstantMaterial
    prop_names = 'diffpotential'
    prop_values = '8.854e-12'
  []

  [gas_species_0]
    type = ADHeavySpecies
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    block = 0
  []
  [Ar_species]
    type = ADHeavySpecies
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  []
  [gas_species_1]
    type = ADHeavySpecies
    heavy_species_name = Ar2p
    heavy_species_mass = 13.28e-26
    heavy_species_charge = 1.0
    block = 0
  []
  [gas_species_2]
    type = ADHeavySpecies
    heavy_species_name = Ars
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  []

[]

[Reactions]
  [Argon]
    species = 'em Arp Ar2p Ars'
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
                 em + Ar -> em + Ars              : EEDF [-11.5]   (reaction2)
                 em + Ar -> em + em + Arp         : EEDF [-15.76]  (reaction3)
                 em + Ars -> em + Ar              : EEDF [11.5]    (reaction4)
                 em + Ars -> em + em + Arp        : EEDF [-4.3]    (reaction5)
                 Ar2p + em -> Ars + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
                 Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tg * exp(-15130/Tg)}
                 Ars + Ars -> Ar2p + em           : 3.6132e8
                 Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
                 Ars + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
                 Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tg/300)^(-0.4)}'
  []
[]
