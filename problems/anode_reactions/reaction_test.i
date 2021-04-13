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
  #end_time = 1e-1
  #end_time = 1e10
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  #resid_vs_jac_scaling_param = 0.5
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  #solve_type = newton
  solve_type = pjfnk
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu NONZERO 1.e-10 superlu_dist'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_package -snes_stol'
  #petsc_options_value = 'lu NONZERO 1.e-10 mumps 0'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  #petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-9
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
  [out_01]
    type = Exodus
  [../]
[]

[Debug]
  show_var_residual_norms = true
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
    #charged_particle = 'Arp Ar2p'
    #Neutrals = 'Ar* H2O OH'
    charged_particle = 'Arp'
    Neutrals = 'Ar*'
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
    charged_particle = 'em_aq OHm_aq H3Op_aq'
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
  #[H2O]
  #  block = 0
  #  #initial_condition = 0.0367321
  #[]
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
  [./em_aq]
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
  #[./Ar2p]
  #  block = 0
  #  initial_condition = -20.693147
  #[../]
  #[OH]
  #  block = 0
  #  initial_condition = -20
  #[]
  #[./OH_aq]
  #  block = 1
  #  #initial_condition = -20
  #  initial_condition = -5
  #[../]

  [./Ar*]
    block = 0
    initial_condition = -25
  [../]

  [./mean_en]
    block = 0
    initial_condition = -20
    # scaling = 1e-1
  [../]

  [./OHm_aq]
    block = 1
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    #initial_condition = -9.210340
    initial_condition = -14
  [../]
  [H3Op_aq]
    block = 1
    initial_condition = -14
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
  [./ADCurrent_OHm]
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
[]

[AuxKernels]
  [./H2O_aq_density]
    type = DensityMoles
    variable = H2O_aq_density
    density_log = H2O_aq
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
  [./ADCurrent_OHm]
    block = 1
    type = ADCurrent
    #potential = potential_liq
    potential = potential
    density_log = OHm_aq
    variable = ADCurrent_OHm
    art_diff = false
    position_units = ${dom1Scale}
  [../]
[]

[InterfaceKernels]
  [./em_advection]
    type = ADInterfaceAdvection
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
    neighbor_var = em
    variable = em_aq
    boundary = water_left
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
  [../]

  # Now we test the henry interface condition
  #[OH_diff]
  #  type = InterfaceDiffusionTest
  #  variable = OH_aq
  #  neighbor_var = OH
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
  #[OH_henry]
  #  type = InterfaceReactionTest
  #  variable = OH_aq
  #  neighbor_var = OH
  #  #kf = 1
  #  #kb = 1
  #  h = 6.48e3
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #  boundary = 'water_left'
  #[]
[]

[BCs]
  [h3op_test]
    type = WaterCircuitBC
    variable = H3Op_aq
    input_current = gamma_e 
    output_current = gamma_w
    potential = potential
    current_carrier = 'em_aq OHm_aq'
    water_surface_area = 5.02e-7 # Formerly 3.14e-6
    position_units = ${dom1Scale}
    boundary = 'right'
  []
  # H2O evaporation boundary condition
  #[H2O_interface]
  #  type = DirichletBC
  #  preset = true
  #  variable = H2O
  #  value = 0.367321
  #  boundary = 'gas_right'
  #[]

  [./Arex_physical_diffusion]
    type = ADHagelaarIonDiffusionBC
    variable = Ar*
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./Ar2p_physical_diffusion]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = Ar2p
  #  boundary = 'left gas_right'
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_physical_advection]
  #  type = ADHagelaarIonAdvectionBC
  #  variable = Ar2p
  #  boundary = 'left gas_right'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./potential_left]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_dirichlet_right]
    type = DirichletBC
    preset = true
    variable = potential
    boundary = right
    value = 0
  [../]
  [./em_bc]
    type = ADHagelaarElectronBC
    variable = em
    boundary = 'left gas_right'
    potential = potential
    mean_en = mean_en
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion_bc]
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'left gas_right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_advection_bc]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'left gas_right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_bc]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'left gas_right'
    potential = potential
    em = em
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    #ip = 'Arp Ar2p'
    ip = 'Arp'
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
  [./OHm_physical]
    type = ADDCIonBC
    variable = OHm_aq
    boundary = 'right'
    #potential = potential_liq
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./em_aq_diffusion_bc]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = em_aq
  #  boundary = 'right'
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[./em_aq_advection_bc]
  #  type = ADHagelaarIonAdvectionBC
  #  variable = em_aq
  #  boundary = 'right'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  #
  #[./OHm_aq_diffusion_bc]
  #  type = ADHagelaarIonDiffusionBC
  #  variable = OHm_aq
  #  boundary = 'right'
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[./OHm_aq_advection_bc]
  #  type = ADHagelaarIonAdvectionBC
  #  variable = OHm_aq
  #  boundary = 'right'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[H3Op_physical]
  #  type = ADDCIonBC
  #  variable = H3Op_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[]
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
    value = -1.5
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
  [./se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'left gas_right'
  [../]
 [./GasBasics]
   type = ADGasElectronMoments
   em = em
   mean_en = mean_en
   #property_tables_file = '/home/shane/projects/zapdos/problems/argon_cm_test/electron_mobility_diffusion.txt'
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

  [./OHm_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OHm_aq
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]
  [./H3O+_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H3Op_aq
    heavy_species_mass = 3e-26
    # ^ just estimated...whatever
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
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
[]

[Postprocessors]
  [gamma_e]
    type = ADElectronFluxIntegral
    variable = em
    r = 0
    potential = potential
    mean_energy = mean_en
    position_units = ${dom0Scale}
    boundary = 'gas_right'
  []
  [gamma_w]
    type = ADTotalFluxIntegral
    variable = em_aq
    ions = 'em_aq OHm_aq'
    potential = potential
    position_units = ${dom1Scale}
    boundary = 'right'
  []
[]

[Reactions]
  [./Argon]
    #species = 'em Arp Ar2p Ar* OH H2O'
    species = 'em Arp Ar*'
    #aux_species = 'Ar H2O'
    reaction_coefficient_format = 'townsend'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'argon_chemistry_rates'
    equation_constants = 'Tgas'
    equation_values = '300'
    equation_variables = 'e_temp'
    potential = 'potential'
    use_log = true
    position_units = ${dom0Scale}
    use_ad = true
    block = 0

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (reaction1)
                 em + Ar -> em + Ar*              : EEDF [-11.5]   (reaction2)
                 em + Ar -> em + em + Arp         : EEDF [-15.76]  (reaction3)
                 em + Ar* -> em + Ar              : EEDF [11.5]    (reaction4)
                 em + Ar* -> em + em + Arp        : EEDF [-4.3]    (reaction5)'
                 #Ar2p + em -> Ar* + Ar            : {5.1187e11*(e_temp*11600/300)^(-0.67)}
                 #Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12/Tgas*exp(-15130/Tgas)}
                 #Ar* + Ar* -> Ar2p + em           : 3.6132e8
                 #Arp + em + em -> Ar + em         : {3.17314235e9*(e_temp)^(-4.5)}
                 #Ar* + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
                 #Arp + Ar + Ar -> Ar2p + Ar       : {81595.089*(Tgas/300)^(-0.4)}
                 #Ar* + H2O -> Ar + OH + H         : 2.89056e8'
                 #Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}'
  [../]


  [./liquid_phase_reactions]
    species = 'em_aq OHm_aq H3Op_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'em_aq + em_aq -> H2 + OHm_aq + OHm_aq       : 3.0703e8
                 em_aq + H3Op_aq -> H_aq + H2O_aq                 : 2.3e7
                 H3Op_aq + OHm_aq -> H_aq + OH_aq + H2O_aq     : 6e7'
  [../]
[]
