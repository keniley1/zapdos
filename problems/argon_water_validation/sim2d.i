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
    #file = 'mesh01.msh'
    file = 'water_2d.msh'
  [../]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
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
  #compute_scaling_once = false
  #resid_vs_jac_scaling_param = 0.5
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
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 2e-13
    growth_factor = 1.4
    optimal_iterations = 10
  [../]
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  [out_002]
    type = Exodus
    interval = 10
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
    
    # Following Necip's experimental parameters
    # Set to 651 kOhm. May be varied in increments of 160 kOhm.
    #ballast_resist = 3.31e5 
    #ballast_resist = 4.91e5
    ballast_resist = 6.51e5  # 651 kOhm
    #ballast_resist = 8.11e5 
    #ballast_resist = 9.71e5
    e = 1.6e-19
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[DriftDiffusionAction]
  [./Water]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    #charged_particle = 'em_aq OHm_aq'
    #Neutrals = 'OH_aq'
    #charged_particle = 'em_aq H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
    #Neutrals = 'H_aq H2O2_aq OH_aq O2_aq O_aq H2_aq HO2_aq O3_aq HO3_aq'
    charged_particle = 'em_aq OHm_aq H3Op_aq'
    #charged_particle = 'em_aq'
    #Is_potential_unique = false
    potential = potential
    using_offset = true
    offset = 40
    use_ad = true
    order = SECOND
    position_units = ${dom1Scale}
    block = 0
  [../]
  [./Salt]
    # Missing Na+, Cl-, NO2-, NO2_2-, NO3-, NO3_2-
    charged_particle = 'Nap_aq Clm_aq'
    First_DriftDiffusionAction_in_block = false
    potential = potential
    using_offset = true
    #offset = -2.3026
    offset = 40
    use_ad = true
    order = SECOND 
    position_units = ${dom1Scale}
    block = 0
  [../]
[]

[Variables]
  [potential]
    block = 0
    order = SECOND
  []
  ###########################
  # Water electrolyte
  ###########################
  [Nap_aq]
    block = 0
    order = SECOND
    # 100 mM
    #initial_condition = 4.60517
    # 10 mM
    #initial_condition = 2.3026

    # 20 mM 
    initial_condition = 2.99575
  [] 
  [Clm_aq]
    block = 0
    order = SECOND
    # 100 mM
    #initial_condition = 4.60517
    # 10 mM
    #initial_condition = 2.3026

    # 20 mM 
    initial_condition = 2.99575
  [] 

  # Water Species
  # e- starts at some small value
  # (10^15 isn't exactly SMALL, but water density is O(10^28) so...
  # small enough!)
  [./em_aq]
    block = 0
    order = SECOND
    #initial_condition = -24
    #initial_condition = -21
    # scaling = 1e-5
    #initial_condition = -14
    #initial_condition = -24
    #initial_condition = -20
    initial_condition = -20
  [../]

  # pH = 7 -- equal parts OH- and H3O+
  [./OHm_aq]
    block = 0
    order = SECOND
    # scaling = 1e-5
    #initial_condition = -24
    #initial_condition = -21
    #initial_condition = -9.210340
    #initial_condition = -20
    initial_condition = -16.1175
  [../]

  [H3Op_aq]
    block = 0
    order = SECOND
    initial_condition = -16.1175
  []

  #[./OH_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./H3Op_aq]
  #  block = 0
  #  #initial_condition = -20
  #  initial_condition = -14
  #[../]
  #[./O2m_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./HO2m_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./H2Op_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./O3m_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./H_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./H2O2_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./O2_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./O_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./H2_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./HO2_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./HO3_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./O3_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
  #[./Om_aq]
  #  block = 0
  #  initial_condition = -20
  #[../]
[]

[AuxVariables]
  [bc_test]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
    block = 0
  []
  [./H2O_aq]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 10.92252
    block = 0
  [../]
  [./H2O_aq_density]
    order = CONSTANT
    family = MONOMIAL
    block = 0
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
  [./rholiq]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./Efield_x]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./Efield_y]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
[]

[AuxKernels]
  [test]
    type = FunctionAux
    variable = bc_test
    function = electron_bc_gaussian
    execute_on = 'initial timestep_end'
  []
  [./H2O_aq_density]
    type = DensityMoles
    variable = H2O_aq_density
    density_log = H2O_aq
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_l]
    type = Position
    variable = x
    position_units = ${dom1Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_nl]
    type = Position
    variable = x_node
    position_units = ${dom1Scale}
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [rholiq_calc]
    type = ChargeDensity
    variable = rholiq
    #charged = 'em_aq H3Op_aq OHm_aq O2m_aq Om_aq HO2m_aq H2Op_aq O3m_aq'
    #charged = 'em_aq OHm_aq Nap_aq Clm_aq'
    charged = 'em_aq'
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  []
  [./Efield_x]
    type = Efield
    component = 0
    #potential = potential_liq
    potential = potential
    variable = Efield_x
    position_units = ${dom1Scale}
    block = 0
  [../]
  [./Efield_y]
    type = Efield
    component = 1
    #potential = potential_liq
    potential = potential
    variable = Efield_y
    position_units = ${dom1Scale}
    block = 0
  [../]
[]

[BCs]
  [em_top]
    type = FunctionFlux
    variable = em_aq
    function = electron_bc_gaussian
    position_units = ${dom0Scale}
    boundary = top
  []
  #[potential_neumann_left]
  #  type = NeumannBC
  #  variable = potential
  #  boundary = top
  #  #value = -6.36942675
  #  value = -1e-8
  #  # divided by 1000 to be in kV/m
  #[]
  [potential_neumann_bc]
    type = FunctionFlux
    variable = potential
    function = potential_bc_gaussian
    position_units = ${dom0Scale}
    boundary = top
  []
  [ground]
    type = DirichletBC
    variable = potential
    value = 0
    boundary = ground
  []
  #[./Nap_aq_physical]
  #  type = ADDCIonBC
  #  variable = Nap_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
  #[./Clm_aq_physical]
  #  type = ADDCIonBC
  #  variable = Clm_aq
  #  boundary = 'right'
  #  potential = potential
  #  position_units = ${dom1Scale}
  #[../]
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
[]

[Functions]
  [electron_bc_gaussian]
    type = ParsedFunction
    vars = 'sigma mu'
    vals = '1e-5 0'
    #value = 'log(1/(sigma*sqrt(2*3.14159)) * exp(-0.5*(x-mu)^2/(sigma^2)))'
    value = '0.05 * exp(-0.5*(x-mu)^2/(sigma^2)) * tanh(1e9*t)'
  []
  [potential_bc_gaussian]
    type = ParsedFunction
    vars = 'sigma mu'
    vals = '1e-5 0'
    #value = 'log(1/(sigma*sqrt(2*3.14159)) * exp(-0.5*(x-mu)^2/(sigma^2)))'
    value = '1e-8 * exp(-0.5*(x-mu)^2/(2*sigma^2)) * tanh(1e9*t)'
  []
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
    value = '-0.8 * (1.001e-3 - y)'
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
    block = 0
  [../]
  [./Clm_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Clm_aq
    heavy_species_mass = 5.887e-26
    heavy_species_charge = -1.0
    diffusivity = 2e-9
    block = 0
  [../]
  [H2Op_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2Op
    heavy_species_mass = 2.9907e-26
    heavy_species_charge = 1
    diffusivity = 2.3e-5
  []

  [./electron_data]
    type = ADGenericConstantMaterial
    prop_names = 'diffem_aq muem_aq Tem_aq'
    prop_values = '4.5e-9 0.000173913 300'
    block = 0
  [../]
  [./electron_sign]
    # I increased the electron mass by a factor of 10 here
    # Not sure what it's really supposed to be
    type = GenericConstantMaterial
    prop_names = 'sgnem_aq massem_aq'
    prop_values = '-1 9.11e-30'
    block = 0
  [../]

  [./N_A_mat]
    type = GenericConstantMaterial
    prop_names = 'N_A e diffpotential diffpotential_liq T_gas p_gas'
    prop_values = '6.022e23 1.602e-19 7.0832e-10 7.0832e-10 300 1.013e5'
    block = 0
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
    block = 0
  [../]

  [./O_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O_aq
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 0
  [../]

  [./O3_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O3_aq
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = 0
    diffusivity = 5e-9
    block = 0
  [../]

  [./O3m_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O3m_aq
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 0
  [../]

  [./OH_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = OH_aq
    heavy_species_mass = 2.82431e-26
    heavy_species_charge = 0
    diffusivity = 4.5e-9
    block = 0
  [../]

  [./HO2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO2_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 0
  [../]

  [./HO3_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO3_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 0
  [../]

  [./H2O2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H2O2_aq
    heavy_species_mass = 5.64840e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 0
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
    block = 0
  [../]

  [./H_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H_aq
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 0
    diffusivity = 5.0e-9
    block = 0
  [../]

  [./H+_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Hp_aq
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 0
  [../]

  [./H3O+_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = H3Op_aq
    heavy_species_mass = 3e-26
    # ^ just estimated...whatever
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    # (Is this really the same as H3O+? I don't understand)
    block = 0
  [../]

  [./HO2-_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = HO2m_aq
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 0
  [../]

  [./Om_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Om_aq
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 0
  [../]

  [./O2m_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2m_aq
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = -1
    diffusivity = 5e-9 
    block = 0
  [../]

  [./O2_aq_mat]
    type = ADHeavySpeciesMaterial
    heavy_species_name = O2_aq
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 2e-9
    block = 0
  [../]
[]

[Reactions]
  # Rate coefficients are in m^3 s^-1
  # Taken from Wei Tian's thesis
  # Note the difference in values (1 L = 1000 m^-3)
  [water2]
    species = 'em_aq OHm_aq H3Op_aq'
    #species = 'em_aq'
    aux_species = 'H2O_aq'
    use_log = true
    position_units = ${dom1Scale}
    track_rates = false
    reaction_coefficient_format = 'rate'
    block = 0

    reactions = 'em_aq + em_aq -> H2_aq + OHm_aq + OHm_aq   : 5.5e6
                 em_aq + H3Op_aq -> H_aq + H2O_aq           : 2.5e7
                 H3Op_aq + OHm_aq -> H_aq + OH_aq + H2O_aq  : 6e7'
    #reactions = 'em_aq + H2O_aq -> H_aq + OHm_aq            : 1.9e-2
    #             #em_aq + H2Op_aq -> H_aq + OH_aq            : 6e-8 
    #             em_aq + H2Op_aq -> H_aq + OH_aq            : 6e8 
    #             # The units on this next one make no sense and have no consistency 
    #             # across literature. This is the Buxton value (halved)
    #             #em_aq + H_aq + H2O_aq -> H2_aq + OHm_aq    : 2.5e7
    #             em_aq + OH_aq -> OHm_aq                    : 3e7
    #             # where does this one come from???
    #             #em_aq + em_aq + H2O_aq -> OHm_aq + OHm_aq  : 2.2e4
    #             em_aq + H3Op_aq -> H_aq + H2O_aq                 : 2.3e7
    #             em_aq + H2O2_aq -> OH_aq + OHm_aq                : 1.1e7 
    #             #em_aq + HO2m_aq + H2O_aq -> OH_aq + OHm_aq + OHm_aq     : 3.5e3
    #             em_aq + HO2m_aq -> OH_aq + OHm_aq + OHm_aq   : 3.5e6
    #             em_aq + O2_aq -> O2m_aq                        : 1.9e7
    #             # Next one is approximated by analogy (probably wrong...)
    #             em_aq + O_aq -> Om_aq                          : 1.9e7
    #             # This one is listed as 1e10 in Chens work. Completely different.
    #             # I am going with this value because I have seen it in multiple references.
    #             H_aq + OH_aq -> H2O_aq                        : 7e6
    #             H_aq + OHm_aq -> em_aq + H2O_aq                  : 2.2e4
    #             H_aq + H2O2_aq -> OH_aq + H2O_aq                 : 9e4
    #             H_aq + O2_aq -> HO2_aq                        : 2.1e7
    #             H_aq + HO2_aq -> H2O2_aq                      : 1e7
    #             O_aq + H2O_aq -> OH_aq + OH_aq                   : 1.3e1
    #             O_aq + O2_aq -> O3_aq                         : 3e6
    #             OH_aq + OH_aq -> H2O2_aq    : 5.5e6
    #             OH_aq + Om_aq -> HO2m_aq    : 2e7
    #             OH_aq + OHm_aq -> Om_aq + H2O_aq   : 1.3e7
    #             OH_aq + HO2_aq -> H2O_aq + O2_aq   : 6e6
    #             OH_aq + O2m_aq -> OHm_aq + O2_aq     : 8e6
    #             Om_aq + H2O_aq -> OHm_aq + OH_aq     : 1.8e3
    #             Om_aq + H2O2_aq -> O2m_aq + H2O_aq   : 5e5
    #             Om_aq + HO2m_aq -> O2m_aq + OHm_aq   : 4e5
    #             Om_aq + O2_aq -> O3m_aq           : 3.6e6
    #             #Om_aq + O2m_aq + H2O_aq -> OHm_aq + OHm_aq + O2_aq   : 6e2
    #             Om_aq + O2m_aq -> OHm_aq + OHm_aq + O2_aq   : 6e5
    #             OH_aq + H2O2_aq -> H2O_aq + HO2_aq     : 2.7e4
    #             OH_aq + HO2m_aq -> OHm_aq + HO2_aq     : 7.5e6
    #             # What the fuck
    #             #H2Op_aq + H2O_aq -> OHm_aq + HO2_aq    : 6
    #             H2Op_aq + H2O_aq -> H3Op_aq + OH_aq    : 6
    #             H3Op_aq + OHm_aq -> H_aq + OH_aq + H2O_aq     : 6e7
    #             HO2_aq + H2O_aq -> H3Op_aq + O2m_aq        : 2
    #             H3Op_aq + O2m_aq -> HO2_aq + H2O_aq        : 6e-2
    #             ##################################
    #             # Additional reactions from Chen
    #             #
    #             # Some of these are from: 
    #             # Elliot, A. John and McCracken, David R. "Computer modelling 
    #             # of the radiolysis in an aqueous lithium salt blanket: 
    #             # Suppression of radiolysis by addition of hydrogen." 
    #             # Fusion Engineering and Design 13 (1990) 21-27
    #             # doi: 10.1016/0920-3796(90)90028-5 
    #             # 
    #             # Note the reactions with H2O are often given with wrong
    #             # (or confusing) rate coefficients. 
    #             # e.g. k = 1.3e10 / [H2O] - 
    #             # this means that the rate coefficient is essentially
    #             # for a two body reaction since H2O is already included
    #             #################################
    #             O_aq + O_aq -> O2_aq          : 2.8e7
    #             #em_aq + O2m_aq + H2O_aq -> HO2m_aq + OHm_aq   : 1.3e4
    #             em_aq + O2m_aq -> HO2m_aq + OHm_aq   : 1.3e7
    #             em_aq + HO2_aq -> HO2m_aq     : 2e7
    #             #em_aq + Om_aq + H2O_aq -> OHm_aq + OHm_aq     : 2.2e4
    #             # This one is listed with conflicting units in literature. 
    #             # (Three body reaction with a two body reaction rate coefficient.)
    #             em_aq + Om_aq -> OHm_aq + OHm_aq       : 2.2e7
    #             #em_aq + O3m_aq + H2O_aq -> O2_aq + OHm_aq + OHm_aq   : 1.6e4 
    #             em_aq + O3m_aq -> O2_aq + OHm_aq + OHm_aq : 1.6e7
    #             em_aq + O3_aq -> O3m_aq     : 3.6e7
    #             H_aq + Om_aq -> OHm_aq      : 1.1e7
    #             H_aq + HO2m_aq -> OHm_aq + OH_aq   : 9e7 
    #             H_aq + O3m_aq -> OHm_aq + O2_aq    : 1e7
    #             H_aq + O2m_aq -> HO2m_aq        : 1.8e7
    #             # Include HO3_aq or no?
    #             H_aq + O3_aq -> HO3_aq          : 3.8e7 
    #             OH_aq + O3m_aq ->  O3_aq + OHm_aq    : 2.6e6
    #             #OH_aq + O3m_aq -> O2_aq + O2_aq + Hp     : 6e6
    #             OH_aq + O3_aq -> HO2_aq + O2_aq          : 1.1e5
    #             HO2_aq + O2m_aq -> HO2m_aq + O2_aq       : 8e4
    #             HO2_aq + HO2_aq -> H2O2_aq + O2_aq       : 7e2
    #             HO2_aq + Om_aq -> O2_aq + OHm_aq         : 6e6
    #             HO2_aq + H2O2_aq -> OH_aq + O2_aq + H2O_aq    : 5e-4
    #             HO2_aq + HO2m_aq -> OHm_aq + OH_aq + O2_aq    : 5e-4
    #             O2m_aq + O2m_aq -> H2O2_aq + O2_aq + OHm_aq + OHm_aq : 1e-1
    #             Om_aq + O2m_aq -> OHm_aq + OHm_aq + O2_aq       : 6e5
    #             O2m_aq + H2O2_aq -> OH_aq + O2_aq + OHm_aq            : 1.3e-4
    #             O2m_aq + HO2m_aq -> Om_aq + O2_aq + OHm_aq            : 1.3e-4
    #             #O2m_aq + O3m_aq + H2O_aq -> OHm_aq + OHm_aq + O2_aq + O2_aq : 1e1 
    #             O2m_aq + O3m_aq -> OHm_aq + OHm_aq + O2_aq + O2_aq   : 1e1
    #             O2m_aq + O3_aq -> O3m_aq + O2_aq      : 1.5e6
    #             Om_aq + Om_aq -> OHm_aq + HO2m_aq    : 1e6 
    #             O3m_aq + H_aq -> O2_aq + OH_aq       : 9e7
    #             O_aq + OHm_aq -> HO2m_aq          : 1.1e2
    #             O_aq + H2O2_aq -> OH_aq + HO2_aq     : 1.6e5
    #             O_aq + HO2m_aq -> OH_aq + O2m_aq     : 5.3e6
    #             O3_aq + H2O2_aq -> OH_aq + HO2_aq + O2_aq   : 3e9
    #             #HO3_aq -> O2_aq + OH_aq           : 1e2
    #             O3m_aq + H3Op_aq -> O2_aq + OH_aq + H2O_aq      : 9e7
    #             HO3_aq -> O2_aq + OH_aq          : 1.1e5
    #             H2O2_aq -> OH_aq + OH_aq               : 4.4e-9
    #             HO2m_aq -> Om_aq + OH_aq         : 1e-5'
  []

  #[water3]
  #  name = 'test'
  #  species = 'em_aq H2_aq OHm_aq H_aq Om_aq OH_aq'
  #  aux_species = 'H2O_aq'
  #  use_log = true
  #  position_units = ${dom1Scale}
  #  track_rates = true
  #  reaction_coefficient_format = 'rate'
  #  block = 0
  #
  #  reactions = 'em_aq + em_aq -> H2_aq + OHm_aq + OHm_aq   : 5.5e6 
  #               em_aq + H_aq -> H2_aq + OHm_aq             : 2.5e7
  #               H_aq + H2O_aq -> H2_aq + OH_aq                   : 1e-2
  #               H_aq + H_aq -> H2_aq                          : 7.5e6
  #               H2_aq + H2O2_aq -> H_aq + OH_aq + H2O_aq            : 6e3
  #               OH_aq + H2_aq -> H_aq + H2O_aq     : 4.2e4
  #               Om_aq + H2_aq -> OHm_aq + H_aq       : 8e7'
  #[]
[]
