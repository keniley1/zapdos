dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  #offset = 20
  offset = 30
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  [./file]
    type = FileMeshGenerator
    file = 'air_plasma_dry_mesh.msh'
  [../]
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = file
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
  automatic_scaling = true
  compute_scaling_once = false
  end_time = 1e1
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  # petsc_options = '-snes_test_display'
  solve_type = NEWTON
  #solve_type = PJFNK
  line_search = 'basic'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  dtmin = 1e-18
  l_max_its = 20
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    # dt = 1.1
    growth_factor = 1.2
   optimal_iterations = 30
  [../]
  #[./TimeIntegrator]
  #  #type = LStableDirk2
  #  type = BDF2
  #[../]
[]

[Outputs]
  perf_graph = true
  # print_linear_residuals = false
  [./out]
    type = Exodus
    #execute_on = 'final'
  [../]
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
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
  [../]
[]

[DriftDiffusion]
  [./N2_Neutrals]
    species = 'N2A3 N2B3 N2D N2P N'
    mass = '2 2 2 2 2'
    charge = '0 0 0 0 0'
    boundaries = 'left right'
    #add_bcs = false
    skip_advection = true
    potential = potential
  [../]
  [./N_Ions]
    species = 'N+ N2+ N3+ N4+'
    mass = '2 2 2 2'
    charge = '1 1 1 1'
    boundaries = 'left right'
    potential = potential
  [../]
  [./O_Neutrals]
    species = 'O2a1 O2b1 O24_ev O O1D O1S O3'
    mass = '1 1 1 1 1 1 1' 
    charge = '0 0 0 0 0 0 0' 
    boundaries = 'left right'
    potential = potential
    skip_advection = true
  [../]
  [./O_positive_ions]
    species = 'O+ O2+ O4+'
    mass = '1 1 1'
    charge = '1 1 1'
    boundaries = 'left right'
    potential = potential
  [../]
  [./O_negative_ions]
    species = 'O- O2- O3- O4-'
    mass = '1 1 1 1'
    charge = '-1 -1 -1 -1'
    boundaries = 'left right'
    potential = potential
  [../]
  [./NxOy_Species]
    #species = 'NO N2O NO2 NO3 N2O5 H H2 OH HO2 H2O2 HNO HNO2 HNO3'
    #mass = '1 1 1 1 1 1 1 1 1 1 1 1 1' 
    #charge = '0 0 0 0 0 0 0 0 0 0 0 0 0' 
    species = 'NO N2O NO2 NO3 N2O5'
    mass = '1 1 1 1 1' 
    charge = '0 0 0 0 0' 
    boundaries = 'left right'
    potential = potential
    skip_advection = true
  [../]
  [./NyOz_positive_ions]
    #species = 'NO+ N2O+ NO2+ O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+'
    #mass = '1 1 1 1 1 1 1 1 1 1'
    #charge = '1 1 1 1 1 1 1 1 1 1'
    species = 'NO+ N2O+ NO2+ O2pN2'
    mass = '1 1 1 1'
    charge = '1 1 1 1'
    boundaries = 'left right'
    potential = potential
  [../]
  [./NyOz_negative_ions]
    #species = 'NO- N2O- NO2- NO3- H- OH-'
    #mass = '2 2 2 2 2 2'
    #charge = '-1 -1 -1 -1 -1 -1'
    species = 'NO- N2O- NO2- NO3-'
    mass = '2 2 2 2'
    charge = '-1 -1 -1 -1'
    boundaries = 'left right'
    potential = potential
  [../]
[]

[Kernels]
  [./em_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = em
    block = 0
  [../]
  [./em_advection]
    type = ADEFieldAdvection
    variable = em
    potential = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = ADCoeffDiffusion
    variable = em
    mean_en = mean_en
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    block = 0
  [../]


  [./potential_diffusion_dom1]
    type = ADCoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = em
    block = 0
  [../]

  [./mean_en_time_deriv]
    #type = ElectronTimeDerivative
    type = ADTimeDerivativeLog
    variable = mean_en
    block = 0
  [../]
  [./mean_en_advection]
    type = ADEFieldAdvection
    variable = mean_en
    potential = potential
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_diffusion]
    type = ADCoeffDiffusion
    variable = mean_en
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_joule_heating]
    #type = JouleHeating
    type = ADJouleHeating
    variable = mean_en
    potential = potential
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_log_stabilization]
    type = LogStabilizationMoles
    variable = mean_en
    block = 0
    #offset = 15
  [../]
[]

[Variables]
  [./potential]
  [../]

  [./mean_en]
    initial_condition = -22.5187
  [../]

  [./em]
    initial_condition = -22.5187
  [../]

  ########################################
  # Main air gases: N2, O2, H2O (humid)
  ########################################
  #[./N2]
  #  initial_condition = 3.48094
  #[../]

  #[./O2]
  #  # Air is N2-O2 at ratio of 4:1
  #  initial_condition = -11.49773
  #[../]

  #[./H2O]
  #  initial_condition = -14.6930
  #[../]

  ##############################
  # Other Nitrogen species
  #############################
  #[./N2v1]
  #  #initial_condition = 2.55878e10
  #  initial_condition = -30.7895
  #[../]

  [./N2A3]
    initial_condition = -30
  [../]

  [./N2B3]
    initial_condition = -30
  [../]

  [./N]
    initial_condition = -26.2055
  [../]

  [./N2D]
    initial_condition = -34.7430
  [../]

  [./N2P]
    initial_condition = -30
  [../]

  [./N+]
    initial_condition = -22.5689
  [../]

  [./N2+]
    initial_condition = -30
  [../]

  [./N3+]
    initial_condition = -30
  [../]

  [./N4+]
    initial_condition = -30
  [../]

  ######################################
  # Oxygen species
  ######################################
  #[./O2v1]
  #  initial_condition = -30.1169
  #[../]

  [./O2a1]
    initial_condition = -35.6115
  [../]

  [./O2b1]
    initial_condition = -36.4883
  [../]

  [./O24_ev]
    initial_condition = -30
  [../]

  [./O]
    initial_condition = -26.7471
  [../]

  [./O1D]
    initial_condition = -30.4032
  [../]

  [./O1S]
    initial_condition = -37.8725
  [../]

  [./O3]
    initial_condition = -30
  [../]

  [./O+]
    initial_condition = -28.9925
  [../]

  [./O2+]
    initial_condition = -26.6322
  [../]

  [./O4+]
    initial_condition = -30
  [../]

  [./O-]
    initial_condition = -30
  [../]

  [./O2-]
    initial_condition = -31.6695
  [../]

  [./O3-]
    initial_condition = -30
  [../]

  [./O4-]
    initial_condition = -35.4103
  [../]

  #####################################
  # Nx-Oy species
  #####################################
  [./NO]
    initial_condition = -28.9390
  [../]

  [./NO+]
    initial_condition = -26.7122
  [../]

  [./NO-]
    initial_condition = -30
  [../]

  [./O2pN2]
    initial_condition = -34.1682
  [../]

  [./N2O]
    initial_condition = -30
  [../]

  [./NO2]
    initial_condition = -30
  [../]

  [./NO3]
    initial_condition = -30
  [../]

  [./N2O5]
    initial_condition = -30
  [../]

  [./N2O+]
    initial_condition = -30
  [../]

  [./NO2+]
    initial_condition = -33.3563
  [../]

  [./N2O-]
    initial_condition = -30
  [../]

  [./NO2-]
    initial_condition = -30
  [../]

  [./NO3-]
    initial_condition = -30
  [../]

  ##################################
  # Hydrogen-based species
  # (For humid simulations)
  ##################################
  #[./H+]
  #  initial_condition = -30
  #[../]
  #
  #[./H2+]
  #  initial_condition = -30
  #[../]
  #
  #[./H3+]
  #  initial_condition = -30
  #[../]
  #
  #[./OH+]
  #  initial_condition = -30
  #[../]
  #
  #[./H2O+]
  #  initial_condition = -28.4108
  #[../]
  #
  #[./H3O+]
  #  initial_condition = -33.8704
  #[../]
  #
  #[./H-]
  #  initial_condition = -30
  #[../]
  #
  #[./OH-]
  #  initial_condition = -30
  #[../]
  #
  #[./H]
  #  initial_condition = -30
  #[../]
  #
  #[./H2]
  #  initial_condition = -29.5567
  #[../]
  #
  #[./OH]
  #  initial_condition = -33.5219
  #[../]
  #
  #[./HO2]
  #  initial_condition = -30
  #[../]
  #
  #[./H2O2]
  #  initial_condition = -30
  #[../]
  #
  #[./HNO]
  #  initial_condition = -30
  #[../]
  #
  #[./HNO2]
  #  initial_condition = -30
  #[../]
  #
  #[./HNO3]
  #  initial_condition = -30
  #[../]
[]

[AuxVariables]
  #[./H2O]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_condition = -0.982269  # 1% H2O
  #  block = 1
  #[../]
  [./N2]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.48094   # 80% N2, 20% O2
    #initial_condition = 3.19325   # 60% N2, 40% O2
    #initial_condition = 2.78779   # 40% N2, 60% O2
    #initial_condition = 2.09464   # 20% N2, 80% O2
  [../]
  [./O2]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 2.09464   # 80% N2, 20% O2
    #initial_condition = 2.78779   # 60% N2, 40% O2
    #initial_condition = 3.19325   # 40% N2, 60% O2
    #initial_condition = 3.48094   # 20% N2, 80% O2
  [../]
  [./NEUTRAL]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = -12.3991 
  [../]
  [./e_temp]
    block = 0
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
    block = 0
  [../]
  [./em_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Efield]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Current_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2A3_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2B3_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2D_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2P_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Np_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N3p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N4p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2a1_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2b1_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O24_ev_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O1D_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O1S_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O3_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Op_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O4p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Om_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2m_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O3m_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O4m_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2O_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO2_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO3_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2O5_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NOp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2Op_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO2p_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./O2pN2_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NOm_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./N2Om_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO2m_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./NO3m_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]


  ###########################
  # Water Species
  ###########################
  #[./H_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H2_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./OH_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./HO2_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H2O2_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./HNO_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./HNO2_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./HNO3_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H2O_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./Hp_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H2p_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H3p_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./OHp_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H2Op_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./H3Op_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./Hm_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  #[./OHm_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
[]

[AuxKernels]
  [./e_temp]
    type = ElectronTemperature
    variable = e_temp
    electron_density = em
    mean_en = mean_en
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./x_g]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./em_lin]
    type = DensityMoles
    variable = em_lin
    density_log = em
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
  [./N2A3_lin]
    type = DensityMoles
    variable = N2A3_lin
    density_log = N2A3
    block = 0
  [../]
  [./N2B3_lin]
    type = DensityMoles
    variable = N2B3_lin
    density_log = N2B3
    block = 0
  [../]
  [./N2D_lin]
    type = DensityMoles
    variable = N2D_lin
    density_log = N2D
    block = 0
  [../]
  [./N2P_lin]
    type = DensityMoles
    variable = N2P_lin
    density_log = N2P
    block = 0
  [../]
  [./N_lin]
    type = DensityMoles
    variable = N_lin
    density_log = N 
    block = 0
  [../]
  [./N2_lin]
    type = DensityMoles
    variable = N2_lin
    density_log = N2
    block = 0
  [../]
  [./Np_lin]
    type = DensityMoles
    variable = Np_lin
    density_log = N+
    block = 0
  [../]
  [./N2p_lin]
    type = DensityMoles
    variable = N2p_lin
    density_log = N2+
    block = 0
  [../]
  [./N3p_lin]
    type = DensityMoles
    variable = N3p_lin
    density_log = N3+
    block = 0
  [../]
  [./N4p_lin]
    type = DensityMoles
    variable = N4p_lin
    density_log = N4+
    block = 0
  [../]
  [./O2_lin]
    type = DensityMoles
    variable = O2_lin
    density_log = O2
    block = 0
  [../]
  [./O2a1_lin]
    type = DensityMoles
    variable = O2a1_lin
    density_log = O2a1
    block = 0
  [../]
  [./O2b1_lin]
    type = DensityMoles
    variable = O2b1_lin
    density_log = O2b1
    block = 0
  [../]
  [./O24_ev_lin]
    type = DensityMoles
    variable = O24_ev_lin
    density_log = O24_ev
    block = 0
  [../]
  [./O_lin]
    type = DensityMoles
    variable = O_lin
    density_log = O
    block = 0
  [../]
  [./O1D_lin]
    type = DensityMoles
    variable = O1D_lin
    density_log = O1D
    block = 0
  [../]
  [./O1S_lin]
    type = DensityMoles
    variable = O1S_lin
    density_log = O1S
    block = 0
  [../]
  [./O3_lin]
    type = DensityMoles
    variable = O3_lin
    density_log = O3
    block = 0
  [../]
  [./Op_lin]
    type = DensityMoles
    variable = Op_lin
    density_log = O+
    block = 0
  [../]
  [./O2p_lin]
    type = DensityMoles
    variable = O2p_lin
    density_log = O2+
    block = 0
  [../]
  [./O4p_lin]
    type = DensityMoles
    variable = O4p_lin
    density_log = O4+
    block = 0
  [../]
  [./Om_lin]
    type = DensityMoles
    variable = Om_lin
    density_log = O-
    block = 0
  [../]
  [./O2m_lin]
    type = DensityMoles
    variable = O2m_lin
    density_log = O2-
    block = 0
  [../]
  [./O3m_lin]
    type = DensityMoles
    variable = O3m_lin
    density_log = O3-
    block = 0
  [../]
  [./O4m_lin]
    type = DensityMoles
    variable = O4m_lin
    density_log = O4-
    block = 0
  [../]
  [./NO_lin]
    type = DensityMoles
    variable = 'NO_lin'
    density_log = 'NO'
    block = 0
  [../]
  [./N2O_lin]
    type = DensityMoles
    variable = N2O_lin
    density_log = N2O
    block = 0
  [../]
  [./NO2_lin]
    type = DensityMoles
    variable = NO2_lin
    density_log = NO2
    block = 0
  [../]
  [./NO3_lin]
    type = DensityMoles
    variable = NO3_lin
    density_log = NO3
    block = 0
  [../]
  [./N2O5_lin]
    type = DensityMoles
    variable = N2O5_lin
    density_log = N2O5
    block = 0
  [../]
  [./NOp_lin]
    type = DensityMoles
    variable = NOp_lin
    density_log = NO+
    block = 0
  [../]
  [./N2Op_lin]
    type = DensityMoles
    variable = N2Op_lin
    density_log = N2O+
    block = 0
  [../]
  [./NO2p_lin]
    type = DensityMoles
    variable = NO2p_lin
    density_log = NO2+
    block = 0
  [../]
  [./O2pN2_lin]
    type = DensityMoles
    variable = O2pN2_lin
    density_log = O2pN2
    block = 0
  [../]
  [./NOm_lin]
    type = DensityMoles
    variable = NOm_lin
    density_log = NO-
    block = 0
  [../]
  [./N2Om_lin]
    type = DensityMoles
    variable = N2Om_lin
    density_log = N2O-
    block = 0
  [../]
  [./NO2m_lin]
    type = DensityMoles
    variable = NO2m_lin
    density_log = NO2-
    block = 0
  [../]
  [./NO3m_lin]
    type = DensityMoles
    variable = NO3m_lin
    density_log = NO3- 
    block = 0
  [../]

  ###########################
  # Water species
  ###########################
  #[./H_lin]
  #  type = DensityMoles
  #  variable = H_lin
  #  density_log = H
  #  block = 0
  #[../]
  #[./H2_lin]
  #  type = DensityMoles
  #  variable = H2_lin
  #  density_log = H2
  #  block = 0
  #[../]
  #[./OH_lin]
  #  type = DensityMoles
  #  variable = OH_lin
  #  density_log = OH
  #  block = 0
  #[../]
  #[./HO2_lin]
  #  type = DensityMoles
  #  variable = HO2_lin
  #  density_log = HO2
  #  block = 0
  #[../]
  #[./H2O2_lin]
  #  type = DensityMoles
  #  variable = H2O2_lin
  #  density_log = H2O2
  #  block = 0
  #[../]
  #[./HNO_lin]
  #  type = DensityMoles
  #  variable = HNO_lin
  #  density_log = HNO
  #  block = 0
  #[../]
  #[./HNO2_lin]
  #  type = DensityMoles
  #  variable = HNO2_lin
  #  density_log = HNO2
  #  block = 0
  #[../]
  #[./HNO3_lin]
  #  type = DensityMoles
  #  variable = HNO3_lin
  #  density_log = HNO3
  #  block = 0
  #[../]
  #[./H2O_lin]
  #  type = DensityMoles
  #  variable = H2O_lin
  #  density_log = H2O
  #  block = 0
  #[../]
  #[./Hp_lin]
  #  type = DensityMoles
  #  variable = Hp_lin
  #  density_log = H+ 
  #  block = 0
  #[../]
  #[./H2p_lin]
  #  type = DensityMoles
  #  variable = H2p_lin
  #  density_log = H2+
  #  block = 0
  #[../]
  #[./H3p_lin]
  #  type = DensityMoles
  #  variable = H3p_lin
  #  density_log = H3+
  #  block = 0
  #[../]
  #[./OHp_lin]
  #  type = DensityMoles
  #  variable = OHp_lin
  #  density_log = OH+
  #  block = 0
  #[../]
  #[./H2Op_lin]
  #  type = DensityMoles
  #  variable = H2Op_lin
  #  density_log = H2O+
  #  block = 0
  #[../]
  #[./H3Op_lin]
  #  type = DensityMoles
  #  variable = H3Op_lin
  #  density_log = H3O+
  #  block = 0
  #[../]
  #[./Hm_lin]
  #  type = DensityMoles
  #  variable = Hm_lin
  #  density_log = H-
  #  block = 0
  #[../]
  #[./OHm_lin]
  #  type = DensityMoles
  #  variable = OHm_lin
  #  density_log = OH-
  #  block = 0
  #[../]
[]

[BCs]
  [./mean_en_physical_right]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'right'
    potential = potential
    em = em
    # Use this with HUMID simulation
    #ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH-'

    # Use this with DRY simulation
    ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2'
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    # Use this with HUMID simulation
    #ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH-'

    # Use this with DRY simulation
    ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./secondary_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    # Use this with HUMID simulation
    #ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH-'

    # Use this with DRY simulation
    ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2'
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./potential_left]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    # Use this with HUMID simulation
    #ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH-'

    # Use this with DRY simulation
    ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2'
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_dirichlet_right]
    type = ADDirichletBC
    variable = potential
    boundary = right
    value = 0
  [../]
  [./em_physical_right]
    type = ADHagelaarElectronBC
    variable = em
    boundary = 'right'
    potential = potential
    mean_en = mean_en
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
    # Use this with HUMID simulation
    #ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH-'

    # Use this with DRY simulation
    ip = 'N+ N2+ N3+ N4+ O+ O2+ O4+ O- O2- O3- O4- NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2'
    mean_en = mean_en
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
  [./potential_bc_func]
    type = ParsedFunction
    #value = 1.25
    value = 4.0
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-1.25 * (1.0001e-3 - x)'
    value = '-4 * (1.0001e-3 - x)'
  [../]
[]

[Materials]
  [./electron_moments]
    type = ADGasElectronMoments
    block = 0
    em = em
    mean_en = mean_en
    property_tables_file = 'air_plasma_dry/files_N2-80_O2-20/electron_mobility_diffusion.txt'
  [../]

  [./gas_constants]
    type = GenericConstantMaterial
    block = 0
    prop_names = ' e         N_A      diffpotential k_boltz eps  se_coeff se_energy T_gas massem   p_gas  n_gas'
    prop_values = '1.6e-19 6.022e23 8.85e-12      1.38e-23 8.854e-12 0.05     3.        300   9.11e-31 1.01e5 40.4915'
  [../]

  [./gas_species_N2]
    type = HeavySpeciesMaterial
    heavy_species_name = N2 
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5 
  [../]

  [./gas_species_N2v1]
    type = HeavySpeciesMaterial
    heavy_species_name = N2v1
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  [./gas_species_N2A3]
    type = HeavySpeciesMaterial
    heavy_species_name = N2A3
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 5.9e-5
  [../]

  [./gas_species_N2B3]
    type = HeavySpeciesMaterial
    heavy_species_name = N2B3
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 5.9e-5
  [../]

  [./gas_species_N]
    type = HeavySpeciesMaterial
    heavy_species_name = N
    heavy_species_mass = 2.325922e-26
    heavy_species_charge = 0
    diffusivity = 5.7e-5
  [../]

  [./gas_species_N2D]
    type = HeavySpeciesMaterial
    heavy_species_name = N2D
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 5.9e-5
  [../]

  [./gas_species_N2P]
    type = HeavySpeciesMaterial
    heavy_species_name = N2P
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 0
    diffusivity = 5.9e-5
  [../]

  [./gas_species_Np]
    type = HeavySpeciesMaterial
    heavy_species_name = N+
    heavy_species_mass = 2.325922e-26
    heavy_species_charge = 1
    diffusivity = 5.7e-5
  [../]

  [./gas_species_N2p]
    type = HeavySpeciesMaterial
    heavy_species_name = N2+
    heavy_species_mass = 4.651843e-26
    heavy_species_charge = 1
    diffusivity = 5.9e-5
  [../]

  [./gas_species_N3p]
    type = HeavySpeciesMaterial
    heavy_species_name = N3+
    heavy_species_mass = 6.977765e-26
    heavy_species_charge = 1
    diffusivity = 4.1e-5
  [../]

  [./gas_species_N4p]
    type = HeavySpeciesMaterial
    heavy_species_name = N4+
    heavy_species_mass = 9.303686e-26
    heavy_species_charge = 1
    diffusivity = 4.1e-5
  [../]

  [./gas_species_O2]
    type = HeavySpeciesMaterial
    heavy_species_name = O2
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  [./gas_species_O2v1]
    type = HeavySpeciesMaterial
    heavy_species_name = O2v1
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  #[./gas_species_O2v2]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = O2v2
  #  heavy_species_mass = 5.313650e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2.1e-5
  #[../]

  #[./gas_species_O2v3]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = O2v3
  #  heavy_species_mass = 5.313650e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2.1e-5
  #[../]

  #[./gas_species_O2v4]
    #type = HeavySpeciesMaterial
    #heavy_species_name = O2v4
    #heavy_species_mass = 5.313650e-26
    #heavy_species_charge = 0
    #diffusivity = 2.1e-5
  #[../]

  [./gas_species_O2a1]
    type = HeavySpeciesMaterial
    heavy_species_name = O2a1
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  [./gas_species_O2b1]
    type = HeavySpeciesMaterial
    heavy_species_name = O2b1
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  [./gas_species_O24_ev]
    type = HeavySpeciesMaterial
    heavy_species_name = O24_ev
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 0
    diffusivity = 2.1e-5
  [../]

  [./gas_species_O]
    type = HeavySpeciesMaterial
    heavy_species_name = O
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 0
    diffusivity = 6e-5
  [../]

  [./gas_species_O1D]
    type = HeavySpeciesMaterial
    heavy_species_name = O1D
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 0
    diffusivity = 6e-5
  [../]

  [./gas_species_O1S]
    type = HeavySpeciesMaterial
    heavy_species_name = O1S
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 0
    diffusivity = 6e-5
  [../]

  #########
  # Diffusion coefficient unknown
  [./gas_species_O3]
    type = HeavySpeciesMaterial
    heavy_species_name = O3
    heavy_species_mass = 7.970475e-26
    heavy_species_charge = 0
    diffusivity = 2e-5 
  [../]

  [./gas_species_Op]
    type = HeavySpeciesMaterial
    heavy_species_name = O+
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = 1
    diffusivity = 5.8e-5
  [../]

  [./gas_species_O2p]
    type = HeavySpeciesMaterial
    heavy_species_name = O2+
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = 1
    diffusivity = 5.6e-5
  [../]

  [./gas_species_O4p]
    type = HeavySpeciesMaterial
    heavy_species_name = O4+
    heavy_species_mass = 1.062730e-25
    heavy_species_charge = 1
    diffusivity = 4.1e-5
  [../]

  [./gas_species_Om]
    type = HeavySpeciesMaterial
    heavy_species_name = O-
    heavy_species_mass = 2.656825e-26
    heavy_species_charge = -1
    diffusivity = 7.0e-5
  [../]

  [./gas_species_O2m]
    type = HeavySpeciesMaterial
    heavy_species_name = O2-
    heavy_species_mass = 5.313650e-26
    heavy_species_charge = -1
    diffusivity = 5.6e-5
  [../]

  ###########
  # Diffusion coefficient unknown
  [./gas_species_O3m]
    type = HeavySpeciesMaterial
    heavy_species_name = O3- 
    heavy_species_mass = 7.970475e-26
    heavy_species_charge = -1
    diffusivity = 5.6e-5
  [../]

  ##########
  # Diffusion coefficient unknown
  [./gas_species_O4m]
    type = HeavySpeciesMaterial
    heavy_species_name = O4- 
    heavy_species_mass = 1.062730e-25
    heavy_species_charge = -1
    diffusivity = 5.6e-5
  [../]

  [./gas_species_NO]
    type = HeavySpeciesMaterial
    heavy_species_name = NO 
    heavy_species_mass = 4.982747e-26
    heavy_species_charge = 0
    diffusivity = 2e-5
  [../]

  [./gas_species_NOp] # nitrogen monoxide
    type = HeavySpeciesMaterial
    heavy_species_name = NO+ 
    heavy_species_mass = 4.982747e-26
    heavy_species_charge = 1
    diffusivity = 5.6e-5
  [../]

  [./gas_species_NOm] # nitrogen monoxide
    type = HeavySpeciesMaterial
    heavy_species_name = NO- 
    heavy_species_mass = 4.982747e-26
    heavy_species_charge = -1
    diffusivity = 5.6e-5
  [../]

  [./gas_species_O2pN2]
    type = HeavySpeciesMaterial
    heavy_species_name = O2pN2
    heavy_species_mass = 9.965493e-26
    heavy_species_charge = 1
    diffusivity = 0.5e-5
  [../]

  # Additional nitrogen-oxygen species

  [./gas_species_N2O]
    type = HeavySpeciesMaterial
    heavy_species_name = N2O
    heavy_species_mass = 7.308668e-26
    heavy_species_charge = 0
    diffusivity = 1.6e-5
  [../]

  [./gas_species_NO2]
    type = HeavySpeciesMaterial
    heavy_species_name = NO2
    heavy_species_mass = 7.639572e-26
    heavy_species_charge = 0
    diffusivity = 1.7e-5
  [../]

  [./gas_species_NO3]
    type = HeavySpeciesMaterial
    heavy_species_name = NO3
    heavy_species_mass = 1.029640e-25
    heavy_species_charge = 0
    diffusivity = 0.9e-5
  [../]

  [./gas_species_N2O5] # dinitrogen pentoxide
    type = HeavySpeciesMaterial
    heavy_species_name = N2O5
    heavy_species_mass = 1.793597e-25
    heavy_species_charge = 0
    diffusivity = 1e-5
  [../]

  [./gas_species_N2Op]
    type = HeavySpeciesMaterial
    heavy_species_name = N2O+
    heavy_species_mass = 7.308668e-26
    heavy_species_charge = 1
    diffusivity = 4.8e-5
  [../]

  [./gas_species_NO2p]
    type = HeavySpeciesMaterial
    heavy_species_name = NO2+
    heavy_species_mass = 7.639572e-26
    heavy_species_charge = 1
    diffusivity = 4.6e-5
  [../]

  [./gas_species_N2Om]
    type = HeavySpeciesMaterial
    heavy_species_name = N2O-
    heavy_species_mass = 7.308668e-26
    heavy_species_charge = -1
    diffusivity = 4.8e-5
  [../]

  [./gas_species_NO2m]
    type = HeavySpeciesMaterial
    heavy_species_name = NO2-
    heavy_species_mass = 7.639572e-26
    heavy_species_charge = -1
    diffusivity = 4.6e-5
  [../]

  [./gas_species_NO3m]
    type = HeavySpeciesMaterial
    heavy_species_name = NO3-
    heavy_species_mass = 1.029640e-25
    heavy_species_charge = -1
    diffusivity = 0.9e-5 
  [../]

  # Hydrogen species
  [./gas_species_Hp]
    type = HeavySpeciesMaterial
    heavy_species_name = H+
    heavy_species_mass = 1.673597e-27
    heavy_species_charge = 1
    diffusivity = 8.8e-5
  [../]

  [./gas_species_H2p]
    type = HeavySpeciesMaterial
    heavy_species_name = H2+
    heavy_species_mass = 3.347526e-27
    heavy_species_charge = 1
    diffusivity = 7e-5
  [../]

  [./gas_species_H3p]
    type = HeavySpeciesMaterial
    heavy_species_name = H3+
    heavy_species_mass = 5.021289e-27
    heavy_species_charge = 1
    diffusivity = 9e-5
  [../]

  [./gas_species_OHp]
    type = HeavySpeciesMaterial
    heavy_species_name = OH+
    heavy_species_mass = 2.824311e-26
    heavy_species_charge = 1
    diffusivity = 4e-5
  [../]

  #[./gas_species_H2Op]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H2O+
  #  heavy_species_mass = 2.988e-26 
  #  heavy_species_charge = 1
  #  diffusivity = 5.9e-5
  #[../]
  #
  #[./gas_species_H3Op]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H3O+
  #  heavy_species_mass = 3.158951e-26 
  #  heavy_species_charge = 1
  #  diffusivity = 6e-5
  #[../]
  #
  #[./gas_species_Hm]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H-
  #  heavy_species_mass = 1.673597e-27
  #  heavy_species_charge = -1
  #  diffusivity = 8.8e-5
  #[../]
  #
  #[./gas_species_OHm]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = OH-
  #  heavy_species_mass = 2.824311e-26
  #  heavy_species_charge = -1
  #  diffusivity = 7e-5
  #[../]
  #
  #[./gas_species_H]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H
  #  heavy_species_mass = 1.673597e-27
  #  heavy_species_charge = 0
  #  diffusivity = 8.8e-5
  #[../]
  #
  #[./gas_species_H2]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H2
  #  heavy_species_mass = 3.347526e-27
  #  heavy_species_charge = 0
  #  diffusivity = 7.8e-5
  #[../]
  #
  #[./gas_species_OH]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = OH
  #  heavy_species_mass = 2.824311e-26
  #  heavy_species_charge = 0
  #  diffusivity = 4e-5
  #[../]
  #
  #[./gas_species_HO2]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = HO2
  #  heavy_species_mass = 5.481069e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2e-5
  #[../]
  #
  #[./gas_species_H2O2]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H2O2
  #  heavy_species_mass = 5.64829e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2e-5
  #[../]
  #
  #[./gas_species_HNO] # Nitroxyl
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = HNO
  #  heavy_species_mass = 5.150116e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2.1e-5
  #[../]
  #
  #[./gas_species_HNO2] # Nitrous acid
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = HNO2
  #  heavy_species_mass = 7.8087e-26
  #  heavy_species_charge = 0
  #  diffusivity = 2.1e-5
  #[../]
  #
  #[./gas_species_HNO3] # Nitric acid
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = HNO3
  #  heavy_species_mass = 1.04633e-25
  #  heavy_species_charge = 0
  #  diffusivity = 2.1e-5
  #[../]
  #
  #[./gas_species_H2O]
  #  type = HeavySpeciesMaterial
  #  heavy_species_name = H2O
  #  heavy_species_mass = 2.988e-26
  #  heavy_species_charge = 0 
  #  diffusivity = 2.3e-5 
  #[../]
[]

[Reactions]
  # This argon reaction network based on a ZDPlasKin example:
  # zdplaskin.laplace.univ-tlse.fr
  # Example: "Micro-cathode sustained discharged in Ar"

  [./DryAir]
    name = 'dry'
    #species = 'N2 N2v1 N2A3 N2B3 N N2D N2P N+ N2+ N3+ N4+ O2 O2v1 O2a1 O2b1 O24_ev O O1D O1S O3 O+ O2+ O4+ O- O2- O3- O4- NO N2O NO2 NO3 N2O5 NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 H+ H2+ H3+ OH+ H2O+ H3O+ H- OH- H H2 OH HO2 H2O2 HNO HNO2 HNO3 H2O e'
    species = 'N2A3 N2B3 N N2D N2P N+ N2+ N3+ N4+ O2a1 O2b1 O24_ev O O1D O1S O3 O+ O2+ O4+ O- O2- O3- O4- NO N2O NO2 NO3 N2O5 NO+ N2O+ NO2+ NO- N2O- NO2- NO3- O2pN2 em'
    aux_species = 'N2 O2'
    reaction_coefficient_format = 'townsend'
    gas_species = 'N2'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'air_plasma_dry/files_N2-80_O2-20'
    equation_constants = 'Tgas Te TeffN TeffN2 TeffN3 TeffN4 TionN TionN2 TionN3 TionN4'
    equation_values = '300 30933.33 778.875 659.1565 587.3252 539.438 1018.383 1018.383 1018.383 1018.383'  # Te = 2/3 * mean_energy * 11600 
    potential = 'potential'
    use_log = true
    position_units = ${dom0Scale}
    #track_rates = true
    use_ad = true
    convert_to_moles = true
    convert_to_meters = 1e-2
    block = 0


    reactions = 'em + N2 -> em + N2                    : EEDF [elastic] (C1_N2_Elastic)
                 em + N2O -> em + N2O                  : EEDF [elastic] (C74_N2O_Elastic)
                 em + N2A3 -> em + N2A3                : EEDF [elastic] (C42_N2A3_Effective_momentum)
                 #em + N2 -> em + N2v1                 : EEDF (C2_N2_Excitation_0.0000_eV)
                 #em + N2 -> em + N2v1                 : EEDF (C3_N2_Excitation_0.29_eV) 
                 #em + N2v1 -> em + N2                 : EEDF [0.29] (C50_N2v1_De-excitation_0.29_eV)
                 #em + O2 -> em + O2v1                 : EEDF [-0.19] (C27_O2_Excitation_0.19_eV)
                 #em + O2 -> em + O2v1                 : EEDF [-0.19] (C28_O2_Excitation_0.19_eV)
                 #em + O2 -> em + O2v2                 : EEDF (C29_O2_Excitation_0.38_eV)
                 #em + O2 -> em + O2v2                 : EEDF (C30_O2_Excitation_0.38_eV)
                 #em + O2 -> em + O2v3                 : EEDF (C31_O2_Excitation_0.57_eV)
                 #em + O2 -> em + O2v4                 : EEDF (C32_O2_Excitation_0.75_eV)
                 #em + O2v1 -> em + O2                 : EEDF [0.19] (C54_O2v1_De-excitation_0.19_eV)
                 #em + O2v2 -> em + O2                 : EEDF (C55_O2v2_De-excitation_0.38_eV)
                 #em + O2v3 -> em + O2                 : EEDF (C56_O2v3_De-excitation_0.57_eV)
                 #em + O2v4 -> em + O2                 : EEDF (C57_O2v4_De-excitation_0.75_eV)
                 ######
                 # Vibrational-translational relaxation (Capitelli2000, page 105)
                 ######
                 #N2v1 + N2 -> N2 + N2               : {7.8e-12 * Tgas * exp(-218.0 / Tgas^(1./3.) + 690. / Tgas) / (1.0 - exp(-0.290*11605/Tgas))}
                 #N2 + N2 -> N2v1 + N2               : {7.8e-12 * Tgas * exp(-218.0 / Tgas^(1./3.) + 690. / Tgas) / (1.0 - exp(-0.290*11605/Tgas)) * exp(-0.29*11605/Tgas)}
                 #N2v1 + N -> N2 + N                 : {4.0e-16 * (Tgas / 300.0)^0.5}
                 #N2 + N -> N2v1 + N                 : {4.0e-16 * (Tgas / 300.0)^0.5 * exp(-0.29*11605.0)}
                 #N2v1 + O -> N2 + O                 : {1.20e-13 * exp( - 27.6 / Tgas^(1.0/3.0))}
                 #N2 + O -> N2v1 + O                 : {1.20e-13 * exp( - 27.6 / Tgas^(1.0/3.0)) * exp(-0.28*11605.0)}
                 #O2v1 + O2 -> O2 + O2               : {1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) )}
                 #O2v2 + O2 -> O2v1 + O2             : {2.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) )}
                 #O2v3 + O2 -> O2v2 + O2             : {3.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) )}
                 #O2v4 + O2 -> O2v3 + O2             : {4.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) )}
                 #O2 + O2 -> O2v1 + O2               : {1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) ) * exp(-0.19 * 11605)}
                 #O2v1 + O2 -> O2v2 + O2             : {2.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) ) * exp(-0.19 * 11605)}
                 #O2v2 + O2 -> O2v3 + O2             : {3.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) ) * exp(-0.19 * 11605)}
                 #O2v3 + O2 -> O2v4 + O2             : {4.0 * 1.35e-12 * Tgas * exp( - 137.9 / Tgas^(1.0/3.0) ) / ( 1.0 - exp(-0.19*11605.0) ) * exp(-0.19 * 11605)}
                 #O2v1 + O -> O2 + O                 : {4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v2 + O -> O2v1 + O               : {2.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v3 + O -> O2v2 + O               : {3.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v4 + O -> O2v3 + O               : {4.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2 + O -> O2v1 + O                 : {4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v1 + O -> O2v2 + O               : {2.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v2 + O -> O2v3 + O               : {3.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #O2v3 + O -> O2v4 + O               : {4.0 * 4.5e-15 * Tgas * exp(-0.19 * 11605)}
                 #####
                 # excitation of electronic levels by electron impact (Bolsig+)
                 # Note that CRANE will need to be modified to allow duplicate reactions here...
                 #####
                 em + N2 -> em + N2A3                 : EEDF [-6.17] (C11_N2_Excitation_6.17_eV)
                 em + N2 -> em + N2A3                 : EEDF [-7.0] (C12_N2_Excitation_7.00_eV)
                 em + N2 -> em + N2A3                 : EEDF [-7.8] (C13_N2_Excitation_7.80_eV)
                 em + N2 -> em + N2B3                 : EEDF [-7.35] (C14_N2_Excitation_7.35_eV)
                 em + N2 -> em + N2B3                 : EEDF [-7.36] (C15_N2_Excitation_7.36_eV)
                 em + N2 -> em + N2B3                 : EEDF [-8.16] (C16_N2_Excitation_8.16_eV)
                 em + N2 -> em + N + N2D              : EEDF [-13.0] (C23_N2_Excitation_13.00_eV)
                 em + O2 -> em + O2a1                 : EEDF [-0.98] (C33_O2_Excitation_0.98_eV)
                 em + O2 -> em + O2b1                 : EEDF [-1.63] (C34_O2_Excitation_1.63_eV)
                 em + O2 -> em + O24_ev               : EEDF [-4.5] (C35_O2_Excitation_4.50_eV)
                 em + O2 -> em + O + O                : EEDF [-6.0] (C36_O2_Excitation_6.00_eV)
                 em + O2 -> em + O + O1D              : EEDF [-8.4] (C37_O2_Excitation_8.40_eV)
                 em + O2 -> em + O + O1S              : EEDF [-9.97] (C38_O2_Excitation_9.97_eV)
                 em + O2a1 -> em + O + O              : EEDF [-6.0] (C51_O2a1_Excitation_6.00_eV)
                 em + O -> em + O1D                   : EEDF [-1.97] (C56_O_Excitation_1.97_eV)
                 em + O -> em + O1S                   : EEDF [-4.19] (C57_O_Excitation_4.19_eV)
                 ####
                 # de-excitation of electronic levels by electron impact (BOLSIG+)
                 ####
                 em + N2A3 -> em + N2                 : EEDF [6.17] (C41_N2A3_De-excitation_6.17_eV)
                 em + O2a1 -> em + O2                : EEDF [0.98] (C48_O2a1_De-excitation_0.98_eV)
                 ####
                 # ionization by electron impact (BOLSIG+)
                 # note the missing data section in the inp file (4 extra reactions - not shown here)
                 ####
                 em + N -> em + em + N+                : EEDF [-14.55] (C54_N_Ionization_14.55_eV)
                 em + O -> em + em + O+                : EEDF [-13.62] (C58_O_Ionization_13.62_eV)
                 em + N2 -> em + em + N2+              : EEDF [-15.6] (C24_N2_Ionization_15.60_eV)
                 em + N2A3 -> em + em + N2+            : EEDF [-10.79] (C43_N2A3_Ionization_10.79_eV)
                 em + O2 -> em + em + O2+              : EEDF [-12.06] (C39_O2_Ionization_12.06_eV)
                 em + O2a1 -> em + em + O2+            : EEDF [-11.0] (C52_O2a1_Ionization_11.00_eV)
                 em + NO -> em + em + NO+              : EEDF [-9.26] (C68_NO_Ionization_9.26_eV)
                 ####
                 # electron-ion recombination
                 ####
                 em + N2+ -> N + N                   : {1.8e-7 * (300/Te)^0.39 * 0.5}
                 em + N2+ -> N + N2D                 : {1.8e-7 * (300/Te)^0.39 * 0.45}
                 em + N2+ -> N + N2P                 : {1.8e-7 * (300/Te)^0.39 * 0.05}
                 em + O2+ -> O + O                   : {2.7e-7 * (300/Te)^0.7 * 0.55}
                 em + O2+ -> O + O1D                 : {2.7e-7 * (300/Te)^0.7 * 0.40}
                 em + O2+ -> O + O1S                 : {2.7e-7 * (300/Te)^0.7 * 0.05}
                 em + NO+ -> O + N                   : {4.2e-7 * (300/Te)^0.85 * 0.2}
                 em + NO+ -> O + N2D                 : {4.2e-7 * (300/Te)^0.85 * 0.8}
                 em + N3+ -> N2 + N                  : {2.0e-7 * (300/Te)^0.5}
                 em + N4+ -> N2 + N2                 : {2.3e-6 * (300/Te)^0.53}
                 em + O4+ -> O2 + O2                 : {1.4e-6 * (300/Te)^0.5}
                 em + O2pN2 -> O2 + N2               : {1.3e-6 * (300/Te)^0.5}
                 em + em + N+ -> N + em                    : {7.0e-20 * (300/Te)^4.5}
                 em + em + O+ -> O + em                    : {7.0e-20 * (300/Te)^4.5}
                 em + N+ + NEUTRAL -> N + NEUTRAL    : {6.0e-27 * (300/Te)^1.5}
                 em + O+ + NEUTRAL -> O + NEUTRAL    : {6.0e-27 * (300/Te)^1.5}
                 ####
                 # electron attachment
                 # em + O2 + O2 -> O2- + O2              : EEDF (
                 # ^ Not sure what to do with this one. Never seen anything like it before.
                 ####
                 em + O2 -> O- + O                     : EEDF (C25_O2_Attachment)
                 em + NO -> O- + N                     : EEDF (C59_NO_Attachment)
                 em + O3 -> O- + O2                    : EEDF (C69_O3_Attachment)
                 em + O3 -> O2- + O                    : EEDF (C70_O3_Attachment)
                 em + O + O2 -> O- + O2                : 1e-31
                 em + O + O2 -> O2- + O                : 1e-31
                 em + O3 + NEUTRAL -> O3- + NEUTRAL    : 1e-31
                 em + NO + NEUTRAL -> NO- + NEUTRAL    : 1e-31
                 em + O2 + N2 -> O2- + N2              : {1.1e-31 * (300/Te)^2 * exp(-70/Tgas) * exp(1500*(Te-Tgas)/(Te*Tgas))}
                 ####
                 # electron detachment
                 ####
                 O-  + O -> O2  + em                   : 1.4e-10
                 O-  + N -> NO  + em                   : 2.6e-10
                 O-  + O2 -> O3  + em                  : 5.0e-15
                 O-  + O2a1 -> O3  + em                : 3.0e-10
                 O-  + O2b1 -> O   + O2 + em           : 6.9e-10
                 O-  + N2A3 -> O   + N2 + em           : 2.2e-9
                 O-  + N2B3 -> O   + N2 + em           : 1.9e-9
                 O-  + O3 -> O2  + O2 + em             : 3.0e-10
                 O2- + O -> O3  + em                   : 1.5e-10
                 O2- + O2 -> O2  + O2 + em             : {2.7e-10 * (TeffN2/300)^0.5 * exp(-5590/TeffN2)}
                 O2- + O2a1 -> O2  + O2 + em           : 2.0e-10
                 O2- + O2b1 -> O2  + O2 + em           : 3.6e-10
                 O2- + N2 -> O2  + N2 + em             : {1.9e-12 * (TeffN2/300)^0.5 * exp(-4990/TeffN2)}
                 O2- + N2A3 -> O2  + N2 + em           : 2.1e-9
                 O2- + N2B3 -> O2  + N2 + em           : 2.5e-9
                 O3- + O -> O2  + O2 + em              : 3.0e-10
                 ####
                 # Detachment for O3- NO- N2O- NO2- NO3- has to be verified (from inp source file)
                 ####
                 O3- + N -> NO + O2 + em               : 5e-10
                 O3- + N2A3 -> O3 + N2 + em            : 2.1e-9
                 NO- + N2A3 -> NO + N2 + em            : 2.1e-9
                 O3- + N2B3 -> O3 + N2 + em            : 2.5e-9
                 NO- + N2B3 -> NO + N2 + em            : 2.5e-9
                 ####
                 # optical transitions and predissociation (Capitelli2000, page 157)
                 ####
                 N2A3 -> N2                           : 0.5 
                 N2B3 -> N2A3                         : 1.34e5
                 O2a1 -> O2                           : 2.6e-4
                 O2b1 -> O2a1                         : 1.5e-3
                 O2b1 -> O2                           : 8.5e-2
                 O24_ev -> O2                         : 11
                 ####
                 # quenching and excitation of N2 (Capitelli2000, page 159)
                 ####
                 N2A3 + O -> NO + N2D                 : 7e-12
                 N2A3 + O -> NO + O1S                 : 2.1e-11
                 N2A3 + N -> N2 + N                   : 2e-12
                 N2A3 + N -> N2 + N2P                 : {4.0e-11 * (300/Tgas)^0.667}
                 N2A3 + O2 -> N2 + O + O1D            : {2.1e-12 * (Tgas/300)^0.55}
                 N2A3 + O2 -> N2 + O2a1               : {2.0e-13 * (Tgas/300)^0.55}
                 N2A3 + O2 -> N2 + O2b1               : {2.0e-13 * (Tgas/300)^0.55}
                 N2A3 + N2 -> N2 + N2                 : 3e-16 
                 N2A3 + NO -> N2 + NO                 : 6.9e-11
                 N2A3 + N2A3 -> N2 + N2B3             : 3e-10
                 N2B3 + N2 -> N2A3 + N2               : 3e-11
                 N2B3 + N2 -> N2 + N2                 : 2e-12
                 N2B3 + O2 -> N2 + O + O              : 3e-10
                 N2B3 + NO -> N2A3 + NO               : 2.4e-10
                 N + N + N2 -> N2A3 + N2              : 1.7e-33
                 N + N + O2 -> N2A3 + O2              : 1.7e-33
                 N + N + NO -> N2A3 + NO              : 1.7e-33
                 N + N + N -> N2A3 + N                : 1e-32
                 N + N + O -> N2A3 + O                : 1e-32
                 N + N + N2 -> N2B3 + N2              : 2.4e-33
                 N + N + O2 -> N2B3 + O2              : 2.4e-33
                 N + N + NO -> N2B3 + NO              : 2.4e-33
                 N + N + N -> N2B3 + N                : 1.4e-32
                 N + N + O -> N2B3 + O                : 1.4e-32
                 ####
                 # deactivation of N metastables (Capitelli2000, page 161)
                 ####
                 N2D + O -> N + O1D                   : 4e-13
                 N2D + O2 -> NO + O                   : 5.2e-12
                 N2D + NO -> N2 + O                   : 1.8e-10
                 N2D + N2 -> N + N2                   : {1.0e-13 * exp(-510/Tgas)}
                 N2P + N -> N + N                     : 1.8e-12
                 N2P + O -> N + O                     : 1.0e-12
                 N2P + N -> N2D + N                   : 6.0e-13
                 N2P + N2 -> N + N2                   : 6.0e-14
                 N2P + N2D -> N2+ + em                 : 1.0e-13
                 N2P + O2 -> NO + O                   : 2.6e-12
                 N2P + NO -> N2A3 + O                 : 3.0e-11
                 ####
                 # quenching and excitation of O2 (Capitelli2000, page 160)
                 ####
                 O2a1 + O -> O2 + O                   : 7.0e-16
                 O2a1 + N -> NO + O                   : {2.0e-14 * exp(-600/Tgas)}
                 O2a1 + O2 -> O2 + O2                 : {3.8e-18 * exp(-205/Tgas)}
                 O2a1 + N2 -> O2 + N2                 : 3.0e-21
                 O2a1 + NO -> O2 + NO                 : 2.5e-11
                 O2a1 + O3 -> O2 + O2 + O1D           : {5.2e-11 * exp(-2840/Tgas)}
                 O2a1 + O2a1 -> O2 + O2b1             : {7.0e-28 * Tgas^3.8 * exp(700/Tgas)}
                 O + O3 -> O2 + O2a1                  : {1.0e-11 * exp(-2300/Tgas)}
                 O2b1 + O -> O2a1 + O                 : 8.1e-14
                 O2b1 + O -> O2 + O1D                 : {3.4e-11 * (300/Tgas)^0.1 * exp(-4200/Tgas)}
                 O2b1 + O2 -> O2a1 + O2               : {4.3e-22 * Tgas^2.4 * exp(-281/Tgas)}
                 O2b1 + N2 -> O2a1 + N2               : {1.7e-15 * (Tgas/300)}
                 O2b1 + NO -> O2a1 + NO               : 6.0e-14
                 O2b1 + O3 -> O2 + O2 + O             : 2.2e-11
                 O24_ev + O -> O2 + O1S               : 9.0e-12
                 O24_ev + O2 -> O2b1 + O2b1           : 3.0e-13
                 O24_ev + N2 -> O2b1 + N2             : 9.0e-15
                 ####
                 # deactivation of O metastables (Capitelli2000, page 161)
                 ####
                 O1D + O -> O + O                     : 8.0e-12
                 O1D + O2 -> O + O2                   : {6.4e-12 * exp(67/Tgas)}
                 O1D + O2 -> O + O2a1                 : 1.0e-12
                 O1D + O2 -> O + O2b1                 : {2.6e-11 * exp(67/Tgas)}
                 O1D + N2 -> O + N2                   : 2.3e-11
                 O1D + O3 -> O2 + O + O               : 1.2e-10
                 O1D + O3 -> O2 + O2                  : 1.2e-10
                 O1D + NO -> O2 + N                   : 1.7e-10
                 O1S + O -> O1D + O                   : 5.0e-11 * exp(-300/Tgas)}
                 O1S + N -> O + N                     : 1.0e-12
                 O1S + O2 -> O1D + O2                 : {1.3e-12 * exp(-850/Tgas)}
                 O1S + O2 -> O + O + O                : {3.0e-12 * exp(-850/Tgas)}
                 O1S + N2 -> O + N2                   : 1.0e-17
                 O1S + O2a1 -> O + O24_ev             : 1.1e-10
                 O1S + O2a1 -> O1D + O2b1             : 2.9e-11
                 O1S + O2a1 -> O + O + O              : 3.2e-11
                 O1S + NO -> O + NO                   : 2.9e-10
                 O1S + NO -> O1D + NO                 : 5.1e-10
                 O1S + O3 -> O2 + O2                  : 2.9e-10
                 O1S + O3 -> O2 + O + O1D             : 2.9e-10
                 ####
                 # bimolecular nitrogen-oxygen reactions (Capitelli2000, page 168)
                 # Two missing reactions: 
                 # N + O3 -> NO + O2 : < 2.0e-16
                 # O + N2O5 -> product : < 3.0e-16
                 ####
                 N + NO -> O + N2                     : {1.8e-11 * (Tgas/300.0)^0.5}
                 N + O2 -> O + NO                     : {3.2e-12 * (Tgas/300.0) * exp(-3150/Tgas)}
                 O + N2 -> N + NO                     : {3.0e-10 * exp(-38370/Tgas)}
                 O + NO -> N + O2                     : {7.5e-12 * (Tgas/300.0) * exp(-19500/Tgas)}
                 NO + NO -> N2 + O2                   : {5.1e-13 * exp(-33660/Tgas)}
                 O2 + O2 -> O + O3                    : {2.0e-11 * exp(-49800/Tgas)}
                 N + N -> N2+ + em                     : {2.7e-11 * exp(-6.74e4/Tgas)}
                 N + O -> NO+ + em                     : {1.6e-12 * (Tgas/300)^0.5 * (0.19+8.6*Tgas) * exp(-32000/Tgas)}
                 ####
                 # dissociation of nitrogen-oxygen molecules (Capitelli2000, page 169)
                 ####
                 N2 + N2 -> N + N + N2                : {5.4e-8 * (1.0 - exp(-3354/Tgas)) * exp(-113200/Tgas)}
                 O2 + N2 -> N + N + O2                : {5.4e-8 * (1.0 - exp(-3354/Tgas)) * exp(-113200/Tgas)}
                 NO + N2 -> N + N + NO                : {5.4e-8 * (1.0 - exp(-3354/Tgas)) * exp(-113200/Tgas)}
                 O + N2 -> N + N + O                  : {5.4e-8 * (1.0 - exp(-3354/Tgas)) * exp(-113200/Tgas) * 6.6}
                 N + N2 -> N + N + N                  : {5.4e-8 * (1.0 - exp(-3354/Tgas)) * exp(-113200/Tgas) * 6.6}
                 O2 + N2 -> O + O + N2                : {6.1e-9 * (1.0 - exp(-2240/Tgas)) * exp(-59380/Tgas)}
                 O2 + O2 -> O + O + O2                : {6.1e-9 * (1.0 - exp(-2240/Tgas)) * exp(-59380/Tgas) * 5.9}
                 O2 + O -> O + O + O                  : {6.1e-9 * (1.0 - exp(-2240/Tgas)) * exp(-59380/Tgas) * 21}
                 O2 + N -> O + O + N                  : {6.1e-9 * (1.0 - exp(-2240/Tgas)) * exp(-59380/Tgas)}
                 O2 + NO -> O + O + NO                : {6.1e-9 * (1.0 - exp(-2240/Tgas)) * exp(-59380/Tgas)}
                 NO + N2 -> N + O + N2                : {8.7e-9 * exp(-75994/Tgas)}
                 NO + O2 -> N + O + O2                : {8.7e-9 * exp(-75994/Tgas)}
                 NO + O -> N + O + O                  : {8.7e-9 * exp(-75994/Tgas) * 20}
                 NO + N -> N + O + N                  : {8.7e-9 * exp(-75994/Tgas) * 20}
                 NO + NO -> N + O + NO                : {8.7e-9 * exp(-75994/Tgas) * 20}
                 O3 + N2 -> O2 + O + N2               : {6.6e-10 * exp(-11600/Tgas)}
                 O3 + O2 -> O2 + O + O2               : {6.6e-10 * exp(-11600/Tgas) * 0.38}
                 O3 + N -> O2 + O + N                 : {6.6e-10 * exp(-11600/Tgas) * 6.3*exp(170/Tgas)}
                 O3 + O -> O2 + O + O                 : {6.6e-10 * exp(-11600/Tgas) * 6.3*exp(170/Tgas)}
                 ####
                 # recombination of nitrogen-oxygen molecules (Capitelli, page 170)
                 # note "max" rate coefficients in the source file.
                 # Do not know how to implement something similar in CRANE.
                 ####
                 N + N + N2 -> N2 + N2                : {8.3e-34 * exp(500/Tgas)}
                 N + N + O2 -> N2 + O2                : {1.8e-33 * exp(435/Tgas)}
                 N + N + NO -> N2 + NO                : {1.8e-33 * exp(435/Tgas)}
                 N + N + N -> N2 + N                  : {1.8e-33 * exp(435/Tgas) * 3}
                 N + N + O -> N2 + O                  : {1.8e-33 * exp(435/Tgas) * 3}
                 O + O + N2 -> O2 + N2                : {2.8e-34 * exp(720/Tgas)}
                 O + O + O2 -> O2 + O2                : {4.0e-33 * (300/Tgas)^0.41}
                 O + O + N -> O2 + N                  : {4.0e-33 * (300/Tgas)^0.41 * 0.8}
                 O + O + O -> O2 + O                  : {4.0e-33 * (300/Tgas)^0.41 * 3.6}
                 O + O + NO -> O2 + NO                : {4.0e-33 * (300/Tgas)^0.41 * 0.17}
                 N + O + N2 -> NO + N2                : {1.0e-32 * (300/Tgas)^0.5}
                 N + O + O2 -> NO + O2                : {1.0e-32 * (300/Tgas)^0.5}
                 N + O + N -> NO + N                  : {1.8e-31 * (300/Tgas)}
                 N + O + O -> NO + O                  : {1.8e-31 * (300/Tgas)}
                 N + O + NO -> NO + NO                : {1.8e-31 * (300/Tgas)}
                 O + O2 + N2 -> O3 + N2               : {5.8e-34 * (300/Tgas)^2.8}
                 O + O2 + O2 -> O3 + O2               : {7.6e-34 * (300/Tgas)^1.9}
                 O + O2 + NO -> O3 + NO               : {7.6e-34 * (300/Tgas)^1.9}
                 O + O2 + N -> O3 + N                 : {3.9e-33 * (300/Tgas)^1.9}
                 O + O2 + O -> O3 + O                 : {3.9e-33 * (300/Tgas)^1.9}
                 ####
                 # positive ion reactions (Capitelli, 179)
                 ####
                 N+ + O   -> N + O+                              : 1.0e-12
                 N+ + O2  -> O2+ + N                              : 2.8e-10
                 N+ + O2  -> NO+ + O                              : 2.5e-10
                 N+ + O2  -> O+ + NO                              : 2.8e-11
                 N+ + O3  -> NO+ + O2                             : 5.0e-10
                 N+ + NO  -> NO+ + N                              : 8.0e-10
                 N+ + NO  -> N2+ + O                              : 3.0e-12
                 N+ + NO  -> O+ + N2                              : 1.0e-12
                 O+ + N2 -> NO+ + N                               : {( 1.5 - 2.0e-3*TeffN + 9.6e-7*TeffN^2 ) * 1.0e-12}
                 O+ + O2 -> O2+ + O                               : {2.0e-11 * (300/TeffN)^0.5}
                 O+ + O3 -> O2+ + O2                              : 1.0e-10
                 O+ + NO -> NO+ + O                               : 2.4e-11
                 O+ + NO -> O2+ + N                               : 3.0e-12
                 O+ + N2D -> N+ + O                               : 1.3e-10
                 N2+ + O2 -> O2+ + N2                             : {6.0e-11 * (300/TeffN2)^0.5}
                 N2+ + O  -> NO+ + N                              : {1.3e-10 * (300/TeffN2)^0.5}
                 N2+ + O3 -> O2+ + O + N2                         : 1.0e-10
                 N2+ + N  -> N+ + N2                              : {7.2e-13 * (TeffN2/300)}
                 N2+ + NO -> NO+ + N2                             : 3.3e-10
                 O2+ + N2 -> NO+ + NO                             : 1.0e-17
                 O2+ + N  -> NO+ + O                              : 1.2e-10
                 O2+ + NO -> NO+ + O2                             : 6.3e-10
                 N3+ + O2 -> O2+ + N + N2                         : 2.3e-11
                 N3+ + N  -> N2+ + N2                             : 6.6e-11
                 N3+ + NO -> NO+ + N + N2                         : 7.0e-11
                 N4+ + N2 -> N2+ + N2 + N2                        : 1.0e-10
                 N4+ + O2 -> O2+ + N2 + N2                        : 2.5e-10
                 N4+ + O  -> O+ + N2 + N2                         : 2.5e-10
                 N4+ + N  -> N+ + N2 + N2                         : 1.0e-11
                 N4+ + NO -> NO+ + N2 + N2                        : 4.0e-10
                 O4+ + N2 -> O2pN2 + O2                           : {4.6e-12 * (TeffN4/300)^2.5 * exp(-2650/TeffN4)}
                 O4+ + O2 -> O2+ + O2 + O2                        : {3.3e-6  * (300/TeffN4)^4   * exp(-5030/TeffN4)}
                 O4+ + O2a1 -> O2+ + O2 + O2                      : 1.0e-10
                 O4+ + O2b1 -> O2+ + O2 + O2                      : 1.0e-10
                 O4+ + O      -> O2+ + O3                         : 3.0e-10
                 O4+ + NO     -> NO+ + O2 + O2                    : 1.0e-10
                 O2pN2 + N2 -> O2+ + N2 + N2                      : {1.1e-6 * (300/TeffN4)^5.3 * exp(-2360/TeffN4)}
                 O2pN2 + O2 -> O4+ + N2                           : 1.0e-9
                 N+ + N2 + N2 -> N3+ + N2                         : {1.7e-29 * (300/TeffN)^2.1}
                 N+ + O + NEUTRAL -> NO+ + NEUTRAL        : 1.0e-29
                 N+ + N + NEUTRAL -> N2+ + NEUTRAL        : 1.0e-29
                 O+ + N2 + NEUTRAL -> NO+ + N + NEUTRAL   : {6.0e-29 * (300/TeffN)^2}
                 O+ + O  + NEUTRAL -> O2+ + NEUTRAL       : 1.0e-29
                 O+ + N  + NEUTRAL -> NO+ + NEUTRAL       : 1.0e-29
                 N2+ + N2 + N2 -> N4+ + N2                        : {5.2e-29 * (300/TeffN2)^2.2}
                 N2+ + N  + N2 -> N3+ + N2                        : {9.0e-30 * exp(400/TeffN2)}
                 O2+ + O2 + O2 -> O4+ + O2                        : {2.4e-30 * (300/TeffN2)^3.2}
                 O2+ + N2 + N2 -> O2pN2 + N2                      : {9.0e-31 * (300/TeffN2)^2}
                 ####
                 # negative ion reactions (Capitelli, 182-183)
                 # NOTE missing reactions: 
                 # O2^- + N2O    => O3^-  + N2                       ! < 1.0d-12
                 #O3^- + N2     => NO2^- + NO                       ! < 5.0d-14
                 ####
                 O-   + O2a1 -> O2-  + O                          : 1.0e-10
                 O-   + O3 -> O3-  + O                            : 8.0e-10
                 O2-  + O -> O-   + O2                            : 3.3e-10
                 O2-  + O3 -> O3-  + O2                           : 3.5e-10
                 O3-  + O -> O2-  + O2                            : 1.0e-11
                 NO-  + O2 -> O2-  + NO                           : 5.0e-10
                 O4- + N2 -> O2- + O2 + N2                        : {1.0e-10 * exp(-1044/TeffN4)}
                 O4- + O2 -> O2- + O2 + O2                        : {1.0e-10 * exp(-1044/TeffN4)}
                 O4- + O -> O3-  + O2                             : 4.0e-10
                 O4- + O -> O-   + O2  + O2                       : 3.0e-10
                 O4- + O2a1 -> O2-  + O2  + O2                    : 1.0e-10
                 O4- + O2b1 -> O2-  + O2  + O2                    : 1.0e-10
                 O-  + O2 + NEUTRAL -> O3- + NEUTRAL      : {1.1e-30 * (300./TeffN)}
                 O2- + O2 + NEUTRAL -> O4- + NEUTRAL      : {3.5e-31 * (300./TeffN2)}
                 ####
                 # ion-ion recombination (Kossyi1992)
                 ####
                 O- + N+ -> O + N                                 : {2e-7 * (300/TionN)^0.5}
                 O- + N2+ -> O + N2                               : {2e-7 * (300/TionN)^0.5}
                 O- + O+ -> O + O                                 : {2e-7 * (300/TionN)^0.5}
                 O- + O2+ -> O + O2                               : {2e-7 * (300/TionN)^0.5}
                 O- + NO+ -> O + NO                               : {2e-7 * (300/TionN)^0.5}
                 O2- + N+ -> O2 + N                               : {2e-7 * (300/TionN)^0.5}
                 O2- + N2+ -> O2 + N2                              : {2e-7 * (300/TionN)^0.5}
                 O2- + O+ -> O2 + O                                 : {2e-7 * (300/TionN)^0.5}
                 O2- + O2+ -> O2 + O2                               : {2e-7 * (300/TionN)^0.5}
                 O2- + NO+ -> O2 + NO                               : {2e-7 * (300/TionN)^0.5}
                 O3- + N+ -> O3 + N                                 : {2e-7 * (300/TionN)^0.5}
                 O3- + N2+ -> O3 + N2                               : {2e-7 * (300/TionN)^0.5}
                 O3- + O+ -> O3 + O                                 : {2e-7 * (300/TionN)^0.5}
                 O3- + O2+ -> O3 + O2                               : {2e-7 * (300/TionN)^0.5}
                 O3- + NO+ -> O3 + NO                               : {2e-7 * (300/TionN)^0.5}
                 NO- + N+ -> NO + N                                 : {2e-7 * (300/TionN)^0.5}
                 NO- + N2+ -> NO + N2                               : {2e-7 * (300/TionN)^0.5}
                 NO- + O+ -> NO + O                                 : {2e-7 * (300/TionN)^0.5}
                 NO- + O2+ -> NO + O2                               : {2e-7 * (300/TionN)^0.5}
                 NO- + NO+ -> NO + NO                               : {2e-7 * (300/TionN)^0.5}
                 O- + N2+ -> O + N + N                              : 1e-7
                 O- + N3+ -> O + N + N2                             : 1e-7
                 O- + N4+ -> O + N2 + N2                            : 1e-7
                 O- + O2+ -> O + O + O                              : 1e-7
                 O- + O4+ -> O + O2 + O2                            : 1e-7
                 O- + NO+ -> O + N + O                              : 1e-7
                 O- + O2pN2 -> O + O2 + N2                          : 1e-7
                 O2- + N2+ -> O2 + N + N                              : 1e-7
                 O2- + N3+ -> O2 + N + N2                             : 1e-7
                 O2- + N4+ -> O2 + N2 + N2                            : 1e-7
                 O2- + O2+ -> O2 + O + O                              : 1e-7
                 O2- + O4+ -> O2 + O2 + O2                            : 1e-7
                 O2- + NO+ -> O2 + N + O                              : 1e-7
                 O2- + O2pN2 -> O2 + O2 + N2                          : 1e-7
                 O3- + N2+ -> O3 + N + N                              : 1e-7
                 O3- + N3+ -> O3 + N + N2                             : 1e-7
                 O3- + N4+ -> O3 + N2 + N2                            : 1e-7
                 O3- + O2+ -> O3 + O + O                              : 1e-7
                 O3- + O4+ -> O3 + O2 + O2                            : 1e-7
                 O3- + NO+ -> O3 + N + O                              : 1e-7
                 O3- + O2pN2 -> O3 + O2 + N2                          : 1e-7
                 NO- + N2+ -> NO + N + N                              : 1e-7
                 NO- + N3+ -> NO + N + N2                             : 1e-7
                 NO- + N4+ -> NO + N2 + N2                            : 1e-7
                 NO- + O2+ -> NO + O + O                              : 1e-7
                 NO- + O4+ -> NO + O2 + O2                            : 1e-7
                 NO- + NO+ -> NO + N + O                              : 1e-7
                 NO- + O2pN2 -> NO + O2 + N2                          : 1e-7
                 O4- + N+ -> O2 + O2 + N                              : 1e-7 
                 O4- + N2+ -> O2 + O2 + N2                            : 1e-7 
                 O4- + O+ -> O2 + O2 + O                              : 1e-7 
                 O4- + O2+ -> O2 + O2 + O2                            : 1e-7 
                 O4- + NO+ -> O2 + O2 + NO                            : 1e-7 
                 O4- + N3+ -> O2 + O2 + N2 + N                        : 1e-7 
                 O4- + N4+ -> O2 + O2 + N2 + N2                       : 1e-7 
                 O4- + O4+ -> O2 + O2 + O2 + O2                       : 1e-7 
                 O4- + O2pN2 -> O2 + O2 + O2 + N2                     : 1e-7 
                 O- + N+ + NEUTRAL -> O + N + NEUTRAL                 : {2e-25 * (300/TionN)^2.5}
                 O- + N2+ + NEUTRAL -> O + N2 + NEUTRAL               : {2e-25 * (300/TionN)^2.5}
                 O- + O+ + NEUTRAL -> O + O + NEUTRAL                 : {2e-25 * (300/TionN)^2.5}
                 O- + O2+ + NEUTRAL -> O + O2 + NEUTRAL               : {2e-25 * (300/TionN)^2.5}
                 O- + NO+ + NEUTRAL -> O + NO + NEUTRAL               : {2e-25 * (300/TionN)^2.5}
                 O2- + N+ + NEUTRAL -> O2 + N + NEUTRAL               : {2e-25 * (300/TionN)^2.5}
                 O2- + N2+ + NEUTRAL -> O2 + N2 + NEUTRAL             : {2e-25 * (300/TionN)^2.5}
                 O2- + O+ + NEUTRAL -> O2 + O + NEUTRAL               : {2e-25 * (300/TionN)^2.5}
                 O2- + O2+ + NEUTRAL -> O2 + O2 + NEUTRAL             : {2e-25 * (300/TionN)^2.5}
                 O2- + NO+ + NEUTRAL -> O2 + NO + NEUTRAL             : {2e-25 * (300/TionN)^2.5}
                 O- + N+ + NEUTRAL -> NO + NEUTRAL                    : {2e-25 * (300/TionN)^2.5}
                 O- + O+ + NEUTRAL -> O2 + NEUTRAL                    : {2e-25 * (300/TionN)^2.5}
                 O- + O2+ + NEUTRAL -> O3 + NEUTRAL                   : {2e-25 * (300/TionN)^2.5}
                 O2- + O+ + NEUTRAL -> O3 + NEUTRAL                   : {2e-25 * (300/TionN)^2.5}
                 ####
                 # Three-body recombination of O3- NO- N2O- NO2- NO3- has to be verified
                 ####
                 O3- + N+ + NEUTRAL -> O3 + N + NEUTRAL               : {2e-25 * (300/TionN2)^2.5}
                 O3- + N2+ + NEUTRAL -> O3 + N2 + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 O3- + O+ + NEUTRAL -> O3 + O + NEUTRAL               : {2e-25 * (300/TionN2)^2.5}
                 O3- + O2+ + NEUTRAL -> O3 + O2 + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 O3- + NO+ + NEUTRAL -> O3 + NO + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO- + N+ + NEUTRAL -> NO + N + NEUTRAL               : {2e-25 * (300/TionN2)^2.5}
                 NO- + N2+ + NEUTRAL -> NO + N2 + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO- + O+ + NEUTRAL -> NO + O + NEUTRAL               : {2e-25 * (300/TionN2)^2.5}
                 NO- + O2+ + NEUTRAL -> NO + O2 + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO- + NO+ + NEUTRAL -> NO + NO + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 em + N2O -> em + em + N2O+            : EEDF (C81_N2O_Ionization_12.89_eV)
                 em + N2O+ -> N2 + O                 : {2.0e-7 * (300/Te)^0.5}
                 em + NO2+ -> NO + O                 : {2.0e-7 * (300/Te)^0.5}
                 em + N2O -> NO- + N                   : EEDF (C73_N2O_Attachment)
                 em + NO2 -> O- + NO                   : 1e-11
                 em + N2O + NEUTRAL -> N2O- + NEUTRAL  : 1e-31
                 O-  + NO -> NO2 + em                  : 2.6e-10
                 O-  + N2 -> N2O + em                  : 5.0e-13
                 NO- + N -> N2O + em                   : 5e-10 
                 N2O- + N -> NO + N2 + em              : 5e-10
                 NO2- + N -> NO + NO + em              : 5e-10
                 NO3- + N -> NO + NO2 + em             : 5e-10
                 NO- + O -> NO2 + em                   : 1.5e-10
                 N2O- + O -> NO + NO + em              : 1.5e-10
                 NO2- + O -> NO + O2 + em              : 1.5e-10
                 NO3- + O -> NO + O3 + em              : 1.5e-10
                 N2O- + N2A3 -> N2O + N2 + em          : 2.1e-9
                 NO2- + N2A3 -> NO2 + N2 + em          : 2.1e-9
                 NO3- + N2A3 -> NO3 + N2 + em          : 2.1e-9
                 N2O- + N2B3 -> N2O + N2 + em          : 2.5e-9
                 NO2- + N2B3 -> NO2 + N2 + em          : 2.5e-9
                 NO3- + N2B3 -> NO3 + N2 + em          : 2.5e-9
                 N2A3 + O2 -> N2O + O                 : {2.0e-14 * (Tgas/300)^0.55}
                 N2A3 + N2O -> N2 + N + NO            : 1.0e-11
                 N2A3 + NO2 -> N2 + O + NO            : 1.0e-12
                 N2D + N2O -> NO + N2                 : 3.5e-12
                 O1D + N2O -> NO + NO                 : 7.2e-11
                 O1D + N2O -> O2 + N2                 : 4.4e-11
                 O1S + N2O -> O + N2O                 : 6.3e-12
                 O1S + N2O -> O1D + N2O               : 3.1e-12
                 N + NO2 -> O + O + N2                : 9.1e-13
                 N + NO2 -> O + N2O                   : 3.0e-12
                 N + NO2 -> N2 + O2                   : 7.0e-13
                 N + NO2 -> NO + NO                   : 2.3e-12
                 O2- + N -> NO2 + em                   : 5.0e-10
                 O + NO -> NO2                        : 4.2e-18
                 O + N2O -> N2 + O2                   : {8.3e-12 * exp(-14000/Tgas)}
                 O + N2O -> NO + NO                   : {1.5e-10 * exp(-14090/Tgas)}
                 O + NO2 -> NO + O2                   : {9.1e-12 * (Tgas/300)^0.18}
                 O + NO3 -> O2 + NO2                  : 1.0e-11
                 N2 + O2 -> O + N2O                   : {2.5e-10 * exp(-50390/Tgas)}
                 NO + NO -> N + NO2                   : {3.3e-16 * (300/Tgas)^0.5 * exp(-39200/Tgas)}
                 NO + NO -> O + N2O                   : {2.2e-12 * exp(-32100/Tgas)}
                 NO + O2 -> O + NO2                   : {2.8e-12 * exp(-23400/Tgas)}
                 NO + O3 -> O2 + NO2                  : {2.5e-13 * exp(-765/Tgas)}
                 NO + N2O -> N2 + NO2                 : {4.6e-10 * exp(-25170/Tgas)}
                 NO + NO3 -> NO2 + NO2                : 1.7e-11
                 O2 + NO2 -> NO + O3                  : {2.8e-12 * exp(-25400/Tgas)}
                 NO2 + NO2 -> NO + NO + O2            : {3.3e-12 * exp(-13500/Tgas)}
                 NO2 + NO2 -> NO + NO3                : {4.5e-10 * exp(-18500/Tgas)}
                 NO2 + O3 -> O2 + NO3                 : {1.2e-13 * exp(-2450/Tgas)}
                 NO2 + NO3 -> NO + NO2 + O2           : {2.3e-13 * exp(-1600/Tgas)}
                 NO3 + O2 -> NO2 + O3                 : {1.5e-12 * exp(-15020/Tgas)}
                 NO3 + NO3 -> O2 + NO2 + NO2          : {4.3e-12 * exp(-3850/Tgas)}
                 N2O + N2 -> N2 + O + N2              : {1.2e-8 * (300/Tgas) * exp(-29000/Tgas)}
                 N2O + O2 -> N2 + O + O2              : {1.2e-8 * (300/Tgas) * exp(-29000/Tgas)}
                 N2O + NO -> N2 + O + NO              : {1.2e-8 * (300/Tgas) * exp(-29000/Tgas) * 2}
                 N2O + N2O -> N2 + O + N2O            : {1.2e-8 * (300/Tgas) * exp(-29000/Tgas) * 4}
                 NO2 + N2 -> N2 + O + N2              : {6.8e-6 * (300/Tgas)^2 * exp(-36180/Tgas)}
                 NO2 + O2 -> N2 + O + O2              : {6.8e-6 * (300/Tgas)^2 * exp(-36180/Tgas) * 0.78}
                 NO2 + NO -> N2 + O + NO              : {6.8e-6 * (300/Tgas)^2 * exp(-36180/Tgas) * 7.8}
                 NO2 + NO2 -> N2 + O + NO2            : {6.8e-6 * (300/Tgas)^2 * exp(-36180/Tgas) * 5.9}
                 NO3 + N2 -> NO2 + O + N2             : {3.1e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + O2 -> NO2 + O + O2             : {3.1e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + NO -> NO2 + O + NO             : {3.1e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + N -> NO2 + O + N               : {3.1e-5 * (300/Tgas)^2 * exp(-25000/Tgas) * 10}
                 NO3 + O -> NO2 + O + O               : {3.1e-5 * (300/Tgas)^2 * exp(-25000/Tgas) * 10}
                 NO3 + N2 -> NO + O2 + N2             : {6.2e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + O2 -> NO + O2 + O2             : {6.2e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + NO -> NO + O2 + NO             : {6.2e-5 * (300/Tgas)^2 * exp(-25000/Tgas)}
                 NO3 + N -> NO + O2 + N               : {6.2e-5 * (300/Tgas)^2 * exp(-25000/Tgas) * 12}
                 NO3 + O -> NO + O2 + O               : {6.2e-5 * (300/Tgas)^2 * exp(-25000/Tgas) * 12}
                 N2O5 + NEUTRAL -> NO2 + NO3 + NEUTRAL : {2.1e-11 * (300/Tgas)^4.4 * exp(-11080/Tgas)}
                 O + N2 + NEUTRAL -> N2O + NEUTRAL    : {3.9e-35 * exp(-10400/Tgas)}
                 O + NO + N2 -> NO2 + N2              : {1.2e-31 * (300/Tgas)^1.8}
                 O + NO + O2 -> NO2 + O2              : {1.2e-31 * (300/Tgas)^1.8 * 0.78}
                 O + NO + NO -> NO2 + NO              : {1.2e-31 * (300/Tgas)^1.8 * 0.78}
                 O + NO2 + N2 -> NO3 + N2             : {8.9e-32 * (300/Tgas)^2}
                 O + NO2 + O2 -> NO3 + O2             : {8.9e-32 * (300/Tgas)^2}
                 O + NO2 + N -> NO3 + N               : {8.9e-32 * (300/Tgas)^2 * 13}
                 O + NO2 + O -> NO3 + O               : {8.9e-32 * (300/Tgas)^2 * 13}
                 O + NO2 + NO -> NO3 + NO             : {8.9e-32 * (300/Tgas)^2 * 2.4}
                 NO2 + NO3 + NEUTRAL -> N2O5 + NEUTRAL : {3.7e-30 * (300/Tgas)^4.1}
                 N+ + N2O -> NO+ + N2                             : 5.5e-10
                 O+ + N2O -> NO+ + NO                             : 2.3e-10
                 O+ + N2O -> N2O+ + O                             : 2.2e-10
                 O+ + N2O -> O2+ + N2                             : 2.0e-11
                 O+ + NO2 -> NO2+ + O                             : 1.6e-9
                 N2+ + N2O -> N2O+ + N2                           : 5.0e-10
                 N2+ + N2O -> NO+ + N + N2                        : 4.0e-10
                 N3+ + NO -> N2O+ + N2                            : 7.0e-11
                 NO2+ + NO -> NO+ + NO2                           : 2.9e-10
                 N2O+ + NO -> NO+ + N2O                           : 2.9e-10
                 O-   + NO2 -> NO2- + O                           : 1.2e-9
                 O-   + N2O -> NO-  + NO                          : 2.0e-10
                 O-   + N2O -> N2O- + O                           : 2.0e-12
                 NO-  + NO2 -> NO2- + NO                          : 7.4e-10
                 NO-  + N2O -> NO2- + N2                          : 2.8e-14
                 NO2- + O3 -> NO3- + O2                           : 1.8e-11
                 NO2- + NO2 -> NO3- + NO                          : 4.0e-12
                 NO2- + NO3 -> NO3- + NO2                         : 5.0e-10
                 NO2- + N2O5 -> NO3- + NO2 + NO2                  : 7.0e-10
                 NO3- + NO -> NO2- + NO2                          : 3.0e-15
                 O- + N2O+ -> O + N2O                             : {2e-7 * (300/TionN)^0.5}
                 O- + NO2+ -> O + NO2                             : {2e-7 * (300/TionN)^0.5}
                 O2- + N2O+ -> O2 + N2O                             : {2e-7 * (300/TionN)^0.5}
                 O2- + NO2+ -> O2 + NO2                             : {2e-7 * (300/TionN)^0.5}
                 O3- + N2O+ -> O3 + N2O                             : {2e-7 * (300/TionN)^0.5}
                 O3- + NO2+ -> O3 + NO2                             : {2e-7 * (300/TionN)^0.5}
                 NO- + N2O+ -> NO + N2O                             : {2e-7 * (300/TionN)^0.5}
                 NO- + NO2+ -> NO + NO2                             : {2e-7 * (300/TionN)^0.5}
                 N2O- + N+ -> N2O + N                                 : {2e-7 * (300/TionN)^0.5}
                 N2O- + N2+ -> N2O + N2                               : {2e-7 * (300/TionN)^0.5}
                 N2O- + O+ -> N2O + O                                 : {2e-7 * (300/TionN)^0.5}
                 N2O- + O2+ -> N2O + O2                               : {2e-7 * (300/TionN)^0.5}
                 N2O- + NO+ -> N2O + NO                               : {2e-7 * (300/TionN)^0.5}
                 N2O- + N2O+ -> N2O + N2O                             : {2e-7 * (300/TionN)^0.5}
                 N2O- + NO2+ -> N2O + NO2                             : {2e-7 * (300/TionN)^0.5}
                 NO2- + N+ -> NO2 + N                                 : {2e-7 * (300/TionN)^0.5}
                 NO2- + N2+ -> NO2 + N2                               : {2e-7 * (300/TionN)^0.5}
                 NO2- + O+ -> NO2 + O                                 : {2e-7 * (300/TionN)^0.5}
                 NO2- + O2+ -> NO2 + O2                               : {2e-7 * (300/TionN)^0.5}
                 NO2- + NO+ -> NO2 + NO                               : {2e-7 * (300/TionN)^0.5}
                 NO2- + N2O+ -> NO2 + N2O                             : {2e-7 * (300/TionN)^0.5}
                 NO2- + NO2+ -> NO2 + NO2                             : {2e-7 * (300/TionN)^0.5}
                 NO3- + N+ -> NO3 + N                                 : {2e-7 * (300/TionN)^0.5}
                 NO3- + N2+ -> NO3 + N2                               : {2e-7 * (300/TionN)^0.5}
                 NO3- + O+ -> NO3 + O                                 : {2e-7 * (300/TionN)^0.5}
                 NO3- + O2+ -> NO3 + O2                               : {2e-7 * (300/TionN)^0.5}
                 NO3- + NO+ -> NO3 + NO                               : {2e-7 * (300/TionN)^0.5}
                 NO3- + N2O+ -> NO3 + N2O                             : {2e-7 * (300/TionN)^0.5}
                 NO3- + NO2+ -> NO3 + NO2                             : {2e-7 * (300/TionN)^0.5}
                 O- + N2O+ -> O + N2 + O                            : 1e-7
                 O- + NO2+ -> O + N + O2                            : 1e-7
                 O2- + N2O+ -> O2 + N2 + O                            : 1e-7
                 O2- + NO2+ -> O2 + N + O2                            : 1e-7
                 O3- + N2O+ -> O3 + N2 + O                            : 1e-7
                 O3- + NO2+ -> O3 + N + O2                            : 1e-7
                 NO- + N2O+ -> NO + N2 + O                            : 1e-7
                 NO- + NO2+ -> NO + N + O2                            : 1e-7
                 N2O- + N2+ -> N2O + N + N                            : 1e-7
                 N2O- + N3+ -> N2O + N + N2                           : 1e-7
                 N2O- + N4+ -> N2O + N2 + N2                          : 1e-7
                 N2O- + O2+ -> N2O + O + O                            : 1e-7
                 N2O- + O4+ -> N2O + O2 + O2                          : 1e-7
                 N2O- + NO+ -> N2O + N + O                            : 1e-7
                 N2O- + N2O+ -> N2O + N2 + O                          : 1e-7
                 N2O- + NO2+ -> N2O + N + O2                          : 1e-7
                 N2O- + O2pN2 -> N2O + O2 + N2                        : 1e-7
                 NO2- + N2+ -> NO2 + N + N                            : 1e-7
                 NO2- + N3+ -> NO2 + N + N2                           : 1e-7
                 NO2- + N4+ -> NO2 + N2 + N2                          : 1e-7
                 NO2- + O2+ -> NO2 + O + O                            : 1e-7
                 NO2- + O4+ -> NO2 + O2 + O2                          : 1e-7
                 NO2- + NO+ -> NO2 + N + O                            : 1e-7
                 NO2- + N2O+ -> NO2 + N2 + O                          : 1e-7
                 NO2- + NO2+ -> NO2 + N + O2                          : 1e-7
                 NO2- + O2pN2 -> NO2 + O2 + N2                        : 1e-7
                 NO3- + N2+ -> NO3 + N + N                            : 1e-7
                 NO3- + N3+ -> NO3 + N + N2                           : 1e-7
                 NO3- + N4+ -> NO3 + N2 + N2                          : 1e-7
                 NO3- + O2+ -> NO3 + O + O                            : 1e-7
                 NO3- + O4+ -> NO3 + O2 + O2                          : 1e-7
                 NO3- + NO+ -> NO3 + N + O                            : 1e-7
                 NO3- + N2O+ -> NO3 + N2 + O                          : 1e-7
                 NO3- + NO2+ -> NO3 + N + O2                          : 1e-7
                 NO3- + O2pN2 -> NO3 + O2 + N2                        : 1e-7
                 O4- + N2O+ -> O2 + O2 + N2O                          : 1e-7 
                 O4- + NO2+ -> O2 + O2 + NO2                          : 1e-7 
                 O- + N2+ + NEUTRAL -> N2O + NEUTRAL                  : {2e-25 * (300/TionN)^2.5}
                 O3- + N2O+ + NEUTRAL -> O3 + N2O + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 O3- + NO2+ + NEUTRAL -> O3 + NO2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO- + N2O+ + NEUTRAL -> NO + N2O + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO- + NO2+ + NEUTRAL -> NO + NO2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 N2O- + N+ + NEUTRAL -> N2O + N + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 N2O- + N2+ + NEUTRAL -> N2O + N2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 N2O- + O+ + NEUTRAL -> N2O + O + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 N2O- + O2+ + NEUTRAL -> N2O + O2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 N2O- + NO+ + NEUTRAL -> N2O + NO + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 N2O- + N2O+ + NEUTRAL -> N2O + N2O + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 N2O- + NO2+ + NEUTRAL -> N2O + NO2 + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 NO2- + N+ + NEUTRAL -> NO2 + N + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO2- + N2+ + NEUTRAL -> NO2 + N2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO2- + O+ + NEUTRAL -> NO2 + O + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO2- + O2+ + NEUTRAL -> NO2 + O2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO2- + NO+ + NEUTRAL -> NO2 + NO + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO2- + N2O+ + NEUTRAL -> NO2 + N2O + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 NO2- + NO2+ + NEUTRAL -> NO2 + NO2 + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 NO3- + N+ + NEUTRAL -> NO3 + N + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO3- + N2+ + NEUTRAL -> NO3 + N2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO3- + O+ + NEUTRAL -> NO3 + O + NEUTRAL             : {2e-25 * (300/TionN2)^2.5}
                 NO3- + O2+ + NEUTRAL -> NO3 + O2 + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO3- + NO+ + NEUTRAL -> NO3 + NO + NEUTRAL           : {2e-25 * (300/TionN2)^2.5}
                 NO3- + N2O+ + NEUTRAL -> NO3 + N2O + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 NO3- + NO2+ + NEUTRAL -> NO3 + NO2 + NEUTRAL         : {2e-25 * (300/TionN2)^2.5}
                 O2+ + NO2 -> NO+ + O3                            : 1.0e-11
                 O2+ + NO2 -> NO2+ + O2                           : 6.6e-10
                 N3+ + O2 -> NO2+ + N2                            : 4.4e-11
                 O2-  + NO2 -> NO2- + O2                          : 7.0e-10
                 O2-  + NO3 -> NO3- + O2                          : 5.0e-10
                 O3-  + NO -> NO3- + O                            : 1.0e-11
                 O3-  + NO -> NO2- + O2                           : 2.6e-12
                 O3-  + NO2 -> NO2- + O3                          : 7.0e-11
                 O3-  + NO2 -> NO3- + O2                          : 2.0e-11
                 O3-  + NO3 -> NO3- + O3                          : 5.0e-10
                 O-  + NO + NEUTRAL -> NO2- + NEUTRAL     : 1.0e-29
                 O- + NO+ + NEUTRAL -> NO2 + NEUTRAL                  : {2e-25 * (300/TionN)^2.5}
                 O2- + N+ + NEUTRAL -> NO2 + NEUTRAL                  : {2e-25 * (300/TionN)^2.5}
                 O4- + NO -> NO3- + O2                            : 2.5e-10
                 O2- + NO+ + NEUTRAL -> NO3 + NEUTRAL                 : {2e-25 * (300/TionN)^2.5}'
  [../]
[]
