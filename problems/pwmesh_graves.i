#dom0Scale=1.0e-3
dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  offset = 40
  potential_units = kV
  use_moles = true
[]

[Mesh]
  type = FileMesh
  file = 'pwmesh_04.msh'
  #file = 'pwmesh_exodus_mesh.e'
  #parallel_type = DISTRIBUTED
[]

[MeshModifiers]
  [./interface]
    type = SideSetsBetweenSubdomains
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    #boundary = 'interface'
    # depends_on = 'box'
  [../]
  [./interface_again]
    type = SideSetsBetweenSubdomains
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
    #boundary = 'interface'
    # depends_on = 'box'
  [../]
  # [./box]
  #   type = SubdomainBoundingBox
  #   bottom_left = '0.55 0 0'
  #   top_right = '1.1 1. 0'
  #   block_id = 1
  # [../]
[]

[Problem]
  coord_type = RZ
  type = FEProblem
[]

[Preconditioning]
  #active = fsp_schur
  active = smp

  [./smp]
    type = SMP
    #full = true
    off_diag_row = '   em  Arp  em  Ar* potential potential potential  OHm   emliq H2    emliq emliq mean_en mean_en   em        Arp'
    off_diag_column = 'Arp em   Ar* em  em        Arp       emliq      emliq O2m   emliq H2O2  HO2m  em      potential potential potential'
    petsc_options_iname = '-snes_linesearch_type -pc_type -pc_factor_shift_type'
    petsc_options_value = 'basic lu NONZERO'
  [../]

  [./fsp_schur]
    type = FSP
    topsplit = 'pc'
    #full = 'true'
    off_diag_row = '   em  em  Ar* potential potential potential potential potential OHm   emliq H2    emliq emliq emliq emliq'
    off_diag_column = 'Arp Ar* em  em        Arp       emliq     Hp        OHm       emliq O2m   emliq H2O2  OH    HO2   HO2m'
    [./pc]
      splitting = 'p c'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
      petsc_options_value = 'hypre boomeramg 0.4'
    [../]
    [./c]
      vars = 'em Arp Ar* mean_en emliq OHm Om O2m O3m HO2m Hp O O2_1 O3 H H2 HO2 HO3 OH H2O2'
      petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      petsc_options_value = 'asm 1 ilu NONZERO 1e-10'
    [../]
  [../]

  [./fsp]
    type = FSP
    topsplit = 'pc'
    #full = 'true'
    [./pc]
      splitting = 'p c'
      splitting_type = additive
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type -pc_hypre_type'
      petsc_options_value = 'hypre boomeramg'
    [../]
    [./c]
      vars = 'em Arp Ar* mean_en emliq OHm Om O2m O3m HO2m Hp O O2_1 O3 H H2 HO2 HO3 OH H2O2'
      petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      petsc_options_value = 'asm 2 ilu NONZERO 1e-10'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  end_time = 1e-1
  automatic_scaling = true
  compute_scaling_once = false
  verbose = true

  #line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  #solve_type = newton
  solve_type = pjfnk
  #petsc_options_iname = '-snes_linesearch_type -pc_type -pc_factor_shift_type'
  #petsc_options_value = 'basic lu NONZERO'
  #petsc_options_iname = '-snes_linesearch_type -pc_type -pc_gamg_sym_graph -pc_factor_shift_type'
  #petsc_options_value = 'basic gamg true NONZERO'
  #petsc_options_iname = '-snes_linesearch_type -pc_type -sub_pc_type -pc_asm_overlap -pc_factor_shift_type'
  #petsc_options_value = 'basic asm ilu 4 NONZERO'
  #petsc_options_iname = '-snes_linesearch_type -pc_type'
  #petsc_options_value = 'basic asm'
  nl_rel_tol = 1e-5
  l_tol = 1e-3
  #nl_abs_tol = 7.6e-5
  dtmin = 1e-16
  l_max_its = 20
  #num_steps = 1
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-12
    growth_factor = 1.2
    optimal_iterations = 15
  [../]
[]

[Outputs]
  #checkpoint = true
  perf_graph = true
  exodus = true
  #nemesis = true
[]

[Debug]
  show_var_residual_norms = true
[]

[UserObjects]
  [./data_provider]
    type = ProvideMobility
    #electrode_area = 5.02e-7 # Formerly 3.14e-6
    #ballast_resist = 1e6
    e = 1.6e-19
    #electrode_area = 5.654867e-4
    #electrode_area = 5.654867e-3
    electrode_area = 1.0
    ballast_resist = 2.2e5
  [../]
[]

[Kernels]
  [./em_time_deriv]
    type = ElectronTimeDerivative
    variable = em
    block = 0
  [../]
  [./em_advection]
    type = EFieldAdvectionElectrons
    variable = em
    potential = potential
    mean_en = mean_en
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_diffusion]
    type = CoeffDiffusionElectrons
    variable = em
    mean_en = mean_en
    initial_condition = -20
    block = 0
    position_units = ${dom0Scale}
  [../]
  #[./em_dd]
  #  type = DriftDiffusionElectrons
  #  variable = em
  #  potential = potential
  #  mean_en = mean_en
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    block = 0
    #offset = 25
  [../]

  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential
    block = 1
    position_units = ${dom1Scale}
  [../]

  [./Arp_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = Arp
    block = 0
  [../]
  [./em_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = em
    block = 0
  [../]

  [./Arp_time_deriv]
    type = ElectronTimeDerivative
    variable = Arp
    block = 0
  [../]
  [./Arp_advection]
    type = EFieldAdvection
    variable = Arp
    potential = potential
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Arp_diffusion]
    type = CoeffDiffusion
    variable = Arp
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_log_stabilization]
    type = LogStabilizationMoles
    variable = Arp
    block = 0
    #offset = 31
  [../]
  #[./Arp_advection_stabilization]
  #  type = EFieldArtDiff
  #  variable = Arp
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]

  [./Arex_time_deriv]
    type = ElectronTimeDerivative
    variable = Ar*
    block = 0
  [../]
  [./Arex_diffusion]
    type = CoeffDiffusion
    variable = Ar*
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex_log_stabilization]
    type = LogStabilizationMoles
    variable = Ar*
    block = 0
    #offset = 25
  [../]

  [./mean_en_time_deriv]
    type = ElectronTimeDerivative
    variable = mean_en
    block = 0
  [../]
  [./mean_en_advection]
    type = EFieldAdvectionEnergy
    variable = mean_en
    potential = potential
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_diffusion]
    type = CoeffDiffusionEnergy
    variable = mean_en
    em = em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_joule_heating]
    type = JouleHeating
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
    #offset = 21
    #offset = 15
  [../]


  [./emliq_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = emliq
    block = 1
  [../]
  [./OHm_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = OHm
    block = 1
  [../]
  [./Hp_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = Hp
    block = 1
  [../]
  [./H2Op_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = H2Op
    block = 1
  [../]
  [./H3Op_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = H3Op
    block = 1
  [../]
  [./HO2m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = HO2m
    block = 1
  [../]
  [./O2m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = O2m
    block = 1
  [../]
  [./O3m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential
    charged = O3m
    block = 1
  [../]


  [./emliq_time_deriv]
    type = ElectronTimeDerivative
    variable = emliq
    block = 1
  [../]
  [./emliq_advection]
    type = EFieldAdvection
    variable = emliq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./emliq_diffusion]
    type = CoeffDiffusion
    variable = emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./emliq_log_stabilization]
    type = LogStabilizationMoles
    variable = emliq
    block = 1
    offset = 28
    #offset = 21
    #offset = 15
  [../]

  [./OHm_time_deriv]
    type = ElectronTimeDerivative
    variable = OHm
    block = 1
  [../]
  [./OHm_advection]
    type = EFieldAdvection
    variable = OHm
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./OHm_diffusion]
    type = CoeffDiffusion
    variable = OHm
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./OHm_log_stabilization]
    type = LogStabilizationMoles
    variable = OHm
    #offset = 28
    #offset = 31
    offset = 11
    block = 1
  [../]

  [./Om_time_deriv]
    type = ElectronTimeDerivative
    variable = Om
    block = 1
  [../]
  [./Om_advection]
    type = EFieldAdvection
    variable = Om
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Om_diffusion]
    type = CoeffDiffusion
    variable = Om
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Om_log_stabilization]
    type = LogStabilizationMoles
    variable = Om
    block = 1
    offset = 31
  [../]


  [./O2m_time_deriv]
    type = ElectronTimeDerivative
    variable = O2m
    block = 1
  [../]
  [./O2m_advection]
    type = EFieldAdvection
    variable = O2m
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O2m_diffusion]
    type = CoeffDiffusion
    variable = O2m
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O2m_log_stabilization]
    type = LogStabilizationMoles
    variable = O2m
    block = 1
    offset = 31
  [../]


  [./O3m_time_deriv]
    type = ElectronTimeDerivative
    variable = O3m
    block = 1
  [../]
  [./O3m_advection]
    type = EFieldAdvection
    variable = O3m
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O3m_diffusion]
    type = CoeffDiffusion
    variable = O3m
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O3m_log_stabilization]
    type = LogStabilizationMoles
    variable = O3m
    block = 1
    offset = 31
  [../]


  [./HO2m_time_deriv]
    type = ElectronTimeDerivative
    variable = HO2m
    block = 1
  [../]
  [./HO2m_advection]
    type = EFieldAdvection
    variable = HO2m
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO2m_diffusion]
    type = CoeffDiffusion
    variable = HO2m
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO2m_log_stabilization]
    type = LogStabilizationMoles
    variable = HO2m
    block = 1
    offset = 31
  [../]


  [./Hp_time_deriv]
    type = ElectronTimeDerivative
    variable = Hp
    block = 1
  [../]
  [./Hp_advection]
    type = EFieldAdvection
    variable = Hp
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Hp_diffusion]
    type = CoeffDiffusion
    variable = Hp
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Hp_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = Hp
    block = 1
  [../]

  [./O_time_deriv]
    type = ElectronTimeDerivative
    variable = O
    block = 1
  [../]
  [./O_advection]
    type = EFieldAdvection
    variable = O
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O_diffusion]
    type = CoeffDiffusion
    variable = O
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = O
    block = 1
  [../]

  [./O3_time_deriv]
    type = ElectronTimeDerivative
    variable = O3
    block = 1
  [../]
  [./O3_advection]
    type = EFieldAdvection
    variable = O3
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O3_diffusion]
    type = CoeffDiffusion
    variable = O3
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O3_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = O3
    block = 1
  [../]


  [./H_time_deriv]
    type = ElectronTimeDerivative
    variable = H
    block = 1
  [../]
  [./H_advection]
    type = EFieldAdvection
    variable = H
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H_diffusion]
    type = CoeffDiffusion
    variable = H
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = H
    block = 1
  [../]


  [./H2_time_deriv]
    type = ElectronTimeDerivative
    variable = H2
    block = 1
  [../]
  [./H2_advection]
    type = EFieldAdvection
    variable = H2
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2_diffusion]
    type = CoeffDiffusion
    variable = H2
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = H2
    block = 1
  [../]


  [./HO2_time_deriv]
    type = ElectronTimeDerivative
    variable = HO2
    block = 1
  [../]
  [./HO2_advection]
    type = EFieldAdvection
    variable = HO2
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO2_diffusion]
    type = CoeffDiffusion
    variable = HO2
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO2_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = HO2
    block = 1
  [../]


  [./HO3_time_deriv]
    type = ElectronTimeDerivative
    variable = HO3
    block = 1
  [../]
  [./HO3_advection]
    type = EFieldAdvection
    variable = HO3
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO3_diffusion]
    type = CoeffDiffusion
    variable = HO3
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./HO3_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = HO3
    block = 1
  [../]


  [./OH_time_deriv]
    type = ElectronTimeDerivative
    variable = OH
    block = 1
  [../]
  [./OH_advection]
    type = EFieldAdvection
    variable = OH
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./OH_diffusion]
    type = CoeffDiffusion
    variable = OH
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./OH_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = OH
    block = 1
  [../]


  [./H2O2_time_deriv]
    type = ElectronTimeDerivative
    variable = H2O2
    block = 1
  [../]
  [./H2O2_advection]
    type = EFieldAdvection
    variable = H2O2
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2O2_diffusion]
    type = CoeffDiffusion
    variable = H2O2
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2O2_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = H2O2
    block = 1
  [../]

  [./H2Op_time_deriv]
    type = ElectronTimeDerivative
    variable = H2Op
    block = 1
  [../]
  [./H2Op_advection]
    type = EFieldAdvection
    variable = H2Op
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2Op_diffusion]
    type = CoeffDiffusion
    variable = H2Op
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2Op_log_stabilization]
    offset = 31
    type = LogStabilizationMoles
    variable = H2Op
    block = 1
  [../]

  [./H3Op_time_deriv]
    type = ElectronTimeDerivative
    variable = H3Op
    block = 1
  [../]
  [./H3Op_advection]
    type = EFieldAdvection
    variable = H3Op
    #potential = potential_liq
    potential = potential
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H3Op_diffusion]
    type = CoeffDiffusion
    variable = H3Op
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H3Op_log_stabilization]
    offset = 11
    type = LogStabilizationMoles
    variable = H3Op
    block = 1
  [../]

[]

[Variables]
  [./potential]
  [../]

  [./em]
    initial_condition = -40
    block = 0
  [../]

  [./Arp]
    initial_condition = -40
    block = 0
  [../]

  [./Ar*]
    initial_condition = -40
    block = 0
  [../]

  [./mean_en]
    initial_condition = -40
    block = 0
  [../]

  [./emliq]
    block = 1
    initial_condition = -28
  [../]

  [./OHm]
    block = 1
    initial_condition = -9.2103
  [../]
  [./Om]
    block = 1
    initial_condition = -31
  [../]

  [./O2m]
    block = 1
    initial_condition = -31
  [../]

  [./O3m]
    block = 1
    initial_condition = -31
  [../]

  [./HO2m]
    block = 1
    initial_condition = -31
  [../]

  [./Hp]
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
    initial_condition = -31
  [../]

  [./H2]
    block = 1
    initial_condition = -31
  [../]

  [./HO2]
    block = 1
    initial_condition = -31
  [../]

  [./HO3]
    block = 1
    initial_condition = -31
  [../]

  [./OH]
    block = 1
    initial_condition = -31
  [../]

  [./H2O2]
    block = 1
    initial_condition = -31
  [../]

  [./H3Op]
    block = 1
    initial_condition = -9.2103
  [../]

  [./H2Op]
    block = 1
    initial_condition = -31
  [../]

  #[./Nap]
  #  block = 1
  #  initial_condition = -9.2103
  #[../]
  #
  #[./Clm]
  #  block = 1
  #  initial_condition = -9.2103
  #[../]



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
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./em_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Arp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  #[./Ar2p_lin]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[../]
  [./Arex_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./emliq_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./OHm_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]

  [./H2O]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 10.92252
    block = 1
  [../]
  [./O2]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = -0.609203
    block = 1
  [../]
  
  [./Hp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./H2Op_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./H3Op_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
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
  [./em_lin]
    type = DensityMoles
    convert_moles = true
    variable = em_lin
    density_log = em
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./Arp_lin]
    type = DensityMoles
    convert_moles = true
    variable = Arp_lin
    density_log = Arp
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  #[./Ar2p_lin]
  #  type = DensityMoles
  #  convert_moles = true
  #  variable = Ar2p_lin
  #  density_log = Ar2p
  #  execute_on = 'initial timestep_end'
  #  block = 0
  #[../]
  [./Arex_lin]
    type = DensityMoles
    convert_moles = true
    variable = Arex_lin
    density_log = Ar*
    execute_on = 'initial timestep_end'
    block = 0
  [../]
  [./emliq_lin]
    type = DensityMoles
    convert_moles = true
    variable = emliq_lin
    density_log = emliq
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./OHm_lin]
    type = DensityMoles
    convert_moles = true
    variable = OHm_lin
    density_log = OHm
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./H2Op_lin]
    type = DensityMoles
    convert_moles = true
    variable = H2Op_lin
    density_log = H2Op
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./H3Op_lin]
    type = DensityMoles
    convert_moles = true
    variable = H3Op_lin
    density_log = H3Op
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./Hp_lin]
    type = DensityMoles
    convert_moles = true
    variable = Hp_lin
    density_log = Hp
    execute_on = 'initial timestep_end'
    block = 1
  [../]
[]

[Postprocessors]
  [./electrode_flux]
    type = SideCurrent
    mobility = 'muem'
    variable = em
    potential = potential
    mean_en = mean_en
    Arp = Arp
    #Ar2p = Ar2p
    r = 0
    position_units = ${dom0Scale}
    boundary = 'electrode'
    execute_on = 'nonlinear linear'
  [../]
[]

[BCs]
  [./potential_bottom]
    type = DirichletBC
    variable = potential
    boundary = bottom
    value = 0
  [../]
  #[./potential_electrode]
  #  type = NeumannCircuitVoltageMoles_KV
  #  variable = potential
  #  boundary = electrode
  #  function = potential_bc_func
  #  ip = Arp
  #  data_provider = data_provider
  #  em = em
  #  mean_en = mean_en
  #  r = 0.0
  #  position_units = ${dom0Scale}
  #[../]
  #[./potential_electrode]
  #  type = DirichletBC
  #  variable = potential
  #  boundary = electrode
  #  value = -1.0
  #[../]
  [./potential_electrode]
    type = CircuitDirichletPotential
    variable = potential
    boundary = electrode
    surface_potential = potential_bc_func
    surface = 'cathode'
    resist = 2.5e5
    position_units = ${dom0Scale}
    #A = 1e-3
    #A = 1.519e-3
    current = electrode_flux
  [../]
  [./em_physical_bottom]
    type = HagelaarElectronBC
    variable = em
    #boundary = bottom
    boundary = 'master0_interface'
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_bottom_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    #boundary = bottom
    boundary = 'master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_bottom_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    #boundary = bottom
    boundary = 'master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_bottom]
    type = HagelaarEnergyBC
    variable = mean_en
    #boundary = bottom
    boundary = 'master0_interface'
    potential = potential
    em = em
    ip = Arp
    #ip2 = Ar2p
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./em_physical_electrode]
    type = HagelaarElectronBC
    variable = em
    boundary = electrode
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = SecondaryElectronBC
    variable = em
    boundary = electrode
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./sec_electrons_left_02]
  #  type = SecondaryElectronBC
  #  variable = em
  #  boundary = electrode
  #  potential = potential
  #  ip = Ar2p
  #  mean_en = mean_en
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./Arp_physical_inlet_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = electrode
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_inlet_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = electrode
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./Ar2p_physical_inlet_diffusion]
  #  type = HagelaarIonDiffusionBC
  #  variable = Ar2p
  #  boundary = electrode
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_physical_inlet_advection]
  #  type = HagelaarIonAdvectionBC
  #  variable = Ar2p
  #  boundary = electrode
  #  potential = potential
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./mean_en_physical_inlet]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = electrode
    potential = potential
    em = em
    ip = Arp
    #ip2 = Ar2p
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./Arp_right]
    type = DriftDiffusionDoNothingBC
    variable = Arp
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'right'
  [../]
  [./Arp_top]
    type = DriftDiffusionDoNothingBC
    variable = Arp
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'top'
  [../]
  #[./Ar2p_right]
  #  type = DriftDiffusionDoNothingBC
  #  variable = Ar2p
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  boundary = 'right'
  #[../]
  #[./Ar2p_top]
  #  type = DriftDiffusionDoNothingBC
  #  variable = Ar2p
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  boundary = 'top'
  #[../]

  [./Arex_right]
    type = DriftDiffusionDoNothingBC
    variable = Ar*
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'right'
  [../]
  [./Arex_top]
    type = DriftDiffusionDoNothingBC
    variable = Ar*
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'top'
  [../]
  [./Arex_electrode]
    type = DriftDiffusionDoNothingBC
    variable = Ar*
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'electrode'
  [../]
  [./Arex_bottom]
    type = DriftDiffusionDoNothingBC
    variable = Ar*
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'master0_interface'
  [../]

  [./em_right]
    type = DriftDiffusionDoNothingElectronBC
    variable = em
    potential = potential
    mean_en = mean_en
    position_units = ${dom0Scale}
    boundary = 'right'
  [../]
  [./em_top]
    type = DriftDiffusionDoNothingElectronBC
    variable = em
    potential = potential
    mean_en = mean_en
    position_units = ${dom0Scale}
    boundary = 'top'
  [../]

  [./mean_en_right]
    type = DoNothingEnergyBC
    variable = mean_en
    em = em
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'right'
  [../]
  [./mean_en_top]
    type = DoNothingEnergyBC
    variable = mean_en
    em = em
    potential = potential
    position_units = ${dom0Scale}
    boundary = 'top'
  [../]

  #[./potential_right]
  #  type = DoNothingPotentialBC
  #  variable = potential
  #  position_units = ${dom0Scale}
  #  boundary = 'right'
  #[../]
  #[./potential_right]
  #  type = DirichletBC
  #  variable = potential
  #  value = 0.0
  #  boundary = 'right'
  #[../]
  #[./potential_top]
  #  type = DoNothingPotentialBC
  #  variable = potential
  #  position_units = ${dom0Scale}
  #  boundary = 'top'
  #[../]


  [./emliq_right]
    type = DCIonBC
    variable = emliq
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]

  [./OHm_physical]
    type = DCIonBC
    variable = OHm
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2Op_physical]
    type = DCIonBC
    variable = H2Op
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H3Op_physical]
    type = DCIonBC
    variable = H3Op
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Hp_interface]
    type = DriftDiffusionDoNothingBC
    variable = Hp
    boundary = 'master1_interface'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2Op_interface]
    type = DriftDiffusionDoNothingBC
    variable = H2Op
    boundary = 'master1_interface'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H3Op_interface]
    type = DriftDiffusionDoNothingBC
    variable = H3Op
    boundary = 'master1_interface'
    potential = potential
    position_units = ${dom1Scale}
  [../]

  [./Om_physical]
    type = DCIonBC
    variable = Om
    #boundary = 'right'
    #potential = potential_liq
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O2m_physical]
    type = DCIonBC
    variable = O2m
    #boundary = 'right'
    #potential = potential_liq
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O3m_physical]
    type = DCIonBC
    variable = O3m
    #boundary = 'right'
    #potential = potential_liq
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./HO2m_physical]
    type = DCIonBC
    variable = HO2m
    #boundary = 'right'
    #potential = potential_liq
    boundary = 'bottom'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Hp_physical]
    type = DCIonBC
    variable = Hp
    #boundary = 'right'
    boundary = 'bottom'
    #potential = potential_liq
    potential = potential
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
[]

[ICs]
  #[./em_ic]
  #  type = FunctionIC
  #  variable = em
  #  function = em_ic_func
  #  block = 0
  #[../]
  #[./Arp_ic]
  #  type = FunctionIC
  #  variable = Arp
  #  function = charged_ic_func
  #  block = 0
  #[../]
  #[./Ar2p_ic]
  #  type = FunctionIC
  #  variable = Ar2p
  #  function = charged_ic_func
  #  block = 0
  #[../]
  #[./mean_en_ic]
  #  type = ConstantIC
  #  variable = mean_en
  #  #value = -23
  #  value = -24
  #  block = 0
  #[../]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
    #block = 0
  [../]
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    #value = '-2.0'
    #value = '-1.9*tanh(1e9*t) - 0.1'
    value = '-0.7*tanh(1e9*t) - 0.1'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-1.25 * (y - 1.000e-3)'
    # Below is a test of a gaussian initial condition
    vars = 'sigmax sigmay A x0 y0'
    vals = '4e-4 0.5e-3 -0.1 0 2.5e-3'
    value = 'A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2))'
    #vars = 'sigmax sigmay A x0 y0'
    #vals = '4e-1 0.9e0 -1.0 0 2.5e0'
    #value = 'A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2)) - 0.3 * (y - 1.000e-3)'
  [../]
  [./em_ic_func]
    type = ParsedFunction
    value = '-21'
  [../]
  [./charged_ic_func]
    type = ParsedFunction
    # Another gaussian initial condition, positioned just below the electrode tip
    #vars = 'sigmax sigmay A x0 y0'
    #vals = '2e-5 5e-5 4.5e14 0 8.75e-4'
    #value = 'if(y<=0.2e-4&z<=1.001e-3&z>=9.45e-4,log(A*exp(-((y-y0)^2 / (2*sigma^2) + (z-z0)^2/(2*sigma)^2))/6.022e23)+10,-30)'
    #value = 'if(log(A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2))/6.022e23)>-24,log(A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2))/6.022e23),-24)'
    #value = '-24'
    value = '-21.6931471805599452'
  [../]
[]

[Materials]
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
   property_tables_file = dc_2d_test/electron_moments.txt
   #property_tables_file = electron_moments_test.txt
   # property_tables_file = dc_stankov/electron_moments.txt
   position_units = ${dom0Scale}
   block = 0
 [../]

  [./gas_species_0]
    type = HeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
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
  [./gas_species_2]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0
    block = 0
  [../]
  [./Ar_species]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  [../]

 #[./water_block]
 #  type = Water
 #  block = 1
 #  potential = potential
 #[../]
  [./H2Op_mat]
    type = HeavySpeciesMod
    heavy_species_name = H2Op
    heavy_species_mass = 2.992e-26
    heavy_species_charge = 1
    diffusivity = 5e-9
    block = 1
  [../]

  [./H3Op_mat]
    type = HeavySpeciesMod
    heavy_species_name = H3Op
    heavy_species_mass = 2.992e-26
    heavy_species_charge = 1
    diffusivity = 9.3e-9
    block = 1
  [../]

  [./OHm_mat]
    type = HeavySpeciesMod
    heavy_species_name = OHm
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.27e-9
    block = 1
  [../]

  [./O_mat]
    type = HeavySpeciesMod
    heavy_species_name = O
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 5e-9
    block = 1
  [../]

  #[./O2_1_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = O2_1
  #  heavy_species_mass = 5.31365e-26
  #  heavy_species_charge = 0
  #  diffusivity = 1.97e-9
  #  block = 1
  #[../]

  [./O3_mat]
    type = HeavySpeciesMod
    heavy_species_name = O3
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = 0
    diffusivity = 5.0e-9
    block = 1
  [../]

  [./OH_mat]
    type = HeavySpeciesMod
    heavy_species_name = OH
    heavy_species_mass = 2.82431e-26
    heavy_species_charge = 0
    diffusivity = 4.5e-9
    block = 1
  [../]

  [./HO2_mat]
    type = HeavySpeciesMod
    heavy_species_name = HO2
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./HO3_mat]
    type = HeavySpeciesMod
    heavy_species_name = HO3
    heavy_species_mass = 8.13785e-26
    heavy_species_charge = 0
    diffusivity = 1.00e-9
    block = 1
  [../]

  [./H2O2_mat]
    type = HeavySpeciesMod
    heavy_species_name = H2O2
    heavy_species_mass = 5.64840e-26
    heavy_species_charge = 0
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./H2_mat]
    type = HeavySpeciesMod
    heavy_species_name = H2
    heavy_species_mass = 3.34752e-26
    heavy_species_charge = 0
    diffusivity = 4.50e-9
    block = 1
  [../]

  [./H_mat]
    type = HeavySpeciesMod
    heavy_species_name = H
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 0
    diffusivity = 4.50e-9
    block = 1
  [../]

  [./Hp_mat]
    type = HeavySpeciesMod
    heavy_species_name = Hp
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 1
    #diffusivity = 9.31e-9
    diffusivity = 4.50e-9
    block = 1
  [../]

  [./HO2m_mat]
    type = HeavySpeciesMod
    heavy_species_name = HO2m
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = -1
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./Om_mat]
    type = HeavySpeciesMod
    heavy_species_name = Om
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = -1
    diffusivity = 5.00e-9
    block = 1
  [../]

  [./O2m_mat]
    type = HeavySpeciesMod
    heavy_species_name = O2m
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = -1
    diffusivity = 5e-9
    block = 1
  [../]

  [./O3m_mat]
    type = HeavySpeciesMod
    heavy_species_name = O3m
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = -1
    diffusivity = 5.0e-9
    block = 1
  [../]

  [./O2_mat]
    type = HeavySpeciesMod
    heavy_species_name = O2
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 2e-9
    block = 1
  [../]

  #[./Nap_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = Nap
  #  heavy_species_mass = 3.8175407e-26
  #  heavy_species_charge = 1
  #  diffusivity = 2.0e-9
  #  block = 1
  #[../]
  #
  #[./Clm_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = Clm
  #  heavy_species_mass = 5.8871086e-26
  #  heavy_species_charge = -1
  #  diffusivity = 2.0e-9
  #  block = 1
  #[../]

  #[./NO2m_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = NO2m
  #  heavy_species_mass = 7.65e-26
  #  heavy_species_charge = -1
  #  diffusivity = 5.0e-9
  #  block = 1
  #[../]
  #
  #[./NO2_2m_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = NO2_2m
  #  heavy_species_mass = 7.65e-26
  #  heavy_species_charge = -2
  #  diffusivity = 5e-9
  #  block = 1
  #[../]
  #
  #[./NO3m_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = NO3m
  #  heavy_species_mass = 1.03e-25
  #  heavy_species_charge = -1
  #  diffusivity = 5e-9
  #  block = 1
  #[../]
  #
  #[./NO3_2m_mat]
  #  type = HeavySpeciesMod
  #  heavy_species_name = NO3_2m
  #  heavy_species_mass = 1.03e-25
  #  heavy_species_charge = -2
  #  diffusivity = 5e-9
  #  block = 1
  #[../]
  
  [./electron_data]
    type = GenericConstantMaterial
    prop_names = 'diffemliq muemliq sgnemliq'
    prop_values = '4.5e-9 0.000173913 -1'
    block = 1
  [../]

  [./N_A_mat]
    type = GenericConstantMaterial
    prop_names = 'N_A e diffpotential diffpotential_liq'
    prop_values = '6.022e23 1.602e-19 7.0832e-10 7.0832e-10'
    block = 1
  [../]

[]

[Reactions]
  [./gas_phase_reactions]
    species = 'em Arp Ar*'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    electron_energy = 'mean_en'
    electron_density = 'em'
    file_location = 'dc_2d_test'
    potential = 'potential'
    equation_variables = 'e_temp'
    position_units = ${dom0Scale}
    use_log = true
    block = 0


    reactions = 'em + Ar -> em + Ar : EEDF [elastic]
                 em + Ar -> em + Ar*  : EEDF [-11.5]
                 em + Ar -> em + em + Arp  : EEDF [-1.576e1]
                 em + Ar* -> em + Ar : EEDF [11.5]
                 em + Ar* -> em + em + Arp : EEDF [-4.3]
                 Arp + em + em -> em + Ar : {3.17314e9 * (e_temp)^(-4.5)}
                 Ar* + Ar + Ar -> Ar + Ar + Ar : 5.077028e3'
                 #Ar2p + em -> Ar* + Ar : {5.1187e11 * (e_temp*11600/300)^(-0.67)}
                 #Ar2p + Ar -> Arp + Ar + Ar : {3.649332e12 / 300 * exp(-15130/300)}
                 #Ar* + Ar* -> Ar2p + em : 3.6132e8
                 #Arp + Ar + Ar -> Ar2p + Ar : 81595.089'
  [../]

  [./liquid_phase_reactions]
    # removed Od1 and its two reactions:
    #
    # Od1 + H2O -> H2O2      : 1.8e7
    # Od1 + H2O -> OH + OH   : 2.3e-13
    #
    species = 'emliq H OHm H2Op OH H2 Om H3Op H2O2 HO2m O2 O2m O HO2 O3 O3m'
    #NO2m NO2_2m NO3m NO3_2m'
    aux_species = 'H2O O2'
    use_log = true
    position_units = ${dom1Scale}
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'emliq + H2O -> H + OHm : 1.9e-2
                 emliq + H2Op -> H + OH : 6e8
                 emliq + emliq -> H2 + OHm + OHm : 3.07029658079e8
                 emliq + H + H2O -> H2 + OHm     : 2.5e4
                 emliq + OH -> OHm               : 3e7
                 emliq + Om + H2O -> OHm + OHm   : 2.2e7
                 emliq + H3Op -> H + H2O         : 2.3e7
                 emliq + H2O2 -> OH + OHm        : 1.1e7
                 emliq + HO2m + H2O -> OH + OHm + OHm : 3.5e3
                 emliq + O2 -> O2m                    : 1.9e7
                 emliq + O -> Om                      : 1.9e7
                 H + H2O -> H2 + OH                   : 1e-2
                 H + H -> H2                          : 7.5e6
                 H + OH -> H2O                        : 7e6
                 H + OHm -> H2O + em                  : 2.2e4
                 H + H2O2 -> H2O + OH                 : 9e4
                 H2 + H2O2 -> H2O + OH + H            : 6e3
                 H + O2 -> HO2                        : 2.1e7
                 H + HO2 -> H2O2                      : 1e7
                 O + H2O -> OH + OH                   : 1.3e1
                 O + O2 -> O3                         : 3e6
                 OH + OH -> H2O2                      : 5.5e6
                 OH + Om -> HO2m                      : 2e8
                 OH + H2 -> H + H2O                   : 4.2e4
                 OH + OHm -> Om + H2O                 : 1.3e7
                 OH + HO2 -> O2 + H2O                 : 6.0e6
                 OH + O2m -> O2 + OHm                 : 8e6
                 Om + H2O -> OHm + OH                 : 1.8e3
                 Om + H2 -> OHm + H                   : 8e4
                 Om + H2O2 -> O2m + H2O               : 5e5
                 Om + HO2m -> O2m + OHm               : 4e5
                 Om + O2 -> O3m                       : 3.6e6
                 Om + O2m + H2O -> OHm + OHm + O2     : 6e2
                 OH + H2O2 -> H2O + HO2               : 2.7e4
                 OH + HO2m -> OHm + HO2               : 7.5e6
                 H2Op + H2O -> H3Op + OH              : 6.0
                 H3Op + OHm -> H2O + H + OH           : 6e7
                 HO2 + H2O -> H3Op + O2m              : 2.0
                 H3Op + O2m -> H2O + HO2              : 6.0e-2
                 H2O -> Hp + OHm                      : 1.3963e-3
                 Hp + OHm -> H2O                      : 1.3973e8'
                 #emliq + NO2m -> NO2_2m               : 5.2e6
                 #emliq + NO3m -> NO3_2m               : 7e6
  [../]
[]
