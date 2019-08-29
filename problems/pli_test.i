#dom0Scale=1e-3
#dom1Scale=1e-1
dom0Scale=1.0
dom1Scale=1.0

[GlobalParams]
  offset = 20
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  type = FileMesh
  file = 'test_mesh.msh'
[]

[MeshModifiers]
  [./interface]
    type = SideSetsBetweenSubdomains
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    # depends_on = 'box'
  [../]
  [./interface_again]
    type = SideSetsBetweenSubdomains
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
    # depends_on = 'box'
  [../]
  [./left]
    type = SideSetsFromNormals
    normals = '-1 0 0'
    new_boundary = 'left'
  [../]
  [./right]
    type = SideSetsFromNormals
    normals = '1 0 0'
    new_boundary = 'right'
  [../]
[]

[Problem]
  type = FEProblem
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
#    ksp_norm = none
  [../]
[]

[Executioner]
  type = Transient
  end_time = 1e-2
  # end_time = 10
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  # petsc_options = '-snes_test_display'
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -ksp_type -snes_linesearch_minlambda'
  petsc_options_value = 'lu NONZERO 1.e-10 fgmres 1e-3'
  # petsc_options_iname = '-snes_type'
  # petsc_options_value = 'test'
  nl_rel_tol = 1e-4
  nl_abs_tol = 7.6e-5
  dtmin = 1e-12
  l_max_its = 20
#  dtmax = 1e-4
  steady_state_check = true
  steady_state_tol = 1e-4
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    # dt = 1.1
    growth_factor = 1.2
    optimal_iterations = 15
  [../]
[]

[Outputs]
  # print_perf_log = true
  print_linear_residuals = false
  [./out]
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
    e = 1.6e-19
  [../]
[]

[Kernels]
  #[./OHm_g_time_deriv]
  #  type = ElectronTimeDerivative
  #  variable = OHm_g
  #  block = 0
  #[../]
  #[./OHm_g_advection]
  #  type = EFieldAdvection
  #  variable = OHm_g
  #  potential = potential
  #  block = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[./OHm_g_diffusion]
  #  type = CoeffDiffusion
  #  variable = OHm_g
  #  block = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[./OHm_g_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = OHm_g
  #  block = 0
  #[../]


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
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    block = 0
  [../]

  [./emliq_time_deriv]
    type = ElectronTimeDerivative
    variable = emliq
    block = 1
  [../]
  [./emliq_advection]
    type = EFieldAdvection
    variable = emliq
    potential = potential_liq
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
  [../]

  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential_liq
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
  [./emliq_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = emliq
    block = 1
  [../]
  [./OHm_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = OHm
    block = 1
  [../]
  [./O2m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = O2m
    block = 1
  [../]
  [./O3m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = O3m
    block = 1
  [../]
  [./HO2m_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = HO2m
    block = 1
  [../]
  [./Hp_charge_source]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = H+
    block = 1
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
  [../]
  # [./Arp_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = Arp
  #   potential = potential
  #   block = 0
  # [../]

  [./OHm_time_deriv]
    type = ElectronTimeDerivative
    variable = OHm
    block = 1
  [../]
  [./OHm_advection]
    type = EFieldAdvection
    variable = OHm
    potential = potential_liq
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
    block = 1
  [../]
  # [./OHm_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = OHm
  #   potential = potential
  #   block = 1
  # [../]

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
    offset = 15
  [../]
  # [./mean_en_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = mean_en
  #   potential = potential
  #   block = 0
  # [../]


  [./Om_time_deriv]
    type = ElectronTimeDerivative
    variable = Om
    block = 1
  [../]
  [./Om_advection]
    type = EFieldAdvection
    variable = Om
    potential = potential_liq
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
  [../]


  [./O2m_time_deriv]
    type = ElectronTimeDerivative
    variable = O2m
    block = 1
  [../]
  [./O2m_advection]
    type = EFieldAdvection
    variable = O2m
    potential = potential_liq
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
  [../]


  [./O3m_time_deriv]
    type = ElectronTimeDerivative
    variable = O3m
    block = 1
  [../]
  [./O3m_advection]
    type = EFieldAdvection
    variable = O3m
    potential = potential_liq
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
  [../]


  [./HO2m_time_deriv]
    type = ElectronTimeDerivative
    variable = HO2m
    block = 1
  [../]
  [./HO2m_advection]
    type = EFieldAdvection
    variable = HO2m
    potential = potential_liq
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
  [../]


  [./Hp_time_deriv]
    type = ElectronTimeDerivative
    variable = H+
    block = 1
  [../]
  [./Hp_advection]
    type = EFieldAdvection
    variable = H+
    potential = potential_liq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Hp_diffusion]
    type = CoeffDiffusion
    variable = H+
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Hp_log_stabilization]
    type = LogStabilizationMoles
    variable = H+
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
    potential = potential_liq
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
    type = LogStabilizationMoles
    variable = O
    block = 1
  [../]


 # [./Od1_time_deriv]
 #   type = ElectronTimeDerivative
 #   variable = Od1
 #   block = 1
 # [../]
 # [./Od1_advection]
 #   type = EFieldAdvection
 #   variable = Od1
 #   potential = potential
 #   block = 1
 #   position_units = ${dom1Scale}
 # [../]
 # [./Od1_diffusion]
 #   type = CoeffDiffusion
 #   variable = Od1
 #   block = 1
 #   position_units = ${dom1Scale}
 # [../]
 # [./Od1_log_stabilization]
 #   type = LogStabilizationMoles
 #   variable = Od1
 #   block = 1
 # [../]


  [./O2_1_time_deriv]
    type = ElectronTimeDerivative
    variable = O2_1
    block = 1
  [../]
  [./O2_1_advection]
    type = EFieldAdvection
    variable = O2_1
    potential = potential_liq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O2_1_diffusion]
    type = CoeffDiffusion
    variable = O2_1
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./O2_1_log_stabilization]
    type = LogStabilizationMoles
    variable = O2_1
    block = 1
  [../]


 # [./O2_time_deriv]
 #   type = ElectronTimeDerivative
 #   variable = O2
 #   block = 1
 # [../]
 # [./O2_advection]
 #   type = EFieldAdvection
 #   variable = O2
 #   potential = potential
 #   block = 1
 #   position_units = ${dom1Scale}
 # [../]
 # [./O2_diffusion]
 #   type = CoeffDiffusion
 #   variable = O2
 #   block = 1
 #   position_units = ${dom1Scale}
 # [../]
 # [./O2_log_stabilization]
 #   type = LogStabilizationMoles
 #   variable = O2
 #   block = 1
 # [../]


  [./O3_time_deriv]
    type = ElectronTimeDerivative
    variable = O3
    block = 1
  [../]
  [./O3_advection]
    type = EFieldAdvection
    variable = O3
    potential = potential_liq
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
    potential = potential_liq
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
    potential = potential_liq
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
    potential = potential_liq
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
    potential = potential_liq
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
    potential = potential_liq
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
    potential = potential_liq
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
    type = LogStabilizationMoles
    variable = H2O2
    block = 1
  [../]
[]

[Variables]
  [./potential]
    block = 0
  [../]
  [./potential_liq]
    block = 1
  [../]

  [./em]
    block = 0
  [../]
  [./emliq]
    block = 1
    # scaling = 1e-5
  [../]

  [./Arp]
    block = 0
  [../]

  [./mean_en]
    block = 0
    # scaling = 1e-1
  [../]

 # [./OHm_g]
 #   block = 0
 #   initial_condition = -31
 # [../]
  [./OHm]
    block = 1
    # scaling = 1e-5
  [../]
  [./Om]
    block = 1
    initial_condition = -22
  [../]

  [./O2m]
    block = 1
    initial_condition = -22
  [../]

  [./O3m]
    block = 1
    initial_condition = -22
  [../]


  [./HO2m]
    block = 1
    initial_condition = -22
  [../]


  [./H+]
    block = 1
    initial_condition = -22
  [../]


  [./O]
    block = 1
    initial_condition = -22
  [../]


  #[./Od1]
  #  block = 1
  #  initial_condition = -22
  #[../]


  [./O2_1]
    block = 1
    initial_condition = -22
  [../]


 # [./O2]
 #   block = 1
 #   initial_condition = -0.609203
 # [../]


  [./O3]
    block = 1
    initial_condition = -22
  [../]


  [./H]
    block = 1
    initial_condition = -22
  [../]


  [./H2]
    block = 1
    initial_condition = -22
  [../]


  [./HO2]
    block = 1
    initial_condition = -22
  [../]


  [./HO3]
    block = 1
    initial_condition = -22
  [../]


  [./OH]
    block = 1
    initial_condition = -22
  [../]


  [./H2O2]
    block = 1
    initial_condition = -22
  [../]
[]

[AuxVariables]
  # BACKGROUND WATER - assuming it won't actually be affected by reactions
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
  [./H2O_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
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
  [./rholiq]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./em_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./emliq_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Arp_lin]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./OHm_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
 # [./OHm_g_lin]
 #   block = 0
 #   order = CONSTANT
 #   family = MONOMIAL
 # [../]
  [./Efield_gas]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Efield_liq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Current_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Current_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./Current_Arp]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./Current_OHm]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tot_gas_current]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./tot_liq_current]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tot_flux_OHm]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./EFieldAdvAux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./DiffusiveFlux_em]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./EFieldAdvAux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./DiffusiveFlux_emliq]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./PowerDep_em]
   order = CONSTANT
   family = MONOMIAL
   block = 0
  [../]
  [./PowerDep_Arp]
   order = CONSTANT
   family = MONOMIAL
   block = 0
  [../]
  [./ProcRate_el]
   order = CONSTANT
   family = MONOMIAL
   block = 0
  [../]
  [./ProcRate_ex]
   order = CONSTANT
   family = MONOMIAL
   block = 0
  [../]
  [./ProcRate_iz]
   order = CONSTANT
   family = MONOMIAL
   block = 0
  [../]
  [./Ar]
    order = CONSTANT
    family = MONOMIAL
    block = 0
    initial_condition = 3.704261
  [../]

  ### LIQUID DENSITIES
  
  [./Om_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O2m_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O3m_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./HO2m_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Hp_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

 # [./Od1_lin]
 #   block = 1
 #   order = CONSTANT
 #   family = MONOMIAL
 # [../]

  [./O2_1_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O2_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O3_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./H_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./H2_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./HO2_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./HO3_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./OH_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./H2O2_lin]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./PowerDep_em]
    type = PowerDep
    density_log = em
    potential = potential
    art_diff = false
    potential_units = kV
    variable = PowerDep_em
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./PowerDep_Arp]
    type = PowerDep
    density_log = Arp
    potential = potential
    art_diff = false
    potential_units = kV
    variable = PowerDep_Arp
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./e_temp]
    type = ElectronTemperature
    variable = e_temp
    electron_density = em
    mean_en = mean_en
    block = 0
  [../]
  [./x_g]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./x_l]
    type = Position
    variable = x
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./x_nl]
    type = Position
    variable = x_node
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./rho]
    type = ParsedAux
    variable = rho
    args = 'em_lin Arp_lin'
    function = 'Arp_lin - em_lin'
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./rholiq]
    type = ParsedAux
    variable = rholiq
    args = 'emliq_lin OHm_lin' # H3Op_lin OHm_lin'
    function = '-emliq_lin - OHm_lin' # 'H3Op_lin - em_lin - OHm_lin'
    execute_on = 'timestep_end'
    block = 1
  [../]
  [./tot_gas_current]
    type = ParsedAux
    variable = tot_gas_current
    args = 'Current_em Current_Arp'
    function = 'Current_em + Current_Arp'
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./tot_liq_current]
    type = ParsedAux
    variable = tot_liq_current
    args = 'Current_emliq Current_OHm' # Current_H3Op Current_OHm'
    function = 'Current_emliq + Current_OHm' # + Current_H3Op + Current_OHm'
    execute_on = 'timestep_end'
    block = 1
  [../]
  [./em_lin]
    type = DensityMoles
    convert_moles = true
    variable = em_lin
    density_log = em
    block = 0
  [../]
  [./emliq_lin]
    type = DensityMoles
    convert_moles = true
    variable = emliq_lin
    density_log = emliq
    block = 1
  [../]
  [./Arp_lin]
    type = DensityMoles
    convert_moles = true
    variable = Arp_lin
    density_log = Arp
    block = 0
  [../]
  [./OHm_lin]
    type = DensityMoles
    convert_moles = true
    variable = OHm_lin
    density_log = OHm
    block = 1
  [../]
 # [./OHm_g_lin]
 #   type = DensityMoles
 #   convert_moles = true
 #   variable = OHm_g_lin
 #   density_log = OHm_g
 #   block = 0
 # [../]
  [./Efield_g]
    type = Efield
    component = 0
    potential = potential
    variable = Efield_gas
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Efield_l]
    type = Efield
    component = 0
    potential = potential_liq
    variable = Efield_liq
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./Current_em]
    type = Current
    potential = potential
    density_log = em
    variable = Current_em
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Current_emliq]
    type = Current
    potential = potential_liq
    density_log = emliq
    variable = Current_emliq
    art_diff = false
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./Current_Arp]
    type = Current
    potential = potential
    density_log = Arp
    variable = Current_Arp
    art_diff = false
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./Current_OHm]
    block = 1
    type = Current
    potential = potential_liq
    density_log = OHm
    variable = Current_OHm
    art_diff = false
    position_units = ${dom1Scale}
  [../]
  [./tot_flux_OHm]
    block = 1
    type = TotalFlux
    potential = potential_liq
    density_log = OHm
    variable = tot_flux_OHm
  [../]
  [./EFieldAdvAux_em]
    type = EFieldAdvAux
    potential = potential
    density_log = em
    variable = EFieldAdvAux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./DiffusiveFlux_em]
    type = DiffusiveFlux
    density_log = em
    variable = DiffusiveFlux_em
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./EFieldAdvAux_emliq]
    type = EFieldAdvAux
    potential = potential_liq
    density_log = emliq
    variable = EFieldAdvAux_emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./DiffusiveFlux_emliq]
    type = DiffusiveFlux
    density_log = emliq
    variable = DiffusiveFlux_emliq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./H2O_lin]
    type = DensityMoles
    convert_moles = true
    variable = H2O_lin
    density_log = H2O
    block = 1 
  [../]
  [./Om_lin]
    type = DensityMoles
    convert_moles = true
    variable = Om_lin
    density_log = Om
    block = 1 
  [../]
  [./O2m_lin]
    type = DensityMoles
    convert_moles = true
    variable = O2m_lin
    density_log = O2m
    block = 1 
  [../]
  [./O3m_lin]
    type = DensityMoles
    convert_moles = true
    variable = O3m_lin
    density_log = O3m
    block = 1 
  [../]
  [./HO2m_lin]
    type = DensityMoles
    convert_moles = true
    variable = HO2m_lin
    density_log = HO2m
    block = 1 
  [../]
  [./Hp_lin]
    type = DensityMoles
    convert_moles = true
    variable = Hp_lin
    density_log = H+
    block = 1 
  [../]
  [./O_lin]
    type = DensityMoles
    convert_moles = true
    variable = O_lin
    density_log = O
    block = 1 
  [../]
 # [./Od1_lin]
 #   type = DensityMoles
 #   convert_moles = true
 #   variable = Od1_lin
 #   density_log = Od1
 #   block = 1 
 # [../]
  [./O2_1_lin]
    type = DensityMoles
    convert_moles = true
    variable = O2_1_lin
    density_log = O2_1
    block = 1 
  [../]
  [./O2_lin]
    type = DensityMoles
    convert_moles = true
    variable = O2_lin
    density_log = O2
    block = 1 
  [../]
  [./O3_lin]
    type = DensityMoles
    convert_moles = true
    variable = O3_lin
    density_log = O3
    block = 1 
  [../]
  [./H_lin]
    type = DensityMoles
    convert_moles = true
    variable = H_lin
    density_log = H
    block = 1 
  [../]
  [./H2_lin]
    type = DensityMoles
    convert_moles = true
    variable = H2_lin
    density_log = H2
    block = 1 
  [../]
  [./HO2_lin]
    type = DensityMoles
    convert_moles = true
    variable = HO2_lin
    density_log = HO2
    block = 1 
  [../]
  [./HO3_lin]
    type = DensityMoles
    convert_moles = true
    variable = HO3_lin
    density_log = HO3
    block = 1 
  [../]
  [./OH_lin]
    type = DensityMoles
    convert_moles = true
    variable = OH_lin
    density_log = OH
    block = 1 
  [../]
  [./H2O2_lin]
    type = DensityMoles
    convert_moles = true
    variable = H2O2_lin
    density_log = H2O2
    block = 1 
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

  [./water_interface]
    type = InterfaceFluxConservation
    neighbor_var = potential_liq
    variable = potential
    boundary = 'master0_interface'
    region_name = 'diffpotential'
    neighbor_region_name = 'diffpotential_liq'
    position_units = ${dom0Scale}
    neighbor_position_units = ${dom1Scale}
  [../]
[]

[BCs]
  [./potential_left]
    type = NeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = left
    function = potential_bc_func
    ip = Arp
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_dirichlet_right]
    type = DirichletBC
    variable = potential_liq
    boundary = right
    value = 0
  [../]
  [./potential_interface]
    type = MatchedValueBC
    variable = potential
    v = potential_liq
    boundary = 'master0_interface'
  [../]
  [./em_physical_right]
    type = HagelaarElectronBC
    variable = em
    boundary = 'master0_interface'
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0.99
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = 'master0_interface'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_right_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = 'master0_interface'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_right]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master0_interface'
    potential = potential
    em = em
    ip = Arp
    r = 0.99
    position_units = ${dom0Scale}
  [../]
  [./em_physical_left]
    type = HagelaarElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    type = SecondaryElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    ip = Arp
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_diffusion]
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = 'left'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_physical_left_advection]
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = 'left'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    ip = Arp
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./emliq_right]
    type = DCIonBC
    variable = emliq
    boundary = right
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./OHm_physical]
    type = DCIonBC
    variable = OHm
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Om_physical]
    type = DCIonBC
    variable = Om
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O2m_physical]
    type = DCIonBC
    variable = O2m
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O3m_physical]
    type = DCIonBC
    variable = O3m
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./HO2m_physical]
    type = DCIonBC
    variable = HO2m
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./Hp_physical]
    type = DCIonBC
    variable = H+
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./O_physical]
    type = DCIonBC
    variable = O
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
 # [./Od1_physical]
 #   type = DCIonBC
 #   variable = Od1
 #   boundary = 'right'
 #   potential = potential
 #   position_units = ${dom1Scale}
 # [../]
  [./O2_1_physical]
    type = DCIonBC
    variable = O2_1
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
 # [./O2_physical]
 #   type = DCIonBC
 #   variable = O2
 #   boundary = 'right'
 #   potential = potential
 #   position_units = ${dom1Scale}
 # [../]
  [./O3_physical]
    type = DCIonBC
    variable = O3
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H_physical]
    type = DCIonBC
    variable = H
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2_physical]
    type = DCIonBC
    variable = H2
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./HO2_physical]
    type = DCIonBC
    variable = HO2
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./HO3_physical]
    type = DCIonBC
    variable = HO3
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./OH_physical]
    type = DCIonBC
    variable = OH
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  [./H2O2_physical]
    type = DCIonBC
    variable = H2O2
    boundary = 'right'
    potential = potential
    position_units = ${dom1Scale}
  [../]
  #[./OHm_physical_right_diffusion]
  #  type = HagelaarIonDiffusionBC
  #  variable = OHm
  #  boundary = 'master1_interface'
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  #[./OHm_physical_right_advection]
  #  type = HagelaarIonAdvectionBC
  #  variable = OHm
  #  boundary = 'master1_interface'
  #  potential = potential
  #  r = 0
  #  position_units = ${dom1Scale}
  #[../]
  
[]

[ICs]
  [./em_ic]
    type = ConstantIC
    variable = em
    value = -21
    block = 0
  [../]
  [./emliq_ic]
    type = ConstantIC
    variable = emliq
    value = -21
    block = 1
  [../]
  [./Arp_ic]
    type = ConstantIC
    variable = Arp
    value = -21
    block = 0
  [../]
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -20
    block = 0
  [../]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
    block = 0
  [../]
  [./potential_liq_ic]
    type = FunctionIC
    variable = potential_liq
    function = potential_ic_func
    block = 1
  [../]
  [./OHm_ic]
    type = ConstantIC
    variable = OHm
    value = -15.6
    block = 1
  [../]
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    # value = '1.25*tanh(1e6*t)'
    value = 1.25
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-1.25 * (1.0001e-3 - x)'
  [../]
[]

[Materials]
 [./GasBasics]
   type = GasBase
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   # user_p_gas = 1.01325e5
   user_p_gas = 101325
   em = em
   potential = potential
   mean_en = mean_en
   user_se_coeff = 0.05
   property_tables_file = cpc_test/e_vals_test.txt
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
 #[./water_block]
 #  type = Water
 #  block = 1
 #  potential = potential
 #[../]
  [./Ar_species]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  [../]

  # WATER SPECIES MATERIALS
  [./O_mat]
    type = HeavySpeciesMod
    heavy_species_name = O
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 2.00e-9
    block = 1
  [../]

 # [./Od1_mat]
 #   type = HeavySpeciesMod
 #   heavy_species_name = Od1
 #   heavy_species_mass = 2.6566962e-26
 #   heavy_species_charge = 0
 #   diffusivity = 2.00e-9
 #   block = 1
 # [../]

  [./O2_1_mat]
    type = HeavySpeciesMod
    heavy_species_name = O2_1
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 1.97e-9
    block = 1
  [../]

  [./O3_mat]
    type = HeavySpeciesMod
    heavy_species_name = O3
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = 0
    diffusivity = 1.75e-9
    block = 1
  [../]

  [./OH_mat]
    type = HeavySpeciesMod
    heavy_species_name = OH
    heavy_species_mass = 2.82431e-26
    heavy_species_charge = 0
    diffusivity = 2.30e-9
    block = 1
  [../]

  [./HO2_mat]
    type = HeavySpeciesMod
    heavy_species_name = HO2
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = 0
    diffusivity = 1.00e-9
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
    diffusivity = 1.00e-9
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

  [./H+_mat]
    type = HeavySpeciesMod
    heavy_species_name = H+
    heavy_species_mass = 1.67376e-26
    heavy_species_charge = 1
    diffusivity = 9.31e-9
    block = 1
  [../]

  ##### THIS ONE IS DECLARED IN WATER.C
  [./OHm_mat]
    type = HeavySpeciesMod
    heavy_species_name = OHm 
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.26e-9
    block = 1
  [../]

  [./HO2m_mat]
    type = HeavySpeciesMod
    heavy_species_name = HO2m
    heavy_species_mass = 5.481026e-26
    heavy_species_charge = -1
    diffusivity = 1.00e-9
    block = 1
  [../]

  [./Om_mat]
    type = HeavySpeciesMod
    heavy_species_name = Om
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = -1
    diffusivity = 2.00e-9
    block = 1
  [../]

  [./O2m_mat]
    type = HeavySpeciesMod
    heavy_species_name = O2m
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = -1
    diffusivity = 1.97e-9
    block = 1
  [../]

  [./O3m_mat]
    type = HeavySpeciesMod
    heavy_species_name = O3m
    heavy_species_mass = 7.97047e-26
    heavy_species_charge = -1
    diffusivity = 1.75e-9
    block = 1
  [../]

  [./O2_mat]
    type = HeavySpeciesMod
    heavy_species_name = O2
    heavy_species_mass = 5.31365e-26
    heavy_species_charge = 0
    diffusivity = 1.97e-9
    block = 1
  [../]

  [./emliq_mat]
    type = HeavySpeciesMod
    heavy_species_name = emliq
    heavy_species_mass = 9.11e-31
    heavy_species_charge = -1
    diffusivity = 1.00e-5
    block = 1
  [../]

  [./N_A_mat]
    type = GenericConstantMaterial
    prop_names = 'N_A e diffpotential diffpotential_liq'
    prop_values = '6.022e23 1.602e-19 7.0832e-10 7.0832e-10'
    #prop_values = '6.022e23 1.602e-19 8.854e-12 7.0832e-10'
    #prop_names = 'N_A e'
    #prop_values = '6.022e23 1.602e-19'
    block = 1
  [../]
[]

[Reactions]
  [./gas_phase_reactions]
    species = 'em Arp'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    electron_energy = 'mean_en'
    electron_density = 'em'
    file_location = 'cpc_test'
    potential = 'potential'
    position_units = ${dom0Scale}
    use_log = true
    block = 0

    reactions = 'em + Ar -> em + Ar : EEDF [elastic]
                 em + Ar -> em + Ar*  : EEDF [-11.5]
                 em + Ar -> em + em + Arp  : EEDF [-1.576e1]'
  [../]

  [./liquid_phase_reactions]
    # removed Od1 and its two reactions:
    #
    # Od1 + H2O -> H2O2      : 1.8e7
    # Od1 + H2O -> OH + OH   : 2.3e-13
    #
    species = 'emliq OHm Om O2m O3m HO2m H+ O O2_1 O3 H H2 HO2 HO3 OH H2O2'
    aux_species = 'H2O O2'
    use_log = true
    position_units = ${dom1Scale}
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'H2O -> H+ + OHm  : 1.4e-3
                 H+ + OHm -> H2O  : 1.4e8
                 H2O2 -> H+ + HO2m : 1.12e-1
                 H+ + HO2m -> H2O2 : 5e7
                 HO2 -> O2m + H+   : 1.35e6
                 O2m + H+ -> HO2   : 5e7
                 OH -> Om + H+     : 1.26e-1
                 Om + H+ -> OH     : 1e8
                 OH + OHm -> H2O + Om : 1.3e7
                 H2O + Om -> OH + OHm : 1.7e3
                 H2O2 + OHm -> HO2m + H2O : 1.3e7
                 HO2m + H2O -> H2O2 + OHm : 5.8e4
                 emliq + H2O -> H + OHm : 1.9e-2
                 H + OHm -> H2O + emliq : 2.2e4
                 OH + OHm -> Om + H2O   : 1.3e7
                 Om + H2O -> OH + OHm   : 1.03e5
                 HO2 + OHm -> O2m + H2O : 5.0e7
                 O2m + H2O -> HO2 + OHm : 18.5767e-3
                 Om + O2 -> O3m         : 3.6e6
                 O3m -> Om + O2         : 3.3e3
                 OH + OH -> H2O2        : 3.6e6
                 H2O2 -> OH + OH        : 2.3e-7
                 H -> emliq + H+        : 3.9
                 emliq + H+ -> H        : 2.3e7
                 O + O2 -> O3           : 4.0e6
                 O + O -> O2            : 2.8e7
                 O2_1 + H2O -> O2 + H2O : 4.9
                 O2_1 + OH -> O2 + OH   : 2.2
                 emliq + OH -> OHm      : 3.0e7
                 emliq + H2O2 -> OH + OHm : 1.1e7
                 emliq + HO2 -> HO2m    : 2.0e7
                 emliq + O2 -> O2m      : 1.9e7
                 emliq + HO2m -> Om + OHm : 3.5e7
                 emliq + O3 -> O3m      : 3.6e7
                 H + H2O -> H2 + OH     : 1.1e7
                 H + Om -> OHm          : 1.0e7
                 H + HO2m -> OHm + OH   : 9.0e7
                 H + O3m -> OHm + O2    : 1.0e7
                 H + H -> H2            : 7.8e6
                 H + OH -> H2O          : 7.0e6
                 H + O2 -> HO2          : 2.1e7
                 H + H2O2 -> OH + H2O   : 9.0e4
                 H + HO2 -> H2O2        : 1.8e7
                 H + O2m -> HO2m        : 1.8e7
                 H + O3 -> HO3          : 3.8e7
                 OH + HO2 -> O2 + H2O   : 6.0e6
                 OH + O2m -> O2 + OHm   : 8.2e6
                 OH + H2 -> H + H2O     : 4.3e4
                 OH + H2O2 -> HO2 + H2O : 2.7e4
                 OH + Om -> HO2m        : 2.5e7
                 OH + HO2m -> HO2 + OHm : 6.0e6
                 OH + O3m -> O3 + OHm   : 2.6e6
                 OH + O3m -> O2 + O2 + H+ : 6.0e6
                 OH + O3 -> HO2 + O2    : 1.1e5
                 HO2 + O2m -> HO2m + O2 : 8.0e4
                 HO2 + HO2 -> H2O2 + O2 : 7.0e2
                 HO2 + Om -> O2 + OHm   : 6.0e6
                 HO2 + H2O2 -> OH + O2 + H2O : 0.5e-3
                 HO2 + HO2m -> OHm + OH + O2 : 0.5e-3
                 O2m + H2O2 -> OH + O2 + OHm : 0.13e-3
                 O2m + HO2m -> Om + O2 + OHm : 0.13e-3
                 O2m + O3 -> O3m + O2        : 1.5e6
                 Om + H2 -> H + OHm          : 8.0e4
                 Om + H2O2 -> H2O + O2m      : 5.0e5
                 Om + HO2m -> OHm + O2m      : 4.0e5
                 Om + O3m -> O2m + O2m       : 7.0e5
                 Om + O3 -> O2m + O2         : 5.0e6
                 O3m + H -> O2 + OH          : 9.0e6
                 HO3 -> O2 + OH              : 1.0e5
                 O + OHm -> HO2m             : 1.1e2
                 O + H2O2 -> OH + HO2        : 1.6e2
                 O + HO2m -> OH + O2m        : 5.3e6
                 O3 + H2O2 -> OH + HO2 + O2  : 3.0e6
                 emliq + O2m -> HO2m + OHm   : 1.3e7
                 emliq + H -> H2 + OHm       : 2.5e7
                 emliq + Om -> OHm + OHm     : 2.2e7
                 emliq + O3m -> O2 + OH + OH : 1.6e7
                 Om + O2m -> OHm + OHm + O2  : 6.0e4
                 O2m + O3m -> OHm + OHm + O2 + O2 : 1.0e1
                 Om + Om -> OHm + HO2m       : 1.0e6
                 emliq + emliq -> H2 + OHm + OHm : 5.5e6
                 O2m + O2m -> H2O2 + O2 + OHm + OHm : 1.0e-1'

    # emliq + emliq + H2O + H2O : H2 + OHm + OHm : {5.5e9/(H2O^2)}
    # O2m + O2m + H2O + H2O -> H2O2 + O2 + OHm + OHm : {1.0e2 / H2O^2}
    # emliq + O2m + H2O -> HO2m + OHm : {1.3e10 / 55410.26}
    # emliq + H + H2O -> H2 + OHm : {2.5e10 / 55410.26}
    # emliq + Om + H2O -> OHm + OHm : {2.2e10/ 55410.26}
    # emliq + O3m + H2O -> O2 + OH + OH : {1.6e10 / 55410.26}
    # Om + O2m + H2O -> OHm + OHm + O2  : {6.0e8 / 55410.26}
    # O2m + O3m + H2O -> OHm + OHm + O2 + O2 : {1.0e4 / 55410.26}
    # Om + Om + H2O -> OHm + HO2m  : {1.0e9 / 55410.26}'



#    reactions = 'emliq -> OHm : 1069.6
#                 emliq + emliq -> OHm + OHm : 3.136e8'
                 #emliq + emliq + OHm -> emliq + emliq : 1'
  [../]
[]

