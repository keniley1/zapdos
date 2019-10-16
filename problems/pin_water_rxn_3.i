#dom0Scale=1.0e-3
dom0Scale=1.0
dom1Scale=1.0

#[Adaptivity]
#    marker = errorfrac
#    max_h_level = 2
#
#    [./Indicators]
#        [./error]
#            type = GradientJumpIndicator
#            variable = Ar*
#        [../]
#    [../]
#
#    [./Markers]
#        [./errorfrac]
#            type = ErrorFractionMarker
#            indicator = error
#            refine = 0.55
#            coarsen = 0.1
#            outputs = none
#        [../]
#    [../]
#[]
[Adaptivity]
  initial_steps = 0
  cycles_per_step = 1
  marker = combo
  initial_marker = combo
  max_h_level = 2
  [Indicators/indicator1]
    type = GradientJumpIndicator
    variable = em
  []
  [Indicators/indicator2]
    type = GradientJumpIndicator
    variable = Ar*
  []
  [Indicators/indicator3]
    type = GradientJumpIndicator
    variable = emliq
  []
  [Markers]
    [./combo]
      type = ComboMarker
      markers = 'marker1 marker2 marker3'
    [../]
    [./marker1]
      type = ErrorFractionMarker
      indicator = indicator1
      coarsen = 0.1
      refine = 0.7
    [../]
    [./marker2]
      type = ErrorFractionMarker
      indicator = indicator2
      coarsen = 0.1
      refine = 0.7
    [../]
    [./marker3]
      type = ErrorFractionMarker
      indicator = indicator3
      coarsen = 0.1
      refine = 0.7
    [../]
  []
[]

[GlobalParams]
  offset = 24
  potential_units = kV
  use_moles = true
[]

[Mesh]
  type = FileMesh
  #file = 'pin_water_large03.msh'
  #file = 'pin_water_large02.msh'
  file = 'pin_water_rxn_3.msh'
  #construct_side_list_from_node_list=true 
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
  active = smp

  [./smp]
    type = SMP
    full = true
  [../]

  [./fsp_schur]
    type = FSP
    topsplit = 'pc'
    full = 'true'
    [./pc]
      splitting = 'p c'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type'
      petsc_options_value = 'hypre'
      #petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
      #petsc_options_value = 'hypre boomeramg 0.25'
    [../]
    [./c]
      vars = 'em Arp mean_en'
      #petsc_options_iname = '-pc_type -sub_pc_type -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      #petsc_options_value = 'bjacobi ilu 100 NONZERO 1e-10'
      #petsc_options_iname = '-pc_type -sub_pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      #petsc_options_value = 'bjacobi ilu super_dist 100 NONZERO 1e-10'
      #petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_minlambda -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      #petsc_options_value = 'bjacobi ilu 1e-3 100 NONZERO 1e-10'

      petsc_options_iname = '-pc_factor_mat_solver_package -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -snes_linesearch_minlambda'
      petsc_options_value = 'mumps asm 4 ilu NONZERO 1e-10 1e-3'
    [../]
  [../]

  [./fsp]
    type = FSP
    topsplit = 'pc'
    full = 'true'
    [./pc]
      splitting = 'p c'
      splitting_type = additive
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
      petsc_options_value = 'hypre boomeramg 0.25'
      #petsc_options_iname = '-pc_type -pc_hypre_type'
      #petsc_options_value = 'hypre boomeramg'
      #petsc_options_iname = '-pc_type -pc_gamg_agg_nsmooths -pc_gamg_threshold -mg_levels_ksp-type -mg_levels_pc_type -mg_levels_ksp_max_it'
      #petsc_options_value = 'gamg 1 0.02 richardson jacobi 2'
    [../]
    [./c]
      vars = 'em Arp mean_en'
      petsc_options_iname = '-pc_type -sub_pc_type -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      petsc_options_value = 'bjacobi ilu 100 NONZERO 1e-10'
      #petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_minlambda -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      #petsc_options_value = 'bjacobi ilu 1e-3 100 NONZERO 1e-10'

      #petsc_options_iname = '-pc_factor_mat_solver_package -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -snes_linesearch_minlambda'
      #petsc_options_value = 'mumps asm 4 ilu NONZERO 1e-10 1e-3'
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
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor -pc_svd_monitor -snes_mf_operator'
  #solve_type = newton 
  solve_type = pjfnk
  petsc_options_iname = '-snes_linesearch_type -pc_type -pc_factor_shift_type'
  petsc_options_value = 'basic lu NONZERO'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  #petsc_options_value = 'lu NONZERO 1e-10'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  #petsc_options_value = 'svd NONZERO 1e-10'


  #petsc_options_iname = '-pc_type -sub_pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
  #petsc_options_value = 'bjacobi ilu mumps 100 NONZERO 1e-10'
  #petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap'
  #petsc_options_value = 'gmres asm lu 100 NONZERO 8'
  #petsc_options_iname = '-sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
  #petsc_options_value = 'NONZERO 1e-10'
  #petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_minlambda -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
  #petsc_options_value = 'bjacobi ilu 1e-3 NONZERO 1e-10'
  #petsc_options_iname = '-pc_factor_mat_solver_package -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -snes_linesearch_minlambda'
  #petsc_options_value = 'mumps asm 4 ilu NONZERO 1e-10 1e-3'
  nl_rel_tol = 1e-5
  #nl_abs_tol = 7.6e-5
  dtmin = 1e-16
  l_max_its = 20
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-12
    growth_factor = 1.2
    optimal_iterations = 15
  [../]
[]

[Outputs]
  checkpoint = true
  perf_graph = true
  exodus = true
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
  #[./em_advection_stabilization]
  #  type = EFieldArtDiff
  #  variable = em
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]

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
  #[./Ar2p_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential
  #  charged = Ar2p
  #  block = 0
  #[../]
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
    offset = 31
  [../]
  #[./Arp_advection_stabilization]
  #  type = EFieldArtDiff
  #  variable = Arp
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]

  #[./Ar2p_time_deriv]
  #  type = ElectronTimeDerivative
  #  variable = Ar2p
  #  block = 0
  #[../]
  #[./Ar2p_advection]
  #  type = EFieldAdvection
  #  variable = Ar2p
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  block = 0
  #[../]
  #[./Ar2p_diffusion]
  #  type = CoeffDiffusion
  #  variable = Ar2p
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Ar2p
  #  block = 0
  #  offset = 24
  #  #offset = 25
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
    offset = 15
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
    block = 1
  [../]
[]

[Variables]
  [./potential]
    #order=SECOND
    #family=LAGRANGE
    #initial_from_file_var = potential
    #block = 0
  [../]

  [./em]
    #order=SECOND
    #family=LAGRANGE
    #initial_from_file_var = em
    initial_condition = -21
    block = 0
  [../]
  
  [./Arp]
    #order=SECOND
    #family=LAGRANGE
    #initial_from_file_var = Arp
    #initial_condition = -31
    initial_condition = -24
    block = 0
  [../]

  #[./Ar2p]
  #  #order=SECOND
  #  #family=LAGRANGE
  #  #initial_from_file_var = Arp
  #  initial_condition = -21
  #  block = 0
  #[../]
  
  [./Ar*]
    #order=SECOND
    #family=LAGRANGE
    #initial_from_file_var = Arp
    initial_condition = -24
    block = 0
  [../]

  [./mean_en]
    #order=SECOND
    #family=LAGRANGE
    #initial_from_file_var = mean_en
    block = 0
  [../]

  [./emliq]
    block = 1
    initial_condition = -20
    # scaling = 1e-5
  [../]

  [./OHm]
    block = 1
    initial_condition = -20
    # scaling = 1e-5
  [../]
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
  #[./Ar2p_physical_bottom_diffusion]
  #  type = HagelaarIonDiffusionBC
  #  variable = Ar2p
  #  boundary = bottom 
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_physical_bottom_advection]
  #  type = HagelaarIonAdvectionBC
  #  variable = Ar2p
  #  boundary = bottom
  #  potential = potential
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
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
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    #value = -23
    value = -20
    block = 0
  [../]
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
    value = '-0.75*tanh(1e9*t) - 0.05'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-1.25 * (y - 1.000e-3)'
    # Below is a test of a gaussian initial condition
    #vars = 'sigmax sigmay A x0 y0'
    #vals = '4e-4 0.5e-3 -2.0 0 2.5e-3'
    #value = 'A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2))'
    vars = 'sigmax sigmay A x0 y0'
    vals = '4e-1 0.9e0 -1.0 0 2.5e0'
    value = 'A*exp(-((x-x0)^2 / (2*sigmax^2) + (y-y0)^2/(2*sigmay)^2)) - 0.3 * (y - 1.000e-3)'
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

 [./water_block]
   type = Water
   block = 1
   potential = potential
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
    species = 'emliq OHm'
    use_log = true
    position_units = ${dom1Scale}
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'emliq -> OHm : 1069.6
                 emliq + emliq -> OHm + OHm : 3.136e8'
  [../]
[]
