dom0Scale=1.0
dom1Scale=1.0

#dom0Scale=1e-3
#dom1Scale=1e-7
[GlobalParams]
  offset = 20
  potential_units = kV
  use_moles = true
[]

#[Mesh]
#  type = GeneratedMesh
#  dim = 1
#  #nx = 17792
#  nx = 71168
#  xmin = 0
#  xmax = 0.0011
#[]


[Mesh]
  # type = GeneratedMesh
  # nx = 2
  # xmax = 1.1
  # dim = 1
  type = FileMesh
  file = 'liquidNew2.msh'
[]

[MeshModifiers]
  [./subdomain0]
    type = ParsedSubdomainMeshModifier
    combinatorial_geometry = 'x < 1e-3'
    #combinatorial_geometry = 'x <= 1'
    block_id = 0 
  [../]
  [./subdomain1]
    type = ParsedSubdomainMeshModifier
    combinatorial_geometry = 'x >= 1e-3'
    #combinatorial_geometry = 'x > 1'
    block_id = 1
    depends_on = subdomain0
  [../]

  [./interface]
    type = SideSetsBetweenSubdomains
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    depends_on = 'subdomain0 subdomain1' 
  [../]
  [./interface_again]
    type = SideSetsBetweenSubdomains
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
    depends_on = 'subdomain0 subdomain1'
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
  # kernel_coverage_check = false
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
  #end_time = 5e-9
  #num_steps = 10
  end_time = 1e-1
  #line_search = 'none'
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = newton 
  petsc_options_iname = '-pc_type -snes_linesearch_minlambda -ksp_type'
  petsc_options_value = 'lu 1e-3 fgmres'
  #petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_minlambda -ksp_type'
  #petsc_options_value = 'bjacobi lu 1e-3 gmres'
  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -snes_linesearch_minlambda'
  #petsc_options_value = 'asm 2 ilu NONZERO 1e-10 1e-3'
  nl_rel_tol = 1e-4
  nl_abs_tol = 7.6e-5
  dtmin = 1e-14
  l_max_its = 20
  #dt = 1e-11
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
  perf_graph = true
  # print_linear_residuals = false
  #exodus = false
  exodus = true
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
    # electrode_area = 1.1
    # ballast_resist = 1.1
    # e = 1.1
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
    block = 0
    position_units = ${dom0Scale}
  [../]
  #[./em_ionization]
  #  type = ElectronsFromIonization
  #  variable = em
  #  potential = potential
  #  mean_en = mean_en
  #  em = em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./em_log_stabilization]
    type = LogStabilizationMoles
    variable = em
    block = 0
  [../]
  # [./em_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = em
  #   potential = potential
  #   block = 0
  # [../]

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
#  [./emliq_reactant_first_order_rxn]
#    type = ReactantFirstOrderRxn
#    variable = emliq
#    block = 1
#  [../]
#  [./emliq_water_bi_sink]
#    type = ReactantAARxn
#    variable = emliq
#    block = 1
#  [../]
  [./emliq_log_stabilization]
    type = LogStabilizationMoles
    variable = emliq
    block = 1
  [../]
  # [./emliq_advection_stabilization]
  #   type = EFieldArtDiff
  #   variable = emliq
  #   potential = potential
  #   block = 1
  # [../]

  [./potential_diffusion_dom1]
    type = CoeffDiffusionLin
    variable = potential
    block = 0
    position_units = ${dom0Scale}
  [../]
  #[./potential_diffusion_dom2]
  #  type = CoeffDiffusionLin
  #  variable = potential
  #  block = 1
  #  position_units = ${dom1Scale}
  #[../]
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
  #[./emliq_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential
  #  charged = emliq
  #  block = 1
  #[../]
  #[./OHm_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential
  #  charged = OHm
  #  block = 1
  #[../]

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
  #[./Arp_ionization]
  #  type = IonsFromIonization
  #  variable = Arp
  #  potential = potential
  #  em = em
  #  mean_en = mean_en
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
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
#  [./OHm_product_first_order_rxn]
#    type = ProductFirstOrderRxn
#    variable = OHm
#    v = emliq
#    block = 1
#  [../]
#  [./OHm_product_aabb_rxn]
#    type = ProductAABBRxn
#    variable = OHm
#    v = emliq
#    block = 1
#  [../]

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
  #[./mean_en_ionization]
  #  type = ElectronEnergyLossFromIonization
  #  variable = mean_en
  #  potential = potential
  #  em = em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./mean_en_elastic]
  #  type = ElectronEnergyLossFromElastic
  #  variable = mean_en
  #  potential = potential
  #  em = em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./mean_en_excitation]
  #  type = ElectronEnergyLossFromExcitation
  #  variable = mean_en
  #  potential = potential
  #  em = em
  #  block = 0
  #  position_units = ${dom0Scale}
  #[../]
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


  [./potential_liq_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential_liq
    block = 1
    position_units = ${dom1Scale}
  [../]
  [./emliq_charge_source_02]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = emliq
    block = 1
  [../]
  [./OHm_charge_source_02]
    type = ChargeSourceMoles_KV
    variable = potential_liq
    charged = OHm
    block = 1
  [../]

[]

[Variables]
  [./potential]
    block = 0
    #scaling = 1e6
  [../]
  [./potential_liq]
    block = 1
    #scaling = 1e6
  [../]
  [./em]
    block = 0
    #scaling = 1e3
  [../]
  [./emliq]
    block = 1
    #scaling = 1e5
  [../]

  [./Arp]
    block = 0
    #scaling = 1e4
  [../]

  [./mean_en]
    block = 0
    #scaling = 1e3
  [../]

  [./OHm]
    block = 1
    #initial_condition = -31
    #scaling = 2
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
  [./Efield]
    order = CONSTANT
    family = MONOMIAL
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
  [./ProcRate_el]
    type = ProcRate
    em = em
    potential = potential
    proc = el
    variable = ProcRate_el
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ProcRate_ex]
    type = ProcRate
    em = em
    potential = potential
    proc = ex
    variable = ProcRate_ex
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./ProcRate_iz]
    type = ProcRate
    em = em
    potential = potential
    proc = iz
    variable = ProcRate_iz
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
    potential = potential_liq
    variable = Efield
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

  #[./water_interface]
  #  type = InterfaceFluxConservation
  #  neighbor_var = potential
  #  variable = potential_liq
  #  boundary = 'master1_interface'
  #  region_name = 'diffpotential_liq'
  #  neighbor_region_name = 'diffpotential'
  #  position_units = ${dom1Scale}
  #  neighbor_position_units = ${dom0Scale}
  #[../]

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
  #[./potential_liq_interface]
  #  type = MatchedValueBC
  #  variable = potential
  #  v = potential_liq
  #  boundary = master0_interface
  #[../]
  [./potential_liq_interface]
    type = MatchedValueBC
    variable = potential_liq
    v = potential
    boundary = master1_interface
  [../]
  [./potential_liq_right]
    type = DirichletBC
    variable = potential_liq
    boundary = right
    value = 0
  [../]
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
#  [./potential_dirichlet_right]
#    type = DirichletBC
#    variable = potential
#    boundary = right
#    value = 0
#  [../]
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
  #[./em_physical_right]
  #  type = MatchedValueLogBC
  #  variable = em
  #  boundary = 'master0_interface'
  #  v = emliq
  #  H = 1e3
  #[../]
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
    potential = potential_liq
    position_units = ${dom1Scale}
  [../]
  [./OHm_physical]
    type = DCIonBC
    variable = OHm
    boundary = 'right'
    potential = potential_liq
    position_units = ${dom1Scale}
  [../]
  
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
  [./OHm_ic]
    type = ConstantIC
    variable = OHm
    value = -15.6
    block = 1
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
[]

[Functions]
  [./potential_bc_func]
    type = ParsedFunction
    #value = '1.25*tanh(1e6*t)'
    value = 1.25
    #value = '1.25*tanh(1e9*t)'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-1.25 * (1.0001e-3 - x)'
  [../]
[]

[Materials]
#  [./gas_block]
#    type = Gas
#    interp_trans_coeffs = true
#    interp_elastic_coeff = true
#    ramp_trans_coeffs = false
#    em = em
#    potential = potential
#    ip = Arp
#    mean_en = mean_en
#    user_se_coeff = .05
#    block = 0
#    property_tables_file = td_argon_mean_en.txt
#  [../]
 # [./water_block]
 #   type = Water
 #   block = 1
 #   potential = potential
 # [../]

 [./GasBasics]
   type = GasBase
   interp_elastic_coeff = true
   interp_trans_coeffs = true
   ramp_trans_coeffs = false
   # user_p_gas = 1.01325e5
   user_p_gas = 1.01e5
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
  [./Ar_species]
    type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 0
  [../]

 # [./water_block]
 #   type = Water
 #   block = 1
 #   potential = potential
 # [../]


  [./OHm_mat]
    type = HeavySpeciesMod
    heavy_species_name = OHm 
    heavy_species_mass = 2.82420e-26
    heavy_species_charge = -1
    diffusivity = 5.26e-9
    block = 1
  [../]

  [./O_mat]
    type = HeavySpeciesMod
    heavy_species_name = O
    heavy_species_mass = 2.6566962e-26
    heavy_species_charge = 0
    diffusivity = 2.00e-9
    block = 1
  [../]

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
    species = 'emliq OHm'
    use_log = true
    position_units = ${dom1Scale}
    block = 1
    reaction_coefficient_format = 'rate'
    reactions = 'emliq -> OHm : 1069.6
                 emliq + emliq -> OHm + OHm : 3.136e8'
  [../]
[]
