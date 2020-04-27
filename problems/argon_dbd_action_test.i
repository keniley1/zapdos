dom0Scale=1
dom1Scale=1
dom2Scale=1

[GlobalParams]
  #offset = 20
  #offset = 30
  offset = 40
  # offset = 0
  potential_units = kV
  use_moles = true
  # potential_units = V
[]

[Mesh]
  #final_generator = dielectric_left_mesh
  [./file]
    type = FileMeshGenerator
    #file = 'argon_csv_bc.msh'
    file = 'argon_dbd_mesh.msh'
  [../]
  [./dielectric_left]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master01_interface'
    input = file
  [../]
  [./plasma_left]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master10_interface'
    input = dielectric_left
  [../]
  [./plasma_right]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '2'
    new_boundary = 'master12_interface'
    input = plasma_left
  [../]
  [./dielectric_right]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '2'
    paired_block = '1'
    new_boundary = 'master21_interface'
    input = plasma_right
  [../]
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = dielectric_right
  [../]
  [./right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  [../]

  #[./dielectric_left_mesh]
  #  type = MeshSideSetGenerator
  #  input = right
  #  boundaries = 'master10_interface'
  #  block_id = 101
  #  block_name = dielectric_left_surface
  #  depends_on = 'plasma_left'
  #[../]
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
  #automatic_scaling = true
  #compute_scaling_once = false
  #end_time = 5e-3
  end_time = 10e-5
  #end_time = 48e-5
  #num_steps = 100
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  #solve_type = NEWTON
  solve_type = PJFNK
  line_search = 'basic'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol -ksp_gmres_restart -pc_factor_mat_ordering_type'
  #petsc_options_value = 'lu NONZERO 1.e-10 0 40 nd'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol -ksp_gmres_restart -snes_linesearch_minlambda'
  petsc_options_value = 'lu NONZERO 1.e-10 0 100 1e-3'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  dtmin = 1e-18
  dtmax = 5e-7
  l_max_its = 100
  nl_max_its = 25
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-8
    growth_factor = 1.2
   optimal_iterations = 30
  [../]
  #[./TimeIntegrator]
  #  type = BDF2
  #[../]
[]

[Outputs]
  perf_graph = true
  [./out]
    type = Exodus
    #execute_on = 'final'
  [../]
  #[./sigma_output]
  #  type = CSV
  #  show = 'sigma_left'
  #[../]
  #[./dof_map]
  #  type = DOFMap
  #[../]
[]

[Debug]
  show_var_residual_norms = true
  #show_actions = true
[]

#[ScalarKernels]
#  [./sigma_left_dt]
#    type = ODETimeDerivative
#    variable = sigma_left
#  [../]
#  [./sigma_left_rhs]
#    type = BoundaryFlux
#    current = left_current
#    variable = sigma_left
#  [../]
#
#  [./sigma_right_dt]
#    type = ODETimeDerivative
#    variable = sigma_right
#  [../]
#  [./sigma_right_rhs]
#    type = BoundaryFlux
#    current = right_current
#    variable = sigma_right
#  [../]
#[]

[Variables]
  #[./sigma_left]
  #  order = FIRST
  #  family = SCALAR
  #  initial_condition = 0
  #[../]
  #[./sigma_right]
  #  order = FIRST
  #  family = SCALAR
  #  initial_condition = 0
  #[../]

  #[./Ar]
  #  initial_condition = 3.70109
  #  block = 1
  #[../]
  [./potential_dom0]
    block = 0
  [../]
  [./potential_dom1]
    block = 1
  [../]
  [./potential_dom2]
    block = 2
  [../]

[]

[DriftDiffusionAction]
  [./Plasma]
    #charged_particle = 'Arp Ar2p'
    charged_particle = 'Arp'
    electrons = 'em'
    Neutrals = 'Ar*'
    potential = potential_dom1
    mean_energy = mean_en
    Is_potential_unique = false
    using_offset = true
    offset = 40
    position_units = ${dom1Scale}
    #Additional_Outputs = 'ElectronTemperature EField'
    Additional_Outputs = 'ElectronTemperature'
    block = 1
  [../]
[]

[Kernels]
  #[./em_time_deriv]
  #  #type = ADTimeDerivativeLog
  #  type = ElectronTimeDerivative
  #  variable = em
  #  block = 1
  #[../]
  #[./em_advection]
  #  #type = ADEFieldAdvection
  #  type = EFieldAdvectionElectrons
  #  variable = em
  #  mean_en = mean_en
  #  potential = potential_dom1
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./em_diffusion]
  #  #type = ADCoeffDiffusion
  #  type = CoeffDiffusionElectrons
  #  variable = em
  #  mean_en = mean_en
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./em_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = em
  #  #offset = 40
  #  block = 1
  #[../]

  [./potential_diffusion_dom0]
    type = CoeffDiffusionLin
    variable = potential_dom0
    block = 0
    position_units = ${dom0Scale}
  [../]
  #[./potential_dom1_diffusion1_block]
  #  type = CoeffDiffusionLin
  #  variable = potential_dom1
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential_dom2
    block = 2
    position_units = ${dom0Scale}
  [../]

  #[./Arp_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential_dom1
  #  charged = Arp
  #  block = 1
  #[../]
  #[./em_charge_source]
  #  type = ChargeSourceMoles_KV
  #  variable = potential_dom1
  #  charged = em
  #  block = 1
  #[../]

  #[./Arp_time_deriv]
  #  #type = ADTimeDerivativeLog
  #  type = ElectronTimeDerivative
  #  variable = Arp
  #  block = 1
  #[../]
  #[./Arp_advection]
  #  #type = ADEFieldAdvection
  #  type = EFieldAdvection
  #  variable = Arp
  #  potential = potential_dom1
  #  position_units = ${dom0Scale}
  #[../]
  #[./Arp_diffusion]
  #  #type = ADCoeffDiffusion
  #  type = CoeffDiffusion
  #  variable = Arp
  #  position_units = ${dom0Scale}
  #[../]
  #[./Arp_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Arp
  #  block = 1
  #[../]

  #[./Arex_time_deriv]
  #  #type = ADTimeDerivativeLog
  #  type = ElectronTimeDerivative
  #  variable = Ar*
  #  block = 1
  #[../]
  #[./Arex_diffusion]
  #  #type = ADCoeffDiffusion
  #  type = CoeffDiffusion
  #  variable = Ar*
  #  position_units = ${dom0Scale}
  #[../]
  #[./Arex_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = Ar*
  #  block = 1
  #[../]

  #[./mean_en_time_deriv]
  #  type = ElectronTimeDerivative
  #  variable = mean_en
  #  block = 1
  #[../]
  #[./mean_en_advection]
  #  type = EFieldAdvectionEnergy
  #  variable = mean_en
  #  potential = potential_dom1
  #  em = em
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./mean_en_diffusion]
  #  type = CoeffDiffusionEnergy
  #  variable = mean_en
  #  em = em
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./mean_en_joule_heating]
  #  type = JouleHeating
  #  variable = mean_en
  #  potential = potential_dom1
  #  em = em
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./mean_en_log_stabilization]
  #  type = LogStabilizationMoles
  #  variable = mean_en
  #  block = 1
  #  #offset = 25
  #  #offset = 15
  #  #offset = 40
  #[../]
[]

[InterfaceKernels] 
  [./potential_left]
    type = InterfaceDiffusion
    neighbor_var = potential_dom0
    variable = potential_dom1
    sigma = sigma_left
    sigma_test = sigma_test
    boundary = master10_interface
    position_units = ${dom0Scale}
    neighbor_position_units = ${dom0Scale}
  [../]

  [./potential_right]
    type = InterfaceDiffusion
    neighbor_var = potential_dom2
    variable = potential_dom1
    sigma = sigma_right
    sigma_test = sigma_test
    boundary = master12_interface
    position_units = ${dom0Scale}
    neighbor_position_units = ${dom0Scale}
  [../]

  #[./potential_left2]
  #  type = PenaltyInterfaceDiffusion
  #  penalty = 1e6
  #  variable = potential_dom0
  #  neighbor_var = potential_dom1
  #  boundary = 'master01_interface'
  #[../]
  #[./potential_right2]
  #  type = PenaltyInterfaceDiffusion
  #  penalty = 1e6
  #  variable = potential_dom2
  #  neighbor_var = potential_dom1
  #  boundary = 'master21_interface'
  #[../]
[]

[AuxVariables]
  [./sigma_test]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./Ar_density]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ar]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.70109
  [../]
  [./x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./rho]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  [./Efield]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./Current_em]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  #[./Current_Arp]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
  #[./tot_gas_current]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 1
  #[../]
[]

[AuxKernels]
  [./sigma_test]
    type = HagelaarChargedFlux
    variable = sigma_test
    potential = potential_dom1
    r_electron = 0
    r_ion = 0
    mean_en = mean_en
    em = em
    ions = 'Arp'
    se_coeff = 0.005
    position_units = ${dom0Scale}
    execute_on = 'LINEAR NONLINEAR'
    boundary = 'master10_interface master12_interface'
  [../]
  [./Ar_lin]
    type = DensityMoles
    variable = Ar_density
    density_log = Ar
    block = 1
  [../]
  #[./e_temp]
  #  type = ElectronTemperature
  #  variable = e_temp
  #  electron_density = em
  #  mean_en = mean_en
  #  block = 1
  #[../]
  [./x_g0]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./x_g1]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./x_g2]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 2
  [../]
  #[./rho]
  #  type = ParsedAux
  #  variable = rho
  #  #args = 'em_lin Arp_lin Ar2p_lin'
  #  #function = 'Arp_lin + Ar2p_lin - em_lin'
  #  args = 'em_lin Arp_lin'
  #  function = 'Arp_lin - em_lin'
  #  execute_on = 'timestep_end'
  #  block = 1
  #[../]
  #[./tot_gas_current]
  #  type = ParsedAux
  #  variable = tot_gas_current
  #  #args = 'Current_em Current_Arp Current_Ar2p'
  #  #function = 'Current_em + Current_Arp + Current_Ar2p'
  #  args = 'Current_em Current_Arp'
  #  function = 'Current_em + Current_Arp'
  #  execute_on = 'timestep_end'
  #  block = 1
  #[../]
  [./Efield_d0]
    type = Efield
    component = 0
    potential = potential_dom0
    #potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./Efield_d1]
    type = Efield
    component = 0
    potential = potential_dom1
    variable = Efield
    position_units = ${dom0Scale}
    block = 1
  [../]
  [./Efield_d2]
    type = Efield
    component = 0
    potential = potential_dom2
    #potential = potential
    variable = Efield
    position_units = ${dom0Scale}
    block = 2
  [../]
  #[./Current_em]
  #  type = Current
  #  potential = potential
  #  density_log = em
  #  variable = Current_em
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
  #[./Current_Arp]
  #  type = Current
  #  potential = potential
  #  density_log = Arp
  #  variable = Current_Arp
  #  art_diff = false
  #  block = 1
  #  position_units = ${dom0Scale}
  #[../]
[]

[BCs]
  [./match_phi_left]
    type = MatchedValueBC
    variable = potential_dom0
    v = potential_dom1
    boundary = 'master01_interface'
  [../]
  [./match_phi_right]
    type = MatchedValueBC
    variable = potential_dom2
    v = potential_dom1
    boundary = 'master21_interface'
  [../]

  [./mean_en_physical_right]
    #type = ADHagelaarEnergyBC
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master12_interface'
    potential = potential_dom1
    em = em
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./mean_en_physical_left]
    #type = ADHagelaarEnergyBC
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'master10_interface'
    potential = potential_dom1
    em = em
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./secondary_energy_left]
    #type = ADSecondaryElectronEnergyBC
    type = SecondaryElectronEnergyBC
    variable = mean_en
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential_dom1
    em = em
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    se_coeff = 0.01
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./secondary_energy_right]
    #type = ADSecondaryElectronEnergyBC
    type = SecondaryElectronEnergyBC
    variable = mean_en
    #boundary = 'left'
    boundary = 'master12_interface'
    potential = potential_dom1
    em = em
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    se_coeff = 0.01
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./em_physical_right]
    #type = ADHagelaarElectronBC
    type = HagelaarElectronBC
    variable = em
    #boundary = 'right'
    boundary = 'master12_interface'
    potential = potential_dom1
    mean_en = mean_en
    #r = 0.99
    r = 0.0
    position_units = ${dom0Scale}
  [../]
  [./em_physical_left]
    #type = ADHagelaarElectronBC
    type = HagelaarElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential_dom1
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_left]
    #type = ADSecondaryElectronBC
    type = SecondaryElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master10_interface'
    potential = potential_dom1
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    se_coeff = 0.01
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./sec_electrons_right]
    #type = ADSecondaryElectronBC
    type = SecondaryElectronBC
    variable = em
    #boundary = 'left'
    boundary = 'master12_interface'
    potential = potential_dom1
    #ip = 'Arp Ar2p'
    ip = 'Arp'
    se_coeff = 0.01
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]

  [./potential_left]
    type = FunctionDirichletBC
    variable = potential_dom0
    #variable = potential
    function = potential_input
    boundary = 'left'
  [../]

  [./potential_dirichlet_right]
    type = DirichletBC
    variable = potential_dom2
    #variable = potential
    boundary = right
    value = 0
  [../]

  [./Arp_advection_left]
    #type = ADHagelaarIonAdvectionBC
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = master10_interface
    potential = potential_dom1
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_advection_right]
    #type = ADHagelaarIonAdvectionBC
    type = HagelaarIonAdvectionBC
    variable = Arp
    boundary = master12_interface
    potential = potential_dom1
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion_left]
    #type = ADHagelaarIonDiffusionBC
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = master10_interface
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_diffusion_right]
    #type = ADHagelaarIonDiffusionBC
    type = HagelaarIonDiffusionBC
    variable = Arp
    boundary = master12_interface
    r = 0
    position_units = ${dom0Scale}
  [../]

  #[./Ar2p_advection_left]
  #  #type = ADHagelaarIonAdvectionBC
  #  type = HagelaarIonAdvectionBC
  #  variable = Ar2p
  #  boundary = master10_interface
  #  potential = potential_dom1
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_advection_right]
  #  #type = ADHagelaarIonAdvectionBC
  #  type = HagelaarIonAdvectionBC
  #  variable = Ar2p
  #  boundary = master12_interface
  #  potential = potential_dom1
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_diffusion_left]
  #  #type = ADHagelaarIonDiffusionBC
  #  type = HagelaarIonDiffusionBC
  #  variable = Ar2p
  #  boundary = master10_interface
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  #[./Ar2p_diffusion_right]
  #  #type = ADHagelaarIonDiffusionBC
  #  type = HagelaarIonDiffusionBC
  #  variable = Ar2p
  #  boundary = master12_interface
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]

  [./Arex_diffusion_left]
    #type = ADHagelaarIonDiffusionBC
    type = HagelaarIonDiffusionBC
    variable = Ar*
    boundary = master10_interface
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arex_diffusion_right]
    #type = ADHagelaarIonDiffusionBC
    type = HagelaarIonDiffusionBC
    variable = Ar*
    boundary = master12_interface 
    r = 0
    position_units = ${dom0Scale}
  [../]
[]

[ICs]
  [./potential_dom0_ic]
    type = FunctionIC
    variable = potential_dom0
    function = potential_ic_func
  [../]
  [./potential_dom1_ic]
    type = FunctionIC
    variable = potential_dom1
    function = potential_ic_func
  [../]
  [./potential_dom2_ic]
    type = FunctionIC
    variable = potential_dom2
    function = potential_ic_func
  [../]

  [./em_ic]
    type = ConstantIC
    variable = em
    value = -20
    block = 1
  [../]
  [./Arp_ic]
    type = ConstantIC
    variable = Arp
    value = -30
    block = 1
  [../]
  #[./Ar2p_ic]
  #  type = ConstantIC
  #  variable = Ar2p
  #  value = -30
  #  block = 1
  #[../]
  [./Arex_ic]
    type = ConstantIC
    variable = Ar*
    value = -20
    block = 1
  [../]
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -20
    block = 1
  [../]
[]

[Functions]
  #[./potential_input]
  #  type = PiecewiseLinear
  #  #data_file = 'voltage_data_01.txt'
  #  data_file = 'voltage_data_temp01.txt'
  #  format = columns
  #[../]
  [./potential_input]
    type = ParsedFunction
    vars = 'f0'
    vals = '50e3'
    value = '-0.75*sin(2*3.1415926*f0*t)'
  [../]
  [./potential_bc_func]
    type = ParsedFunction
    # value = '1.25*tanh(1e6*t)'
    value = 0.8
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '-0.1 * (1.0001e-4 - x)'
  [../]
[]

[Materials]
  #[./electron_moments_ad]
  #  type = ADGasElectronMoments
  #  block = 1
  #  em = em
  #  mean_en = mean_en
  #  property_tables_file = 'argon_cm_test/electron_mobility_diffusion.txt'
  #[../]
  #[./gas_constants]
  #  type = GenericConstantMaterial
  #  block = 1
  #  prop_names = ' e         N_A      diffpotential k_boltz eps  se_coeff se_energy T_gas massem   p_gas  n_gas'
  #  prop_values = '1.6e-19 6.022e23 8.85e-12      1.38e-23 8.854e-12 0.05     3.        300   9.11e-31 1.01e5 40.4915'
  #[../]
  [./Test]
    type = HagelaarSurfaceCharge
    potential = potential_dom1
    r_electron = 0
    r_ion = 0
    mean_en = mean_en
    em = em
    ions = 'Arp'
    se_coeff = 0.005
    position_units = ${dom0Scale}
    boundary = 'master10_interface master12_interface'
  [../]


  [./electron_moments]
    type = GasElectronMoments
    block = 1
    em = em
    mean_en = mean_en
    user_se_coeff = 0.005
    potential = potential_dom1
    user_p_gas = 101325
    property_tables_file = 'argon_cm_test/electron_mobility_diffusion.txt'
    interp_trans_coeffs = true
    interp_elastic_coeff = false
    ramp_trans_coeffs = false
    position_units = ${dom0Scale}
  [../]


  [./dielectric_left_side]
    type = GenericConstantMaterial
    block = 0
    #prop_names = 'diffpotential'
    prop_names = 'diffpotential_dom0'
    prop_values = '8.85e-11'
  [../]
  [./gas_phase]
    type = GenericConstantMaterial
    block = 1
    #prop_names = 'diffpotential'
    prop_names = 'diffpotential_dom1'
    prop_values = '8.85e-12'
  [../]
  [./dielectric_right_side]
    type = GenericConstantMaterial
    block = 2
    #prop_names = 'diffpotential'
    prop_names = 'diffpotential_dom2'
    prop_values = '8.85e-11'
  [../]

  [./gas_species_0]
    #type = ADHeavySpeciesMaterial
    type = HeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    block = 1
  [../]

  [./gas_species_10]
    #type = ADHeavySpeciesMaterial
    type = HeavySpeciesMaterial
    heavy_species_name = Ar2p
    heavy_species_mass = 13.28e-26
    heavy_species_charge = 1.0
    block = 1
  [../]

  [./gas_species_2]
    #type = ADHeavySpeciesMaterial
    type = HeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]

  [./gas_species_3]
    #type = ADHeavySpeciesMaterial
    type = HeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    block = 1
  [../]
[]

#[Postprocessors]
#  #[./left_current]
#  #  type = SideCurrent
#  #  variable = em
#  #  mobility = muem
#  #  potential = potential_dom1
#  #  r = 0
#  #  mean_en = mean_en
#  #  #ions = 'Arp Ar2p'
#  #  ions = 'Arp'
#  #  Arp = Arp
#  #  #se_coeff = 0.01
#  #  se_coeff = 0.0
#  #  position_units = ${dom1Scale}
#  #  execute_on = 'LINEAR NONLINEAR'
#  #  boundary = 'master10_interface'
#  #[../]
#  #
#  #[./right_current]
#  #  type = SideCurrent
#  #  variable = em
#  #  mobility = muem
#  #  potential = potential_dom1
#  #  r = 0
#  #  mean_en = mean_en
#  #  #ions = 'Arp Ar2p'
#  #  ions = 'Arp'
#  #  Arp = Arp
#  #  se_coeff = 0.0
#  #  position_units = ${dom1Scale}
#  #  execute_on = 'LINEAR NONLINEAR'
#  #  boundary = 'master12_interface'
#  #[../]
#
#  [./left_current]
#    type = SurfaceCharge
#    variable = em
#    potential = potential_dom1
#    r_electron = 0
#    r_ion = 0
#    mean_en = mean_en
#    #ions = 'Arp Ar2p'
#    ions = 'Arp'
#    se_coeff = 0.005
#    position_units = ${dom0Scale}
#    execute_on = 'LINEAR NONLINEAR'
#    boundary = 'master10_interface'
#  [../]
#
#  [./right_current]
#    type = SurfaceCharge
#    variable = em
#    potential = potential_dom1
#    r_electron = 0
#    r_ion = 0
#    mean_en = mean_en
#    #ions = 'Arp Ar2p'
#    ions = 'Arp'
#    se_coeff = 0.005
#    position_units = ${dom0Scale}
#    execute_on = 'LINEAR NONLINEAR'
#    boundary = 'master12_interface'
#  [../]
#[]

[Reactions]
  # This argon reaction network based on a ZDPlasKin example:
  # zdplaskin.laplace.univ-tlse.fr
  # Example: "Micro-cathode sustained discharged in Ar"
  [./Argon]
    #species = 'em Arp Ar2p Ar*'
    species = 'em Arp Ar*'
    aux_species = 'Ar'
    #reaction_coefficient_format = 'townsend'
    reaction_coefficient_format = 'rate'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    file_location = 'argon_cm_test/'
    #equation_constants = 'Tgas'
    #equation_values = '300'
    #equation_variables = 'e_temp'
    potential = 'potential_dom1'
    use_log = true
    position_units = ${dom1Scale}
    #use_ad = true
    #convert_to_moles = true
    #convert_to_meters = 1e-2
    block = 1

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
                 em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
                 em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ars_De-excitation_11.50_eV)
                 em + Ar* -> em + em + Arp        : EEDF [-4.43] (C5_Ars_Ionization_4.43_eV)
                 Ar* + Ar -> Ar + Ar              : 1807
                 Ar* + Ar* -> em + Ar + Arp       : 3.3734e8'

    #reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Elastic)
    #             em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
    #             em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
    #             em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ars_De-excitation_11.50_eV)
    #             em + Ar* -> em + em + Arp        : EEDF [-4.43] (C5_Ars_Ionization_4.43_eV)
    #             Arp + em + em -> Ar + em         : {3.17314235e9 * (e_temp/11600)^(-4.5)}
    #             Ar* + Ar + Ar -> Ar + Ar + Ar    : 5077.02776
    #             Ar2p + em -> Ar* + Ar            : {5.1187e11 * (e_temp/300)^(-0.67)}
    #             Ar2p + Ar -> Arp + Ar + Ar       : {3.649332e12 / Tgas * exp(-15130/Tgas)}
    #             Ar* + Ar* -> Ar2p + em           : 3.6132e8
    #             Ar* + Ar -> Ar + Ar              : 1807
    #             Arp + Ar + Ar -> Ar2p + Ar       : {81595.089 * (Tgas/300)^(-0.4)}'
  [../]

[]

