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
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
    input = geo
  [../]
  [./interface2]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = '1'
    paired_block = '0'
    new_boundary = 'master1_interface'
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

  [./fsp]
    type = FSP
    topsplit = 'pc'
    #full = 'true'
    [./pc]
      splitting = 'p c'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    [../]
    [./p]
      vars = 'potential'
      petsc_options_iname = '-pc_type -pc_hypre_type'
      petsc_options_value = 'hypre boomeramg'
    [../]
    [./c]
      vars = 'em Arp mean_en emliq OH- O- O2- O3- HO2- H+ O O2_1 O3 H H2 HO2 HO3 OH H2O2'
      petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
      petsc_options_value = 'asm 2 ilu NONZERO 1e-10'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  #end_time = 1e-1
  #end_time = 1e10
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  #resid_vs_jac_scaling_param = 1
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  #solve_type = pjfnk
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-12
  dtmin = 1e-18
  l_max_its = 20
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    growth_factor = 1.4
    optimal_iterations = 10
  [../]
[]

[Outputs]
  # perf_graph = true
  #print_densityear_residuals = false
  #[./m5_kV_1um_02]
  #[./m4_kV_1um_02]
  #[./m3_kV_1um_02]
  #[./m2_kV_1um_02]
  #[./m1_kV_1um_02]
  [m1_kV_1um_02]
    type = Exodus
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Variables]
  [H2O]
    block = 0
    #initial_condition = 0.0367321
  []
  [emliq]
    block = 1
  []
  [potential]
    block = 0
    initial_condition = 0
  []
[]

[Kernels]
  [H2O_dt]
    type = ADTimeDerivativeLog
    variable = H2O
    block = 0
  []
  [H2O_diff]
    type = ADCoeffDiffusion
    variable = H2O
    block = 0
    position_units = ${dom0Scale}
  []


  [emliq_dt]
    type = ADTimeDerivativeLog
    variable = emliq
    block = 1
  []

  [potential_dt]
    type = ADTimeDerivativeLog
    variable = potential
    block = 0
  []
[]

[BCs]
  [H2O_right]
    type = DirichletBC
    variable = H2O
    value = 0.38
    boundary = 'master0_interface'
  []
  [H2O_left]
    type = ADHagelaarIonDiffusionBC
    variable = H2O
    r = 0
    position_units = ${dom0Scale}
    boundary = 'master0_interface'
  []
  #[H2O_left]
  #  type = DirichletBC
  #  variable = H2O
  #  value = -3
  #  boundary = 'left'
  #[]
  #[H2O_left]
  #  type = ADDriftDiffusionOpenBC
  #  variable = H2O
  #  potential = potential
  #  position_units = ${dom0Scale}
  #  boundary = 'left'
  #[]
[]

[ICs]
  [H2O_ic]
    type = FunctionIC
    variable = H2O
    function = water_ic_func
  []
[]

[Functions]
  [./water_ic_func]
    type = ParsedFunction
    value = 'log(8.6949e23/6.022e23)'
    #value = '-2.15'
  [../]
[]

[Materials]
  [./se_coefficient]
    type = GenericConstantMaterial
    prop_names = 'se_coeff'
    prop_values = '0.01'
    boundary = 'left master0_interface'
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
[]
