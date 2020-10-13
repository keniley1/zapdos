[Mesh]
  #type = GeneratedMesh
  #dim = 2
  #nx = 20
  #ny = 20
  [file]
    type = FileMeshGenerator
    file = 'I_shape.msh'
  []
  #construct_side_list_from_node_list = true
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
  [v]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [u_dt]
    type = TimeDerivative
    variable = u
  []
  [diff_u]
    type = MatDiffusion
    variable = u
    diffusivity = D_u
  []

  #[v_dt]
  #  type = TimeDerivative
  #  variable = v
  #[]
  [diff_v]
    type = MatDiffusion
    variable = v
    diffusivity = D_v
  []

  [force_u]
    type = CoupledForce
    variable = u
    v = v
  []
  
  [force_v]
    type = BodyForce
    variable = v
    value = 1
  []
[]

[BCs]
  [bottom_u]
    type = DirichletBC
    variable = u
    boundary = 'bottom bottom_left bottom_right top top_left top_right right_electrode left_electrode'
    #boundary = 'bottom_left'
    value = 0
  []

  [bottom_v]
    type = DirichletBC
    variable = v
    #boundary = 'bottom bottom_left bottom_right top top_left top_right right_electrode'
    boundary = 'bottom bottom_left bottom_right top_left top_right right_electrode top'
    value = 0
  []
  [left_electrode]
    type = FunctionDirichletBC
    function = 'varying_value'
    variable = v
    boundary = 'left_electrode'
    #boundary = 'top'
  []
[]

[Materials]
  [diffusivity_u]
    type = GenericConstantMaterial
    prop_names = D_u
    prop_values = 0.1
  []

  [diffusivity_v]
    type = GenericConstantMaterial
    prop_names = D_v
    prop_values = 10
  []
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

[Functions]
  [gaussian]
    type = ParsedFunction
    vars = 'A    x0   y0   sigma_x   sigma_y'
    vals = '1.0  0    0    0.5       0.5'
    value = 'A * exp(-((x-x0)^2/(2.0*sigma_x^2) + (y-y0)^2/(2.0*sigma_y^2)))'
  []

  [varying_value]
    type = ParsedFunction
    vars = 'w'
    vals = '100'
    value = '10*sin(w*t)'
  []
[]

[Executioner]
  #type = Steady
  type = Transient
  end_time = 0.10
  solve_type = 'newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 0.005
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
