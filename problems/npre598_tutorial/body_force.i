[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  #[dt]
  #  type = TimeDerivative
  #  variable = u
  #[]
  [diff]
    type = Diffusion
    variable = u
  []
  [force]
    type = BodyForce
    variable = u
    value = 1
    #function = gaussian
  []
[]

[BCs]
  [bottom]
    type = DirichletBC
    variable = u
    boundary = 'bottom'
    value = 0
  []

  [top]
    type = DirichletBC
    variable = u
    boundary = 'top'
    value = 0
  []

  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 0
  []

  [right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 0
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
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  compute_scaling_once = false
  end_time = 10
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = NEWTON
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -ksp_type -snes_linesearch_minlambda'
  #petsc_options_value = 'lu NONZERO 1.e-10 preonly 1e-3'
  nl_rel_tol = 1e-4
  nl_abs_tol = 7.6e-5
  dtmin = 1e-12
  steady_state_detection = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-2
    # dt = 1.1
    growth_factor = 1.2
   optimal_iterations = 15
  [../]
  #type = Steady
  #solve_type = 'newton'
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
