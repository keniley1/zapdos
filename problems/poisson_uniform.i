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
  [diff]
    type = Diffusion
    variable = u
  []
  [force]
    type = BodyForce
    variable = u
    value = 1
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
  type = Steady
  solve_type = 'newton'
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
