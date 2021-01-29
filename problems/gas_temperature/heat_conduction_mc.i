[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 100.0 # Length of test chamber
  ymax =  10.0 # Test chamber radius
[]

[MeshModifiers]
  [block_A]
    type = BoundingBoxNodeSet
    new_boundary = block_A
    top_right = '100 10 0'
    bottom_left = '0 5 0'
  [../]
  [block_B]
    type = BoundingBoxNodeSet
    new_boundary = block_B
    top_right = '100 5 0'
    bottom_left = '0 0 0'
  [../]
[]

[Variables]
  [./temperature]
  [../]
[]

[Kernels]
  [./heat_conduction]
    type = HeatConduction
    variable = temperature
  [../]
[]

[BCs]
  [./inlet_temperature]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 200.0 # (K)
  [../]
  [./outlet_temperature]
    type = NeumannBC
    variable = temperature
    boundary = right
    value = 2.0 # (K)
  [../]
[]

[Materials]
  [./block_A]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10 # K: (W/m*K) from wikipedia @296K
    block = block_A
  [../]
  [./block_B]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 5 # K: (W/m*K) from wikipedia @296K
    block = block_B
  [../]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = X
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
