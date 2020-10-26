# This simulation is modeled after the COMSOL plasma moduleexample problem,
# "Dielectric Barrier Discharge".
#
# The domain is 1D and consists of 3 blocks, each 0.1 mm thick.
# The plasma discharge region is set from  0.1 mm <= x <= 0.2 mm
# and is bounded on either side by a 0.1 mm thick dielectric region.
#
# Only argon ions, electrons, and a single excited state are included for simplicity.
#
# Surface charge is allowed to accumulate on the dielectric surfaces based on the
# total charged particle flux (\sum_i q_i \Gamma_i - \Gamma_e).
# Two different secondary electron emission coefficients are supplied to the
# left and right boundaries through the GenericConstantMaterial interface.
#
# Potential boundary conditions are included as both InterfaceKernels and normal BCs.
#
# In this case the mesh is not scaled at all. 
dom0Scale=1
dom1Scale=1
dom2Scale=1

[GlobalParams]
  potential_units = kV
  use_moles = true
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'argon_dbd_mesh.msh'
  []
  [dielectric_left]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'b0_right'
    input = file
  []
  [plasma_left]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'b1_left'
    input = dielectric_left
  []
  [plasma_right]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '2'
    new_boundary = 'b1_right'
    input = plasma_left
  []
  [dielectric_right]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '2'
    paired_block = '1'
    new_boundary = 'b2_left'
    input = plasma_right
  []
  [left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = dielectric_right
  []
  [right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  []
[]

[Problem]
  type = FEProblem
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = 'lu NONZERO 1.e-10'
  []
[]

[Outputs]
  # Surface charge is provided from a Material property.
  # Here the option output_material_properties is set to
  # true and the surface_charge material property is named.
  # Now "surface_charge" will appear as an AuxVariable in
  # the exodus file.
  perf_graph = true
  [out]
    type = Exodus
    output_material_properties = true
    show_material_properties = 'surface_charge'
  []
[]

[Materials]
  #########
  # Define some necessary constants.
  # Energy of secondary electrons is set to 1 eV here.
  # e:          fundamental charge
  # N_A:        Avogadro's number
  # k_boltz:    Boltzmann constant
  # eps:        permittivity of free space
  # se_energy:  Secondary electron coefficient - emitted energy (set to 1 eV here)
  # T_gas:      Gas temperature (Kelvin)
  # massem:     9.11e-31 (electron mass)
  # p_gas:      Gas pressure (Pascals)
  #########
  [gas_constants]
    type = GenericConstantMaterial
    block = 1
    prop_names = ' e       N_A      k_boltz  eps         se_energy T_gas  massem   p_gas'
    prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12   2.5       400    9.11e-31 1.01e5'
  []
[]
