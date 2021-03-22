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
# In this case the mesh is scaled such that the domain goes
# from 0 to 3.
dom0Scale=1
dom1Scale=1
dom2Scale=1

[GlobalParams]
  offset = 60
  potential_units = kV
  use_moles = true
[]

[Mesh]
  [./file]
    type = FileMeshGenerator
    file = 'argon_dbd_mesh.msh'
  [../]
  [./dielectric_left]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'b0_right'
    input = file
  [../]
  [./plasma_left]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'b1_left'
    input = dielectric_left
  [../]
  [./plasma_right]
    # plasma master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '1'
    paired_block = '2'
    new_boundary = 'b1_right'
    input = plasma_left
  [../]
  [./dielectric_right]
    # left dielectric master
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = '2'
    paired_block = '1'
    new_boundary = 'b2_left'
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
[]

[Problem]
  type = FEProblem
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = 'lu NONZERO 1.e-10'
  [../]
[]

[Debug]
  #show_var_residual_norms = true
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  compute_scaling_once = false
  end_time = 10e-5
 # num_steps = 50
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = NEWTON
  #solve_type = PJFNK
  line_search = 'basic'
  nl_rel_tol = 1e-6
  #nl_rel_tol = 1e-4
  l_max_its = 100
  nl_max_its = 25

  # Adaptive timesteppers are useful for plasma simulations since the
  # timestep requirements may vary widely over the course of a simulation.
  # In this case timesteps of 0.1-10 ns are required between pulses as the
  # electron current becomes large, but larger timesteps are allowed otherwise.
  #
  # A maximum timestep of 100 ns is set to make sure the voltage pulse (period
  # of 20 microseconds) is finely sampled.
  dtmin = 1e-14
  dtmax = 5e-7
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-9
    growth_factor = 1.2
   optimal_iterations = 12
  [../]
[]

[Outputs]
  # Surface charge is provided from a Material property.
  # Here the option output_material_properties is set to
  # true and the surface_charge material property is named.
  # Now "surface_charge" will appear as an AuxVariable in
  # the exodus file.
  perf_graph = true
  [./out_01]
    type = Exodus
    output_material_properties = true
    show_material_properties = 'surface_charge'
  [../]
[]

[DriftDiffusionAction]
  [plasma]
    electrons = em
    charged_particle = 'Arp'
    Neutrals = 'Ar*'
    potential = potential_dom1
    mean_energy = mean_en
    using_offset = false
    offset = 40
    use_ad = true
    order = SECOND
    position_units = ${dom0Scale}
    Additional_Outputs = 'ElectronTemperature'
    block = 1
  []
[]

[Variables]
  [potential_dom0]
    block = 0
  []
  [potential_dom1]
    block = 1
  []
  [potential_dom2]
    block = 2
  []
[]

[Kernels]
  [./potential_diffusion_dom0]
    type = CoeffDiffusionLin
    variable = potential_dom0
    block = 0
    position_units = ${dom0Scale}
  [../]
  [./potential_diffusion_dom2]
    type = CoeffDiffusionLin
    variable = potential_dom2
    block = 2
    position_units = ${dom2Scale}
  [../]
[]

[InterfaceKernels]
  # At the dielectric interfaces, the potential is required to be continuous
  # in value but discontinuous in slope due to surface charge accumulation.
  #
  # The potential requires two different boundary conditions on each side:
  #    (1) An InterfaceKernel to provide the Neumann boundary condition
  #    (2) A MatchedValueBC to ensure that
  [./potential_left]
    type = ADPotentialSurfaceCharge
    neighbor_var = potential_dom0
    variable = potential_dom1
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom0Scale}
    boundary = b1_left
  [../]

  [./potential_right]
    type = ADPotentialSurfaceCharge
    neighbor_var = potential_dom2
    variable = potential_dom1
    position_units = ${dom1Scale}
    neighbor_position_units = ${dom2Scale}
    boundary = b1_right
  [../]
[]

[AuxVariables]
  # Argon is considered to be a static background
  # Here it is defined as an AuxVariable
  # Units: log(mol m^-3)
  [./Ar]
    block = 1
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 3.70109
  [../]
  [./Arp_density]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./em_density]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ar*_density]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x_node]
  [../]
  [Efield]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [./em_lin]
    type = DensityMoles
    variable = em_density
    density_log = em
    block = 1
  [../]
  [./Arp_lin]
    type = DensityMoles
    variable = Arp_density
    density_log = Arp
    block = 1
  [../]
  [./Ar*_lin]
    type = DensityMoles
    variable = Ar*_density
    density_log = Ar*
    block = 1
  [../]
  [./x0]
    type = Position
    variable = x
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./x1]
    type = Position
    variable = x
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./x2]
    type = Position
    variable = x
    position_units = ${dom2Scale}
    block = 2
  [../]
  [./xn0]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
    block = 0
  [../]
  [./xn1]
    type = Position
    variable = x_node
    position_units = ${dom1Scale}
    block = 1
  [../]
  [./xn2]
    type = Position
    variable = x_node
    position_units = ${dom2Scale}
    block = 2
  [../]
  [efield_d0]
    type = Efield
    component = 0
    potential = potential_dom0
    variable = Efield
    position_units = ${dom0Scale}
    block = 0
  []
  [efield_d1]
    type = Efield
    component = 0
    potential = potential_dom1
    variable = Efield
    position_units = ${dom0Scale}
    block = 1
  []
  [efield_d2]
    type = Efield
    component = 0
    potential = potential_dom2
    variable = Efield
    position_units = ${dom0Scale}
    block = 2
  []
[]

[BCs]
  ######
  # POTENTIAL BOUNDARY CONDITIONS
  ######
  # Interface BCs:
  [./match_phi_left]
    type = MatchedValueBC
    variable = potential_dom0
    v = potential_dom1
    boundary = 'b0_right'
  [../]
  [./match_phi_right]
    type = MatchedValueBC
    variable = potential_dom2
    v = potential_dom1
    boundary = 'b2_left'
  [../]

  # Electrode and ground BCs:
  # Electrode is at x = 0, named 'left'
  [./potential_left]
    type = FunctionDirichletBC
    variable = potential_dom0
    function = potential_input
    #function = potential_input_test
    boundary = 'left'
  [../]
  # Ground is at x = 0.3 mm, named 'right'
  [./potential_dirichlet_right]
    type = DirichletBC
    variable = potential_dom2
    #variable = potential
    boundary = right
    value = 0
  [../]

  ######
  # HAGELAAR BCS
  ######
  [./em_physical_bc]
    type = ADHagelaarElectronBC
    variable = em
    boundary = 'b1_left b1_right'
    potential = potential_dom1
    mean_en = mean_en
    r = 0.0
    position_units = ${dom1Scale}
  [../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    boundary = 'b1_left b1_right'
    potential = potential_dom1
    ip = 'Arp'
    mean_en = mean_en
    r = 0
    position_units = ${dom1Scale}
  [../]
  [./mean_en_physical]
    type = ADHagelaarEnergyBC
    variable = mean_en
    boundary = 'b1_left b1_right'
    potential = potential_dom1
    em = em
    r = 0
    position_units = ${dom1Scale}
  [../]
  [./secondary_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'b1_left b1_right'
    potential = potential_dom1
    em = em
    ip = 'Arp'
    r = 0
    position_units = ${dom1Scale}
  [../]
  [./Arp_bcs]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'b1_left b1_right'
    potential = potential_dom1
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_bcs2]
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'b1_left b1_right'
    r = 0
    position_units = ${dom1Scale}
  [../]
  [./Arex_bcs]
    type = ADHagelaarIonDiffusionBC
    variable = Ar*
    boundary = 'b1_left b1_right'
    r = 0
    position_units = ${dom1Scale}
  [../]
  #########
  #########
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
    value = -41
    block = 1
  [../]
  [./Arp_ic]
    type = ConstantIC
    variable = Arp
    value = -41
    block = 1
  [../]
  [./Arex_ic]
    type = ConstantIC
    variable = Ar*
    value = -45
    block = 1
  [../]
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -42
    block = 1
  [../]
[]

[Functions]
  # Define a sinusoidal voltage pulse with a frequency of 50 kHz
  # Amplitude is set to 750 Volts
  # (Note that here the amplitude is set to 0.75. Potential units are
  # typically included as kV, not V. This option is set in the
  # GlobalParams block.)
  [./potential_input]
    type = ParsedFunction
    vars = 'f0'
    vals = '50e3'
    value = '-0.75*sin(2*3.1415926*f0*t)'
  [../]

  # Set the initial condition to a line from -10 V on the left and
  # 0 on the right.
  # (Poisson solver tends to struggle with a uniformly zero potential IC.)
  [./potential_ic_func]
    type = ParsedFunction
    value = '-0.01 * (3.0001e-4 - x)'
  [../]
[]

[Materials]
  #########
  # Define secondary electron emission coefficients on the left and right
  # dielectrics.
  #########
  [./se_left]
    type = GenericConstantMaterial
    boundary = 'b1_left'
    prop_names = 'se_coeff'
    prop_values = '0.01'
  [../]
  [./se_right]
    type = GenericConstantMaterial
    boundary = 'b1_right'
    prop_names = 'se_coeff'
    prop_values = '1e-6'
  [../]

  #########
  # Define electron transport coefficients
  #########
  [./electron_moments_ad]
    type = ADGasElectronMoments
    block = 1
    em = em
    mean_en = mean_en
    property_tables_file = 'data/electron_mobility_diffusion.txt'
  [../]

  #########
  # Define some necessary constants.
  # Energy of secondary electrons is set to 1 eV here.
  # e: fundamental charge
  # N_A: Avogadro's number
  # k_boltz: Boltzmann constant
  # eps: permittivity of free space
  # se_energy: Secondary electron coefficient - emitted energy (set to 1 eV here)
  # T_gas: Gas temperature (Kelvin)
  # massem: 9.11e-31 (electron mass)
  # p_gas: Gas pressure (Pascals)
  #########
  [./gas_constants]
    type = GenericConstantMaterial
    block = 1
    prop_names = ' e       N_A      k_boltz  eps         se_energy T_gas  massem   p_gas'
    prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12   2.5       400    9.11e-31 1.01e5'
  [../]


  [test]
    #type = SurfaceCharge
    type = SurfaceChargeNew
    include_electrons = true
    potential = potential_dom1
    em = 'em'
    mean_en = mean_en
    #ions = 'Arp'
    species = 'Arp'
    position_units = ${dom1Scale}
    boundary = 'b1_left b1_right'
  []
  #[test]
  #  type = ADGenericConstantMaterial
  #  prop_names = 'surface_charge'
  #  prop_values = '0'
  #  boundary = 'b1_left b1_right'
  #[]

  ######
  # Diffusion coefficients for the potential in each region must be defined.
  # Here a GenericConstantMaterial is used for simplicity.
  # For all variables, potential included, diffusivity is defined as:
  # "diff" + [variable name]
  ######
  [./dielectric_left_side]
    type = GenericConstantMaterial
    prop_names = 'diffpotential_dom0'
    prop_values = '8.85e-11'
    block = 0
  [../]
  [./gas_phase]
    type = GenericConstantMaterial
    prop_names = 'diffpotential_dom1'
    prop_values = '8.85e-12'
    block = 1
  [../]
  [./dielectric_right_side]
    type = GenericConstantMaterial
    prop_names = 'diffpotential_dom2'
    prop_values = '8.85e-11'
    block = 2
  [../]

  ######
  # HeavySpeciesMaterial defines charge, mass, transport coefficients, and
  # temperature for each species.
  #
  # Transport coefficients and temperature are defined as ADMaterialProperties.
  # Although they currently (as of June 16, 2020) remain constant, future
  # implementations may include mixture-averaged diffusion coefficients and
  # effective ion temperatures with nonlinear dependence on other variables.
  ######
  [./gas_species_0]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Arp
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 1.0
    diffusivity = 1.6897e-5
    block = 1
  [../]
  [./gas_species_2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    diffusivity = 1.6897e-5
    block = 1
  [../]
  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    diffusivity = 1.6897e-5
    block = 1
  [../]
[]

[Reactions]
  # Here a list of plasma reactions are included to account for
  # electron-impact processes during the discharge.
  # All cross sections were taken from the COMSOL Plasma
  # Module example problem, "Dielectric Barrier Discharge".
  # (The example files are freely available.)
  #
  # Cross sections are typically found from LXCat.
  #
  # Rate coefficients and electron transport coefficients
  # (mobility, diffusivity) are computed from the cross sections
  # through an external Boltzmann solver. (In this case, BOLSIG+
  # was used assuming a Maxwellian electron distribution function.)
  # The tabulated rate coefficients are stored in the dbd_data
  # directory, as specified by the file_location input parameter.

  [./Argon]
    species = 'em Arp Ar*'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    file_location = 'data'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    potential = 'potential_dom1'
    use_log = true
    use_ad = true
    position_units = ${dom1Scale}
    block = 1

    # Note that rate coefficients are in molar units.
    # Two-body reaction units:  m^3 mol^-1 s^-1
    #
    # EEDF - rate coefficients taken from tabulated files located in 'file_location',
    #        computed from cross section data with Bolsig+.
    #
    # Energy changes are included in square brackets, with the units in eV. Note that the
    # sign is with respect to electrons; e.g. a negative sign indicates that electrons
    # lose energy in the reaction.
    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
                 em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
                 em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ars_De-excitation_11.50_eV)
                 em + Ar* -> em + em + Arp        : EEDF [-4.43] (C5_Ars_Ionization_4.43_eV)
                 Ar* + Ar -> Ar + Ar              : 1807
                 Ar* + Ar* -> em + Ar + Arp       : 3.3734e8'
  [../]
[]
