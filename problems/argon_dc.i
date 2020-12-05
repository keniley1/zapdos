dom0Scale=1

[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'dc_mesh.msh'
  []
  [left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = file
  []
  [right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  []
[]

[GlobalParams]
  #offset = 40
  potential_units = kV
  use_moles = true
[]

[Problem]
  type = FEProblem
  # kernel_coverage_check = false
[]

[Debug]
  #show_var_residual_norms = true
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
  compute_scaling_once = false
  # Setting an end time at 1000 seconds
  # In reality the simulation will reach steady state long before this
  end_time = 1000
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = NEWTON
  line_search = 'basic'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_max_its = 100
  nl_max_its = 20

  # Adaptive timesteppers are useful for plasma simulations since the
  # timestep requirements may vary widely over the course of a simulation.
  # In this case timesteps of 0.1-10 ns are required between pulses as the
  # electron current becomes large, but larger timesteps are allowed otherwise.
  #
  # A maximum timestep of 100 ns is set to make sure the voltage pulse (period
  # of 20 microseconds) is finely sampled.
  dtmin = 1e-14
  dtmax = 1
  steady_state_detection = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    growth_factor = 1.2
   optimal_iterations = 10
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
  [../]
[]

[UserObjects]
  [./data_provider]
    type = ProvideMobility
    electrode_area = 0.005
    ballast_resist = 10000
    e = 1.6e-19
  [../]
[]

[Variables]
  [potential]
  []
[]

[DriftDiffusionAction]
  [plasma]
    electrons = em
    charged_particle = 'Arp'
    Neutrals = 'Ar*'
    potential = potential
    mean_energy = mean_en
    using_offset = false
    offset = 30
    use_ad = true
    position_units = ${dom0Scale}
    Additional_Outputs = 'ElectronTemperature'
    block = 0
  []
[]

[AuxVariables]
  # Argon is considered to be a static background
  # Here it is defined as an AuxVariable
  # Units: log(mol m^-3)
  [./Ar]
    block = 0
    order = CONSTANT
    family = MONOMIAL
    initial_condition = -5.23352
    # n_Ar = exp(-5.23352) * 6.022e23
  [../]
[]

[BCs]
  [./potential_right]
    type = ADNeumannCircuitVoltageMoles_KV
    variable = potential
    boundary = 'right'
    function = potential_input
    ip = 'Arp'
    data_provider = data_provider
    em = em
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]

  # Ground is at x = 0 m, named 'left'
  [./potential_dirichlet_left]
    type = DirichletBC
    variable = potential
    boundary = 'left'
    value = 0
  [../]

  ######
  # HAGELAAR BCS
  ######
  #[./em_physical_bc]
  #  type = ADHagelaarElectronBC
  #  variable = em
  #  boundary = 'left right'
  #  potential = potential
  #  mean_en = mean_en
  #  r = 0.0
  #  position_units = ${dom0Scale}
  #[../]
  [./sec_electrons_left]
    type = ADSecondaryElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    ip = 'Arp'
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  [../]
  #[./mean_en_physical]
  #  type = ADHagelaarEnergyBC
  #  variable = mean_en
  #  boundary = 'left right'
  #  potential = potential
  #  em = em
  #  r = 0
  #  position_units = ${dom0Scale}
  #[../]
  [./secondary_energy_left]
    type = ADSecondaryElectronEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = em
    ip = 'Arp'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_bcs]
    type = ADHagelaarIonAdvectionBC
    variable = Arp
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  [../]
  [./Arp_bcs2]
    type = ADHagelaarIonDiffusionBC
    variable = Arp
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  [Arp_bc]
    type = IonBC
    variable = Arp
    boundary = 'left right'
    potential = potential
    r = 0
    position_units = ${dom0Scale}
  []
  [e_bc]
    type = ElectronBC
    variable = em
    boundary = 'left right'
    potential = potential
    mean_en = mean_en
    r = 0
    position_units = ${dom0Scale}
  []
  [mean_en_bc]
    type = ElectronEnergyBC
    variable = mean_en
    boundary = 'left right'
    potential = potential
    em = em
    r = 0
    position_units = ${dom0Scale}
  []
  [./Arex_bcs]
    type = ADHagelaarIonDiffusionBC
    variable = Ar*
    boundary = 'left right'
    r = 0
    position_units = ${dom0Scale}
  [../]
  #########
  #########
[]

[ICs]
  [./potential_ic]
    type = FunctionIC
    variable = potential
    function = potential_ic_func
  [../]

  [./em_ic]
    type = ConstantIC
    variable = em
    value = -41
    block = 0
  [../]
  [./Arp_ic]
    type = ConstantIC
    variable = Arp
    value = -41
    block = 0
  [../]
  [./Arex_ic]
    type = ConstantIC
    variable = Ar*
    value = -45
    block = 0
  [../]
  [./mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -41
    block = 0
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
    #vars = 'f0'
    #vals = '50e3'
    #value = '-0.75*sin(2*3.1415926*f0*t)'
    value = '0.2'
  [../]

  # Set the initial condition to a line from -10 V on the left and
  # 0 on the right.
  # (Poisson solver tends to struggle with a uniformly zero potential IC.)
  [./potential_ic_func]
    type = ParsedFunction
    #value = '-0.01 * (3.0001e-4 - x)'
    #value = '0.2 * (x - 0.4)'
    value = '0.5 * x'
  [../]
[]

[Materials]
  #########
  # Define secondary electron emission coefficients on the left and right
  # dielectrics.
  #########
  [./se_left]
    type = GenericConstantMaterial
    boundary = 'left'
    prop_names = 'se_coeff'
    prop_values = '0.35'
  [../]
  [./se_right]
    type = GenericConstantMaterial
    boundary = 'right'
    prop_names = 'se_coeff'
    prop_values = '0'
  [../]

  #########
  # Define electron transport coefficients
  #########
  [./electron_moments_ad]
    type = ADGasElectronMoments
    block = 0
    em = em
    mean_en = mean_en
    property_tables_file = 'dc_data/electron_mobility_diffusion.txt'
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
    block = 0
    prop_names = ' e       N_A      k_boltz  eps         se_energy T_gas  massem   p_gas'
    prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12   3.        300    9.11e-31 13.3'
    #prop_names = ' e       N_A      k_boltz  eps         se_energy T_gas  p_gas'
    #prop_values = '1.6e-19 6.022e23 1.38e-23 8.854e-12   3.        300  13.3'
  [../]


  ######
  # 'Diffusion coefficients' for the potential in each region must be defined.
  # Here a GenericConstantMaterial is used for simplicity.
  # For all variables, potential included, diffusivity is defined as:
  # "diff" + [variable name]
  ######
  [./gas_phase]
    type = GenericConstantMaterial
    prop_names = 'diffpotential'
    prop_values = '8.85e-12'
    block = 0
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
    diffusivity = 0.04
    block = 0
  [../]
  [./gas_species_2]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    diffusivity = 0.04
    block = 0
  [../]
  [./gas_species_3]
    type = ADHeavySpeciesMaterial
    heavy_species_name = Ar*
    heavy_species_mass = 6.64e-26
    heavy_species_charge = 0.0
    diffusivity = 0.04
    block = 0
  [../]
[]

[Reactions]
  active = 'all'
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
  [all]
    species = 'em Arp Ar*'
    aux_species = 'Ar'
    reaction_coefficient_format = 'townsend'
    file_location = 'dc_data'
    gas_species = 'Ar'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    potential = 'potential'
    use_log = true
    use_ad = true
    position_units = ${dom0Scale}
    block = 0

    reactions = 'em + Ar -> em + Ar               : EEDF [elastic] (C1_Ar_Elastic)
                 em + Ar -> em + Ar*              : EEDF [-11.5] (C2_Ar_Excitation_11.50_eV)
                 em + Ar -> em + em + Arp         : EEDF [-15.76] (C3_Ar_Ionization_15.80_eV)
                 em + Ar* -> em + Ar              : EEDF [11.5] (C4_Ars_De-excitation_11.50_eV)
                 em + Ar* -> em + em + Arp        : EEDF [-4.43] (C5_Ars_Ionization_4.43_eV)
                 Ar* + Ar -> Ar + Ar              : 1807
                 Ar* + Ar* -> em + Ar + Arp       : 3.3734e8'
  []
[]
