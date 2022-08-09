#This tutorial is of an electronegative oxygen discharge
#model. This is explained in Lieberman, page 259.



#A uniform scaling factor of the mesh.
#E.g if set to 1.0, there is not scaling
# and if set to 0.010, there mesh is scaled by a cm
# See Lieberman pg 364 and figure 10.6
dom0Scale=45e-3

[GlobalParams]
  #Scales the potential by V or kV
  potential_units = kV
  #Converts density from #/m^3 to moles/m^3
  use_moles = true
[]

[Mesh]
  #Mesh is define by a Gmsh file
  [./geo]
    type = FileMeshGenerator
    file = 'Lymberopoulos_paper.msh'
  [../]
  #Renames all sides with the specified normal
  #For 1D, this is used to rename the end points of the mesh
  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = geo
  [../]
  [./right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  [../]
  uniform_refine = 1
[]

#Defines the problem type, such as FE, eigen value problem, etc.
[Problem]
  type = FEProblem
[]

[DriftDiffusionAction]
  [./Plasma]
    #User define name for electrons (usually 'em')
    electrons = 'em'
    #User define name for ions
    charged_particle = 'O- O2+'
    #User define name for potential (usually 'potential')
    potential = potential
    #Defines if this potential exist in only one block/material (set 'true' for single gases)
    Is_potential_unique = true
    #User define name for the electron mean energy density (usually 'mean_en')
    mean_energy = mean_en
    #The position scaling for the mesh, define at top of input file
    position_units = ${dom0Scale}
    #Additional outputs, such as ElectronTemperature, Current, and EField.
    Additional_Outputs = 'ElectronTemperature'
  [../]
[]

[Reactions]
  [./Oxygen]
    #Name of reactant species that are variables
    species = 'em O- O2+'
    #Name of reactant species that are auxvariables
    aux_species = 'O2'
    #Type of coefficient (rate or townsend)
    reaction_coefficient_format = 'rate'
    #Name of background gas
    gas_species = 'O2'
    #Name of the electron mean energy density (usually 'mean_en')
    electron_energy = 'mean_en'
    #Name of the electrons (usually 'em')
    electron_density = 'em'
    #Defines if electrons are tracked
    include_electrons = true
    #Defines directory holding rate text files
    file_location = 'rate_coefficients'
    #Name of name for potential (usually 'potential')
    potential = 'potential'
    #Defines if log form is used (true for Zapdos)
    use_log = true
    #Defines if automatic differentiation is used (true for Zapdos)
    use_ad = true
    #The position scaling for the mesh, define at top of input file
    position_units = ${dom0Scale}
    #Name of material block ('0' for an user undefined block)
    block = 0
    #Inputs of the plasma chemsity

    # These are parameters required equation-based rate coefficients
    equation_variables = 'e_temp'

    #e.g. Reaction : Constant or EEDF dependent [Threshold Energy] (Text file name)
    #     ionization, diss. attachment, recombination pg 365
    #     Multiplied by avogadro's to convert m^3/s -> moles*m^3/s
    reactions = 'em + O2 -> em + em + O2+ : {(2.13*10e-14*exp(-14.5/e_temp))*6.022e23} [-12.06]
                 em + O2 -> O + O- : {(7.89*10e-17*exp(-3.07/e_temp))*6.022e23} [0]
                 O2+ + O- -> O2 + O : {(1.4*10e-13)*6.022e23}'
  [../]
[]

[AuxVariables]
  #Add a scaled position units used for plotting other element AuxVariables
  [./x_node]
  [../]

  #Background gas (e.g Ar)
  [./O2]
  [../]
[]

[AuxKernels]
  #Add at scaled position units used for plotting other element AuxVariables
  [./x_ng]
    type = Position
    variable = x_node
    position_units = ${dom0Scale}
  [../]

  #Background gas number density (e.g. for 1Torr)
  [./O2_val]
    type = FunctionAux
    variable = O2
    function = 'log(1.6e21/6.022e23)'
    execute_on = INITIAL
  [../]
[]

#Currently there is no Action for BC (but one is currently in development)
#Below is the Lymberopulos family of BC
#(For other BC example, please look at Tutorial 04 and Tutorial 06)
[BCs]
#Voltage Boundary Condition Ffor a Power-Ground RF Discharge
  [./potential_left]
    type = FunctionDirichletBC
    variable = potential
    boundary = 'left'
    function = potential_left_bc_func
    preset = false
  [../]
  [./potential_dirichlet_right]
    type = FunctionDirichletBC
    variable = potential
    boundary = 'right'
    #value = 0
    function = potential_right_bc_func
    preset = false
  [../]

  #Boundary conditions for electons
  [./em_physical_diffusion]
    type = SakiyamaElectronDiffusionBC
    variable = em
    mean_en = mean_en
    boundary = 'left right'
    position_units = ${dom0Scale}
  [../]

  #Boundary conditions for ions
    [./O-_physical_advection]
      type = SakiyamaIonAdvectionBC
      variable = 'O-'
      potential = potential
      boundary = 'left right'
      position_units = ${dom0Scale}
    [../]

    [./O2+_physical_advection]
      type = SakiyamaIonAdvectionBC
      variable = 'O2+'
      potential = potential
      boundary = 'left right'
      position_units = ${dom0Scale}
    [../]


  #New Boundary conditions for mean energy, should be the same as in paper
    [./mean_en_physical_diffusion]
      type = SakiyamaEnergyDiffusionBC
      variable = mean_en
      em = em
      boundary = 'left right'
      position_units = ${dom0Scale}
    [../]
  []

#Initial conditions for variables.
#If left undefine, the IC is zero
[ICs]
  [./em_ic]
    type = FunctionIC
    variable = em
    function = density_ic_func
  [../]
  [./O-_ic]
    type = FunctionIC
    variable = O-
    function = ion_density_ic_func
  [../]
  [./O2+_ic]
    type = FunctionIC
    variable = O2+
    function = ion_density_ic_func
  [../]
  [./mean_en_ic]
    type = FunctionIC
    variable = mean_en
    function = energy_density_ic_func
  [../]
  #[./potential_ic]
    #type = FunctionIC
    #variable = potential
    #function = potential_ic_func
  #[../]
[]

#Define function used throughout the input file (e.g. BCs and ICs)
[Functions]
  [./potential_left_bc_func]
    type = ParsedFunction
    value = '-0.100*sin(2*3.1415926*13.56e6*t)'
  [../]
  [./potential_right_bc_func]
    type = ParsedFunction
    value = '0.100*sin(2*3.1415926*13.56e6*t)'
  [../]
  [./potential_ic_func]
    type = ParsedFunction
    value = '0.100 * (25.4e-3 - x)'
  [../]
  [./ion_density_ic_func]
    type = ParsedFunction
    value = 'log((1e14 + 1e16 * (1-x/1)^2 * (x/1)^2)/6.022e23)'
  [../]
  [./density_ic_func]
    type = ParsedFunction
    value = 'log((1e13 + 1e15 * (1-x/1)^2 * (x/1)^2)/6.022e23)'
  [../]
  [./energy_density_ic_func]
    type = ParsedFunction
    value = 'log(3./2.) + log((1e13 + 1e15 * (1-x/1)^2 * (x/1)^2)/6.022e23)'
  [../]
[]

[Materials]
  #The material properties for electrons.
  #Also hold universal constant, such as Avogadro's number, elementary charge, etc.
  [./GasBasics]
    type = GasElectronMoments
    #False means constant electron coeff, defined by user
    interp_trans_coeffs = false
    #Leave as false (CRANE accounts of elastic coeff.)
    interp_elastic_coeff = false
    #Leave as false, unless computational error is due to rapid coeff. changes
    ramp_trans_coeffs = false
    #User difine pressure in pa
    user_p_gas = 6.666
    #Name for electrons (usually 'em')
    em = em
    #Name for potential (usually 'potential')
    potential = potential
    #Name for the electron mean energy density (usually 'mean_en')
    mean_en = mean_en
    #User define electron mobility coeff. (define as 0.0 if not used)
    user_electron_mobility = 30.0 #2.0160e24
    #User define electron diffusion coeff. (define as 0.0 if not used)
    user_electron_diffusion_coeff = 119.8757763975 #5.96e24 119.8757763975
    pressure_dependent_electron_coeff = false
    #Name of text file with electron properties
    property_tables_file = rate_coefficients/oxygen_rxn.txt
  [../]
  #The material properties of the ion
  [./gas_species_O-]
    type = ADHeavySpecies
    heavy_species_name = O-
    heavy_species_mass = 2.6559e-26
    heavy_species_charge = -1.0
    mobility = 9.5322
    diffusivity = 0.2468
  [../]
  [./gas_species_O2+]
    type = ADHeavySpecies
    heavy_species_name = O2+
    heavy_species_mass = 5.3137e-26
    heavy_species_charge = 1.0
    mobility = 4.7644
    diffusivity = 0.1234
  [../]
  #The material properties of the background gas
  [./gas_species_2]
    type = ADHeavySpecies
    heavy_species_name = O2
    heavy_species_mass = 5.3137e-26
    heavy_species_charge = 0.0
  [../]
[]

#Preconditioning options
#Learn more at: https://mooseframework.inl.gov/syntax/Preconditioning/index.html
[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]

  [./fdp]
    type = FDP
    full = true
  [../]
[]

#How to execute the problem.
#Defines type of solve (such as steady or transient),
# solve type (Newton, PJFNK, etc.) and tolerances
[Executioner]
  type = Transient
  end_time = 7.3746e-5
  dt = 1e-9
  dtmin = 1e-14
  scheme = bdf2
  solve_type = NEWTON

  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_linesearch_minlambda'
  petsc_options_value = 'lu NONZERO 1.e-10 1e-3'

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-7
  l_max_its = 20
  line_search = none
  automatic_scaling = true
  compute_scaling_once = false

[]

#Defines the output type of the file (multiple output files can be define per run)
[Outputs]
  perf_graph = true
  file_base = 'oxygen_plasma'
  [./out]
    type = Exodus
  [../]
[]
