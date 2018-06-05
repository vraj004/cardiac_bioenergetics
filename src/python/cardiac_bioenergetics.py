# program to simulate 3D energy metabolism in cardiac cells

import os
from opencmiss.iron import iron

#problem parameters
length=5.0
numberOfElements=2
#user numbers
coordsysUserNumber=1
regionUserNumber=1
linbasisUserNumber=1
genmeshUserNumber=1
meshUserNumber=1
decompositionUserNumber=1
equationsSetUserNumber=1
equationsSetFieldUserNumber=3
geometricFieldUserNumber=1
dependentFieldUserNumber=2
materialFieldUserNumber=4
sourceFieldUserNumber=5
problemUserNumber=1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

#_________________Define coord sys and region________________#
coordsys = iron.CoordinateSystem()
coordsys.CreateStart(coordsysUserNumber)
coordsys.DimensionSet(1)
coordsys.CreateFinish()

region=iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem=coordsys
region.CreateFinish()

#_______________Define bases________________________________#
linbasis=iron.Basis()
linbasis.CreateStart(linbasisUserNumber)
linbasis.type=iron.BasisTypes.LAGRANGE_HERMITE_TP
linbasis.numberOfXi = 1
linbasis.interpolationXi=[iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*1
linbasis.quadratureNumberOfGaussXi = [3]
linbasis.CreateFinish()

#_______________Define mesh and decomposition________________#
genmesh=iron.GeneratedMesh()
genmesh.CreateStart(genmeshUserNumber,region)
genmesh.type=iron.GeneratedMeshTypes.REGULAR
genmesh.basis=[linbasis]
genmesh.extent=[length]
genmesh.numberOfElements=[numberOfElements]
mesh=iron.Mesh()
genmesh.CreateFinish(meshUserNumber,mesh)

decomposition=iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type=iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains=numberOfComputationalNodes
decomposition.CreateFinish()

#________________Define geometric field__________________#
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
genmesh.GeometricParametersCalculate(geometricField)

#____________________Define equations set and dependent field_________#
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
    iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CONSTANT_REAC_DIFF]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
equationsSet.DependentCreateFinish()
dependentField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

#__________________Define materials field___________________________#
# Create the material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
equationsSet.MaterialsCreateFinish()
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1.0)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,1.0)

#___________________Define source field______________________________#
sourceField=iron.Field()
equationsSet.SourceCreateStart(sourceFieldUserNumber,sourceField)
equationsSet.SourceCreateFinish()
sourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
sourceField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,1,10.0)

#___________________Define equations field______________________________#

equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#___________________Define problem and problem control loops______________________________#
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
        iron.ProblemSubtypes.CONSTANT_REAC_DIFF_NO_SPLIT]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.TimesSet(0.0,0.5,0.01)
controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
controlLoop.TimeOutputSet(1)
problem.ControlLoopCreateFinish()

#___________________Define problem solvers and boundary conditions______________________________#
# Create problem solver
problem.SolversCreateStart()
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
# Set no flux boundaries
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,1,1,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,3,1,iron.BoundaryConditionsTypes.FIXED,0.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/reaction_diffusion_1d","FORTRAN")
fields.ElementsExport("./results/reaction_diffusion_1d","FORTRAN")
fields.Finalise()
