!> \file
!> $Id: cardiac_ecc.f90 2014-12-16 ghosh_shourya $
!> \author Vijay Rajagopal
!> \brief Main program file to simulate cardiac bioenergetics using opencmfe_ library routines
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is Opencmfe_
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>
!> Main program
PROGRAM Cardiac_bioenergetics

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE ISO_C_BINDING
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=6

  INTEGER(CMISSIntg), PARAMETER :: ATPFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: ATPMaterialsFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ATPSourceFieldUserNumber=17


  INTEGER(CMISSIntg), PARAMETER :: ADPFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ADPMaterialsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: ADPEquationsSetUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: ADPEquationsSetFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: ADPSourceFieldUserNumber=23

  INTEGER(CMISSIntg), PARAMETER :: AMPFieldUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: AMPMaterialsFieldUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: AMPEquationsSetUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: AMPEquationsSetFieldUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: AMPSourceFieldUserNumber=28

  INTEGER(CMISSIntg), PARAMETER :: PCrFieldUserNumber=29
  INTEGER(CMISSIntg), PARAMETER :: PCrMaterialsFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: PCrEquationsSetUserNumber=31
  INTEGER(CMISSIntg), PARAMETER :: PCrEquationsSetFieldUserNumber=32
  INTEGER(CMISSIntg), PARAMETER :: PCrSourceFieldUserNumber=33

  INTEGER(CMISSIntg), PARAMETER :: CrFieldUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: CrMaterialsFieldUserNumber=35
  INTEGER(CMISSIntg), PARAMETER :: CrEquationsSetUserNumber=36
  INTEGER(CMISSIntg), PARAMETER :: CrEquationsSetFieldUserNumber=37
  INTEGER(CMISSIntg), PARAMETER :: CrSourceFieldUserNumber=38

  INTEGER(CMISSIntg), PARAMETER :: PiFieldUserNumber=39
  INTEGER(CMISSIntg), PARAMETER :: PiMaterialsFieldUserNumber=40
  INTEGER(CMISSIntg), PARAMETER :: PiEquationsSetUserNumber=41
  INTEGER(CMISSIntg), PARAMETER :: PiEquationsSetFieldUserNumber=42
  INTEGER(CMISSIntg), PARAMETER :: PiSourceFieldUserNumber=43


  INTEGER(CMISSIntg), PARAMETER :: OxyFieldUserNumber=69
  INTEGER(CMISSIntg), PARAMETER :: OxyMaterialsFieldUserNumber=70
  INTEGER(CMISSIntg), PARAMETER :: OxyEquationsSetUserNumber=71
  INTEGER(CMISSIntg), PARAMETER :: OxyEquationsSetFieldUserNumber=72
  INTEGER(CMISSIntg), PARAMETER :: OxySourceFieldUserNumber=73


  INTEGER(CMISSIntg), PARAMETER :: DPsiFieldUserNumber=44
  INTEGER(CMISSIntg), PARAMETER :: DPsiMaterialsFieldUserNumber=45
  INTEGER(CMISSIntg), PARAMETER :: DPsiEquationsSetUserNumber=46
  INTEGER(CMISSIntg), PARAMETER :: DPsiEquationsSetFieldUserNumber=47
  INTEGER(CMISSIntg), PARAMETER :: DPsiSourceFieldUserNumber=48

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: BCFieldUserNumber=74

  INTEGER(CMISSIntg),DIMENSION(2) :: BCNODES
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER :: node
  REAL(CMISSDP) :: VALUE, BCVALUE
  INTEGER(CMISSIntg) :: MitocondriaIndex, MyofibrilIndex

  
 
  !CMISS variables
  !setting up fields 
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField,BCField
  TYPE(cmfe_FieldType) :: ATPField,ATPMaterialsField,ATPEquationsSetField,ATPSourceField
  TYPE(cmfe_FieldType) :: ADPField,ADPMaterialsField,ADPEquationsSetField,ADPSourceField
  TYPE(cmfe_FieldType) :: AMPField,AMPMaterialsField,AMPEquationsSetField,AMPSourceField
  TYPE(cmfe_FieldType) :: PCrField,PCrMaterialsField,PCrEquationsSetField,PCrSourceField
  TYPE(cmfe_FieldType) :: CrField,CrMaterialsField,CrEquationsSetField,CrSourceField
  TYPE(cmfe_FieldType) :: PiField,PiMaterialsField,PiEquationsSetField,PiSourceField
  TYPE(cmfe_FieldType) :: OxyField,OxyMaterialsField,OxyEquationsSetField,OxySourceField
  TYPE(cmfe_FieldType) :: DPsiField,DPsiMaterialsField,DPsiEquationsSetField,DPsiSourceField
  TYPE(cmfe_EquationsType) :: ATPEquations, ADPEquations,AMPEquations, PCrEquations,CrEquations, PiEquations
  TYPE(cmfe_EquationsType) :: DPsiEquations, OxyEquations
  TYPE(cmfe_EquationsSetType) :: ATPEquationsSet, ADPEquationsSet,AMPEquationsSet, PCrEquationsSet,CrEquationsSet   
  TYPE(cmfe_EquationsSetType) :: PiEquationsSet, DPsiEquationsSet, Oxyequationsset

  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: MeshElements
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  

!Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES
  INTEGER :: st,i,NUMBER_OF_COORDS,NODES_PER_ELE, it, nod, NUMBER_OF_NODES,NUMBER_OF_ELEMENTS,element
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:) :: RyRDensity, VolMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  REAL(CMISSDP) :: nodex,nodey,nodez,val,ave, zeroval
  LOGICAL :: EXPORT_FIELD=.FALSE.
  INTEGER(CMISSIntg) :: Err,EquationsSetIndex,NODE_NUMBER,CONDITION,CellMLIndex,ELEM_NUMBER,NUMBER_OF_MTOS
  INTEGER(CMISSIntg),DIMENSION(166) :: CELLBOUNDARYNODES
  INTEGER(CMISSIntg) :: SL_BD_MARKER,MITO_BD_MARKER,MITO_REGION_MARKER,WITH_MITO_ELEMENTS,ELEM_LABEL

  REAL(CMISSDP) :: startT,endT,Tstep,ODE_TIME_STEP
  REAL(CMISSDP) :: init_ATP,ATPDiffx,ATPDiffy,ATPDiffz
  REAL(CMISSDP) :: init_ADP,ADPDiffx,ADPDiffy,ADPDiffz
  REAL(CMISSDP) :: init_AMP,AMPDiffx,AMPDiffy,AMPDiffz
  REAL(CMISSDP) :: init_PCr,PCrDiffx,PCrDiffy,PCrDiffz
  REAL(CMISSDP) :: init_Cr,CrDiffx,CrDiffy,CrDiffz
  REAL(CMISSDP) :: init_Pi,PiDiffx,PiDiffy,PiDiffz
  REAL(CMISSDP) :: init_Oxy,OxyDiffx,OxyDiffy,OxyDiffz
  REAL(CMISSDP) :: init_DPsi,DPsiDiffx,DPsiDiffy,DPsiDiffz


  INTEGER(CMISSIntg) :: simplex_order
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,CELLPATH,RyRModel,RYRDENSITYFILE,MITOVOLS
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NodeDomain,ElementDomain


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
    
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!_________________________________________________________________________________________________
  !Read inputs.txt 

  open(unit=9,file='input.txt',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening inputs file',st
    STOP
  ELSE
    PRINT *,'inputs file opened correctly'
    READ(9,*) !file info
    READ(9,*) NODEFILE,ELEMFILE
    READ(9,*) !comment
    READ(9,*) simplex_order
    READ(9,*) !atp
    READ(9,*) init_ATP,ATPDiffx,ATPDiffy,ATPDiffz
    READ(9,*) !adp
    READ(9,*) init_ADP,ADPDiffx,ADPDiffy,ADPDiffz
    READ(9,*) !amp
    READ(9,*) init_AMP,AMPDiffx,AMPDiffy,AMPDiffz
    READ(9,*) !pcr
    READ(9,*) init_PCr,PCrDiffx,PCrDiffy,PCrDiffz
    READ(9,*) !cr
    READ(9,*) init_Cr,CrDiffx,CrDiffy,CrDiffz
    READ(9,*) !pi
    READ(9,*) init_Pi,PiDiffx,PiDiffy,PiDiffz
    READ(9,*) !Oxy
    READ(9,*) init_Oxy,OxyDiffx,OxyDiffy,OxyDiffz
    READ(9,*) !DPsi
    READ(9,*) init_DPsi,DPsiDiffx,DPsiDiffy,DPsiDiffz
    READ(9,*) !Extra
    READ(9,*) MITOVOLS
  ENDIF
  CLOSE(9)
  zeroval= 0.0_CMISSDP
  EXPORT_FIELD=.TRUE.
!_________________________________________________________________________________________________
  !Intialise Opencmfe_
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(cmfe_ERRORS_TRAP_ERROR,Err)
  !get computational nodes for parallel processing
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  write (*,*) "NumberOfComputationalNodes:   ",NumberOfComputationalNodes
  write (*,*) "ComputationalNodeNumber:  ",ComputationalNodeNumber

   !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemUserNumber,2,Err)
  !The coordinate system is 3D by default;set it to be 2D in above command.
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_LabelSet(Region,"Cell",Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

!_________________________________________________________________________________________________
  !Start the creation of a trilinear-simplex basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,cmfe_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
  write (*,*) "Simplex_Order:   ",simplex_order
  CALL cmfe_Basis_CreateFinish(Basis,Err)

!___________________________________________________________________________________________________
   !Time to create a mesh - wohoo!
  !CODE TO READ in nodes and elements

  WRITE(*,*) NODEFILE  
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS+2))

    IF(BOUNDARY_MARKER.EQ.1) THEN
      ALLOCATE(NodeNums(NUMBER_OF_NODES,3))

      DO i = 1,NUMBER_OF_NODES
        READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,3)
      ENDDO
    ELSE
      ALLOCATE(NodeNums(NUMBER_OF_NODES,1))
      DO i = 1,NUMBER_OF_NODES
        READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2)
      ENDDO
    ENDIF
  ENDIF
  CLOSE(10)
  !Read in elements

  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES

    IF(ELE_ATTRIBUTES.EQ.0) THEN
      ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,(NODES_PER_ELE+1)))
      IF(simplex_order==1) THEN
        DO i = 1,NUMBER_OF_ELEMENTS
          READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4)
        ENDDO
      ELSEIF(simplex_order==2) THEN
        DO i = 1,NUMBER_OF_ELEMENTS
          READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5), &
            & ElemMap(i,6),ElemMap(i,7),ElemMap(i,8),ElemMap(i,9),ElemMap(i,10),ElemMap(i,11)
        ENDDO
      ENDIF
    ELSE
      ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,(NODES_PER_ELE+2)))
      IF(simplex_order==1) THEN
        DO i = 1,NUMBER_OF_ELEMENTS
          READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
        ENDDO
      ELSEIF(simplex_order==2) THEN
        DO i = 1,NUMBER_OF_ELEMENTS
          READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5), &
            & ElemMap(i,6),ElemMap(i,7),ElemMap(i,8),ElemMap(i,9), &
            & ElemMap(i,10),ElemMap(i,11),ElemMap(i,12)
        ENDDO
      ENDIF

    ENDIF
  ENDIF 
  CLOSE(11)

 WRITE(*,*) MITOVOLS
!____________________________________________________________________________________________________

!Use array of nodes and elements for setting up mesh.
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  write (*,*) "NUMBER_OF_NODES:   ",NUMBER_OF_NODES
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  write (*,*) "NUMBER_OF_ELEMENTS:   ",NUMBER_OF_ELEMENTS
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL cmfe_MeshElements_Initialise(MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  IF(simplex_order==1) THEN
    DO i = 1,NUMBER_OF_ELEMENTS
      element = ElemMap(i,1)
      CALL cmfe_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
       &   ElemMap(i,4)/),Err)
    ENDDO
  ELSEIF(simplex_order==2) THEN
    DO i = 1,NUMBER_OF_ELEMENTS
      element = ElemMap(i,1)
      CALL cmfe_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
       &  ElemMap(i,4),ElemMap(i,5),ElemMap(i,6),ElemMap(i,7),ElemMap(i,8), &
       &  ElemMap(i,9),ElemMap(i,10),ElemMap(i,11)/),Err)
    ENDDO
  ENDIF
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)
!________________________________________________________________________________________________________
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,cmfe_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  write (*,*) "NumberOfComputationalNodes:  ",NumberOfComputationalNodes

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 3 field components in 1 mesh component
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,2,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)


  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)*0.002
      nodey = NodeCoords(i,2)*0.002
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex, Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
  CALL cmfe_Field_ParameterSetGetNode(RegionUserNumber, GeometricFieldUserNumber, &
  & cmfe_FIELD_U_VARIABLE_TYPE, cmfe_FIELD_VALUES_SET_TYPE, 1, 1, node, 1, val, err)	
     ENDIF
    ENDDO
  !Bottom two lines ensure that the updated nodal values are updated across all cpus
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"mesh","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"mesh","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF

!___________________________________________________________________________________________________________________________


  !################( ATP )#############################################################--->

!Create ATP  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(ATPEquationsSet,Err)
  CALL cmfe_Field_Initialise(ATPEquationsSetField,Err)

  CALL cmfe_EquationsSet_CreateStart(ATPEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ATPEquationsSetFieldUserNumber,ATPEquationsSetField,ATPEquationsSet,Err)

  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(ATPEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(ATPField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(ATPEquationsSet,ATPFieldUserNumber,ATPField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,"ATP Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(ATPEquationsSet,Err)
  !set up node 2 to have initial dependent value as 0.5.

  !Create the equations set material field variables
  !currently the opencmfe_ code assumes storage coeff is always 1.
  !by default 3 comps for diff no source i.e. diff coeffs in 3 directions are set constant    	spatially = 1
  CALL cmfe_Field_Initialise(ATPMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(ATPEquationsSet,ATPMaterialsFieldUserNumber,ATPMaterialsField,Err)
  
  !Set to node based interpolation
  CALL cmfe_Field_ComponentInterpolationSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)

  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(ATPEquationsSet,Err)

  ! First set mitochondrial  diffusivity to be zero everhwhere  
   CALL cmfe_Field_ComponentValuesInitialise(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x

  !Set IMS diffusivity to be 1% of myofibril diffusivity

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.AND.NodeCoords(i,3)>10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*ATPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*ATPDiffy,Err)
       ENDIF
     ENDIF
      
     IF(NodeNums(i,3)==20.OR.NodeNums(i,3)==15) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*ATPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*ATPDiffy,Err)
       ENDIF
     ENDIF    
 
   !Set myofibrils diffusivity to be 100% of inputs.txt values       
      IF(NodeCoords(i,3)==10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,ATPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,ATPDiffy,Err)
       ENDIF
     ENDIF
  ENDDO



  CALL cmfe_Field_ParameterSetUpdateStart(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(ATPSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(ATPEquationsSet,ATPSourceFieldUserNumber,ATPSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"ATP Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(ATPEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(ATPSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)


  !################( ADP )#############################################################--->

!Create ADP  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(ADPEquationsSet,Err)
  CALL cmfe_Field_Initialise(ADPEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(ADPEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ADPEquationsSetFieldUserNumber,ADPEquationsSetField,ADPEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(ADPEquationsSet,Err)

  CALL cmfe_Field_Initialise(ADPField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(ADPEquationsSet,ADPFieldUserNumber,ADPField,Err)
  CALL cmfe_Field_VariableLabelSet(ADPField,cmfe_FIELD_U_VARIABLE_TYPE,"ADP Field",Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(ADPEquationsSet,Err)

  CALL cmfe_Field_Initialise(ADPMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(ADPEquationsSet,ADPMaterialsFieldUserNumber,ADPMaterialsField,Err)
   !Set to node based interpolation
  CALL cmfe_Field_ComponentInterpolationSet(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err) 
 CALL cmfe_EquationsSet_MaterialsCreateFinish(ADPEquationsSet,Err)


  CALL cmfe_Field_ComponentValuesInitialise(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.AND.NodeCoords(i,3)>10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*ADPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*ADPDiffy,Err)
       ENDIF
     ENDIF
        
          IF(NodeNums(i,3)==20.OR.NodeNums(i,3)==15) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*ADPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*ADPDiffy,Err)
       ENDIF
     ENDIF    
     		     
      IF(NodeCoords(i,3)==10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,ADPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,ADPDiffy,Err)
       ENDIF
     ENDIF
  ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(ADPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(ADPSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(ADPEquationsSet,ADPSourceFieldUserNumber,ADPSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(ADPSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"ADP Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(ADPEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(ADPSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set


  !################( AMP )#############################################################--->

!Create AMP  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(AMPEquationsSet,Err)
  CALL cmfe_Field_Initialise(AMPEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(AMPEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & AMPEquationsSetFieldUserNumber,AMPEquationsSetField,AMPEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(AMPEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(AMPField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(AMPEquationsSet,AMPFieldUserNumber,AMPField,Err)
  CALL cmfe_Field_VariableLabelSet(AMPField,cmfe_FIELD_U_VARIABLE_TYPE,"AMP Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(AMPEquationsSet,Err)
  !set up node 2 to have initial dependent value as 0.5.

  !Create the equations set material field variables
  !currently the opencmfe_ code assumes storage coeff is always 1.
  !by default 3 comps for diff no source i.e. diff coeffs in 3 directions are set constant    	spatially = 1
  CALL cmfe_Field_Initialise(AMPMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(AMPEquationsSet,AMPMaterialsFieldUserNumber,AMPMaterialsField,Err)
  CALL cmfe_Field_ComponentInterpolationSet(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(AMPEquationsSet,Err)


  CALL cmfe_Field_ComponentValuesInitialise(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x


    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.AND.NodeCoords(i,3)>10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*AMPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*AMPDiffy,Err)
       ENDIF
     ENDIF
        
      IF(NodeNums(i,3)==20.OR.NodeNums(i,3)==15) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.01_CMISSDP*AMPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,0.01_CMISSDP*AMPDiffy,Err)
       ENDIF
     ENDIF 
        
      IF(NodeCoords(i,3)==10) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,AMPDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,AMPDiffy,Err)
       ENDIF
     ENDIF
  ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(AMPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(AMPSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(AMPEquationsSet,AMPSourceFieldUserNumber,AMPSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(AMPSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"AMP Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(AMPEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(AMPSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)



  !################( PCr )#############################################################--->

!Create PCr  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(PCrEquationsSet,Err)
  CALL cmfe_Field_Initialise(PCrEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(PCrEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & PCrEquationsSetFieldUserNumber,PCrEquationsSetField,PCrEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(PCrEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(PCrField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(PCrEquationsSet,PCrFieldUserNumber,PCrField,Err)
  CALL cmfe_Field_VariableLabelSet(PCrField,cmfe_FIELD_U_VARIABLE_TYPE,"PCr Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(PCrEquationsSet,Err)
  !set up node 2 to have initial dependent value as 0.5.

  !Create the equations set material field variables
  !currently the opencmfe_ code assumes storage coeff is always 1.
  !by default 3 comps for diff no source i.e. diff coeffs in 3 directions are set constant    	spatially = 1
  CALL cmfe_Field_Initialise(PCrMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(PCrEquationsSet,PCrMaterialsFieldUserNumber,PCrMaterialsField,Err)
    CALL cmfe_Field_ComponentInterpolationSet(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(PCrEquationsSet,Err)


  CALL cmfe_Field_ComponentValuesInitialise(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x

  !Set same diffusivity values in IMS and myofibrils

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.OR.NodeNums(i,3)==20) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,PCrDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,PCrDiffy,Err)
       ENDIF
     ENDIF
   ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(PCrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(PCrSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(PCrEquationsSet,PCrSourceFieldUserNumber,PCrSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(PCrSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"PCr Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(PCrEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(PCrSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)



  !################( Cr )#############################################################--->

!Create Cr  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(CrEquationsSet,Err)
  CALL cmfe_Field_Initialise(CrEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(CrEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CrEquationsSetFieldUserNumber,CrEquationsSetField,CrEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(CrEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(CrField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(CrEquationsSet,CrFieldUserNumber,CrField,Err)
  CALL cmfe_Field_VariableLabelSet(CrField,cmfe_FIELD_U_VARIABLE_TYPE,"Cr Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(CrEquationsSet,Err)
  !set up node 2 to have initial dependent value as 0.5.

  !Create the equations set material field variables
  !currently the opencmfe_ code assumes storage coeff is always 1.
  !by default 3 comps for diff no source i.e. diff coeffs in 3 directions are set constant    	spatially = 1
  CALL cmfe_Field_Initialise(CrMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(CrEquationsSet,CrMaterialsFieldUserNumber,CrMaterialsField,Err)
  CALL cmfe_Field_ComponentInterpolationSet(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err) 
 !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(CrEquationsSet,Err)


  CALL cmfe_Field_ComponentValuesInitialise(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.OR.NodeNums(i,3)==20) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,CrDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,CrDiffy,Err)
       ENDIF
     ENDIF
   ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CrMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(CrSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(CrEquationsSet,CrSourceFieldUserNumber,CrSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(CrSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"Cr Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(CrEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(CrSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)
  WRITE(*,'(A)') "00000000000000000000000000000000000000000000000000000000000."



 !################( Pi )#############################################################--->

!Create Pi  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(PiEquationsSet,Err)
  CALL cmfe_Field_Initialise(PiEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(PiEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & PiEquationsSetFieldUserNumber,PiEquationsSetField,PiEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish Pieating the equations set
  CALL cmfe_EquationsSet_CreateFinish(PiEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(PiField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(PiEquationsSet,PiFieldUserNumber,PiField,Err)
  CALL cmfe_Field_VariableLabelSet(PiField,cmfe_FIELD_U_VARIABLE_TYPE,"Pi Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(PiEquationsSet,Err)
  !set up node 2 to have initial dependentll value as 0.5.

  !Create the equations set material field variables
  !currently the opencmfe_ code assumes storage coeff is always 1.
  !by default 3 comps for diff no source i.e. diff coeffs in 3 directions are set constant    	spatially = 1
  CALL cmfe_Field_Initialise(PiMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart 	(PiEquationsSet,PiMaterialsFieldUserNumber,PiMaterialsField,Err)
  CALL cmfe_Field_ComponentInterpolationSet(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(PiEquationsSet,Err)


  CALL cmfe_Field_ComponentValuesInitialise(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,0.0_CMISSDP,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,0.0_CMISSDP,Err) !diff coeff in x

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.OR.NodeNums(i,3)==20) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,PiDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,PiDiffy,Err)
       ENDIF
     ENDIF
   ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(PiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(PiSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(PiEquationsSet,PiSourceFieldUserNumber,PiSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(PiSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"Pi Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(PiEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(PiSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)

!################( Oxy )#############################################################--->

!Create Oxy  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(OxyEquationsSet,Err)
  CALL cmfe_Field_Initialise(OxyEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(OxyEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & OxyEquationsSetFieldUserNumber,OxyEquationsSetField,OxyEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(OxyEquationsSet,Err)

  CALL cmfe_Field_Initialise(OxyField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(OxyEquationsSet,OxyFieldUserNumber,OxyField,Err)
  CALL cmfe_Field_VariableLabelSet(OxyField,cmfe_FIELD_U_VARIABLE_TYPE,"Oxy Field",Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(OxyEquationsSet,Err)


    CALL cmfe_Field_Initialise(OxyMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(OxyEquationsSet,OxyMaterialsFieldUserNumber,OxyMaterialsField,Err)
  CALL cmfe_Field_ComponentInterpolationSet(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(OxyEquationsSet,Err)


  !Set up same diffusivity everywhere in the cell

  CALL cmfe_Field_ComponentValuesInitialise(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,OxyDiffx,Err) !diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,OxyDiffx,Err) !diff coeff in x

    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)<20.OR.NodeNums(i,3)==20) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,OxyDiffx,Err)
         CALL cmfe_Field_ParameterSetUpdateNode(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),2,OxyDiffy,Err)
       ENDIF
     ENDIF
   ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(OxyMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(OxySourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(OxyEquationsSet,OxySourceFieldUserNumber,OxySourceField,Err)
  CALL cmfe_Field_VariableLabelSet(OxySourceField,cmfe_FIELD_U_VARIABLE_TYPE,"Oxy Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(OxyEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(OxySourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set

!################( DPsi )#############################################################--->

!Create DPsi  reaction diffusion with constant source equations_set 
  CALL cmfe_EquationsSet_Initialise(DPsiEquationsSet,Err)
  CALL cmfe_Field_Initialise(DPsiEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(DPsiEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & DPsiEquationsSetFieldUserNumber,DPsiEquationsSetField,DPsiEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(DPsiEquationsSet,Err)

  CALL cmfe_Field_Initialise(DPsiField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(DPsiEquationsSet,DPsiFieldUserNumber,DPsiField,Err)
  CALL cmfe_Field_VariableLabelSet(DPsiField,cmfe_FIELD_U_VARIABLE_TYPE,"DPsi Field",Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(DPsiEquationsSet,Err)

    CALL cmfe_Field_Initialise(DPsiMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(DPsiEquationsSet,DPsiMaterialsFieldUserNumber,DPsiMaterialsField,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_NODE_BASED_INTERPOLATION,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(DPsiEquationsSet,Err)

  !Set up zero diffusivity everywhere in the cell

  CALL cmfe_Field_ComponentValuesInitialise(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 1,zeroval,Err) !diff coeff in x
  CALL cmfe_Field_ComponentValuesInitialise(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
   & cmfe_FIELD_VALUES_SET_TYPE, &
   & 2,zeroval,Err) !diff coeff in y

  CALL cmfe_Field_ParameterSetUpdateStart(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DPsiMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
   CALL cmfe_Field_Initialise(DPsiSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(DPsiEquationsSet,DPsiSourceFieldUserNumber,DPsiSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(DPsiSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"DPsi Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(DPsiEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DPsiSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set

 !#############################################################################--->



  WRITE(*,'(A)') "AAAAAAAAAAAAAAAAAAAAAAA."


!Start to set up CellML Fields

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a toy constant source model from a file
  CALL cmfe_CellML_ModelImport(CellML,"mitochondria.cellml",MitocondriaIndex,Err)
  CALL cmfe_CellML_ModelImport(CellML,"myofibril.cellml",MyofibrilIndex,Err)


  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  !These are effectively parameters that you know won't change in the course of the ode solving for one time step. i.e. fixed before running cellml, known in opencmfe_ and 
  !changed only in opencmfe_ - components of the parameters field
  
  CALL cmfe_CellML_VariableSetAsKnown(CellML,MitocondriaIndex,"general_constants/param",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,MyofibrilIndex,"ATP/param",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,MitocondriaIndex,"dO2_dt/V_VO2",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,MyofibrilIndex,"H_ATP/H_ATP",Err)
    WRITE(*,'(A)') "Time for CellML."

  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> Opencmfe_ field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !dependent field, solve the dae, and then put the result of the dae into the source field.


  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"ATPi/ATPi",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"ATPi/ATPi",cmfe_FIELD_VALUES_SET_TYPE, &
    & ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ADPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"ADPi/ADPi",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"ADPi/ADPi",cmfe_FIELD_VALUES_SET_TYPE, &
    & ADPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,AMPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"AMPi/AMPi",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"AMPi/AMPi",cmfe_FIELD_VALUES_SET_TYPE, &
    & AMPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,PCrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"PCri/PCri",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"PCri/PCri",cmfe_FIELD_VALUES_SET_TYPE, &
    & PCrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"Cri/Cri",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"Cri/Cri",cmfe_FIELD_VALUES_SET_TYPE, &
    & CrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,PiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"Pii/Pii",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"Pii/Pii",cmfe_FIELD_VALUES_SET_TYPE, &
    & PiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,OxyField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"dO2_dt/O2",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"dO2_dt/O2",cmfe_FIELD_VALUES_SET_TYPE, &
    & OxyField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DPsiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"dPsi_dt/dPsi",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MitocondriaIndex,"dPsi_dt/dPsi",cmfe_FIELD_VALUES_SET_TYPE, &
    & DPsiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)


  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPSourceField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MitocondriaIndex,"general_constants/param",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPSourceField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"ATP/param",cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"ATP/ATP",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"ATP/ATP",cmfe_FIELD_VALUES_SET_TYPE, &
    & ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ADPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"ADP/ADP",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"ADP/ADP",cmfe_FIELD_VALUES_SET_TYPE, &
    & ADPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,AMPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"AMP/AMP",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"AMP/AMP",cmfe_FIELD_VALUES_SET_TYPE, &
    & AMPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,PCrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"PCr/PCr",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"PCr/PCr",cmfe_FIELD_VALUES_SET_TYPE, &
    & PCrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"Cr/Cr",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"Cr/Cr",cmfe_FIELD_VALUES_SET_TYPE, &
    & CrField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,PiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & MyofibrilIndex,"Pi/Pi",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,MyofibrilIndex,"Pi/Pi",cmfe_FIELD_VALUES_SET_TYPE, &
    & PiField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Finish the creation of CellML <--> Opencmfe_ field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML, CellMLModelsFieldUserNumber, &
    & CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)
  !The CellMLModelsField is an integer field that stores which model is being used by which node.
  !By default all field parameters have default model value of 1, i.e. the first model. But, this command below is for example purposes
  WRITE(*,'(A)') "BBBBBBBBBBBBBBB."


  !By default all field parameters have default model value of 1, i.e. the first model. 
  !First set up myofibrilar model at all nodes

   CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField, & 
     & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1,2_CMISSIntg,Err)
  !Set up mitochondrial model at IMS nodes
     
    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)>10.OR.NodeNums(i,3)==15) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,1_CMISSIntg,Err)
       ENDIF
     ENDIF
        
  !Set up no ODE model at matrix nodes
 IF(NodeCoords(i,3)==20.AND.NodeNums(i,3)==0) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0_CMISSIntg,Err)
       ENDIF
     ENDIF
    ENDDO

   CALL cmfe_Field_ParameterSetUpdateStart(CellMLModelsField, &
    & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
   CALL cmfe_Field_ParameterSetUpdateFinish(CellMLModelsField, &
    & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)


  WRITE(*,'(A)') "CCCCCCCCCCCCCCCCCCCCCCC."

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML, & 
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)
   WRITE(*,'(A)') "DDDDDDDDDDDDDDD."


  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

  !Start the creation of CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of CellML Intermediate
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)

  WRITE(*,'(A)') "ATPPPPPPPPPPPPPPPPPPPPPPPPPPPP."

  !set initial value of the dependent field/state variable,
  CALL cmfe_Field_ComponentValuesInitialise(ATPField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_ATP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(ADPField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_ADP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(AMPField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_AMP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(PCrField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_PCr,Err)
  CALL cmfe_Field_ComponentValuesInitialise(CrField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_Cr,Err)
  CALL cmfe_Field_ComponentValuesInitialise(PiField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_Pi,Err)
  CALL cmfe_Field_ComponentValuesInitialise(OxyField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_Oxy,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DPsiField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)
 
    DO i = 1,NUMBER_OF_NODES
     IF(NodeCoords(i,3)>10.OR.NodeNums(i,3)==15) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(DPsiField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,init_DPsi,Err)
       ENDIF
     ENDIF
     IF(NodeCoords(i,3)==20.AND.NodeNums(i,3)==0) THEN
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
       IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(ATPField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,8000.0_CMISSDP,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(ADPField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,80.0_CMISSDP,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(PiField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,1800.0_CMISSDP,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(DPsiField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NodeNums(i,1),1,0.0_CMISSDP,Err)
       ENDIF
     ENDIF
    ENDDO

  CALL cmfe_Field_ParameterSetUpdateStart(DPsiField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DPsiField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    

  !Create the equations set equations for ATP
  CALL cmfe_Equations_Initialise(ATPEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(ATPEquationsSet,ATPEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(ATPEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(ATPEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !Set the equations set output
  !CALL cmfe_EquationsOutputTypeSet(ATPEquations,cmfe_EquationsNoOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,cmfe_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(ATPEquationsSet,Err)
  
  !Create the equations set equations for ADP
  CALL cmfe_Equations_Initialise(ADPEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(ADPEquationsSet,ADPEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(ADPEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(ADPEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(ADPEquationsSet,Err)

  !Create the equations set equations for AMP
  CALL cmfe_Equations_Initialise(AMPEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(AMPEquationsSet,AMPEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(AMPEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(AMPEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(AMPEquationsSet,Err)

  !Create the equations set equations for PCr
  CALL cmfe_Equations_Initialise(PCrEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(PCrEquationsSet,PCrEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(PCrEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(PCrEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(PCrEquationsSet,Err)

  !Create the equations set equations for Cr
  CALL cmfe_Equations_Initialise(CrEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(CrEquationsSet,CrEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(CrEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(CrEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(CrEquationsSet,Err)

  !Create the equations set equations for Pi
  CALL cmfe_Equations_Initialise(PiEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(PiEquationsSet,PiEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(PiEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(PiEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(PiEquationsSet,Err)

!Create the equations set equations for Oxy
  CALL cmfe_Equations_Initialise(OxyEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(OxyEquationsSet,OxyEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(OxyEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(OxyEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(OxyEquationsSet,Err)

!Create the equations set equations for DPsi
  CALL cmfe_Equations_Initialise(DPsiEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(DPsiEquationsSet,DPsiEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(DPsiEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(DPsiEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(DPsiEquationsSet,Err)

  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
 CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],Problem,Err)
  !Set the problem to be a strang split reaction diffusion problem
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)


  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,cmfe_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSDP,3000.0_CMISSDP,1.0_CMISSDP,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,100,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Set up the problem solvers for Strang splitting
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_DAESolverTypeSet(Solver,cmfe_SOLVER_DAE_EULER,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,0.01_CMISSDP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_NO_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL cmfe_Solver_DynamicThetaSet(Solver,1.0_CMISSDP,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_TIMING_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LibraryTypeSet(LinearSolver,cmfe_SOLVER_LAPACK_LIBRARY,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,cmfe_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LinearDirectTypeSet(LinearSolver,cmfe_SOLVER_DIRECT_LU,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_Solvercmfe_Library,Err)
  !CALL cmfe_SolverLinearTypeSet(LinearSolver,cmfe_SolverLinearDirectSolveType,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_SolverMUMPSLibrary,Err)
  !CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,10000,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_Solver_No_Output,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)

  !Third solver is another DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_Solver_DAESolverTypeSet(Solver,cmfe_SOLVER_DAE_EULER,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,0.01_CMISSDP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverTimingOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverSolverOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)


  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)


  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquationsSparsityTypeSet(SolverEquations,cmfe_SolverEquationsSparseMatrices,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,cmfe_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,ATPEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,ADPEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,PiEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,AMPEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,PCrEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CrEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,OxyEquationsSet,EquationsSetIndex,Err)
  !CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,DPsiEquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!_________________________________________________________________________________________________________
     !Set up the boundary contions
  WRITE(*,*) 'Set up boundary conditions'  
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CONDITION = cmfe_BOUNDARY_CONDITION_FIXED
  BCVALUE = 0.0_CMISSDP

  DO i = 1,NUMBER_OF_NODES
  IF(NodeNums(i,3)==10) THEN
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
     IF(NodeDomain==ComputationalNodeNumber) THEN
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ADPField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,AMPField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,PCrField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CrField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,PiField, &
     & cmfe_FIELD_DELUDELN_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,BCVALUE,Err)
    ENDIF
    ENDIF
  IF(NodeNums(i,3)==10) THEN
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNums(i,1),1,NodeDomain,Err)
     IF(NodeDomain==ComputationalNodeNumber) THEN
     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,OxyField, &
     & cmfe_FIELD_U_VARIABLE_TYPE, &
     & 1, cmfe_NO_GLOBAL_DERIV, &
     & NodeNums(i,1),1,CONDITION,init_Oxy,Err)
    ENDIF
    ENDIF
    ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!_________________________________________________________________________________________________________
   !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"Diffusion","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"Diffusion","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF

  WRITE(*,'(A)') "Program successfully completed."

  CALL cmfe_Finalise(Err)

  STOP
  END PROGRAM Cardiac_bioenergetics
  
  
  
  
  
  
  
  
  
  















