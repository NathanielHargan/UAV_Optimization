# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024.HF3 replay file
# Internal Version: 2024_05_20-11.43.19 RELr426 191376
# Run by natha on Thu Sep 25 05:12:44 2025
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=116.21265411377, 
    height=163.013122558594)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('UAV_Nathaniel_Hargan.cae')
#: The model database "C:\Desktop\UAV_Mission\UAV_Nathaniel_Hargan.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Drone_Slice']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].setValues(displayedObject=None)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    step='Force Apllied 100N')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, optimizationTasks=OFF, 
    geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['Model-1'].parts['Drone_Slice']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.Job(name='Job-2', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numCpus=1, numGPUs=0)
mdb.jobs['UAV_SLICE_ANGLE_CORRECTION'].submit(consistencyChecking=OFF)
#: The job input file "UAV_SLICE_ANGLE_CORRECTION.inp" has been submitted for analysis.
#: Job UAV_SLICE_ANGLE_CORRECTION: Analysis Input File Processor completed successfully.
#: Job UAV_SLICE_ANGLE_CORRECTION: Abaqus/Standard completed successfully.
#: Job UAV_SLICE_ANGLE_CORRECTION completed successfully. 
session.viewports['Viewport: 1'].setValues(displayedObject=None)
o1 = session.openOdb(name='C:/temp/UAV_SLICE_ANGLE_CORRECTION.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/temp/UAV_SLICE_ANGLE_CORRECTION.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          11
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2992.57, 
    farPlane=4011.83, width=2041.16, height=1167.62, cameraPosition=(412.118, 
    -2632.67, 2258), cameraUpVector=(0.142255, 0.862966, 0.484824), 
    cameraTarget=(508.047, -98.1324, 36.5835))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2932.71, 
    farPlane=4070.51, width=2000.33, height=1144.27, cameraPosition=(279.509, 
    -3189.23, 1378.86), cameraUpVector=(0.145358, 0.678107, 0.720446), 
    cameraTarget=(503.102, -118.885, 3.803))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Mises'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='E', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Max. In-Plane Principal'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='E', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'E11'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Mises'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Max. Principal'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2963.67, 
    farPlane=4046.34, width=2021.45, height=1156.35, cameraPosition=(1255.47, 
    -1410.71, 3068.39), cameraUpVector=(0.704111, 0.668274, -0.240079), 
    cameraTarget=(539.335, -52.8569, 66.5273))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3087.31, 
    farPlane=3917.44, width=2105.78, height=1204.59, cameraPosition=(323.922, 
    620.008, 3392.58), cameraUpVector=(0.95871, -0.0183411, -0.283793), 
    cameraTarget=(503.883, 24.4261, 78.8649))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2741.12, 
    farPlane=4247.76, width=1869.65, height=1069.52, cameraPosition=(-2034.46, 
    -713.212, 2247.61), cameraUpVector=(0.842443, 0.217585, 0.492897), 
    cameraTarget=(415.834, -25.349, 36.1182))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
