# -*- coding: mbcs -*-
# Do not delete the following import lines
# THIS FILE IS IN THE BOUNDARY HOLD DIRECTORY BECAUSE I CAN'T GET RID OF THE STUFF IN MY H DRIVE
from abaqus import *
from abaqusConstants import *
import __main__

def PIDs(name): # 'standard.exe'
    import signal
    import os
    
    output = os.popen('wmic process get name, processid').read()
    PIDsList = []
    for line in output.splitlines():
        if line.startswith(name):
            line = line.replace(name,'')
            line = line.replace(' ','')
            PIDsList.append(line)
    return PIDsList

def test_delete(): # worked
    import os, time
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    JobFiles = ['.com','.dat','.inp','.ipm','.log','.msg','.odb','.prt','.sim']
    ToDelete = 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_7p5_LE_25kPa'
    # delete job files
    for extension in JobFiles:
        while os.path.isfile(Directory+ToDelete+extension) == True:
            #time.sleep(0.1)
            try:
                os.remove(Directory+ToDelete+extension)
            except:
                pass


def AA_Data_Last_Frame():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    # ModelName = 'Comp_attempt - Copy'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    SavePath = 'MacroOutputs/'
    
    #OutputTimes = [1]
    Job = ['ZZ_Change_LE_5p4_P_10p81_ID_1p1_IT_3_IR_1p9_IE_9x_LT_1p86_LE_17500']
    try:
        o3 = session.openOdb(name=Directory+Job[0]+'.odb')
        session.viewports[session.currentViewportName].setValues(displayedObject=o3)
        session.viewports[session.currentViewportName].makeCurrent()
        session.mdbData.summary()
        odb = session.odbs[Directory+Job[0]+'.odb']
    except:
        reply = getWarningReply(message=('Database '+Job[0]+' output processing failed: '), buttons=(YES,NO))
        return
    
    step = odb.steps['Compress']
    num_frames = len(step.frames)
    session.viewports[session.currentViewportName].odbDisplay.setFrame(step=0, frame=num_frames) # set last frame
    
    try:
        # LBU2
        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
        NODE_LABELS=[]
        for i in range(len(NODES)):
            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_LBU2'+'.npy', xy1)
        # LBS2
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_LBS2'+'.npy', xy1)
        # P1U2
        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
        NODE_LABELS=[]
        for i in range(len(NODES)):
            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_P1U2'+'.npy', xy1)
        # P1S2
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_P1S2'+'.npy', xy1)
        # P2U2
        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
        NODE_LABELS=[]
        for i in range(len(NODES)):
            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_P2U2'+'.npy', xy1)
        # P2S2
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_P2S2'+'.npy', xy1)
        # LTU2
        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
        NODE_LABELS=[]
        for i in range(len(NODES)):
            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
        session.viewports[session.currentViewportName].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
        pth = session.paths['Path-1']
        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
            removeDuplicateXYPairs=True, includeAllElements=False)
        numpy.save(Directory+SavePath+Job[0]+'_LTU2'+'.npy', xy1)
    except:
        reply = getWarningReply(message=('Path data recovery failed'), buttons=(YES,NO))

    # close odb
    session.odbs[Directory+Job[0]+'.odb'].close()

def AA_Submit_One_LE_SampleStrainControlled():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # ENTER ALL VALUES AS DOUBLES, NO INTEGERS ALLOWED
    Preload = [0.1081] # [0.025,0.05,0.075,0.1] #,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = '10p81' # ['2p5%','5%','7p5%','10%'] #,'12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = 17500.0 #15974.31 , 16410.0, 6.0
    LayerElasticityN = '17500' 
    BackgroundElasticity = 17500.0
    InclusionElasticity = 17500.0*9.0 # Doesn't converge after 7.5x background
    InclusionElasticityN = '9x'
    LayerThickness = 0.00186
    LayerThicknessN = '1p86'
    SampleThickness = 0.0054 # CANNOT CHANGE WITHOUT CHANGING MESH MANUALLY
    InclusionDepth = 0.0011 # Minimum 60e-06
    InclusionDepthN = '1p1'
    InclusionThickness = 0.003 # Minimum 10e-06, NOT allowed to have inclusion touching bottom of sample, need at least 10e-06 clearance
    InclusionThicknessN = '3'
    InclusionRadius = 0.0019 # Minimum 10e-06
    InclusionRadiusN = '1p9'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Change_LE_5p4'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    JobPrefix = 'ZZ'
    
    # Change Sample
    p = mdb.models[ModelName].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    if SampleThickness-InclusionDepth-InclusionThickness == 0: # If Inclusion touches the bottom of the sample
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
        d[2].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[0].setValues(value=50e-06, ) # Small multiple of mesh size
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    else:
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
        d[2].setValues(value=SampleThickness-InclusionDepth-InclusionThickness, ) # Distance from bottom of sample to bottom of inclusion
        d[0].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    p.features['Datum plane-1'].setValues(offset=0.05) # Change sample hold point
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionElasticity, 0.48), )) # Set Inclusion Material Properties
    mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.48), )) # Set Background Material Properties

    # Change Layer
    p = mdb.models[ModelName].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    d[6].setValues(value=LayerThickness, ) # Layer Thickness
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    if LayerThickness > 0.00125:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.0005, ) # first to second particle layer
    elif LayerThickness > 0.00076:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.00001, ) # first to second particle layer
    else:
        d[0].setValues(value=0.0001, ) # bottom to first particle layer
        d[1].setValues(value=0.0001, ) # first to second particle layer
    if LayerThickness == 0.004:
        d[5].setValues(value=0.0006, )
    else:
        d[5].setValues(value=0.0005, )
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models[ModelName].sketches['__edit__']
    if LayerElasticity >= BackgroundElasticity:
        p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.005)
    else:
        p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.0001) # If layer is soft smaller increments needed
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerElasticity, 0.48), )) # Set Layer Material Properties
    
    # Regenerate Assembly
    a = mdb.models[ModelName].rootAssembly
    a.regenerate()
    # Boundary Conditions, gradient of position with time equation derived in notebook
    # B_t = -Preload[-1]*SampleThickness*LayerElasticity/BackgroundElasticity
    # T_t = -Preload[-1]*(LayerThickness+SampleThickness*LayerElasticity/BackgroundElasticity)
    # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
    B_t = -Preload[-1]*SampleThickness
    #T_t = -Preload[-1]*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-Preload[-1])
    #mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
    mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    # Output Times Calc and Set
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    # Amplitude Tupple
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-numpy.array(Preload))
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*(1-numpy.array(Preload))**(BackgroundElasticity/LayerElasticity) - LayerThickness
    Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*numpy.power((1.0-numpy.array(Preload)),(BackgroundElasticity/LayerElasticity)) - LayerThickness
    AmplitudeTupple = ()
    AmplitudeTupple = AmplitudeTupple + ((0,0,),)
    for i in range(0,len(OutputTimes)):
        AmplitudeTupple = AmplitudeTupple + ((OutputTimes[i],Delta_T[i],),)
    mdb.models[ModelName].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=AmplitudeTupple)
    #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
    OutputTimesTupple = ()
    for i in range(0,len(OutputTimes)):
        OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
    mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        
    # Make Job
    JobName = JobPrefix+'_'+ModelName+'_P_'+PreloadN+'_ID_'+InclusionDepthN+'_IT_'+InclusionThicknessN+'_IR_'+InclusionRadiusN+'_IE_'+InclusionElasticityN+'_LT_'+LayerThicknessN+'_LE_'+LayerElasticityN
    j = mdb.Job(name=JobName, 
        model=ModelName, 
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=50, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)
    j.submit(consistencyChecking=OFF)

def A_Submit_One_Change_Particle_LE_SampleStrainControlled():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # ENTER ALL VALUES AS DOUBLES, NO INTEGERS ALLOWED
    LtoP1 = 0.00015 # min 100 um
    LtoP1N = '0p15'
    P1toP2 = 0.00015 # min 100 um
    P1toP2N = '0p15'
    Preload = [0.10] # [0.025,0.05,0.075,0.1] #,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = '10' # ['2p5%','5%','7p5%','10%'] #,'12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = 15000.0 #15974.31 , 16410.0, 6.0
    LayerElasticityN = '15000' 
    BackgroundElasticity = 15000.0
    InclusionElasticity = 15000.0*10.0 # Doesn't converge after 7.5x background
    InclusionElasticityN = '10x'
    LayerThickness = 0.004
    LayerThicknessN = '4'
    SampleThickness = 0.0054 # CANNOT CHANGE WITHOUT CHANGING MESH MANUALLY
    InclusionDepth = 0.0005 # Minimum 60e-06
    InclusionDepthN = '0p5'
    InclusionThickness = 0.002 # Minimum 10e-06, NOT allowed to have inclusion touching bottom of sample, need at least 10e-06 clearance
    InclusionThicknessN = '2'
    InclusionRadius = 0.001 # Minimum 10e-06
    InclusionRadiusN = '1'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Change_LE_5p4'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    
    # Change Sample
    p = mdb.models[ModelName].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    if SampleThickness-InclusionDepth-InclusionThickness == 0: # If Inclusion touches the bottom of the sample
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
        d[2].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[0].setValues(value=50e-06, ) # Small multiple of mesh size
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    else:
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
        d[2].setValues(value=SampleThickness-InclusionDepth-InclusionThickness, ) # Distance from bottom of sample to bottom of inclusion
        d[0].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    p.features['Datum plane-1'].setValues(offset=0.05) # Change sample hold point
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionElasticity, 0.48), )) # Set Inclusion Material Properties
    mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.48), )) # Set Background Material Properties

    # Change Layer
    p = mdb.models[ModelName].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    d[6].setValues(value=LayerThickness, ) # Layer Thickness
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    if LayerThickness > 0.00125:
        d[0].setValues(value=LtoP1, ) # bottom to first particle layer                 # MESSING WITH THESE, PUT BACK TO 500um each at end
        d[1].setValues(value=P1toP2, ) # first to second particle layer
    elif LayerThickness > 0.00076:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.00001, ) # first to second particle layer
    else:
        d[0].setValues(value=0.0001, ) # bottom to first particle layer
        d[1].setValues(value=0.0001, ) # first to second particle layer
    if LayerThickness == 0.004:
        d[5].setValues(value=0.0006, )
    else:
        d[5].setValues(value=0.0005, )
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models[ModelName].sketches['__edit__']
    if LayerElasticity >= BackgroundElasticity:
        p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.005)
    else:
        p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.0001) # If layer is soft smaller increments needed
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerElasticity, 0.48), )) # Set Layer Material Properties
    
    # Regenerate Assembly
    a = mdb.models[ModelName].rootAssembly
    a.regenerate()
    # Boundary Conditions, gradient of position with time equation derived in notebook
    # B_t = -Preload[-1]*SampleThickness*LayerElasticity/BackgroundElasticity
    # T_t = -Preload[-1]*(LayerThickness+SampleThickness*LayerElasticity/BackgroundElasticity)
    # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
    B_t = -Preload[-1]*SampleThickness
    #T_t = -Preload[-1]*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-Preload[-1])
    #mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
    mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    # Output Times Calc and Set
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    # Amplitude Tupple
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-numpy.array(Preload))
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*(1-numpy.array(Preload))**(BackgroundElasticity/LayerElasticity) - LayerThickness
    Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*numpy.power((1.0-numpy.array(Preload)),(BackgroundElasticity/LayerElasticity)) - LayerThickness
    AmplitudeTupple = ()
    AmplitudeTupple = AmplitudeTupple + ((0,0,),)
    for i in range(0,len(OutputTimes)):
        AmplitudeTupple = AmplitudeTupple + ((OutputTimes[i],Delta_T[i],),)
    mdb.models[ModelName].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=AmplitudeTupple)
    #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
    OutputTimesTupple = ()
    for i in range(0,len(OutputTimes)):
        OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
    mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        
    # Make Job
    JobName = 'W_'+ModelName+'_d1_'+LtoP1N+'_d2_'+P1toP2N+'_P_'+PreloadN+'_D_'+InclusionDepthN+'_T_'+InclusionThicknessN+'_R_'+InclusionRadiusN+'_E_'+InclusionElasticityN+'_L_'+LayerThicknessN+'_LE_'+LayerElasticityN
    j = mdb.Job(name=JobName, 
        model=ModelName, 
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=50, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)
    j.submit(consistencyChecking=OFF)

def A_Submit_Hyper_LE_SampleStrainControlled():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # ENTER ALL VALUES AS DOUBLES, NO INTEGERS ALLOWED
    # Only valid if done elasticity values correctly, only allowed one output preload
    Preload = [0.1081] # Strain in the layer
    PreloadN = '10p81'
    LayerElasticity = 16410.0 # NEEDS TO BE THE EQUIVALENT LINEAR VALUE TO REACH THE SAME DISPLACEMENT AT THE SAME STRAIN
    LayerElasticityN = 'HypMR_DS' 
    BackgroundElasticity = LayerElasticity
    InclusionElasticity = LayerElasticity*6.0 # Doesn't converge after 7.5x background
    InclusionElasticityN = 'HypMR_DS'
    LayerThickness = 0.00105
    LayerThicknessN = '1p1'
    SampleThickness = 0.0054 # CANNOT CHANGE WITHOUT CHANGING MESH MANUALLY
    InclusionDepth = 0.0011 # Minimum 60e-06
    InclusionDepthN = '1p1'
    InclusionThickness = 0.003 # Minimum 10e-06, NOT allowed to have inclusion touching bottom of sample, need at least 10e-06 clearance
    InclusionThicknessN = '3'
    InclusionRadius = 0.0019 # Minimum 10e-06
    InclusionRadiusN = '1p9'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Change_LE_HypE_5p4'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    
    # Change Sample
    p = mdb.models[ModelName].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    if SampleThickness-InclusionDepth-InclusionThickness == 0: # If Inclusion touches the bottom of the sample
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
        d[2].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[0].setValues(value=50e-06, ) # Small multiple of mesh size
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    else:
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
        d[2].setValues(value=SampleThickness-InclusionDepth-InclusionThickness, ) # Distance from bottom of sample to bottom of inclusion
        d[0].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    p.features['Datum plane-1'].setValues(offset=0.05) # Change sample hold point
    p.regenerate()
    p.generateMesh()
    #mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionElasticity, 0.45), )) # Set Inclusion Material Properties
    #mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties

    # Change Layer
    p = mdb.models[ModelName].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    d[6].setValues(value=LayerThickness, ) # Layer Thickness
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    if LayerThickness > 0.00125:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.0005, ) # first to second particle layer
    elif LayerThickness > 0.00076:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.00001, ) # first to second particle layer
    else:
        d[0].setValues(value=0.0001, ) # bottom to first particle layer
        d[1].setValues(value=0.0001, ) # first to second particle layer
    if LayerThickness == 0.004:
        d[5].setValues(value=0.0006, )
    else:
        d[5].setValues(value=0.0005, )
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models[ModelName].sketches['__edit__']
    if LayerElasticity >= BackgroundElasticity:
        p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.005)
    else:
        p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.0001) # If layer is soft smaller increments needed
    p.regenerate()
    p.generateMesh()
    # mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerElasticity, 0.45), )) # Set Layer Material Properties
    
    # Regenerate Assembly
    a = mdb.models[ModelName].rootAssembly
    a.regenerate()
    # Boundary Conditions, gradient of position with time equation derived in notebook
    # B_t = -Preload[-1]*SampleThickness*LayerElasticity/BackgroundElasticity
    # T_t = -Preload[-1]*(LayerThickness+SampleThickness*LayerElasticity/BackgroundElasticity)
    # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
    B_t = -Preload[-1]*SampleThickness
    #T_t = -Preload[-1]*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-Preload[-1])
    #mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
    mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    # Output Times Calc and Set
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    # Amplitude Tupple
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-numpy.array(Preload))
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*(1-numpy.array(Preload))**(BackgroundElasticity/LayerElasticity) - LayerThickness
    Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*numpy.power((1.0-numpy.array(Preload)),(BackgroundElasticity/LayerElasticity)) - LayerThickness
    AmplitudeTupple = ()
    AmplitudeTupple = AmplitudeTupple + ((0,0,),)
    for i in range(0,len(OutputTimes)):
        AmplitudeTupple = AmplitudeTupple + ((OutputTimes[i],Delta_T[i],),)
    mdb.models[ModelName].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=AmplitudeTupple)
    #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
    OutputTimesTupple = ()
    for i in range(0,len(OutputTimes)):
        OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
    mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        
    # Make Job
    JobName = 'V_'+ModelName+'_P_'+PreloadN+'_ID_'+InclusionDepthN+'_IT_'+InclusionThicknessN+'_IR_'+InclusionRadiusN+'_IE_'+InclusionElasticityN+'_LT_'+LayerThicknessN+'_LE_'+LayerElasticityN
    j = mdb.Job(name=JobName, 
        model=ModelName, 
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=50, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)
    j.submit(consistencyChecking=OFF)

def Submit_One_Configuration_SampleStrainControlled():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = 25000.0 # this script can change this
    LayerElasticityN = '25kPa'
    BackgroundElasticity = 25000.0
    InclusionElasticity = 25000.0*5.0 # Doesn't converge after 7.5x background
    InclusionElasticityN = '5x'
    LayerThickness = 0.001
    LayerThicknessN = '1'
    SampleThickness = 0.005
    InclusionDepth = 0.001 # Minimum 60e-06
    InclusionDepthN = '1'
    InclusionThickness = 0.003 # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = '3'
    InclusionRadius = 0.002 # Minimum 10e-06
    InclusionRadiusN = '2'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    
    # Change Sample
    p = mdb.models[ModelName].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    if SampleThickness-InclusionDepth-InclusionThickness == 0: # If Inclusion touches the bottom of the sample
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
        d[2].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[0].setValues(value=50e-06, ) # Small multiple of mesh size
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    else:
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
        d[2].setValues(value=SampleThickness-InclusionDepth-InclusionThickness, ) # Distance from bottom of sample to bottom of inclusion
        d[0].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionElasticity, 0.45), )) # Set Inclusion Material Properties
    mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties

    # Change Layer
    p = mdb.models[ModelName].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    d[6].setValues(value=LayerThickness, ) # Layer Thickness
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    if LayerThickness > 0.00125:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.0005, ) # first to second particle layer
    elif LayerThickness > 0.00076:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.00001, ) # first to second particle layer
    else:
        d[0].setValues(value=0.0001, ) # bottom to first particle layer
        d[1].setValues(value=0.0001, ) # first to second particle layer
    if LayerThickness == 0.004:
        d[5].setValues(value=0.0006, )
    else:
        d[5].setValues(value=0.0005, )
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models[ModelName].sketches['__edit__']
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerElasticity, 0.45), )) # Set Layer Material Properties
    
    # Regenerate Assembly
    a = mdb.models[ModelName].rootAssembly
    a.regenerate()
    # Boundary Conditions, gradient of position with time equation derived in notebook
    # B_t = -Preload[-1]*SampleThickness*LayerElasticity/BackgroundElasticity
    # T_t = -Preload[-1]*(LayerThickness+SampleThickness*LayerElasticity/BackgroundElasticity)
    # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
    B_t = -Preload[-1]*SampleThickness
    T_t = -Preload[-1]*(SampleThickness+LayerThickness*BackgroundElasticity/LayerElasticity)
    mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
    mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    # Output Times Calc and Set
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
    OutputTimesTupple = ()
    for i in range(0,len(OutputTimes)):
        OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
    mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        
    # Make Job
    JobName = 'Submit_One_'+ModelName+'_ID_'+InclusionDepthN+'_IT_'+InclusionThicknessN+'_IR_'+InclusionRadiusN+'_IE_'+InclusionElasticityN+'_LT_'+LayerThicknessN+'_LE_'+LayerElasticityN
    j = mdb.Job(name=JobName, 
        model=ModelName, 
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=50, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)
    j.submit(consistencyChecking=OFF)


def A_Data_Submit_One(): # NOT FINISHED
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = 25000
    LayerElasticityN = '25kPa'
    BackgroundElasticity = 25000
    InclusionElasticity = 50000 # Doesn't converge after 7.5x background
    InclusionElasticityN = '2x'
    LayerThickness = 0.003
    LayerThicknessN = '3'
    SampleThickness = 0.005
    InclusionDepth = 500e-06 # Minimum 60e-06
    InclusionDepthN = '0p5'
    InclusionThickness = 0.003 # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = '3'
    InclusionRadius = 0.003 # Minimum 10e-06
    InclusionRadiusN = '3'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    Job = ['Submit_One_'+ModelName+'_ID_'+InclusionDepth+'_IT_'+InclusionThickness+'_IR_'+InclusionRadius+'_IE_'+InclusionElasticity+'_LT_'+LayerThickness+'_LE_'+LayerElasticity]
    try:
        o3 = session.openOdb(name=Directory+Job[0]+'.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        session.viewports['Viewport: 1'].makeCurrent()
        session.mdbData.summary()
        odb = session.odbs[Directory+Job[0]+'.odb']
    except:
        reply = getWarningReply(message=('Database '+Job[0]+' output processing failed: '), buttons=(YES,NO))
        return

    for p in range(0,len(Preload)):
        try:
            if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                # LBU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                # LBS2
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                # P1U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                # P2U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                # LTU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
            else:
                # make a list of failed times attached to the job and eventually export as txt file or similar
                reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Job[0]+'at time '+str(OutputTimes[p])+' failed due to time mismatch'), buttons=(YES,NO))
        except:
            reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Job[0]+' at time '+str(OutputTimes[p])+' failed due to mon-existence'), buttons=(YES,NO))
        
    # close odb
    session.odbs[Directory+Job[0]+'.odb'].close()
    

def Final_Mesh17_DeleteFilesAfterJob():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FailedFrames = []
    FailedJobs = []
    nJobs = 0 # Counting Variable
    nMaxJobs = 15 # Max number of Jobs running at once, ideally we use one core for each job, but found that it crashes due to RAM, solved with restricted memory percentage per job
    nCoresPerJob = 2
    timeOut = 3*60*60 # Max seconds per job
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['2x','5x','7p5x','10x']
    LayerThickness = [0.0005]
    LayerThicknessN = ['0p5']
    SampleThickness = 0.005
    InclusionDepth = [0.002] # Minimum 60e-06
    InclusionDepthN = ['2']
    InclusionThickness = [500e-06] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['0p5']
    InclusionRadius = [0.003] # Minimum 10e-06
    InclusionRadiusN = ['3']
    #InclusionElasticity = [50000,5*25000,7.5*25000,10*25000] # Doesn't converge after 7.5x background
    #InclusionElasticityN = ['2x','5x','7p5x','10x']
    #LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    #LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    #SampleThickness = 0.005
    #InclusionDepth = [60e-06,500e-06,0.001,0.002] # Minimum 60e-06
    #InclusionDepthN = ['0p05','0p5','1','2']
    #InclusionThickness = [500e-06,0.001,0.002,0.003] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    #InclusionThicknessN = ['0p5','1','2','3']
    #InclusionRadius = [500e-06,0.001,0.002,0.003] # Minimum 10e-06
    #InclusionRadiusN = ['0p5','1','2','3']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    # ModelName+'_'+'ID'+'_'+InclusionDepthN[id]+'IT'+'_'+InclusionThicknessN[it]+'IR'+'_'+InclusionRadiusN[ir]+
    #          'IE'+'_'+InclusionElasticityN[ie]+'_'+'LT'+'_'+LayerThicknessN[lt]+'_'+'LE'+'_'+LayerElasticityN[le]+'_'+'P'+'_'PreloadN[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Jobs = []
    ActiveJobs = []
    CompletedJobs = [] # ['Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_0p5_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_1_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_2_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p05_IT_3_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_0p5_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_1_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_2_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_0p5_IT_3_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_0p5_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_1_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_2_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_1_IT_3_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_0p5_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_1_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_2_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_2_IR_3_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_0p5_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_7p5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_1_IE_10x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_2x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_4_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_3_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_2_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_1_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_5x_LT_0p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_7p5x_LT_10_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_7p5x_LT_7p5_LE_25kPa', 'Mesh17-Final_ID_2_IT_3_IR_2_IE_7p5x_LT_5_LE_25kPa']
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    JobFiles = ['.com','.dat','.inp','.ipm','.log','.msg','.odb','.prt','.sim'] # for some reason it can't delete the ODB
    # Alter and Submit Model
    for n1 in range(0,len(Samples)): # Can control where to start from with n1
        #if (n1 < 20 and ): continue # Skips iterations where conditions are satisfied
        # Change Sample
        p = mdb.models[ModelName].parts['Sample']
        s = p.features['Partition face-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        if SampleThickness-Samples[n1][0]-Samples[n1][1] == 0: # If Inclusion touches the bottom of the sample
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            d[2].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[0].setValues(value=50e-06, ) # Small multiple of mesh size
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            d[2].setValues(value=SampleThickness-Samples[n1][0]-Samples[n1][1], ) # Distance from bottom of sample to bottom of inclusion
            d[0].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        p.features['Partition face-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((Samples[n1][3], 0.45), )) # Set Inclusion Material Properties
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties
        for n2 in range(0,len(Layers)):
            tempJobName = ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1]
            if ( CompletedJobs.count(tempJobName) > 0 ): continue # Skips iterations where conditions are satisfied
            # Change Layer
            p = mdb.models[ModelName].parts['Layer']
            s = p.features['Shell planar-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s1 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
            d[6].setValues(value=Layers[n2][0], ) # Layer Thickness
            p.features['Shell planar-1'].setValues(sketch=s1)
            del mdb.models[ModelName].sketches['__edit__']
            s = p.features['Partition face-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s2 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
            if Layers[n2][0] > 0.00125:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.0005, ) # first to second particle layer
            elif Layers[n2][0] > 0.00076:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.00001, ) # first to second particle layer
            else:
                d[0].setValues(value=0.0001, ) # bottom to first particle layer
                d[1].setValues(value=0.0001, ) # first to second particle layer
            if Layers[n2][0] == 0.004:
                d[5].setValues(value=0.0006, )
            else:
                d[5].setValues(value=0.0005, )
            p.features['Partition face-1'].setValues(sketch=s2)
            del mdb.models[ModelName].sketches['__edit__']
            p.regenerate()
            p.generateMesh()
            mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((Layers[n2][1], 0.45), )) # Set Layer Material Properties
            
            # Regenerate Assembly
            a = mdb.models[ModelName].rootAssembly
            a.regenerate()
            # Boundary Conditions, gradient of position with time equation derived in notebook
            B_t = -Preload[-1]*SampleThickness*Layers[n2][1]/BackgroundElasticity
            T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness*Layers[n2][1]/BackgroundElasticity)
            # OLD, calcs preload as total % compression of both sample and layer #T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)
            # OLD, calcs preload as total % compression of both sample and layer #B_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)*(1-Layers[n2][0]*BackgroundElasticity/(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1]))
            mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
            mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            # Output Times Calc and Set
            OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
            #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
            OutputTimesTupple = ()
            for i in range(0,len(OutputTimes)):
                OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
            mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
            mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
            # Make Job
            PrevPIDs = PIDs('standard.exe')
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
            j = mdb.Job(name=Jobs[-1], 
                model=ModelName, 
                description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
                queue=None, memory=5, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCoresPerJob, 
                numDomains=nCoresPerJob, numGPUs=0)
            j.submit(consistencyChecking=OFF)
            submitting = 1
            LoopTime = 0
            while submitting == 1:
                time.sleep(0.1)
                LoopTime += 0.1
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                NewPIDs = PIDs('standard.exe')
                temp = list(set(NewPIDs).difference(PrevPIDs))
                if len(temp) == 1:
                    submitting = 0
                    ActiveJobs.append([Jobs[-1],LoopTime,temp[0]])
                    nJobs = nJobs + 1
                elif LoopTime > 180:
                    submitting = 0
            #reply = getWarningReply(message=(temp[0]), buttons=(YES,NO))
            
            # have up to nMaxJobs running at the same time, submitting waiting for running before next submission
            while nJobs >= nMaxJobs:
                time.sleep(0.1)
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                ToRemove = []
                ToKill = []
                for Job in ActiveJobs:
                    if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                        ToRemove.append(Job)
                    if Job[1] > timeOut:
                        ToKill.append(Job)
                        ToRemove.append(Job)
                for Job in ToKill:
                    os.kill(int(Job[2]), signal.SIGINT)#kill
                for Job in ToRemove:
                    ActiveJobs.remove(Job)
                    nJobs = nJobs - 1
                    try:
                        o3 = session.openOdb(name=Directory+Job[0]+'.odb')
                        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                        session.viewports['Viewport: 1'].makeCurrent()
                        session.mdbData.summary()
                        odb = session.odbs[Directory+Job[0]+'.odb']
                    except:
                        #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
                        FailedJobs.append(Job[0])
                        continue
                    for p in range(0,len(Preload)):
                        try:
                            if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                                #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                                # LBU2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                                # LBS2
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                                # P1U2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                                # P2U2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                                # LTU2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                            else:
                                # make a list of failed times attached to the job and eventually export as txt file or similar
                                FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
                        except:
                            #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                            FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
                    
                    # close odb
                    session.odbs[Directory+Job[0]+'.odb'].close()
                    # delete job files
                    for extension in JobFiles:
                        while os.path.isfile(Directory+Job[0]+extension) == True:
                            #time.sleep(0.1)
                            try:
                                os.remove(Directory+Job[0]+extension)
                            except:
                                pass
                    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
                    numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
                    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))
                
            
        
    
    # Wait for last jobs to be completed
    while nJobs > 0:
        time.sleep(0.1)
        for i,Job in enumerate(ActiveJobs):
            ActiveJobs[i][1] += 0.1 # increment time counter for remaining jobs
        ToRemove = []
        ToKill = []
        for Job in ActiveJobs:
            if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                ToRemove.append(Job)
            if Job[1] > timeOut:
                ToKill.append(Job)
                ToRemove.append(Job)
        for Job in ToKill:
            os.kill(int(Job[2]), signal.SIGINT) # kill
        for Job in ToRemove:
            ActiveJobs.remove(Job)
            nJobs = nJobs - 1
            try:
                o3 = session.openOdb(name=Directory+Job[0]+'.odb')
                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                session.viewports['Viewport: 1'].makeCurrent()
                session.mdbData.summary()
                odb = session.odbs[Directory+Job[0]+'.odb']
            except:
                #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
                FailedJobs.append(Job[0])
                continue
            for p in range(0,len(Preload)):
                try:
                    if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                        session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                        #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                        # LBU2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                        # LBS2
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                        # P1U2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                        # P2U2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                        # LTU2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                    else:
                        # make a list of failed times attached to the job and eventually export as txt file or similar
                        FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
                except:
                    #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                    FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
                
            # close odb
            session.odbs[Directory+Job[0]+'.odb'].close()
            # delete job files
            for extension in JobFiles:
                while os.path.isfile(Directory+Job[0]+extension) == True:
                    #time.sleep(0.1)
                    try:
                        os.remove(Directory+Job[0]+extension)
                    except:
                        pass
            numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
            numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
            numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))
        
    
    
def A_Final_ChangeLE_SampleStrainControlled_DeleteFilesAfterJob():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # ENTER ALL MATH INTENDED VALUES AS DOUBLES INSTEAD OF INTEGERS ie 1.0 instead of 1
    FailedFrames = []
    FailedJobs = []
    nJobs = 0 # Counting Variable
    nMaxJobs = 10 # Max number of Jobs running at once, ideally we use one core for each job, but found that it crashes due to RAM, solved with restricted memory percentage per job
    nCoresPerJob = 2
    timeOut = 3*60*60 # Max seconds per job
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2] #,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    #LayerElasticity = [0.1*25000.0,0.2*25000.0,0.5*25000.0,2.0*25000.0,5.0*25000.0,10.0*25000.0]
    #LayerElasticityN = ['0p1x','0p2x','0p5x','2x','5x','10x']
    LayerElasticity = [0.2*25000.0,0.5*25000.0,0.666666666666*25000.0,1.5*25000.0,2.0*25000.0,5.0*25000.0]
    LayerElasticityN = ['0p2x','0p5x','0p66x','1p5x','2x','5x']
    BackgroundElasticity = 25000.0
    InclusionElasticity = [10.0*25000.0] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['10x']
    LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = [0.001] # Minimum 60e-06
    InclusionDepthN = ['1']
    InclusionThickness = [0.003] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['3']
    InclusionRadius = [0.003] # Minimum 10e-06
    InclusionRadiusN = ['3']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    # ModelName+'_'+'ID'+'_'+InclusionDepthN[id]+'IT'+'_'+InclusionThicknessN[it]+'IR'+'_'+InclusionRadiusN[ir]+
    #          'IE'+'_'+InclusionElasticityN[ie]+'_'+'LT'+'_'+LayerThicknessN[lt]+'_'+'LE'+'_'+LayerElasticityN[le]+'_'+'P'+'_'PreloadN[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Change_LE'
    Jobs = []
    ActiveJobs = []
    CompletedJobs = [] # ['Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_0p2x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_0p5x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_2x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_5x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_10_LE_10x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_7p5_LE_0p1x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_7p5_LE_0p2x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_7p5_LE_0p5x', 'Change_LE_ID_0p5_IT_3_IR_2_IE_10x_LT_7p5_LE_5x']
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    #JobFiles = ['.com','.dat','.inp','.ipm','.log','.msg','.odb','.prt','.sim'] # to delete
    JobFiles = ['.com','.dat','.inp','.ipm','.log','.odb','.prt','.sim'] # to delete but keep message file
    # Alter and Submit Model
    for n1 in range(0,len(Samples)): # Can control where to start from with n1
        #if (n1 < 20 and ): continue # Skips iterations where conditions are satisfied
        # Change Sample
        p = mdb.models[ModelName].parts['Sample']
        s = p.features['Partition face-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        if SampleThickness-Samples[n1][0]-Samples[n1][1] == 0: # If Inclusion touches the bottom of the sample
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            d[2].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[0].setValues(value=50e-06, ) # Small multiple of mesh size
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            d[2].setValues(value=SampleThickness-Samples[n1][0]-Samples[n1][1], ) # Distance from bottom of sample to bottom of inclusion
            d[0].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        p.features['Partition face-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.features['Datum plane-1'].setValues(offset=0.05) # Change sample hold point
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((Samples[n1][3], 0.45), )) # Set Inclusion Material Properties
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties
        for n2 in range(0,len(Layers)):
            tempJobName = ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1]
            if ( CompletedJobs.count(tempJobName) > 0 ): continue # Skips iterations where conditions are satisfied
            # Change Layer
            p = mdb.models[ModelName].parts['Layer']
            s = p.features['Shell planar-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s1 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
            d[6].setValues(value=Layers[n2][0], ) # Layer Thickness
            p.features['Shell planar-1'].setValues(sketch=s1)
            del mdb.models[ModelName].sketches['__edit__']
            s = p.features['Partition face-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s2 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
            if Layers[n2][0] > 0.00125:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.0005, ) # first to second particle layer
            elif Layers[n2][0] > 0.00076:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.00001, ) # first to second particle layer
            else:
                d[0].setValues(value=0.0001, ) # bottom to first particle layer
                d[1].setValues(value=0.0001, ) # first to second particle layer
            if Layers[n2][0] == 0.004:
                d[5].setValues(value=0.0006, )
            else:
                d[5].setValues(value=0.0005, )
            p.features['Partition face-1'].setValues(sketch=s2)
            del mdb.models[ModelName].sketches['__edit__']
            if Layers[n2][1] >= BackgroundElasticity:
                p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
            else:
                p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
            p.regenerate()
            p.generateMesh()
            mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((Layers[n2][1], 0.45), )) # Set Layer Material Properties
            if Layers[n2][1] >= BackgroundElasticity:
                p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
                mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.005)
            else:
                p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
                mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.0001) # If layer is soft smaller increments needed
            
            # Regenerate Assembly 
            a = mdb.models[ModelName].rootAssembly
            a.regenerate()
            # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
            B_t = -Preload[-1]*SampleThickness
            #T_t = -Preload[-1]*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-Preload[-1])
            #mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
            mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            # Output Times Calc and Set
            OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
            # Amplitude Tupple
            #Delta_T = -numpy.array(Preload)*SampleThickness + Layers[n2][0]*BackgroundElasticity/Layers[n2][1]*numpy.log(1-numpy.array(Preload))
            Delta_T = -numpy.array(Preload)*SampleThickness + Layers[n2][0]*numpy.power((1.0-numpy.array(Preload)),(BackgroundElasticity/Layers[n2][1])) - Layers[n2][0]
            AmplitudeTupple = ()
            AmplitudeTupple = AmplitudeTupple + ((0,0,),)
            for i in range(0,len(OutputTimes)):
                AmplitudeTupple = AmplitudeTupple + ((OutputTimes[i],Delta_T[i],),)
            mdb.models[ModelName].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=AmplitudeTupple)
            #OutputTimes tupple transform
            OutputTimesTupple = ()
            for i in range(0,len(OutputTimes)):
                OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
            mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
            mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
            
            # Make Job
            PrevPIDs = PIDs('standard.exe')
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
            j = mdb.Job(name=Jobs[-1], 
                model=ModelName, 
                description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
                queue=None, memory=5, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCoresPerJob, 
                numDomains=nCoresPerJob, numGPUs=0)
            j.submit(consistencyChecking=OFF)
            submitting = 1
            LoopTime = 0
            while submitting == 1:
                time.sleep(0.1)
                LoopTime += 0.1
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                NewPIDs = PIDs('standard.exe')
                temp = list(set(NewPIDs).difference(PrevPIDs))
                if len(temp) == 1:
                    submitting = 0
                    ActiveJobs.append([Jobs[-1],LoopTime,temp[0]])
                    nJobs = nJobs + 1
                elif LoopTime > 180:
                    submitting = 0
            #reply = getWarningReply(message=(temp[0]), buttons=(YES,NO))
            
            # have up to nMaxJobs running at the same time, submitting waiting for running before next submission
            while nJobs >= nMaxJobs:
                time.sleep(0.1)
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                ToRemove = []
                ToKill = []
                for Job in ActiveJobs:
                    if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                        ToRemove.append(Job)
                    if Job[1] > timeOut:
                        ToKill.append(Job)
                        ToRemove.append(Job)
                for Job in ToKill:
                    os.kill(int(Job[2]), signal.SIGINT)#kill
                for Job in ToRemove:
                    ActiveJobs.remove(Job)
                    nJobs = nJobs - 1
                    try:
                        o3 = session.openOdb(name=Directory+Job[0]+'.odb')
                        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                        session.viewports['Viewport: 1'].makeCurrent()
                        session.mdbData.summary()
                        odb = session.odbs[Directory+Job[0]+'.odb']
                    except:
                        #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
                        FailedJobs.append(Job[0])
                        continue
                    for p in range(0,len(Preload)):
                        try:
                            if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                                #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                                # LBU2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputsCurrentSample/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                                # LBS2
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputsCurrentSample/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                                # P1U2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputsCurrentSample/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                                # P2U2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputsCurrentSample/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                                # LTU2
                                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                                NODE_LABELS=[]
                                for i in range(len(NODES)):
                                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                                pth = session.paths['Path-1']
                                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                                    removeDuplicateXYPairs=True, includeAllElements=False)
                                numpy.save(Directory+'MacroOutputsCurrentSample/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                            else:
                                # make a list of failed times attached to the job and eventually export as txt file or similar
                                FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
                        except:
                            #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                            FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
                    
                    # close odb
                    session.odbs[Directory+Job[0]+'.odb'].close()
                    # delete job files
                    for extension in JobFiles:
                        while os.path.isfile(Directory+Job[0]+extension) == True:
                            #time.sleep(0.1)
                            try:
                                os.remove(Directory+Job[0]+extension)
                            except:
                                pass
                    numpy.save(Directory+'MacroOutputsCurrentSample/'+'JobList'+'.npy', numpy.array(Jobs))
                    numpy.save(Directory+'MacroOutputsCurrentSample/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
                    numpy.save(Directory+'MacroOutputsCurrentSample/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))
                
            
        
    
    # Wait for last jobs to be completed
    while nJobs > 0:
        time.sleep(0.1)
        for i,Job in enumerate(ActiveJobs):
            ActiveJobs[i][1] += 0.1 # increment time counter for remaining jobs
        ToRemove = []
        ToKill = []
        for Job in ActiveJobs:
            if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                ToRemove.append(Job)
            if Job[1] > timeOut:
                ToKill.append(Job)
                ToRemove.append(Job)
        for Job in ToKill:
            os.kill(int(Job[2]), signal.SIGINT) # kill
        for Job in ToRemove:
            ActiveJobs.remove(Job)
            nJobs = nJobs - 1
            try:
                o3 = session.openOdb(name=Directory+Job[0]+'.odb')
                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                session.viewports['Viewport: 1'].makeCurrent()
                session.mdbData.summary()
                odb = session.odbs[Directory+Job[0]+'.odb']
            except:
                #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
                FailedJobs.append(Job[0])
                continue
            for p in range(0,len(Preload)):
                try:
                    if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                        session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                        #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                        # LBU2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                        # LBS2
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                        # P1U2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                        # P2U2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                        # LTU2
                        NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                        NODE_LABELS=[]
                        for i in range(len(NODES)):
                            NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                        session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                        pth = session.paths['Path-1']
                        xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                            projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                            removeDuplicateXYPairs=True, includeAllElements=False)
                        numpy.save(Directory+'MacroOutputs/'+Job[0]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                    else:
                        # make a list of failed times attached to the job and eventually export as txt file or similar
                        FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
                except:
                    #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                    FailedFrames.append('Job '+Job[0]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
                
            # close odb
            session.odbs[Directory+Job[0]+'.odb'].close()
            # delete job files
            for extension in JobFiles:
                while os.path.isfile(Directory+Job[0]+extension) == True:
                    #time.sleep(0.1)
                    try:
                        os.remove(Directory+Job[0]+extension)
                    except:
                        pass
            numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
            numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
            numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))
        
    
 
def Final_MissingFrames():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    from abaqus import getWarningReply, YES, NO
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000,5*25000,7.5*25000,10*25000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['2x','5x','7p5x','10x']
    LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = [60e-06,500e-06,0.001,0.002] # Minimum 60e-06
    InclusionDepthN = ['0p05','0p5','1','2']
    InclusionThickness = [500e-06,0.001,0.002,0.003] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['0p5','1','2','3']
    InclusionRadius = [500e-06,0.001,0.002,0.003] # Minimum 10e-06
    InclusionRadiusN = ['0p5','1','2','3']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # List of files to analyse
    Jobs = []
    for n1 in range(0,len(Samples)):
        for n2 in range(0,len(Layers)):
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
    # Output Data
    FailedFrames = []
    FailedJobs = []
    for n in range(0,len(Jobs)):
        try:
            if n > 0:
                session.odbs[Directory+Jobs[n-1]+'.odb'].close()
        except:
            pass
        try:
            o3 = session.openOdb(name=Directory+Jobs[n]+'.odb')
            session.mdbData.summary()
            odb = session.odbs[Directory+Jobs[n]+'.odb']
        except:
            #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
            FailedJobs.append(Jobs[n])
            continue
        for p in range(0,len(Preload)):
            try:
                if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                    pass
                else:
                    # make a list of failed times attached to the job and eventually export as txt file or similar
                    FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
            except:
                #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
            
        
    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))

def Final_JustData():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    from abaqus import getWarningReply, YES, NO
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000,5*25000,7.5*25000,10*25000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['2x','5x','7p5x','10x']
    LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = [60e-06,500e-06,0.001,0.002] # Minimum 60e-06
    InclusionDepthN = ['0p05','0p5','1','2']
    InclusionThickness = [500e-06,0.001,0.002,0.003] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['0p5','1','2','3']
    InclusionRadius = [500e-06,0.001,0.002,0.003] # Minimum 10e-06
    InclusionRadiusN = ['0p5','1','2','3']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # List of files to analyse
    Jobs = []
    for n1 in range(0,len(Samples)):
        for n2 in range(0,len(Layers)):
            if ( ((n1 < 20) or (n1 == 20 and n2 < 5)) and (Layers[n2][0] != 0.004) ) or n1 > 30: continue # Skips iterations where conditions are satisfied
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
    # Output Data
    FailedFrames = []
    FailedJobs = []
    for n in range(0,len(Jobs)):
        try:
            if n > 0:
                session.odbs[Directory+Jobs[n-1]+'.odb'].close()
        except:
            pass
        try:
            o3 = session.openOdb(name=Directory+Jobs[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            odb = session.odbs[Directory+Jobs[n]+'.odb']
        except:
            #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
            FailedJobs.append(Jobs[n])
            continue
        for p in range(0,len(Preload)):
            try:
                if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                    #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                    # LBU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                    # LBS2
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                    # P1U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                    # P2U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                    # LTU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                else:
                    # make a list of failed times attached to the job and eventually export as txt file or similar
                    FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
            except:
                #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
            
        
    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))

def Mesh17_BoundaryHold_UseSysKill():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nJobs = 0 # Counting Variable
    nMaxJobs = 16 # Max number of Jobs running at once, ideally we use one core for each job, but found that it crashes due to RAM, solved with restricted memory percentage per job
    nCoresPerJob = 2
    timeOut = 3*60*60 # Max seconds per job
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000,5*25000,7.5*25000,10*25000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['2x','5x','7p5x','10x']
    LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = [60e-06,500e-06,0.001,0.002] # Minimum 60e-06
    InclusionDepthN = ['0p05','0p5','1','2']
    InclusionThickness = [500e-06,0.001,0.002,0.003] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['0p5','1','2','3']
    InclusionRadius = [500e-06,0.001,0.002,0.003] # Minimum 10e-06
    InclusionRadiusN = ['0p5','1','2','3']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    # ModelName+'_'+'ID'+'_'+InclusionDepthN[id]+'IT'+'_'+InclusionThicknessN[it]+'IR'+'_'+InclusionRadiusN[ir]+
    #          'IE'+'_'+InclusionElasticityN[ie]+'_'+'LT'+'_'+LayerThicknessN[lt]+'_'+'LE'+'_'+LayerElasticityN[le]+'_'+'P'+'_'PreloadN[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Jobs = []
    ActiveJobs = []
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # Alter and Submit Model
    for n1 in range(0,len(Samples)): # Can control where to start from with n1
        #if (n1 < 20 and ): continue # Skips iterations where conditions are satisfied
        # Change Sample
        p = mdb.models[ModelName].parts['Sample']
        s = p.features['Partition face-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        if SampleThickness-Samples[n1][0]-Samples[n1][1] == 0: # If Inclusion touches the bottom of the sample
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            d[2].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[0].setValues(value=50e-06, ) # Small multiple of mesh size
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            d[2].setValues(value=SampleThickness-Samples[n1][0]-Samples[n1][1], ) # Distance from bottom of sample to bottom of inclusion
            d[0].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        p.features['Partition face-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((Samples[n1][3], 0.45), )) # Set Inclusion Material Properties
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties
        for n2 in range(0,len(Layers)):
            if ( ((n1 < 20) or (n1 == 20 and n2 < 5)) and (Layers[n2][0] != 0.004) ): continue # Skips iterations where conditions are satisfied
            # Change Layer
            p = mdb.models[ModelName].parts['Layer']
            s = p.features['Shell planar-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s1 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
            d[6].setValues(value=Layers[n2][0], ) # Layer Thickness
            p.features['Shell planar-1'].setValues(sketch=s1)
            del mdb.models[ModelName].sketches['__edit__']
            s = p.features['Partition face-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s2 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
            if Layers[n2][0] > 0.00125:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.0005, ) # first to second particle layer
            elif Layers[n2][0] > 0.00076:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.00001, ) # first to second particle layer
            else:
                d[0].setValues(value=0.0001, ) # bottom to first particle layer
                d[1].setValues(value=0.0001, ) # first to second particle layer
            if Layers[n2][0] == 0.004:
                d[5].setValues(value=0.0006, )
            else:
                d[5].setValues(value=0.0005, )
            p.features['Partition face-1'].setValues(sketch=s2)
            del mdb.models[ModelName].sketches['__edit__']
            p.regenerate()
            p.generateMesh()
            mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((Layers[n2][1], 0.45), )) # Set Layer Material Properties
            
            # Regenerate Assembly
            a = mdb.models[ModelName].rootAssembly
            a.regenerate()
            # Boundary Conditions, gradient of position with time equation derived in notebook
            B_t = -Preload[-1]*SampleThickness*Layers[n2][1]/BackgroundElasticity
            T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness*Layers[n2][1]/BackgroundElasticity)
            # OLD, calcs preload as total % compression of both sample and layer #T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)
            # OLD, calcs preload as total % compression of both sample and layer #B_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)*(1-Layers[n2][0]*BackgroundElasticity/(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1]))
            mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
            mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            # Output Times Calc and Set
            OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
            #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
            OutputTimesTupple = ()
            for i in range(0,len(OutputTimes)):
                OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
            mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
            mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
            # Make Job
            PrevPIDs = PIDs('standard.exe')
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
            j = mdb.Job(name=Jobs[-1], 
                model=ModelName, 
                description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
                queue=None, memory=5, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCoresPerJob, 
                numDomains=nCoresPerJob, numGPUs=0)
            j.submit(consistencyChecking=OFF)
            submitting = 1
            LoopTime = 0
            while submitting == 1:
                time.sleep(0.1)
                LoopTime += 0.1
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                NewPIDs = PIDs('standard.exe')
                temp = list(set(NewPIDs).difference(PrevPIDs))
                if len(temp) == 1:
                    submitting = 0
                    ActiveJobs.append([Jobs[-1],LoopTime,temp[0]])
                    nJobs = nJobs + 1
                elif LoopTime > 180:
                    submitting = 0
            #reply = getWarningReply(message=(temp[0]), buttons=(YES,NO))
            
            # have up to nMaxJobs running at the same time, submitting waiting for running before next submission
            while nJobs >= nMaxJobs:
                time.sleep(0.1)
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 0.1 # increment time counter for active jobs
                ToRemove = []
                ToKill = []
                for Job in ActiveJobs:
                    if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                        ToRemove.append(Job)
                    if Job[1] > timeOut:
                        ToKill.append(Job)
                        ToRemove.append(Job)
                for Job in ToKill:
                    os.kill(int(Job[2]), signal.SIGINT)#kill
                for Job in ToRemove:
                    ActiveJobs.remove(Job)
                    nJobs = nJobs - 1
                
            
        
    
    # Wait for last jobs to be completed
    while nJobs > 0:
        time.sleep(0.1)
        for i,Job in enumerate(ActiveJobs):
            ActiveJobs[i][1] += 0.1 # increment time counter for remaining jobs
        ToRemove = []
        ToKill = []
        for Job in ActiveJobs:
            if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                ToRemove.append(Job)
            if Job[1] > timeOut:
                ToKill.append(Job)
                ToRemove.append(Job)
        for Job in ToKill:
            os.kill(int(Job[2]), signal.SIGINT) # kill
        for Job in ToRemove:
            ActiveJobs.remove(Job)
            nJobs = nJobs - 1
        
    
    # Output Data
    FailedFrames = []
    FailedJobs = []
    for n in range(0,len(Jobs)):
        try:
            if n > 0:
                session.odbs[Directory+Jobs[n-1]+'.odb'].close()
        except:
            pass
        try:
            o3 = session.openOdb(name=Directory+Jobs[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            odb = session.odbs[Directory+Jobs[n]+'.odb']
        except:
            #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
            FailedJobs.append(Jobs[n])
            continue
        for p in range(0,len(Preload)):
            try:
                if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                    #reply = getWarningReply(message=('test'), buttons=(YES,NO))
                    # LBU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                    # LBS2
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                    # P1U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                    # P2U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                    # LTU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                else:
                    # make a list of failed times attached to the job and eventually export as txt file or similar
                    FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
            except:
                #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                FailedFrames.append('Job '+Jobs[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
            
        
    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedJobs'+'.npy', numpy.array(FailedJobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))


def OLD_Mesh17_BoundaryHold_NoTimeOut():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    from abaqus import getWarningReply, YES, NO
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nJobs = 0 # Counting Variable
    nMaxJobs = 8 # Max number of Jobs running at once
    timeOut = 7500 # Max seconds per job
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.025,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000,5*25000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['2x','5x']
    LayerThickness = [0.0005,0.001,0.002,0.003,0.004,0.005,0.0075,0.01]
    LayerThicknessN = ['0p5','1','2','3','4','5','7p5','10']
    SampleThickness = 0.005
    InclusionDepth = [500e-06,0.001,0.002] # Minimum 50e-06
    InclusionDepthN = ['0p5','1','2']
    InclusionThickness = [500e-06,0.001,0.002] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['0p5','1','2']
    InclusionRadius = [500e-06,0.001,0.002] # Minimum 10e-06
    InclusionRadiusN = ['0p5','1','2']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    # ModelName+'_'+'ID'+'_'+InclusionDepthN[id]+'IT'+'_'+InclusionThicknessN[it]+'IR'+'_'+InclusionRadiusN[ir]+
    #          'IE'+'_'+InclusionElasticityN[ie]+'_'+'LT'+'_'+LayerThicknessN[lt]+'_'+'LE'+'_'+LayerElasticityN[le]+'_'+'P'+'_'PreloadN[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Jobs = []
    ActiveJobs = []
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # Alter and Submit Model
    for n1 in range(0,len(Samples)):
        # Change Sample
        p = mdb.models[ModelName].parts['Sample']
        s = p.features['Partition face-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        if SampleThickness-Samples[n1][0]-Samples[n1][1] == 0: # If Inclusion touches the bottom of the sample
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            d[2].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[0].setValues(value=50e-06, ) # Small multiple of mesh size
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            d[2].setValues(value=SampleThickness-Samples[n1][0]-Samples[n1][1], ) # Distance from bottom of sample to bottom of inclusion
            d[0].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        p.features['Partition face-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((Samples[n1][3], 0.45), )) # Set Inclusion Material Properties
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties
        for n2 in range(0,len(Layers)):
            # Change Layer
            p = mdb.models[ModelName].parts['Layer']
            s = p.features['Shell planar-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s1 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
            d[6].setValues(value=Layers[n2][0], ) # Layer Thickness
            p.features['Shell planar-1'].setValues(sketch=s1)
            del mdb.models[ModelName].sketches['__edit__']
            s = p.features['Partition face-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s2 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
            if Layers[n2][0] > 0.00125:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.0005, ) # first to second particle layer
            elif Layers[n2][0] > 0.00076:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.00001, ) # first to second particle layer
            else:
                d[0].setValues(value=0.0001, ) # bottom to first particle layer
                d[1].setValues(value=0.0001, ) # first to second particle layer
            p.features['Partition face-1'].setValues(sketch=s2)
            del mdb.models[ModelName].sketches['__edit__']
            p.regenerate()
            p.generateMesh()
            mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((Layers[n2][1], 0.45), )) # Set Layer Material Properties
            
            # Regenerate Assembly
            a = mdb.models[ModelName].rootAssembly
            a.regenerate()
            # Boundary Conditions, gradient of position with time equation derived in notebook
            B_t = -Preload[-1]*SampleThickness*Layers[n2][1]/BackgroundElasticity
            T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness*Layers[n2][1]/BackgroundElasticity)
            # OLD, calcs preload as total % compression of both sample and layer #T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)
            # OLD, calcs preload as total % compression of both sample and layer #B_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)*(1-Layers[n2][0]*BackgroundElasticity/(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1]))
            mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
            mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            # Output Times Calc and Set
            OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
            #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
            OutputTimesTupple = ()
            for i in range(0,len(OutputTimes)):
                OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
            mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
            mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
            # Make Job
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
            ActiveJobs.append([Jobs[-1],20])
            j = mdb.Job(name=Jobs[-1], 
                model=ModelName, 
                description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
                queue=None, memory=90, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
                numDomains=4, numGPUs=0)
            j.submit(consistencyChecking=OFF)
            nJobs = nJobs + 1
            
            time.sleep(60)
            #reply = getWarningReply(message=('got past first wait'), buttons=(YES,NO))
            # have up to nMaxJobs running at the same time, submitting with a minute break in between
            while nJobs >= nMaxJobs:
                time.sleep(60)
                ToRemove = []
                ToKill = []
                for Job in ActiveJobs:
                    if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                        ToRemove.append(Job)
                for Job in ToRemove:
                    ActiveJobs.remove(Job)
                    nJobs = nJobs - 1
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 10 # increment time counter for remaining jobs
                
            
        
    
    # Wait for last jobs to be completed
    while nJobs > 0:
        time.sleep(60)
        ToRemove = []
        ToKill = []
        for Job in ActiveJobs:
            if os.path.isfile(Job[0]+'.cid') == False:
                ToRemove.append(Job)
        for Job in ToRemove:
            ActiveJobs.remove(Job)
            nJobs = nJobs - 1
        for i in range(0,len(ActiveJobs)):
            ActiveJobs[i][1] += 10 # increment time counter for remaining jobs
        
    
    # Output Data
    FailedFrames = []
    for n in range(0,len(Jobs)):
        try:
            o3 = session.openOdb(name=Directory+Jobs[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+Jobs[n]+'.odb']
        except:
            #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
            pass
        for p in range(0,len(Preload)):
            try:
                if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                    reply = getWarningReply(message=(session.odbs[Directory+Jobs[n]+'.odb'].steps['Compress'].frames[p+1].frameValue), buttons=(YES,NO))
                    # LBU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                    # LBS2
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                    # P1U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                    # P2U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                    # LTU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                else:
                    # make a list of failed times attached to the job and eventually export as txt file or similar
                    FailedFrames.append('Job '+Job[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
            except:
                #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                FailedFrames.append('Job '+Job[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
            
        
    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))

def FuncTest(): # all good
    from abaqus import getWarningReply, YES, NO
    temp = PIDs('standard.exe')
    for PID in temp:
        getWarningReply(message=(PID), buttons=(YES,NO))
    

def TestSysKill():
    import signal
    import os
    import time
    from abaqus import getWarningReply, YES, NO
    import numpy
    output = os.popen('wmic process get name, processid').read()
    PIDs = []
    for line in output.splitlines():
        if line.startswith('standard.exe'):
            line = line.replace('standard.exe','')
            line = line.replace(' ','')
            PIDs.append(line)
            #getWarningReply(message=(line), buttons=(YES,NO))
    #reply = getWarningReply(message=(output), buttons=(YES,NO))
    #for n in range(0,10):
      #  time.sleep(10)
     #   if n == 5: os.kill(4952, signal.SIGINT) # Seems to have worked as expected, will check if .023 file disappears 


def OLD_Mesh17_BoundaryHold():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import keyboard
    import signal
    from abaqus import getWarningReply, YES, NO
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def thread_function(name):
        time.sleep(10)
    nJobs = 0 # Counting Variable
    nMaxJobs = 2 # Max number of Jobs running at once
    timeOut = 75 # Max seconds per job
    Preload = [0.05,0.1] #,0.15,0.2,0.25,0.3] # Strain in the layer
    PreloadN = ['5%','10%'] #,'15%','20%','25%','30%']
    LayerElasticity = [25000]
    LayerElasticityN = ['25kPa']
    BackgroundElasticity = 25000
    InclusionElasticity = [50000] # Doesn't converge after 7.5x background
    InclusionElasticityN = ['50kPa']
    LayerThickness = [0.0005]
    LayerThicknessN = ['500um']
    #LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    #LayerThicknessN = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = [750e-06] # Minimum 50e-06
    InclusionDepthN = ['750um']
    InclusionThickness = [0.002] # Minimum 10e-06, allowed to have inclusion touching bottom of sample
    InclusionThicknessN = ['2mm']
    InclusionRadius = [0.001] # Minimum 10e-06
    InclusionRadiusN = ['1mm']
    ListSamples = [InclusionDepth, InclusionThickness, InclusionRadius, InclusionElasticity]
    ListSamplesN = [InclusionDepthN, InclusionThicknessN, InclusionRadiusN, InclusionElasticityN]
    Samples = [list(p) for p in itertools.product(*ListSamples)]
    SamplesN = [list(p) for p in itertools.product(*ListSamplesN)]
    ListLayers = [LayerThickness, LayerElasticity]
    ListLayersN = [LayerThicknessN, LayerElasticityN]
    Layers = [list(p) for p in itertools.product(*ListLayers)]
    LayersN = [list(p) for p in itertools.product(*ListLayersN)]
    # ModelName+'_'+'ID'+'_'+InclusionDepthN[id]+'IT'+'_'+InclusionThicknessN[it]+'IR'+'_'+InclusionRadiusN[ir]+
    #          'IE'+'_'+InclusionElasticityN[ie]+'_'+'LT'+'_'+LayerThicknessN[lt]+'_'+'LE'+'_'+LayerElasticityN[le]+'_'+'P'+'_'PreloadN[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh17-Final'
    Jobs = []
    ActiveJobs = []
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # Alter and Submit Model
    for n1 in range(0,len(Samples)):
        # Change Sample
        p = mdb.models[ModelName].parts['Sample']
        s = p.features['Partition face-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        if SampleThickness-Samples[n1][0]-Samples[n1][1] == 0: # If Inclusion touches the bottom of the sample
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            d[2].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[0].setValues(value=50e-06, ) # Small multiple of mesh size
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            d[2].setValues(value=SampleThickness-Samples[n1][0]-Samples[n1][1], ) # Distance from bottom of sample to bottom of inclusion
            d[0].setValues(value=Samples[n1][1], ) # Inclusion Thickness
            d[1].setValues(value=Samples[n1][2], ) # Inclusion Radius
        p.features['Partition face-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((Samples[n1][3], 0.45), )) # Set Inclusion Material Properties
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties
        for n2 in range(0,len(Layers)):
            # Change Layer
            p = mdb.models[ModelName].parts['Layer']
            s = p.features['Shell planar-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s1 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
            d[6].setValues(value=Layers[n2][0], ) # Layer Thickness
            p.features['Shell planar-1'].setValues(sketch=s1)
            del mdb.models[ModelName].sketches['__edit__']
            s = p.features['Partition face-1'].sketch
            mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
            s2 = mdb.models[ModelName].sketches['__edit__']
            g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
            if Layers[n2][0] > 0.00125:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.0005, ) # first to second particle layer
            elif Layers[n2][0] > 0.00076:
                d[0].setValues(value=0.0005, ) # bottom to first particle layer
                d[1].setValues(value=0.00001, ) # first to second particle layer
            else:
                d[0].setValues(value=0.0001, ) # bottom to first particle layer
                d[1].setValues(value=0.0001, ) # first to second particle layer
            p.features['Partition face-1'].setValues(sketch=s2)
            del mdb.models[ModelName].sketches['__edit__']
            p.regenerate()
            p.generateMesh()
            mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((Layers[n2][1], 0.45), )) # Set Layer Material Properties
            
            # Regenerate Assembly
            a = mdb.models[ModelName].rootAssembly
            a.regenerate()
            # Boundary Conditions, gradient of position with time equation derived in notebook
            B_t = -Preload[-1]*SampleThickness*Layers[n2][1]/BackgroundElasticity
            T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness*Layers[n2][1]/BackgroundElasticity)
            # OLD, calcs preload as total % compression of both sample and layer #T_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)
            # OLD, calcs preload as total % compression of both sample and layer #B_t = -Preload[-1]*(Layers[n2][0]+SampleThickness)*(1-Layers[n2][0]*BackgroundElasticity/(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1]))
            mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
            mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
            # Output Times Calc and Set
            OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
            #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
            OutputTimesTupple = ()
            for i in range(0,len(OutputTimes)):
                OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
            mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
            mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
            # Make Job
            Jobs.append(ModelName+'_ID_'+SamplesN[n1][0]+'_IT_'+SamplesN[n1][1]+'_IR_'+SamplesN[n1][2]+'_IE_'+SamplesN[n1][3]+'_LT_'+LayersN[n2][0]+'_LE_'+LayersN[n2][1])
            ActiveJobs.append([Jobs[-1],20])
            j = mdb.Job(name=Jobs[-1], 
                model=ModelName, 
                description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
                queue=None, memory=90, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=8, 
                numDomains=8, numGPUs=0)
            j.submit(consistencyChecking=OFF)
            nJobs = nJobs + 1
            
            x = threading.Thread(target=thread_function, args=(1,))
            x.start()
            reply = getWarningReply(message=('First Wait Joined'), buttons=(YES,NO))
            x.join()
            #reply = getWarningReply(message=('got past first wait'), buttons=(YES,NO))
            # have up to nMaxJobs running at the same time, submitting with a minute break in between
            while nJobs >= nMaxJobs:
                x = threading.Thread(target=thread_function, args=(1,))
                x.start()
                x.join()
                ToRemove = []
                ToKill = []
                for Job in ActiveJobs:
                    if os.path.isfile(Job[0]+'.cid') == False: # Check if job is either in analysis or running by seeing if license file exists
                        ToRemove.append(Job)
                    if Job[1] > timeOut:
                        ToKill.append(Job)
                        ToRemove.append(Job)
                for Job in ToKill:
                    mdb.jobs[Job[0]].kill()
                    #reply = getWarningReply(message=('Tried to kill '+Jobs[0]), buttons=(YES,NO))
                for Job in ToRemove:
                    ActiveJobs.remove(Job)
                    nJobs = nJobs - 1
                for i,Job in enumerate(ActiveJobs):
                    ActiveJobs[i][1] += 10 # increment time counter for remaining jobs
                
            
        
    
    # Wait for last jobs to be completed
    while nJobs > 0:
        x = threading.Thread(target=thread_function, args=(1,))
        x.start()
        x.join()
        #reply = getWarningReply(message=('Does this work as a delay?'), buttons=(YES,NO)) # this would fix the problem though soooo
        ToRemove = []
        ToKill = []
        for Job in ActiveJobs:
            if os.path.isfile(Job[0]+'.cid') == False:
                ToRemove.append(Job)
            if Job[1] > timeOut:
                ToKill.append(Job)
                ToRemove.append(Job)
        for Job in ToKill:
            mdb.jobs[Job[0]].kill()
            while os.path.isfile(Job[0]+'.cid') == True:
                x = threading.Thread(target=thread_function, args=(1,))
                x.start()
                x.join()
        for Job in ToRemove:
            ActiveJobs.remove(Job)
            nJobs = nJobs - 1
        for i in range(0,len(ActiveJobs)):
            ActiveJobs[i][1] += 10 # increment time counter for remaining jobs
        
    
    # Output Data
    FailedFrames = []
    for n in range(0,len(Jobs)):
        try:
            o3 = session.openOdb(name=Directory+Jobs[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+Jobs[n]+'.odb']
        except:
            #reply = getWarningReply(message=('Database '+Jobs[n]+' output processing failed: '), buttons=(YES,NO))
            pass
        for p in range(0,len(Preload)):
            try:
                if numpy.abs(OutputTimes[p]-odb.steps['Compress'].frames[p+1].frameValue) < 0.01: # If reported time is within 0.01 seconds of listed time we accept the frame
                    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                    reply = getWarningReply(message=(session.odbs[Directory+Jobs[n]+'.odb'].steps['Compress'].frames[p+1].frameValue), buttons=(YES,NO))
                    # LBU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOMSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBU2'+'.npy', xy1)
                    # LBS2
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LBS2'+'.npy', xy1)
                    # P1U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P1U2'+'.npy', xy1)
                    # P2U2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_P2U2'+'.npy', xy1)
                    # LTU2
                    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes
                    NODE_LABELS=[]
                    for i in range(len(NODES)):
                        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOPSHORT'].nodes[i].label)
                    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                    pth = session.paths['Path-1']
                    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                        removeDuplicateXYPairs=True, includeAllElements=False)
                    numpy.save(Directory+'MacroOutputs/'+Jobs[n]+'_'+'P'+'_'+PreloadN[p]+'_LTU2'+'.npy', xy1)
                else:
                    # make a list of failed times attached to the job and eventually export as txt file or similar
                    FailedFrames.append('Job '+Job[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to time mismatch')
            except:
                #reply = getWarningReply(message=('Frame '+PreloadN[p]+' in Job '+Jobs[n]+' output processing failed'+message), buttons=(YES,NO))
                FailedFrames.append('Job '+Job[n]+' Frame '+str(p+1)+' at time '+str(OutputTimes[p])+' failed due to non-existence')
            
        
    numpy.save(Directory+'MacroOutputs/'+'JobList'+'.npy', numpy.array(Jobs))
    numpy.save(Directory+'MacroOutputs/'+'FailedFrames'+'.npy', numpy.array(FailedFrames))

def OLD_FullSetBoundaryHold():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Preload = [0.05,0.1,0.15,0.2,0.25,0.3] # Strain in the layer
    LayerE = [25000]
    SampleE = [25000]
    InclusionE = [50000]
    LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    SampleThickness = 0.005
    InclusionDepth = 0.001
    InclusionThickness = 0.002
    InclusionRadius = 0.001
    OutputTimes = Preload*(LayerThickness*SampleE+SampleThickness*LayerE)/(Preload[-1])
    # JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    # JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    JobNames = ['Sample1_L10','Sample1_L7p5','Sample1_L5','Sample1_L4','Sample1_L3','Sample1_L2','Sample1_L1','Sample1_L0p5']
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'LinearElastic_Compress_Inclusion_30p-Automate'
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    for n in range(0,len(ThicknessRange)):
        p = mdb.models[ModelName].parts['Layer']
        s = p.features['Shell planar-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[6].setValues(value=ThicknessRange[n], )
        p = mdb.models[ModelName].parts['Layer']
        p.features['Shell planar-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        if ThicknessRange[n] > 0.0015:
            p.features['Datum plane-2'].setValues(offset=1e-03)
        else:
            p.features['Datum plane-2'].setValues(offset=7.2e-04)
        if ThicknessRange[n] > 0.00075:
            p.features['Datum plane-1'].setValues(offset=5e-04)
        else:
            p.features['Datum plane-2'].setValues(offset=6e-05)
            p.features['Datum plane-1'].setValues(offset=3e-05)
        p = mdb.models[ModelName].parts['Layer']
        p.regenerate()
        p.generateMesh()
        a = mdb.models[ModelName].rootAssembly
        a.regenerate()
        j = mdb.Job(name=JobNames[n], 
            model=ModelName, 
            description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
            queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
            numDomains=4, numGPUs=0)
        j.submit(consistencyChecking=OFF)
        try: # COULD SUBMIT MORE JOBS AT A TIME AND WAIT FOR EVERY Nth JOB TO COMPLETE INSTEAD OF WAITING FOR EACH, WOULD REQUIRE SUBMITTING SHORTER JOBS FIRST
            j.waitForCompletion(60*10)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' timed out', message
    for n in range(0,len(ThicknessRange)):
        try:
            o3 = session.openOdb(name=Directory+JobNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobNames[n]+'.odb']
        except AbaqusException, message:
                print 'Database '+JobNames[n]+' output processing failed', message
        for p in range(0,len(PreloadNames)):
            try:
                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                # LBU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LBU2'+'.npy', xy1)
                # LBS2
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LBS2'+'.npy', xy1)
                # P1U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_P1U2'+'.npy', xy1)
                # P2U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_P2U2'+'.npy', xy1)
                # LTU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LTU2'+'.npy', xy1)
            except AbaqusException, message:
                print 'Frame '+PreloadNames[p]+' in Job '+JobNames[n]+' output processing failed', message

def OLD_FullSet():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    ThicknessRange = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    # JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    JobNames = ['Sample1_L10','Sample1_L7p5','Sample1_L5','Sample1_L4','Sample1_L3','Sample1_L2','Sample1_L1','Sample1_L0p5']
    PreloadNames = ['10%','15%','20%','25%','30%']
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'LinearElastic_Compress_Inclusion_30p-Automate'
    for n in range(0,len(ThicknessRange)):
        p = mdb.models[ModelName].parts['Layer']
        s = p.features['Shell planar-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[6].setValues(value=ThicknessRange[n], )
        p = mdb.models[ModelName].parts['Layer']
        p.features['Shell planar-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        if ThicknessRange[n] > 0.0015:
            p.features['Datum plane-2'].setValues(offset=1e-03)
        else:
            p.features['Datum plane-2'].setValues(offset=7.2e-04)
        if ThicknessRange[n] > 0.00075:
            p.features['Datum plane-1'].setValues(offset=5e-04)
        else:
            p.features['Datum plane-2'].setValues(offset=6e-05)
            p.features['Datum plane-1'].setValues(offset=3e-05)
        p = mdb.models[ModelName].parts['Layer']
        p.regenerate()
        p.generateMesh()
        a = mdb.models[ModelName].rootAssembly
        a.regenerate()
        j = mdb.Job(name=JobNames[n], 
            model=ModelName, 
            description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
            queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
            numDomains=4, numGPUs=0)
        j.submit(consistencyChecking=OFF)
        try: # COULD SUBMIT MORE JOBS AT A TIME AND WAIT FOR EVERY Nth JOB TO COMPLETE INSTEAD OF WAITING FOR EACH, WOULD REQUIRE SUBMITTING SHORTER JOBS FIRST
            j.waitForCompletion(60*10)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' timed out', message
    for n in range(0,len(ThicknessRange)):
        try:
            o3 = session.openOdb(name=Directory+JobNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobNames[n]+'.odb']
        except AbaqusException, message:
                print 'Database '+JobNames[n]+' output processing failed', message
        for p in range(0,len(PreloadNames)):
            try:
                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                # LBU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LBU2'+'.npy', xy1)
                # LBS2
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LBS2'+'.npy', xy1)
                # P1U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_P1U2'+'.npy', xy1)
                # P2U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_P2U2'+'.npy', xy1)
                # LTU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobNames[n]+'_'+PreloadNames[p]+'_LTU2'+'.npy', xy1)
            except AbaqusException, message:
                print 'Frame '+PreloadNames[p]+' in Job '+JobNames[n]+' output processing failed', message
            

def OLD_FullSetBackup():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    ThicknessRange = [0.01] #[0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    # JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 1%, 10%, 20%, 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    JobNames = ['MacroTest_L10','MacroTest_L7p5','MacroTest_L5','MacroTest_L4','MacroTest_L3','MacroTest_L2','MacroTest_L1','MacroTest_L0p5']
    PreloadNames = ['10%','15%','20%','25%','30%']
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    for n in range(0,len(ThicknessRange)):
        p = mdb.models[ModelName].parts['Layer']
        s = p.features['Shell planar-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[6].setValues(value=ThicknessRange[n], )
        p = mdb.models[ModelName].parts['Layer']
        p.features['Shell planar-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        if ThicknessRange[n] > 0.0015:
            p.features['Datum plane-2'].setValues(offset=1e-03)
            p.features['Datum plane-1'].setValues(offset=5e-04)
        else:
            p.features['Datum plane-2'].setValues(offset=6e-05)
            p.features['Datum plane-1'].setValues(offset=3e-05)
        p = mdb.models[ModelName].parts['Layer']
        p.regenerate()
        p.generateMesh()
        a = mdb.models[ModelName].rootAssembly
        a.regenerate()
        j = mdb.Job(name=JobNames[n], 
            model=ModelName, 
            description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
            queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
            numDomains=4, numGPUs=0)
        j.submit(consistencyChecking=OFF)
        try:
            j.waitForCompletion(60*10)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' timed out', message
    for n in range(0,len(ThicknessRange)):
        try:
            o3 = session.openOdb(name=Directory+JobNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobNames[n]+'.odb']
            # LBU2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LBU2'+'.npy', xy1)
            # LBS2
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='S', outputPosition=NODAL, refinement=(COMPONENT, 'S22'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LBS2'+'.npy', xy1)
            # P1U2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_P1U2'+'.npy', xy1)
            # P2U2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_P2U2'+'.npy', xy1)
            # LTU2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LTU2'+'.npy', xy1)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' output processing failed', message
        
        
def OLD_ExportResults(): 
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    #a = mdb.models[ModelName].rootAssembly
    #session.viewports['Viewport: 1'].setValues(displayedObject=a)
    #a = mdb.models[ModelName].rootAssembly
    #session.viewports['Viewport: 1'].setValues(displayedObject=a)
    o3 = session.openOdb(name=Directory+'MacroTest_L5.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    session.mdbData.summary()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
    odb = session.odbs[Directory+'MacroTest_L5.odb']
    NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
    NODE_LABELS=[]
    for i in range(len(NODES)):
        NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
    #xyp = session.XYPlot('XYPlot-3')
    #chartName = xyp.charts.keys()[0]
    #chart = xyp.charts[chartName]
    pth = session.paths['Path-1']
    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
        removeDuplicateXYPairs=True, includeAllElements=False)
    numpy.save("F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/foo2.npy", xy1)
    #c1 = session.Curve(xyData=xy1)
    #chart.setValues(curvesToPlot=(c1, ), )
    #session.charts[chartName].autoColor(lines=True, symbols=True)
    #session.viewports['Viewport: 1'].setValues(displayedObject=xyp)

def OLD_TestExport():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    ThicknessRange = [0.01] #[0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    # Want to run for each thickness, at preload's of 1%, 10%, 20%, 30%, 40%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    JobNames = ['MacroTest_L10','MacroTest_L7p5','MacroTest_L5','MacroTest_L4','MacroTest_L3','MacroTest_L2','MacroTest_L1','MacroTest_L0p5']
    for n in range(0,len(ThicknessRange)):
        try:
            o3 = session.openOdb(name=Directory+JobNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobNames[n]+'.odb']
            # LBU2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LBU2'+'.npy', xy1)
            # LBS2
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LBS2'+'.npy', xy1)
            # P1U2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_P1U2'+'.npy', xy1)
            # P2U2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_P2U2'+'.npy', xy1)
            # LTU2
            NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
            NODE_LABELS=[]
            for i in range(len(NODES)):
                NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
            session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
            pth = session.paths['Path-1']
            xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                removeDuplicateXYPairs=True, includeAllElements=False)
            numpy.save(Directory+JobNames[n]+'_LTU2'+'.npy', xy1)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' output processing failed', message

def ChangeInclusion():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=ON)
    p1 = mdb.models['MacroRecording'].parts['Sample']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.62852, 
        farPlane=0.9273, width=0.126778, height=0.0740087, 
        viewOffsetX=-0.0574821, viewOffsetY=0.00046139)
    p = mdb.models['MacroRecording'].parts['Sample'] # sets current part, only has to be done once
    p.features['Datum plane-1'].setValues(offset=0.002) # radius of cylindrical inclusion
    p = mdb.models['MacroRecording'].parts['Sample']
    p.regenerate()
    p = mdb.models['MacroRecording'].parts['Sample']
    p.features['Datum plane-2'].setValues(offset=0.003) # distance from sample bottom to bottom of cylindrical inclusion
    p = mdb.models['MacroRecording'].parts['Sample']
    p.regenerate()
    p = mdb.models['MacroRecording'].parts['Sample']
    p.features['Datum plane-3'].setValues(offset=0.0045) # distance from sample bottom to top of cylindrical inclusion
    p = mdb.models['MacroRecording'].parts['Sample']
    p.regenerate() # has to be done once before generating mesh and after all changes
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['MacroRecording'].parts['Sample'] # used to regenerate mesh
    p.generateMesh() # regenerating mesh
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON, mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    mdb.models['MacroRecording'].materials['Inclusion'].elastic.setValues(table=((
        75000.0, 0.45), )) # changing elasticity of inclusion, should progrably regenerate instances before starting sim


def OLD_OverflowTests():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Range = [1,2,3,4,5]
    Directory = 'F:/Rhys/AbaqusModels/LinearElastic_OP_Investigation/'
    JobNames = ['OverflowTests_Rad_1_longer','OverflowTests_Rad_2_longer','OverflowTests_Rad_3_longer','OverflowTests_Rad_4_longer','OverflowTests_Rad_5_longer']
    ModelName = 'LinearElastic_Compress_CentralStick_30p-Automate-DispControlled-DenseMesh'
    for n in range(0,len(Range)):
        p = mdb.models[ModelName].parts['RigidPlateBottom']
        s = p.features['2D Analytic rigid shell-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[2].setValues(value=Range[n], )
        p.features['2D Analytic rigid shell-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        p = mdb.models[ModelName].parts['RigidPlateTop']
        s = p.features['2D Analytic rigid shell-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[1].setValues(value=Range[n], )
        p.features['2D Analytic rigid shell-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        p.regenerate()
        a = mdb.models[ModelName].rootAssembly
        a.regenerate()
        j = mdb.Job(name=JobNames[n], 
            model=ModelName, 
            description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
            queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
            numDomains=4, numGPUs=0)
        j.submit(consistencyChecking=OFF)
        try:
            j.waitForCompletion(60*10)
        except AbaqusException, message:
            print 'Job '+JobNames[n]+' timed out', message                
        

def OLD_DataOnly():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    Preload = [0.05,0.1,0.15,0.2,0.25,0.3] # Strain in the layer
    LayerE = 25000
    SampleE = 25000
    InclusionE = 50000
    LayerThickness = [0.01,0.0005] # [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    LayerThicknessNames = ['10','0p5'] # ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = 0.001
    InclusionThickness = 0.002
    InclusionRadius = 0.001
    PreloadNames = ['5%','10%','15%','20%','25%','30%']
    ## JobName+'_'+'LT'+'_'+LayerThickness[j]+'_'+'P'+'_'Preload[p]+'_'+OutputNames[o]+'.npy'
    ## JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    ## Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    ## Can encode preload data in the load magnitude, ramp function and output instance sampling times
    ## Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    JobName = 'Mesh3Lagrange0p4'
    ModelName = 'Mesh3'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # Output Data
    for n in range(0,len(LayerThickness)):
        try:
            o3 = session.openOdb(name=Directory+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'.odb']
        except AbaqusException, message:
                print 'Database '+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+' output processing failed', message
        for p in range(0,len(Preload)):
            try:
                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                # LBU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LBU2'+'.npy', xy1)
                # LBS2
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LBS2'+'.npy', xy1)
                # P1U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_P1U2'+'.npy', xy1)
                # P2U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_P2U2'+'.npy', xy1)
                # LTU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LTU2'+'.npy', xy1)
            except AbaqusException, message:
                print 'Frame '+PreloadNames[p]+' in Job '+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+' output processing failed', message


def OLD_BoundaryHold():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import os, time
    # RUN COMMAND "showStopButtonInGui" in CAE COMMAND LINE BEFORE BEGINNING MACRO TO HAVE AN ABORT BUTTON
    nJobs = 0 # Counting Variable
    nMaxJobs = 4 # Max number of Jobs running at once
    timeOut = 1800 # Max seconds per job
    InclusionMeshSize = 2.5e-06 # Doesn't change the mesh size, this program just needs to know it
    LayerMeshSize = 25e-06 # Doesn't change the mesh size, this program just needs to know it
    Preload = [0.05,0.1,0.15,0.2,0.25,0.3] # Strain in the layer
    PreloadNames = ['5%','10%','15%','20%','25%','30%']
    LayerE = 25000
    BackgroundE = 25000
    InclusionE = 250000 #50000
    LayerThickness = [0.0005,0.01]
    LayerThicknessNames = ['0p5','10']
    #LayerThickness = [0.01,0.0075,0.005,0.004,0.003,0.002,0.001,0.0005]
    #LayerThicknessNames = ['10','7p5','5','4','3','2','1','0p5']
    SampleThickness = 0.005
    InclusionDepth = 100e-06 #0.001
    InclusionThickness = 0.003 #0.002
    InclusionRadius = 0.003 #0.001
    # JobName+'_'+'LT'+'_'+LayerThickness[j]+'_'+'P'+'_'Preload[p]+'_'+OutputNames[o]+'.npy'
    # JobNames[j]+'_'+PreloadNames[p]+'_'+OutputNames[o]+'.npy'
    # Want to run for each thickness, at preload's of 10%, 15%, 20%, 25% & 30%, then want to run for different sample-inclusion mech contrasts and shapes
    # Can encode preload data in the load magnitude, ramp function and output instance sampling times
    # Need extra for loops for the sample changes (or can do that manually and repeat macro), then another file to export the data if possible.
    JobName = 'Mesh8_test_100umDepth' #'Mesh7-ExtremeParametersLagrange0p4'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Mesh8'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Automate/'
    # Alter and Submit Model
    for n in range(0,len(LayerThickness)):
        # Change Layer
        p = mdb.models[ModelName].parts['Layer']
        s = p.features['Shell planar-1'].sketch
        mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
        s1 = mdb.models[ModelName].sketches['__edit__']
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        d[6].setValues(value=LayerThickness[n], )
        p = mdb.models[ModelName].parts['Layer']
        p.features['Shell planar-1'].setValues(sketch=s1)
        del mdb.models[ModelName].sketches['__edit__']
        if LayerThickness[n] > 0.0015:
            p.features['Datum plane-2'].setValues(offset=1e-03)
        else:
            p.features['Datum plane-2'].setValues(offset=1e-03-LayerMeshSize)
        if LayerThickness[n] > 0.00075:
            p.features['Datum plane-1'].setValues(offset=500e-06)
        else:
            p.features['Datum plane-2'].setValues(offset=500e-06-LayerMeshSize)
            p.features['Datum plane-1'].setValues(offset=500e-06-11*LayerMeshSize)
        p = mdb.models[ModelName].parts['Layer']
        p.regenerate()
        p.generateMesh()
        mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerE, 0.45), ))
        # Change Sample
        p = mdb.models[ModelName].parts['Sample'] # sets current part, only has to be done once
        if SampleThickness == InclusionThickness:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundC']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionC']) # Inclusion section node set
            p.features['Datum plane-1'].setValues(offset=InclusionRadius)
            p.features['Datum plane-2'].setValues(offset=InclusionMeshSize)
            p.features['Datum plane-3'].setValues(offset=2*InclusionMeshSize)
        elif InclusionDepth == 0:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundU']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionU']) # Inclusion section node set
            p.features['Datum plane-1'].setValues(offset=InclusionRadius)
            p.features['Datum plane-2'].setValues(offset=SampleThickness-InclusionThickness-InclusionMeshSize)
            p.features['Datum plane-3'].setValues(offset=SampleThickness-InclusionThickness)
        elif SampleThickness-InclusionDepth-InclusionThickness == 0:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
            p.features['Datum plane-1'].setValues(offset=InclusionRadius)
            p.features['Datum plane-2'].setValues(offset=InclusionThickness)
            p.features['Datum plane-3'].setValues(offset=InclusionThickness+InclusionMeshSize)
        else:
            mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
            mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
            p.features['Datum plane-1'].setValues(offset=InclusionRadius) # radius of cylindrical inclusion
            p.features['Datum plane-2'].setValues(offset=SampleThickness-InclusionDepth-InclusionThickness) # distance from sample bottom to bottom of cylindrical inclusion
            p.features['Datum plane-3'].setValues(offset=SampleThickness-InclusionDepth) # distance from sample bottom to top of cylindrical inclusion
        p.regenerate() # has to be done once before generating mesh and after all changes
        p = mdb.models[ModelName].parts['Sample'] # used to regenerate mesh
        p.generateMesh() # regenerating mesh
        mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundE, 0.45), ))
        mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionE, 0.45), )) 
        # Regenerate Assembly
        a = mdb.models[ModelName].rootAssembly
        a.regenerate()
        # Boundary Conditions
        T_t = -Preload[-1]*(LayerThickness[n]+SampleThickness)
        B_t = -Preload[-1]*(LayerThickness[n]+SampleThickness)*(1-LayerThickness[n]*BackgroundE/(LayerThickness[n]*BackgroundE+SampleThickness*LayerE))
        mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
        mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
        mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
        # Output Times Calc and Set
        OutputTimes = (numpy.array(Preload)*(LayerThickness[n]*BackgroundE+SampleThickness*LayerE)/(Preload[-1]*(LayerThickness[n]+SampleThickness)*BackgroundE)).tolist()
        OutputTimesTupple = ()
        for i in range(0,len(OutputTimes)):
            OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
        mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
        mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        # Make Job
        j = mdb.Job(name=JobName+'_'+'LT'+'_'+LayerThicknessNames[n], 
            model=ModelName, 
            description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
            queue=None, memory=90, memoryUnits=PERCENTAGE, 
            getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
            resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
            numDomains=4, numGPUs=0)
        j.submit(consistencyChecking=OFF)
        nJobs = nJobs + 1
        
        
        total_time = 60
        time.sleep(60)
        # have up to nMaxJobs running at the same time, submitting with a minute break in between
        while nJobs >= nMaxJobs:
            if os.path.isfile(JobName+'_'+'LT'+'_'+LayerThicknessNames[n+1-nMaxJobs]+'.023') == False:
                nJobs = nJobs - 1
                total_time = 0
            if total_time > timeOut:
                mdb.jobs[JobName+'_'+'LT'+'_'+LayerThicknessNames[n+1-nMaxJobs]].kill()
                total_time = 0
            total_time += 30
            time.sleep(30)
        
        #try: # COULD SUBMIT MORE JOBS AT A TIME AND WAIT FOR EVERY Nth JOB TO COMPLETE INSTEAD OF WAITING FOR EACH, WOULD REQUIRE SUBMITTING SHORTER JOBS FIRST
        #    j.waitForCompletion(60*1)
        #except AbaqusException, message:
        #    mdb.jobs[JobName+'_'+'LT'+'_'+LayerThicknessNames[n]].kill()
        #    print 'Job '+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+' timed out', message
    # Wait for last nMaxJobs - 1 jobs to be completed
    CheckDone = 0
    total_time = 0
    while CheckDone == 0:
        for n in range(0,len(LayerThickness)):
            if os.path.isfile(JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'.023') == True:
                CheckDone = 0
                break
            CheckDone = 1
        if total_time > timeOut:
            mdb.jobs[JobName+'_'+'LT'+'_'+LayerThicknessNames[n]].kill()
        total_time += 30
        time.sleep(30)
    
    # Output Data
    for n in range(0,len(LayerThickness)):
        try:
            o3 = session.openOdb(name=Directory+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            session.viewports['Viewport: 1'].makeCurrent()
            session.mdbData.summary()
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.pickingExpression._processNodeSets(selectedNodeSets='LAYER-1.LAYERBOTTOM,')
            odb = session.odbs[Directory+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'.odb']
        except AbaqusException, message:
                print 'Database '+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+' output processing failed', message
        for p in range(0,len(Preload)):
            try:
                session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=p+1 )
                # LBU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERBOTTOM'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LBU2'+'.npy', xy1)
                # LBS2
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S22'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LBS2'+'.npy', xy1)
                # P1U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['PARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_P1U2'+'.npy', xy1)
                # P2U2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['SECONDPARTICLELAYER'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_P2U2'+'.npy', xy1)
                # LTU2
                NODES=odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes
                NODE_LABELS=[]
                for i in range(len(NODES)):
                    NODE_LABELS.append(odb.rootAssembly.instances['LAYER-1'].nodeSets['LAYERTOP'].nodes[i].label)
                session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', NODE_LABELS), ))
                session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
                pth = session.paths['Path-1']
                xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
                numpy.save(Directory+'MacroOutputs/'+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+'_LTU2'+'.npy', xy1)
            except AbaqusException, message:
                print 'Frame '+PreloadNames[p]+' in Job '+JobName+'_'+'LT'+'_'+LayerThicknessNames[n]+'_'+'P'+'_'+PreloadNames[p]+' output processing failed', message
                
        



def OLD_Mesh8Make():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior

    p = mdb.models['Mesh8'].parts['Layer']
    for n in range(4,int(0.0075/0.00005)):
        p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=n*0.00005)



def OLD_PartitionOverInclusion():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior

    p = mdb.models['Mesh9'].parts['Sample']
    f1, e1, d2 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(
        0, 0.005, 0.0))
    s1 = mdb.models['Mesh9'].ConstrainedSketch(name='__profile__', sheetSize=0.015, 
        gridSpacing=0.00025, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(decimalPlaces=4)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    for n in range(1,int(0.003/0.00005)):
        s1.Line(point1=(0.00005*n, -1e-05), point2=(0.00005*n, -4e-05))
    s1.Line(point1=(0, -1e-05), point2=(0.00005*(n+1), -1e-05))
    s1.Line(point1=(0, -4e-05), point2=(0.00005*(n+1), -4e-05))
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1000 ]', ), )
    e, d1 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Mesh9'].sketches['__profile__']
    
def OLD_PartitionMiddleSample():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior

    p = mdb.models['Mesh9'].parts['Sample']
    f1, e1, d2 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f1[72], sketchPlaneSide=SIDE1, origin=(
        0.003, 0.005, 0.0))
    s1 = mdb.models['Mesh9'].ConstrainedSketch(name='__profile__', sheetSize=0.015, 
        gridSpacing=0.00025, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(decimalPlaces=4)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    n=0
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    for n in range(1,int(0.0045/0.00005)+1):
        s1.Line(point1=(0.00005*n, -1e-05), point2=(0.00005*n, -4e-05))
    s1.Line(point1=(0, -1e-05), point2=(0.00005*(n+1), -1e-05))
    s1.Line(point1=(0, -4e-05), point2=(0.00005*(n+1), -4e-05))
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#0:2 #100 ]', ), )
    e, d1 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Mesh9'].sketches['__profile__']
    

def TestIter():
    import itertools
    l1 = [1,2]
    l2 = [4,5]
    l3 = [7,8]
    l = [l1,l2,l3]

    combinations = [list(p) for p in itertools.product(*l)]


def ChangeSample():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p1 = mdb.models['Mesh17-Final'].parts['Sample']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Mesh17-Final'].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models['Mesh17-Final'].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models['Mesh17-Final'].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, 
        upToFeature=p.features['Partition face-1'], filter=COPLANAR_EDGES)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.205535, 
        farPlane=0.21499, width=0.0220649, height=0.0124167, cameraPosition=(
        0.0093098, 0.00220864, 0.210262), cameraTarget=(0.0093098, 0.00220864, 
        0))
    d[2].setValues(value=0.0015, )
    d[0].setValues(value=0.003, )
    d[1].setValues(value=0.00325, )
    s1.unsetPrimaryObject()
    p = mdb.models['Mesh17-Final'].parts['Sample']
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models['Mesh17-Final'].sketches['__edit__']
    p = mdb.models['Mesh17-Final'].parts['Sample']
    p.regenerate()
    p = mdb.models['Mesh17-Final'].parts['Sample']
    p.generateMesh()
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON, mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    mdb.models['Mesh17-Final'].materials['Inclusion'].elastic.setValues(table=((
        250000.0, 0.45), ))
    mdb.models['Mesh17-Final'].materials['Layer'].elastic.setValues(table=((
        25001.0, 0.45), ))
    mdb.models['Mesh17-Final'].materials['Layer'].elastic.setValues(table=((
        25000.0, 0.45), ))


def ChangeLayer():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Mesh17-Final'].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models['Mesh17-Final'].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models['Mesh17-Final'].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, 
        upToFeature=p.features['Shell planar-1'], filter=COPLANAR_EDGES)
    d[6].setValues(value=0.0005, ) # Layer Thickness
    s1.unsetPrimaryObject()
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models['Mesh17-Final'].sketches['__edit__']
    p = mdb.models['Mesh17-Final'].parts['Layer']
    s = p.features['Partition face-1'].sketch
    mdb.models['Mesh17-Final'].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models['Mesh17-Final'].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    s2.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s2, 
        upToFeature=p.features['Partition face-1'], filter=COPLANAR_EDGES)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.205449, 
        farPlane=0.217402, width=0.00531672, height=0.00299189, 
        cameraPosition=(0.00870002, 0.000734773, 0.211426), cameraTarget=(
        0.00870002, 0.000734773, 0))
    d[1].setValues(value=0.0001, ) # first to second particle layer
    d[0].setValues(value=0.0001, ) # bottom to first particle layer
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.205632, 
        farPlane=0.21722, width=0.00328085, height=0.00184624, cameraPosition=(
        0.00792419, 0.000465607, 0.211426), cameraTarget=(0.00792419, 
        0.000465607, 0))
    s2.unsetPrimaryObject()
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models['Mesh17-Final'].sketches['__edit__']
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.regenerate()
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.201685, 
        farPlane=0.203581, width=0.00750934, height=0.00422575, 
        viewOffsetX=-0.0425291, viewOffsetY=-0.00450957)
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.generateMesh()


def KillTest():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import time
    import threading
    import ctypes
    from abaqus import getWarningReply, YES, NO
    def thread_function(name):
        time.sleep(2)
        ctypes.windll.user32.keybd_event(0x0D, 0, 0, 0) #Enter Key Down
        ctypes.windll.user32.keybd_event(0x0D, 0, 0x0002, 0) #Enter Key Up
    x = threading.Thread(target=thread_function, args=(1,), daemon=True)
    x.start()
    reply = getWarningReply(message='hi', buttons=(YES,NO))
    x.join()
    reply = getWarningReply(message=('Skip Job Submission?'), buttons=(YES,NO))
    if reply == NO:
        mdb.jobs['Soft2'].submit(consistencyChecking=OFF)
    time.sleep(1)
    reply = getWarningReply(message=('1 second'), buttons=(YES,NO))
    time.sleep(1)
    reply = getWarningReply(message=('2 seconds'), buttons=(YES,NO))
    time.sleep(2)
    reply = getWarningReply(message=('3 seconds'), buttons=(YES,NO))
    time.sleep(3)
    reply = getWarningReply(message=('Kill Job?'), buttons=(YES,NO))
    if reply == YES:
        mdb.jobs['Soft2'].kill()
        reply = getWarningReply(message=('Wait for its death'), buttons=(YES,NO))
    

def Macro1():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    a = mdb.models['Mesh17-Final'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a = mdb.models['Mesh17-Final'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    o3 = session.openOdb(
        name='F:/Rhys/AbaqusModels/BoundaryHold_Automate/Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_2x_LT_0p5_LE_25kPa.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    session.mdbData.summary()
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
        'Magnitude'), )
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(
        plotState=CONTOURS_ON_DEF)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'), 
        )
    session.pickingExpression._processNodeSets(
        selectedNodeSets='LAYER-1.LAYERBOTTOMSHORT,')
    session.Path(name='Path-1', type=NODE_LIST, expression=(('LAYER-1', (18, 19, 
        329, 1792, 1793, 1794, 1795, 1796, 1797, 1798, 1799, 1800, 1801, 1802, 
        1803, 1804, 1805, 1806, 1807, 1808, 1809, 1810, 1811, 1812, 1813, 1814, 
        1815, 1816, 1817, 1818, 1819, 1820, 1821, 1822, 1823, 1824, 1825, 1826, 
        1827, 1828, 1829, 1830, 1831, 1832, 1833, 1834, 1835, 1836, 1837, 1838, 
        1839, 1840, 1841, 1842, 1843, 1844, 1845, 1846, 1847, 1848, 1849, 1850, 
        1851, 1852, 1853, 1854, 1855, 1856, 1857, 1858, 1859, 1860, 1861, 1862, 
        1863, 1864, 1865, 1866, 1867, 1868, 1869, 1870, 1871, 1872, 1873, 1874, 
        1875, 1876, 1877, 1878, 1879, 1880, 1881, 1882, 1883, 1884, 1885, 1886, 
        1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 
        1899, 1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 
        1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921, 1922, 
        1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1931, 1932, 1933, 1934, 
        1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 
        1947, 1948, 1949, 1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 
        1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 
        1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 
        1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 
        1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 
        2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 
        2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 
        2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 
        2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2054, 
        2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 
        2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 
        2079, 2080, 2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 
        2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100, 2101, 2102, 
        2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 
        2115, 2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2124, 2125, 2126, 
        2127, 2128, 2129, 2130, 2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138, 
        2139, 2140, 2141, 2142, 2143, 2144, 2145, 2146, 2147, 2148, 2149, 2150, 
        2151, 2152, 2153, 2154, 2155, 2156, 2157, 2158, 2159, 2160, 2161, 2162, 
        2163, 2164, 2165, 2166, 2167, 2168, 2169, 2170, 2171, 2172, 2173, 2174, 
        2175, 2176, 2177, 2178, 2179, 2180, 2181, 2182, 2183, 2184, 2185, 2186, 
        2187, 2188, 2189, 2190, 2191, 2192, 2193, 2194, 2195, 2196, 2197, 2198, 
        2199, 2200, 2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210, 
        2211, 2212, 2213, 2214, 2215, 2216, 2217, 2218, 2219, 2220, 2221, 2222, 
        2223, 2224, 2225, 2226, 2227, 2228, 2229, 2230, 2231, 2232, 2233, 2234, 
        2235, 2236, 2237, 2238, 2239, 2240, 2241, 2242, 2243, 2244, 2245, 2246, 
        2247, 2248, 2249, 2250, 2251, 2252, 2253, 2254, 2255, 2256, 2257, 2258, 
        2259, 2260, 2261, 2262, 2263, 2264, 2265, 2266, 2267, 2268, 2269, 2270, 
        2271, 2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 2280, 2281, 2282, 
        2283, 2284, 2285, 2286, 2287, 2288, 2289, 2290, 2291, 2292, 2293, 2294, 
        2295, 2296, 2297, 2298, 2299, 2300, 2301, 2302, 2303, 2304, 2305, 2306, 
        2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 
        2319, 2320, 2321, 2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329, 2330, 
        2331, 2332, 2333, 2334, 2335, 2336, 2337, 2338, 2339, 2340, 2341, 2342, 
        2343, 2344, 2345, 2346, 2347, 2348, 2349, 2350, 2351, 2352, 2353, 2354, 
        2355, 2356, 2357, 2358, 2359, 2360, 2361, 2362, 2363, 2364, 2365, 2366, 
        2367, 2368, 2369, 2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 
        2379, 2380, 2381, 2382, 2383, 2384, 2385, 2386, 2387, 2388, 2389, 2390, 
        2391, 2392, 2393, 2394, 2395, 2396, 2397, 2398, 2399, 2400, 2401, 2402, 
        2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412, 2413, 2414, 
        2415, 2416, 2417, 2418, 2419, 2420, 2421, 2422, 2423, 2424, 2425, 2426, 
        2427, 2428, 2429, 2430, 2431, 2432, 2433, 2434, 2435, 2436, 2437, 2438, 
        2439, 2440, 2441, 2442, 2443, 2444, 2445, 2446, 2447, 2448, 2449, 2450, 
        2451, 2452, 2453, 2454, 2455, 2456, 2457, 2458, 2459, 2460, 2461, 2462, 
        2463, 2464, 2465, 2466, 2467, 2468, 2469, 2470, 2471, 2472, 2473, 2474, 
        2475, 2476, 2477, 2478, 2479, 2480, 2481, 2482, 2483, 2484, 2485, 2486, 
        2487, 2488, 2489, 2490, 2491, 2492, 2493, 2494, 2495, 2496, 2497, 2498, 
        2499, 2500, 2501, 2502, 2503, 2504, 2505, 2506, 2507, 2508, 2509, 2510, 
        2511, 2512, 2513, 2514, 2515, 2516, 2517, 2518, 2519, 2520, 2521, 2522, 
        2523, 2524, 2525, 2526, 2527, 2528, 2529, 2530, 2531, 2532, 2533, 2534, 
        2535, 2536, 2537, 2538, 2539, 2540, 2541, 2542, 2543, 2544, 2545, 2546, 
        2547, 2548, 2549, 2550, 2551, 2552, 2553, 2554, 2555, 2556, 2557, 2558, 
        2559, 2560, 2561, 2562, 2563, 2564, 2565, 2566, 2567, 2568, 2569, 2570, 
        2571, 2572, 2573, 2574, 2575, 2576, 2577, 2578, 2579, 2580, 2581, 2582, 
        2583, 2584, 2585, 2586, 2587, 2588, 2589, 2590, 2591, 2592, 2593, 2594, 
        2595, 2596, 2597, 2598, 2599, 2600, 2601, 2602, 2603, 2604, 2605, 2606, 
        2607, 2608, 2609, 2610, 2611, 2612, 2613, 2614, 2615, 2616, 2617, 2618, 
        2619, 2620, 2621, 2622, 2623, 2624, 2625, 2626, 2627, 2628, 2629, 2630, 
        2631, 2632, 2633, 2634, 2635, 2636, 2637, 2638, 2639, 2640, 2641, 2642, 
        2643, 2644, 2645, 2646, 2647, 2648, 2649, 2650, 2651, 2652, 2653, 2654, 
        2655, 2656, 2657, 2658, 2659, 2660, 2661, 2662, 2663, 2664, 2665, 2666, 
        2667, 2668, 2669, 2670, 2671, 2672, 2673, 2674, 2675, 2676, 2677, 2678, 
        2679, 2680, 2681, 2682, 2683, 2684, 2685, 2686, 2687, 2688, 2689, 2690, 
        2691, 2692, 2693, 2694, 2695, 2696, 2697, 2698, 2699, 2700, 2701, 2702, 
        2703, 2704, 2705, 2706, 2707, 2708, 2709, 2710, 2711, 2712, 2713, 2714, 
        2715, 2716, 2717, 2718, 2719, 2720, 2721, 2722, 2723, 2724, 2725, 2726, 
        2727, 2728, 2729, 2730, 2731, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 
        2739, 2740, 2741, 2742, 2743, 2744, 2745, 2746, 2747, 2748, 2749, 2750, 
        2751, 2752, 2753, 2754, 2755, 2756, 2757, 2758, 2759, 2760, 2761, 2762, 
        2763, 2764, 2765, 2766, 2767, 2768, 2769, 2770, 2771, 2772, 2773, 2774, 
        2775, 2776, 2777, 2778, 2779, 2780, 2781, 2782, 2783, 2784, 2785, 2786, 
        2787, 2788, 2789, 2790, 2791, 2792, 2793, 2794, 2795, 2796, 2797, 2798, 
        2799, 2800, 2801, 2802, 2803, 2804, 2805, 2806, 2807, 2808, 2809, 2810, 
        2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818, 2819, 2820, 2821, 2822, 
        2823, 2824, 2825, 2826, 2827, 2828, 2829, 2830, 2831, 2832, 2833, 2834, 
        2835, 2836, 2837, 2838, 2839, 2840, 2841, 2842, 2843, 2844, 2845, 2846, 
        2847, 2848, 2849, 2850, 2851, 2852, 2853, 2854, 2855, 2856, 2857, 2858, 
        2859, 2860, 2861, 2862, 2863, 2864, 2865, 2866, 2867, 2868, 2869, 2870, 
        2871, 2872, 2873, 2874, 2875, 2876, 2877, 2878, 2879, 2880, 2881, 2882, 
        2883, 2884, 2885, 2886, 2887, 2888, 2889, 2890, 2891, 2892, 2893, 2894, 
        2895, 2896, 2897, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2906, 
        2907, 2908, 2909, 2910, 2911, 2912, 2913, 2914, 2915, 2916, 2917, 2918, 
        2919, 2920, 2921, 2922, 2923, 2924, 2925, 2926, 2927, 2928, 2929, 2930, 
        2931, 2932, 2933, 2934, 2935, 2936, 2937, 2938, 2939, 2940, 2941, 2942, 
        2943, 2944, 2945, 2946, 2947, 2948, 2949, 2950, 2951, 2952, 2953, 2954, 
        2955, 2956, 2957, 2958, 2959, 2960, 2961, 2962, 2963, 2964, 2965, 2966, 
        2967, 2968, 2969, 2970, 2971, 2972, 2973, 2974, 2975, 2976, 2977, 2978, 
        2979, 2980, 2981, 2982, 2983, 2984, 2985, 2986, 2987, 2988, 2989, 2990, 
        2991, 2992, 2993, 2994, 2995, 2996, 2997, 2998, 2999, 3000, 3001, 3002, 
        3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 
        3015, 3016, 3017, 3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 3026, 
        3027, 3028, 3029, 3030, 3031, 3032, 3033, 3034, 3035, 3036, 3037, 3038, 
        3039, 3040, 3041, 3042, 3043, 3044, 3045, 3046, 3047, 3048, 3049, 3050, 
        3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3060, 3061, 3062, 
        3063, 3064, 3065, 3066, 3067, 3068, 3069, 3070, 3071, 3072, 3073, 3074, 
        3075, 3076, 3077, 3078, 3079, 3080, 3081, 3082, 3083, 3084, 3085, 3086, 
        3087, 3088, 3089, 3090, 3091, 3092, 3093, 3094, 3095, 3096, 3097, 3098, 
        3099, 3100, 3101, 3102, 3103, 3104, 3105, 3106, 3107, 3108, 3109, 3110, 
        3111, 3112, 3113, 3114, 3115, 3116, 3117, 3118, 3119, 3120, 3121, 3122, 
        3123, 3124, 3125, 3126, 3127, 3128, 3129, 3130, 3131, 3132, 3133, 3134, 
        3135, 3136, 3137, 3138, 3139, 3140, 3141, 3142, 3143, 3144, 3145, 3146, 
        3147, 3148, 3149, 3150, 3151, 3152, 3153, 3154, 3155, 3156, 3157, 3158, 
        3159, 3160, 3161, 3162, 3163, 3164, 3165, 3166, 3167, 3168, 3169, 3170, 
        3171, 3172, 3173, 3174, 3175, 3176, 3177, 3178, 3179, 3180, 3181, 3182, 
        3183, 3184, 3185, 3186, 3187, 3188, 3189, 3190, 3191, 3192, 3193, 3194, 
        3195, 3196, 3197, 3198, 3199, 3200, 3201, 3202, 3203, 3204, 3205, 3206, 
        3207, 3208, 3209, 3210, 3211, 3212, 3213, 3214, 3215, 3216, 3217, 3218, 
        3219, 3220, 3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230, 
        3231, 3232, 3233, 3234, 3235, 3236, 3237, 3238, 3239, 3240, 3241, 3242, 
        3243, 3244, 3245, 3246, 3247, 3248, 3249, 3250, 3251, 3252, 3253, 3254, 
        3255, 3256, 3257, 3258, 3259, 3260, 3261, 3262, 3263, 3264, 3265, 3266, 
        3267, 3268, 3269, 3270, 3271, 3272, 3273, 3274, 3275, 3276, 3277, 3278, 
        3279, 3280, 3281, 3282, 3283, 3284, 3285, 3286, 3287, 3288, 3289, 3290, 
        3291, 3292, 3293, 3294, 3295, 3296, 3297, 3298, 3299, 3300, 3301, 3302, 
        3303, 3304, 3305, 3306, 3307, 3308, 3309, 3310, 3311, 3312, 3313, 3314, 
        3315, 3316, 3317, 3318, 3319, 3320, 3321, 3322, 3323, 3324, 3325, 3326, 
        3327, 3328, 3329, 3330, 3331, 3332, 3333, 3334, 3335, 3336, 3337, 3338, 
        3339, 3340, 3341, 3342, 3343, 3344, 3345, 3346, 3347, 3348, 3349, 3350, 
        3351, 3352, 3353, 3354, 3355, 3356, 3357, 3358, 3359, 3360, 3361, 3362, 
        3363, 3364, 3365, 3366, 3367, 3368, 3369, 3370, 3371, 3372, 3373, 3374, 
        3375, 3376, 3377, 3378, 3379, 3380, 3381, 3382, 3383, 3384, 3385, 3386, 
        3387, 3388, 3389, 3390, 3391, 3392, 3393, 3394, 3395, 3396, 3397, 3398, 
        3399, 3400, 3401, 3402, 3403, 3404, 3405, 3406, 3407, 3408, 3409, 3410, 
        3411, 3412, 3413, 3414, 3415, 3416, 3417, 3418, 3419, 3420, 3421, 3422, 
        3423, 3424, 3425, 3426, 3427, 3428, 3429, 3430, 3431, 3432, 3433, 3434, 
        3435, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 3446, 
        3447, 3448, 3449, 3450, 3451, 3452, 3453, 3454, 3455, 3456, 3457, 3458, 
        3459, 3460, 3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3470, 
        3471, 3472, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480, 3481, 3482, 
        3483, 3484, 3485, 3486, 3487, 3488, 3489, 3490, 3491, 3492, 3493, 3494, 
        3495, 3496, 3497, 3498, 3499, 3500, 3501, 3502, 3503, 3504, 3505, 3506, 
        3507, 3508, 3509, 3510, 3511, 3512, 3513, 3514, 3515, 3516, 3517, 3518, 
        3519, 3520, 3521, 3522, 3523, 3524, 3525, 3526, 3527, 3528, 3529, 3530, 
        3531, 3532, 3533, 3534, 3535, 3536, 3537, 3538, 3539, 3540, 3541, 3542, 
        3543, 3544, 3545, 3546, 3547, 3548, 3549, 3550, 3551, 3552, 3553, 3554, 
        3555, 3556, 3557, 3558, 3559, 3560, 3561, 3562, 3563, 3564, 3565, 3566, 
        3567, 3568, 3569, 3570, 3571, 3572, 3573, 3574, 3575, 3576, 3577, 3578, 
        3579, 3580, 3581, 3582, 3583, 3584, 3585, 3586, 3587, 3588, 3589, 3590, 
        3591, 3592, 3593, 3594, 3595, 3596, 3597, 3598, 3599, 3600, 3601, 3602, 
        3603, 3604, 3605, 3606, 3607, 3608, 3609, 3610, 3611, 3612, 3613, 3614, 
        3615, 3616, 3617, 3618, 3619, 3620, 3621, 3622, 3623, 3624, 3625, 3626, 
        3627, 3628, 3629, 3630, 3631, 3632, 3633, 3634, 3635, 3636, 3637, 3638, 
        3639, 3640, 3641, 3642, 3643, 3644, 3645, 3646, 3647, 3648, 3649, 3650, 
        3651, 3652, 3653, 3654, 3655, 3656, 3657, 3658, 3659, 3660, 3661, 3662, 
        3663, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 
        3675, 3676, 3677, 3678, 3679, 3680, 3681, 3682, 3683, 3684, 3685, 3686, 
        3687, 3688, 3689, 3690, 3691, 3692, 3693, 3694, 3695, 3696, 3697, 3698, 
        3699, 3700, 3701, 3702, 3703, 3704, 3705, 3706, 3707, 3708, 3709, 3710, 
        3711, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 
        3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731, 3732, 3733, 3734, 
        3735, 3736, 3737, 3738, 3739, 3740, 3741, 3742, 3743, 3744, 3745, 3746, 
        3747, 3748, 3749, 3750, 3751, 3752, 3753, 3754, 3755, 3756, 3757, 3758, 
        3759, 3760, 3761, 3762, 3763, 3764, 3765, 3766, 3767, 3768, 3769, 3770, 
        3771, 3772, 3773, 3774, 3775, 3776, 3777, 3778, 3779, 3780, 3781, 3782, 
        3783, 3784, 3785, 3786, 3787, 3788, 3789, 3790, 3791, 3792, 3793, 3794, 
        3795, 3796, 3797, 3798, 3799, 3800, 3801, 3802, 3803, 3804, 3805, 3806, 
        3807, 3808, 3809, 3810, 3811, 3812, 3813, 3814, 3815, 3816, 3817, 3818, 
        3819, 3820, 3821, 3822, 3823, 3824, 3825, 3826, 3827, 3828, 3829, 3830, 
        3831, 3832, 3833, 3834, 3835, 3836, 3837, 3838, 3839, 3840, 3841, 3842, 
        3843, 3844, 3845, 3846, 3847, 3848, 3849, 3850, 3851, 3852, 3853, 3854, 
        3855, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3863, 3864, 3865, 3866, 
        3867, 3868, 3869, 3870, 3871, 3872, 3873, 3874, 3875, 3876, 3877, 3878, 
        3879, 3880, 3881, 3882, 3883, 3884, 3885, 3886, 3887, 3888, 3889, 3890, 
        3891, 3892, 3893, 3894, 3895, 3896, 3897, 3898, 3899, 3900, 3901, 3902, 
        3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 
        3915, 3916, 3917, 3918, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3926, 
        3927, 3928, 3929, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 
        3939, 3940, 3941, 3942, 3943, 3944, 3945, 3946, 3947, 3948, 3949, 3950, 
        3951, 3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 
        3963, 3964, 3965, 3966, 3967, 3968, 3969, 3970, 3971, 3972, 3973, 3974, 
        3975, 3976, 3977, 3978, 3979, 3980, 3981, 3982, 3983, 3984, 3985, 3986, 
        3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997, 3998, 
        3999, 4000, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009, 4010, 
        4011, 4012, 4013, 4014, 4015, 4016, 4017, 4018, 4019, 4020, 4021, 4022, 
        4023, 4024, 4025, 4026, 4027, 4028, 4029, 4030, 4031, 4032, 4033, 4034, 
        4035, 4036, 4037, 4038, 4039, 4040, 4041, 4042, 4043, 4044, 4045, 4046, 
        4047, 4048, 4049, 4050, 4051, 4052, 4053, 4054, 4055, 4056, 4057, 4058, 
        4059, 4060, 4061, 4062, 4063, 4064, 4065, 4066, 4067, 4068, 4069, 4070, 
        4071, 4072, 4073, 4074, 4075, 4076, 4077, 4078, 4079, 4080, 4081, 4082, 
        4083, 4084, 4085, 4086, 4087, 4088, 4089, 4090, 4091, 4092, 4093, 4094, 
        4095, 4096, 4097, 4098, 4099, 4100, 4101, 4102, 4103, 4104, 4105, 4106, 
        4107, 4108, 4109, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117, 4118, 
        4119, 4120, 4121, 4122, 4123, 4124, 4125, 4126, 4127, 4128, 4129, 4130, 
        4131, 4132, 4133, 4134, 4135, 4136, 4137, 4138, 4139, 4140, 4141, 4142, 
        4143, 4144, 4145, 4146, 4147, 4148, 4149, 4150, 4151, 4152, 4153, 4154, 
        4155, 4156, 4157, 4158, 4159, 4160, 4161, 4162, 4163, 4164, 4165, 4166, 
        4167, 4168, 4169, 4170, 4171, 4172, 4173, 4174, 4175, 4176, 4177, 4178, 
        4179, 4180, 4181, 4182, 4183, 4184, 4185, 4186, 4187, 4188, 4189, 4190, 
        4191, 4192, 4193, 4194, 4195, 4196, 4197, 4198, 4199, 4200, 4201, 4202, 
        4203, 4204, 4205, 4206, 4207, 4208, 4209, 4210, 4211, 4212, 4213, 4214, 
        4215, 4216, 4217, 4218, 4219, 4220, 4221, 4222, 4223, 4224, 4225, 4226, 
        4227, 4228, 4229, 4230, 4231, 4232, 4233, 4234, 4235, 4236, 4237, 4238, 
        4239, 4240, 4241, 4242, 4243, 4244, 4245, 4246, 4247, 4248, 4249, 4250, 
        4251, 4252, 4253, 4254, 4255, 4256, 4257, 4258, 4259, 4260, 4261, 4262, 
        4263, 4264, 4265, 4266, 4267, 4268, 4269, 4270, 4271, 4272, 4273, 4274, 
        4275, 4276, 4277, 4278, 4279, 4280, 4281, 4282, 4283, 4284, 4285, 4286, 
        4287, 4288, 4289, 4290, 4291, 4292, 4293, 4294, 4295, 4296, 4297, 4298, 
        4299, 4300, 4301, 4302, 4303, 4304, 4305, 4306, 4307, 4308, 4309, 4310, 
        4311, 4312, 4313, 4314, 4315, 4316, 4317, 4318, 4319, 4320, 4321, 4322, 
        4323, 4324, 4325, 4326, 4327, 4328, 4329, 4330, 4331, 4332, 4333, 4334, 
        4335, 4336, 4337, 4338, 4339, 4340, 4341, 4342, 4343, 4344, 4345, 4346, 
        4347, 4348, 4349, 4350, 4351, 4352, 4353, 4354, 4355, 4356, 4357, 4358, 
        4359, 4360, 4361, 4362, 4363, 4364, 4365, 4366, 4367, 4368, 4369, 4370, 
        4371, 4372, 4373, 4374, 4375, 4376, 4377, 4378, 4379, 4380, 4381, 4382, 
        4383, 4384, 4385, 4386, 4387, 4388, 4389, 4390, 4391, 4392, 4393, 4394, 
        4395, 4396, 4397, 4398, 4399, 4400, 4401, 4402, 4403, 4404, 4405, 4406, 
        4407, 4408, 4409, 4410, 4411, 4412, 4413, 4414, 4415, 4416, 4417, 4418, 
        4419, 4420, 4421, 4422, 4423, 4424, 4425, 4426, 4427, 4428, 4429, 4430, 
        4431, 4432, 4433, 4434, 4435, 4436, 4437, 4438, 4439, 4440, 4441, 4442, 
        4443, 4444, 4445, 4446, 4447, 4448, 4449, 4450, 4451, 4452, 4453, 4454, 
        4455, 4456, 4457, 4458, 4459, 4460, 4461, 4462, 4463, 4464, 4465, 4466, 
        4467, 4468, 4469, 4470, 4471, 4472, 4473, 4474, 4475, 4476, 4477, 4478, 
        4479, 4480, 4481, 4482, 4483, 4484, 4485, 4486, 4487, 4488, 4489, 4490, 
        4491, 4492, 4493, 4494, 4495, 4496, 4497, 4498, 4499, 4500, 4501, 4502, 
        4503, 4504, 4505, 4506, 4507, 4508, 4509, 4510, 4511, 4512, 4513, 4514, 
        4515, 4516, 4517, 4518, 4519, 4520, 4521, 4522, 4523, 4524, 4525, 4526, 
        4527, 4528, 4529, 4530, 4531, 4532, 4533, 4534, 4535, 4536, 4537, 4538, 
        4539, 4540, 4541, 4542, 4543, 4544, 4545, 4546, 4547, 4548, 4549, 4550, 
        4551, 4552, 4553, 4554, 4555, 4556, 4557, 4558, 4559, 4560, 4561, 4562, 
        4563, 4564, 4565, 4566, 4567, 4568, 4569, 4570, 4571, 4572, 4573, 4574, 
        4575, 4576, 4577, 4578, 4579, 4580, 4581, 4582, 4583, 4584, 4585, 4586, 
        4587, 4588, 4589, 4590, 4591, 4592, 4593, 4594, 4595, 4596, 4597, 4598, 
        4599, 4600, 4601, 4602, 4603, 4604, 4605, 4606, 4607, 4608, 4609, 4610, 
        4611, 4612, 4613, 4614, 4615, 4616, 4617, 4618, 4619, 4620, 4621, 4622, 
        4623, 4624, 4625, 4626, 4627, 4628, 4629, 4630, 4631, 4632, 4633, 4634, 
        4635, 4636, 4637, 4638, 4639, 4640, 4641, 4642, 4643, 4644, 4645, 4646, 
        4647, 4648, 4649, 4650, 4651, 4652, 4653, 4654, 4655, 4656, 4657, 4658, 
        4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666, 4667, 4668, 4669, 4670, 
        4671, 4672, 4673, 4674, 4675, 4676, 4677, 4678, 4679, 4680, 4681, 4682, 
        4683, 4684, 4685, 4686, 4687, 4688, 4689, 4690, 4691, 4692, 4693, 4694, 
        4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 
        4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 
        4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730, 
        4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 
        4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 
        4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766, 
        4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 
        4779, 4780, 4781, 4782, 4783, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 
        4791, 4792, 4793, 4794, 4795, 4796, 4797, 4798, 4799, 4800, 4801, 4802, 
        4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4813, 4814, 
        4815, 4816, 4817, 4818, 4819, 4820, 4821, 4822, 4823, 4824, 4825, 4826, 
        4827, 4828, 4829, 4830, 4831, 4832, 4833, 4834, 4835, 4836, 4837, 4838, 
        4839, 4840, 4841, 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 
        4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 
        4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 
        4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 
        4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 
        4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 
        4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 
        4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 
        4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 
        4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 
        4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969, 4970, 
        4971, 4972, 4973, 4974, 4975, 4976, 4977, 4978, 4979, 4980, 4981, 4982, 
        4983, 4984, 4985, 4986, 4987, 4988, 4989, 4990, 4991, 4992, 4993, 4994, 
        4995, 4996, 4997, 4998, 4999, 5000, 5001, 5002, 5003, 5004, 5005, 5006, 
        5007, 5008, 5009, 5010, 5011, 5012, 5013, 5014, 5015, 5016, 5017, 5018, 
        5019, 5020, 5021, 5022, 5023, 5024, 5025, 5026, 5027, 5028, 5029, 5030, 
        5031, 5032, 5033, 5034, 5035, 5036, 5037, 5038, 5039, 5040, 5041, 5042, 
        5043, 5044, 5045, 5046, 5047, 5048, 5049, 5050, 5051, 5052, 5053, 5054, 
        5055, 5056, 5057, 5058, 5059, 5060, 5061, 5062, 5063, 5064, 5065, 5066, 
        5067, 5068, 5069, 5070, 5071, 5072, 5073, 5074, 5075, 5076, 5077, 5078, 
        5079, 5080, 5081, 5082, 5083, 5084, 5085, 5086, 5087, 5088, 5089, 5090, 
        5091, 5092, 5093, 5094, 5095, 5096, 5097, 5098, 5099, 5100, 5101, 5102, 
        5103, 5104, 5105, 5106, 5107, 5108, 5109, 5110, 5111, 5112, 5113, 5114, 
        5115, 5116, 5117, 5118, 5119, 5120, 5121, 5122, 5123, 5124, 5125, 5126, 
        5127, 5128, 5129, 5130, 5131, 5132, 5133, 5134, 5135, 5136, 5137, 5138, 
        5139, 5140, 5141, 5142, 5143, 5144, 5145, 5146, 5147, 5148, 5149, 5150, 
        5151, 5152, 5153, 5154, 5155, 5156, 5157, 5158, 5159, 5160, 5161, 5162, 
        5163, 5164, 5165, 5166, 5167, 5168, 5169, 5170, 5171, 5172, 5173, 5174, 
        5175, 5176, 5177, 5178, 5179, 5180, 5181, 5182, 5183, 5184, 5185, 5186, 
        5187, 5188, 5189, 5190, 5191, 5192, 5193, 5194, 5195, 5196, 5197, 5198, 
        5199, 5200, 5201, 5202, 5203, 5204, 5205, 5206, 5207, 5208, 5209, 5210, 
        5211, 5212, 5213, 5214, 5215, 5216, 5217, 5218, 5219, 5220, 5221, 5222, 
        5223, 5224, 5225, 5226, 5227, 5228, 5229, 5230, 5231, 5232, 5233, 5234, 
        5235, 5236, 5237, 5238, 5239, 5240, 5241, 5242, 5243, 5244, 5245, 5246, 
        5247, 5248, 5249, 5250, 5251, 5252, 5253, 5254, 5255, 5256, 5257, 5258, 
        5259, 5260, 5261, 5262, 5263, 5264, 5265, 5266, 5267, 5268, 5269, 5270, 
        5271, 5272, 5273, 5274, 5275, 5276, 5277, 5278, 5279, 5280, 5281, 5282, 
        5283, 5284, 5285, 5286, 5287, 5288, 5289, 5290, 5291, 5292, 5293, 5294, 
        5295, 5296, 5297, 5298, 5299, 5300, 5301, 5302, 5303, 5304, 5305, 5306, 
        5307, 5308, 5309, 5310, 5311, 5312, 5313, 5314, 5315, 5316, 5317, 5318, 
        5319, 5320, 5321, 5322, 5323, 5324, 5325, 5326, 5327, 5328, 5329, 5330, 
        5331, 5332, 5333, 5334, 5335, 5336, 5337, 5338, 5339, 5340, 5341, 5342, 
        5343, 5344, 5345, 5346, 5347, 5348, 5349, 5350, 5351, 5352, 5353, 5354, 
        5355, 5356, 5357, 5358, 5359, 5360, 5361, 5362, 5363, 5364, 5365, 5366, 
        5367, 5368, 5369, 5370, 5371, 5372, 5373, 5374, 5375, 5376, 5377, 5378, 
        5379, 5380, 5381, 5382, 5383, 5384, 5385, 5386, 5387, 5388, 5389, 5390, 
        5391, 5392, 5393, 5394, 5395, 5396, 5397, 5398, 5399, 5400, 5401, 5402, 
        5403, 5404, 5405, 5406, 5407, 5408, 5409, 5410, 5411, 5412, 5413, 5414, 
        5415, 5416, 5417, 5418, 5419, 5420, 5421, 5422, 5423, 5424, 5425, 5426, 
        5427, 5428, 5429, 5430, 5431, 5432, 5433, 5434, 5435, 5436, 5437, 5438, 
        5439, 5440, 5441, 5442, 5443, 5444, 5445, 5446, 5447, 5448, 5449, 5450, 
        5451, 5452, 5453, 5454, 5455, 5456, 5457, 5458, 5459, 5460, 5461, 5462, 
        5463, 5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472, 5473, 5474, 
        5475, 5476, 5477, 5478, 5479, 5480, 5481, 5482, 5483, 5484, 5485, 5486, 
        5487, 5488, 5489, 5490, 5491, 5492, 5493, 5494, 5495, 5496, 5497, 5498, 
        5499, 5500, 5501, 5502, 5503, 5504, 5505, 5506, 5507, 5508, 5509, 5510, 
        5511, 5512, 5513, 5514, 5515, 5516, 5517, 5518, 5519, 5520, 5521, 5522, 
        5523, 5524, 5525, 5526, 5527, 5528, 5529, 5530, 5531, 5532, 5533, 5534, 
        5535, 5536, 5537, 5538, 5539, 5540, 5541, 5542, 5543, 5544, 5545, 5546, 
        5547, 5548, 5549, 5550, 5551, 5552, 5553, 5554, 5555, 5556, 5557, 5558, 
        5559, 5560, 5561, 5562, 5563, 5564, 5565, 5566, 5567, 5568, 5569, 5570, 
        5571, 5572, 5573, 5574, 5575, 5576, 5577, 5578, 5579, 5580, 5581, 5582, 
        5583, 5584, 5585, 5586, 5587, 5588, 5589, 5590, 5591, 5592, 5593, 5594, 
        5595, 5596, 5597, 5598, 5599, 5600, 5601, 5602, 5603, 5604, 5605, 5606, 
        5607, 5608, 5609, 5610, 5611, 5612, 5613, 5614, 5615, 5616, 5617, 5618, 
        5619, 5620, 5621, 5622, 5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 
        5631, 5632, 5633, 5634, 5635, 5636, 5637, 5638, 5639, 5640, 5641, 5642, 
        5643, 5644, 5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5653, 5654, 
        5655, 5656, 5657, 5658, 5659, 5660, 5661, 5662, 5663, 5664, 5665, 5666, 
        5667, 5668, 5669, 5670, 5671, 5672, 5673, 5674, 5675, 5676, 5677, 5678, 
        5679, 5680, 5681, 5682, 5683, 5684, 5685, 5686, 5687, 5688, 5689, 5690, 
        5691, 5692, 5693, 5694, 5695, 5696, 5697, 5698, 5699, 5700, 5701, 5702, 
        5703, 5704, 5705, 5706, 5707, 5708, 5709, 5710, 5711, 5712, 5713, 5714, 
        5715, 5716, 5717, 5718, 5719, 5720, 5721, 5722, 5723, 5724, 5725, 5726, 
        5727, 5728, 5729, 5730, 5731, 5732, 5733, 5734, 5735, 5736, 5737, 5738, 
        5739, 5740, 5741, 5742, 5743, 5744, 5745, 5746, 5747, 5748, 5749, 5750, 
        5751, 5752, 5753, 5754, 5755, 5756, 5757, 5758, 5759, 5760, 5761, 5762, 
        5763, 5764, 5765, 5766, 5767, 5768, 5769, 5770, 5771, 5772, 5773, 5774, 
        5775, 5776, 5777, 5778, 5779, 5780, 5781, 5782, 5783, 5784, 5785, 5786, 
        5787, 5788, 5789, 5790, 5791, 5792, 5793, 5794, 5795, 5796, 5797, 5798, 
        5799, 5800, 5801, 5802, 5803, 5804, 5805, 5806, 5807, 5808, 5809, 5810, 
        5811, 5812, 5813, 5814, 5815, 5816, 5817, 5818, 5819, 5820, 5821, 5822, 
        5823, 5824, 5825, 5826, 5827, 5828, 5829, 5830, 5831, 5832, 5833, 5834, 
        5835, 5836, 5837, 5838, 5839, 5840, 5841, 5842, 5843, 5844, 5845, 5846, 
        5847, 5848, 5849, 5850, 5851, 5852, 5853, 5854, 5855, 5856, 5857, 5858, 
        5859, 5860, 5861, 5862, 5863, 5864, 5865, 5866, 5867, 5868, 5869, 5870, 
        5871, 5872, 5873, 5874, 5875, 5876, 5877, 5878, 5879, 5880, 5881, 5882, 
        5883, 5884, 5885, 5886, 5887, 5888, 5889, 5890, 5891, 5892, 5893, 5894, 
        5895, 5896, 5897, 5898, 5899, 5900, 5901, 5902, 5903, 5904, 5905, 5906, 
        5907, 5908, 5909, 5910, 5911, 5912, 5913, 5914, 5915, 5916, 5917, 5918, 
        5919, 5920, 5921, 5922, 5923, 5924, 5925, 5926, 5927, 5928, 5929, 5930, 
        5931, 5932, 5933, 5934, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 
        5943, 5944, 5945, 5946, 5947, 5948, 5949, 5950, 5951, 5952, 5953, 5954, 
        5955, 5956, 5957, 5958, 5959, 5960, 5961, 5962, 5963, 5964, 5965, 5966, 
        5967, 5968, 5969, 5970, 5971, 5972, 5973, 5974, 5975, 5976, 5977, 5978, 
        5979, 5980, 5981, 5982, 5983, 5984, 5985, 5986, 5987, 5988, 5989, 5990, 
        5991, 5992, 5993, 5994, 5995, 5996, 5997, 5998, 5999, 6000, 6001, 6002, 
        6003, 6004, 6005, 6006, 6007, 6008, 6009, 6010, 6011, 6012, 6013, 6014, 
        6015, 6016, 6017, 6018, 6019, 6020, 6021, 6022, 6023, 6024, 6025, 6026, 
        6027, 6028, 6029, 6030, 6031, 6032, 6033, 6034, 6035, 6036, 6037, 6038, 
        6039, 6040, 6041, 6042, 6043, 6044, 6045, 6046, 6047, 6048, 6049, 6050, 
        6051, 6052, 6053, 6054, 6055, 6056, 6057, 6058, 6059, 6060, 6061, 6062, 
        6063, 6064, 6065, 6066, 6067, 6068, 6069, 6070, 6071, 6072, 6073, 6074, 
        6075, 6076, 6077, 6078, 6079, 6080, 6081, 6082, 6083, 6084, 6085, 6086, 
        6087, 6088, 6089, 6090, 6091, 6092, 6093, 6094, 6095, 6096, 6097, 6098, 
        6099, 6100, 6101, 6102, 6103, 6104, 6105, 6106, 6107, 6108, 6109, 6110, 
        6111, 6112, 6113, 6114, 6115, 6116, 6117, 6118, 6119, 6120, 6121, 6122, 
        6123, 6124, 6125, 6126, 6127, 6128, 6129, 6130, 6131, 6132, 6133, 6134, 
        6135, 6136, 6137, 6138, 6139, 6140, 6141, 6142, 6143, 6144, 6145, 6146, 
        6147, 6148, 6149, 6150, 6151, 6152, 6153, 6154, 6155, 6156, 6157, 6158, 
        6159, 6160, 6161, 6162, 6163, 6164, 6165, 6166, 6167, 6168, 6169, 6170, 
        6171, 6172, 6173, 6174, 6175, 6176, 6177, 6178, 6179, 6180, 6181, 6182, 
        6183, 6184, 6185, 6186, 6187, 6188, 6189, 6190, 6191, 6192, 6193, 6194, 
        6195, 6196, 6197, 6198, 6199, 6200, 6201, 6202, 6203, 6204, 6205, 6206, 
        6207, 6208, 6209, 6210, 6211, 6212, 6213, 6214, 6215, 6216, 6217, 6218, 
        6219, 6220, 6221, 6222, 6223, 6224, 6225, 6226, 6227, 6228, 6229, 6230, 
        6231, 6232, 6233, 6234, 6235, 6236, 6237, 6238, 6239, 6240, 6241, 6242, 
        6243, 6244, 6245, 6246, 6247, 6248, 6249, 6250, 6251, 6252, 6253, 6254, 
        6255, 6256, 6257, 6258, 6259, 6260, 6261, 6262, 6263, 6264, 6265, 6266, 
        6267, 6268, 6269, 6270, 6271, 6272, 6273, 6274, 6275, 6276, 6277, 6278, 
        6279, 6280, 6281, 6282, 6283, 6284, 6285, 6286, 6287, 6288, 6289, 6290, 
        6291, 6292, 6293, 6294, 6295, 6296, 6297, 6298, 6299, 6300, 6301, 6302, 
        6303, 6304, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312, 6313, 6314, 
        6315, 6316, 6317, 6318, 6319, 6320, 6321, 6322, 6323, 6324, 6325, 6326, 
        6327, 6328, 6329, 6330, 6331, 6332, 6333, 6334, 6335, 6336, 6337, 6338, 
        6339, 6340, 6341, 6342, 6343, 6344, 6345, 6346, 6347, 6348, 6349, 6350, 
        6351, 6352, 6353, 6354, 6355, 6356, 6357, 6358, 6359, 6360, 6361, 6362, 
        6363, 6364, 6365, 6366, 6367, 6368, 6369, 6370, 6371, 6372, 6373, 6374, 
        6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382, 6383, 6384, 6385, 6386, 
        6387, 6388, 6389, 6390, 6391, 6392, 6393, 6394, 6395, 6396, 6397, 6398, 
        6399, 6400, 6401, 6402, 6403, 6404, 6405, 6406, 6407, 6408, 6409, 6410, 
        6411, 6412, 6413, 6414, 6415, 6416, 6417, 6418, 6419, 6420, 6421, 6422, 
        6423, 6424, 6425, 6426, 6427, 6428, 6429, 6430, 6431, 6432, 6433, 6434, 
        6435, 6436, 6437, 6438, 6439, 6440, 6441, 6442, 6443, 6444, 6445, 6446, 
        6447, 6448, 6449, 6450, 6451, 6452, 6453, 6454, 6455, 6456, 6457, 6458, 
        6459, 6460, 6461, 6462, 6463, 6464, 6465, 6466, 6467, 6468, 6469, 6470, 
        6471, 6472, 6473, 6474, 6475, 6476, 6477, 6478, 6479, 6480, 6481, 6482, 
        6483, 6484, 6485, 6486, 6487, 6488, 6489, 6490, 6491, 6492, 6493, 6494, 
        6495, 6496, 6497, 6498, 6499, 6500, 6501, 6502, 6503, 6504, 6505, 6506, 
        6507, 6508, 6509, 6510, 6511, 6512, 6513, 6514, 6515, 6516, 6517, 6518, 
        6519, 6520, 6521, 6522, 6523, 6524, 6525, 6526, 6527, 6528, 6529, 6530, 
        6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6539, 6540, 6541, 6542, 
        6543, 6544, 6545, 6546, 6547, 6548, 6549, 6550, 6551, 6552, 6553, 6554, 
        6555, 6556, 6557, 6558, 6559, 6560, 6561, 6562, 6563, 6564, 6565, 6566, 
        6567, 6568, 6569, 6570, 6571, 6572, 6573, 6574, 6575, 6576, 6577, 6578, 
        6579, 6580, 6581, 6582, 6583, 6584, 6585, 6586, 6587, 6588, 6589, 6590, 
        6591, 6592, 6593, 6594, 6595, 6596, 6597, 6598, 6599, 6600, 6601, 6602, 
        6603, 6604, 6605, 6606, 6607, 6608, 6609, 6610, 6611, 6612, 6613, 6614, 
        6615, 6616, 6617, 6618, 6619, 6620, 6621, 6622, 6623, 6624, 6625, 6626, 
        6627, 6628, 6629, 6630, 6631, 6632, 6633, 6634, 6635, 6636, 6637, 6638, 
        6639, 6640, 6641, 6642, 6643, 6644, 6645, 6646, 6647, 6648, 6649, 6650, 
        6651, 6652, 6653, 6654, 6655, 6656, 6657, 6658, 6659, 6660, 6661, 6662, 
        6663, 6664, 6665, 6666, 6667, 6668, 6669, 6670, 6671, 6672, 6673, 6674, 
        6675, 6676, 6677, 6678, 6679, 6680, 6681, 6682, 6683, 6684, 6685, 6686, 
        6687, 6688, 6689, 6690, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 
        6699, 6700, 6701, 6702, 6703, 6704, 6705, 6706, 6707, 6708, 6709, 6710, 
        6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6720, 6721, 6722, 
        6723, 6724, 6725, 6726, 6727, 6728, 6729, 6730, 6731, 6732, 6733, 6734, 
        6735, 6736, 6737, 6738, 6739, 6740, 6741, 6742, 6743, 6744, 6745, 6746, 
        6747, 6748, 6749, 6750, 6751, 6752, 6753, 6754, 6755, 6756, 6757, 6758, 
        6759, 6760, 6761, 6762, 6763, 6764, 6765, 6766, 6767, 6768, 6769, 6770, 
        6771, 6772, 6773, 6774, 6775, 6776, 6777, 6778, 6779, 6780, 6781, 6782, 
        6783, 6784, 6785, 6786, 6787, 6788, 6789, 6790, 6791, 6792, 6793, 6794, 
        6795, 6796, 6797, 6798, 6799, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 
        6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 6816, 6817, 6818, 
        6819, 6820, 6821, 6822, 6823, 6824, 6825, 6826, 6827, 6828, 6829, 6830, 
        6831, 6832, 6833, 6834, 6835, 6836, 6837, 6838, 6839, 6840, 6841, 6842, 
        6843, 6844, 6845, 6846, 6847, 6848, 6849, 6850, 6851, 6852, 6853, 6854, 
        6855, 6856, 6857, 6858, 6859, 6860, 6861, 6862, 6863, 6864, 6865, 6866, 
        6867, 6868, 6869, 6870, 6871, 6872, 6873, 6874, 6875, 6876, 6877, 6878, 
        6879, 6880, 6881, 6882, 6883, 6884, 6885, 6886, 6887, 6888, 6889, 6890, 
        6891, 6892, 6893, 6894, 6895, 6896, 6897, 6898, 6899, 6900, 6901, 6902, 
        6903, 6904, 6905, 6906, 6907, 6908, 6909, 6910, 6911, 6912, 6913, 6914, 
        6915, 6916, 6917, 6918, 6919, 6920, 6921, 6922, 6923, 6924, 6925, 6926, 
        6927, 6928, 6929, 6930, 6931, 6932, 6933, 6934, 6935, 6936, 6937, 6938, 
        6939, 6940, 6941, 6942, 6943, 6944, 6945, 6946, 6947, 6948, 6949, 6950, 
        6951, 6952, 6953, 6954, 6955, 6956, 6957, 6958, 6959, 6960, 6961, 6962, 
        6963, 6964, 6965, 6966, 6967, 6968, 6969, 6970, 6971, 6972, 6973, 6974, 
        6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6985, 6986, 
        6987, 6988, 6989, 6990, 6991, 6992, 6993, 6994, 6995, 6996, 6997, 6998, 
        6999, 7000, 7001, 7002, 7003, 7004, 7005, 7006, 7007, 7008, 7009, 7010, 
        7011, 7012, 7013, 7014, 7015, 7016, 7017, 7018, 7019, 7020, 7021, 7022, 
        7023, 7024, 7025, 7026, 7027, 7028, 7029, 7030, 7031, 7032, 7033, 7034, 
        7035, 7036, 7037, 7038, 7039, 7040, 7041, 7042, 7043, 7044, 7045, 7046, 
        7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7057, 7058, 
        7059, 7060, 7061, 7062, 7063, 7064, 7065, 7066, 7067, 7068, 7069, 7070, 
        7071, 7072, 7073, 7074, 7075, 7076, 7077, 7078, 7079, 7080, 7081, 7082, 
        7083, 7084, 7085, 7086, 7087, 7088, 7089, 7090, 7091, 7092, 7093, 7094, 
        7095, 7096, 7097, 7098, 7099, 7100, 7101, 7102, 7103, 7104, 7105, 7106, 
        7107, 7108, 7109, 7110, 7111, 7112, 7113, 7114, 7115, 7116, 7117, 7118, 
        7119, 7120, 7121, 7122, 7123, 7124, 7125, 7126, 7127, 7128, 7129, 7130, 
        7131, 7132, 7133, 7134, 7135, 7136, 7137, 7138, 7139, 7140, 7141, 7142, 
        7143, 7144, 7145, 7146, 7147, 7148, 7149, 7150, 7151, 7152, 7153, 7154, 
        7155, 7156, 7157, 7158, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 
        7167, 7168, 7169, 7170, 7171, 7172, 7173, 7174, 7175, 7176, 7177, 7178, 
        7179, 7180, 7181, 7182, 7183, 7184, 7185, 7186, 7187, 7188, 7189, 7190, 
        7191, 7192, 7193, 7194, 7195, 7196, 7197, 7198, 7199, 7200, 7201, 7202, 
        7203, 7204, 7205, 7206, 7207, 7208, 7209, 7210, 7211, 7212, 7213, 7214, 
        7215, 7216, 7217, 7218, 7219, 7220, 7221, 7222, 7223, 7224, 7225, 7226, 
        7227, 7228, 7229, 7230, 7231, 7232, 7233, 7234, 7235, 7236, 7237, 7238, 
        7239, 7240, 7241, 7242, 7243, 7244, 7245, 7246, 7247, 7248, 7249, 7250, 
        7251, 7252, 7253, 7254, 7255, 7256, 7257, 7258, 7259, 7260, 7261, 7262, 
        7263, 7264, 7265, 7266, 7267, 7268, 7269, 7270, 7271, 7272, 7273, 7274, 
        7275, 7276, 7277, 7278, 7279, 7280, 7281, 7282, 7283, 7284, 7285, 7286, 
        7287, 7288, 7289, 7290, 7291, 7292, 7293, 7294, 7295, 7296, 7297, 7298, 
        7299, 7300, 7301, 7302, 7303, 7304, 7305, 7306, 7307, 7308, 7309, 7310, 
        7311, 7312, 7313, 7314, 7315, 7316, 7317, 7318, 7319, 7320, 7321, 7322, 
        7323, 7324, 7325, 7326, 7327, 7328, 7329, 7330, 7331, 7332, 7333, 7334, 
        7335, 7336, 7337, 7338, 7339, 7340, 7341, 7342, 7343, 7344, 7345, 7346, 
        7347, 7348, 7349, 7350, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 
        7359, 7360, 7361, 7362, 7363, 7364, 7365, 7366, 7367, 7368, 7369, 7370, 
        7371, 7372, 7373, 7374, 7375, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 
        7383, 7384, 7385, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7394, 
        7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7404, 7405, 7406, 
        7407, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7416, 7417, 7418, 
        7419, 7420, 7421, 7422, 7423, 7424, 7425, 7426, 7427, 7428, 7429, 7430, 
        7431, 7432, 7433, 7434, 7435, 7436, 7437, 7438, 7439, 7440, 7441, 7442, 
        7443, 7444, 7445, 7446, 7447, 7448, 7449, 7450, 7451, 7452, 7453, 7454, 
        7455, 7456, 7457, 7458, 7459, 7460, 7461, 7462, 7463, 7464, 7465, 7466, 
        7467, 7468, 7469, 7470, 7471, 7472, 7473, 7474, 7475, 7476, 7477, 7478, 
        7479, 7480, 7481, 7482, 7483, 7484, 7485, 7486, 7487, 7488, 7489, 7490, 
        7491, 7492, 7493, 7494, 7495, 7496, 7497, 7498, 7499, 7500, 7501, 7502, 
        7503, 7504, 7505, 7506, 7507, 7508, 7509, 7510, 7511, 7512, 7513, 7514, 
        7515, 7516, 7517, 7518, 7519, 7520, 7521, 7522, 7523, 7524, 7525, 7526, 
        7527, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7535, 7536, 7537, 7538, 
        7539, 7540, 7541, 7542, 7543, 7544, 7545, 7546, 7547, 7548, 7549, 7550, 
        7551, 7552, 7553, 7554, 7555, 7556, 7557, 7558, 7559, 7560, 7561, 7562, 
        7563, 7564, 7565, 7566, 7567, 7568, 7569, 7570, 7571, 7572, 7573, 7574, 
        7575, 7576, 7577, 7578, 7579, 7580, 7581, 7582, 7583, 7584, 7585, 7586, 
        7587, 7588, 7589, 7590, 7591, 7592, 7593, 7594, 7595, 7596, 7597, 7598, 
        7599, 7600, 7601, 7602, 7603, 7604, 7605, 7606, 7607, 7608, 7609, 7610, 
        7611, 7612, 7613, 7614, 7615, 7616, 7617, 7618, 7619, 7620, 7621, 7622, 
        7623, 7624, 7625, 7626, 7627, 7628, 7629, 7630, 7631, 7632, 7633, 7634, 
        7635, 7636, 7637, 7638, 7639, 7640, 7641, 7642, 7643, 7644, 7645, 7646, 
        7647, 7648, 7649, 7650, 7651, 7652, 7653, 7654, 7655, 7656, 7657, 7658, 
        7659, 7660, 7661, 7662, 7663, 7664, 7665, 7666, 7667, 7668, 7669, 7670, 
        7671, 7672, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7680, 7681, 7682, 
        7683, 7684, 7685, 7686, 7687, 7688, 7689, 7690, 7691, 7692, 7693, 7694, 
        7695, 7696, 7697, 7698, 7699, 7700, 7701, 7702, 7703, 7704, 7705, 7706, 
        7707, 7708, 7709, 7710, 7711, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 
        7719, 7720, 7721, 7722, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 
        7731, 7732, 7733, 7734, 7735, 7736, 7737, 7738, 7739, 7740, 7741, 7742, 
        7743, 7744, 7745, 7746, 7747, 7748, 7749, 7750, 7751, 7752, 7753, 7754, 
        7755, 7756, 7757, 7758, 7759, 7760, 7761, 7762, 7763, 7764, 7765, 7766, 
        7767, 7768, 7769, 7770, 7771, 7772, 7773, 7774, 7775, 7776, 7777, 7778, 
        7779, 7780, 7781, 7782, 7783, 7784, 7785, 7786, 7787, 7788, 7789, 7790, 
        7791, 7792, 7793, 7794, 7795, 7796, 7797, 7798, 7799, 7800, 7801, 7802, 
        7803, 7804, 7805, 7806, 7807, 7808, 7809, 7810, 7811, 7812, 7813, 7814, 
        7815, 7816, 7817, 7818, 7819, 7820, 7821, 7822, 7823, 7824, 7825, 7826, 
        7827, 7828, 7829, 7830, 7831, 7832, 7833, 7834, 7835, 7836, 7837, 7838, 
        7839, 7840, 7841, 7842, 7843, 7844, 7845, 7846, 7847, 7848, 7849, 7850, 
        7851, 7852, 7853, 7854, 7855, 7856, 7857, 7858, 7859, 7860, 7861, 7862, 
        7863, 7864, 7865, 7866, 7867, 7868, 7869, 7870, 7871, 7872, 7873, 7874, 
        7875, 7876, 7877, 7878, 7879, 7880, 7881, 7882, 7883, 7884, 7885, 7886, 
        7887, 7888, 7889, 7890, 7891, 7892, 7893, 7894, 7895, 7896, 7897, 7898, 
        7899, 7900, 7901, 7902, 7903, 7904, 7905, 7906, 7907, 7908, 7909, 7910, 
        7911, 7912, 7913, 7914, 7915, 7916, 7917, 7918, 7919, 7920, 7921, 7922, 
        7923, 7924, 7925, 7926, 7927, 7928, 7929, 7930, 7931, 7932, 7933, 7934, 
        7935, 7936, 7937, 7938, 7939, 7940, 7941, 7942, 7943, 7944, 7945, 7946, 
        7947, 7948, 7949, 7950, 7951, 7952, 7953, 7954, 7955, 7956, 7957, 7958, 
        7959, 7960, 7961, 7962, 7963, 7964, 7965, 7966, 7967, 7968, 7969, 7970, 
        7971, 7972, 7973, 7974, 7975, 7976, 7977, 7978, 7979, 7980, 7981, 7982, 
        7983, 7984, 7985, 7986, 7987, 7988, 7989, 7990, 7991, 7992, 7993, 7994, 
        7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8003, 8004, 8005, 8006, 
        8007, 8008, 8009, 8010, 8011, 8012, 8013, 8014, 8015, 8016, 8017, 8018, 
        8019, 8020, 8021, 8022, 8023, 8024, 8025, 8026, 8027, 8028, 8029, 8030, 
        8031, 8032, 8033, 8034, 8035, 8036, 8037, 8038, 8039, 8040, 8041, 8042, 
        8043, 8044, 8045, 8046, 8047, 8048, 8049, 8050, 8051, 8052, 8053, 8054, 
        8055, 8056, 8057, 8058, 8059, 8060, 8061, 8062, 8063, 8064, 8065, 8066, 
        8067, 8068, 8069, 8070, 8071, 8072, 8073, 8074, 8075, 8076, 8077, 8078, 
        8079, 8080, 8081, 8082, 8083, 8084, 8085, 8086, 8087, 8088, 8089, 8090, 
        8091, 8092, 8093, 8094, 8095, 8096, 8097, 8098, 8099, 8100, 8101, 8102, 
        8103, 8104, 8105, 8106, 8107, 8108, 8109, 8110, 8111, 8112, 8113, 8114, 
        8115, 8116, 8117, 8118, 8119, 8120, 8121, 8122, 8123, 8124, 8125, 8126, 
        8127, 8128, 8129, 8130, 8131, 8132, 8133, 8134, 8135, 8136, 8137, 8138, 
        8139, 8140, 8141, 8142, 8143, 8144, 8145, 8146, 8147, 8148, 8149, 8150, 
        8151, 8152, 8153, 8154, 8155, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 
        8163, 8164, 8165, 8166, 8167, 8168, 8169, 8170, 8171, 8172, 8173, 8174, 
        8175, 8176, 8177, 8178, 8179, 8180, 8181, 8182, 8183, 8184, 8185, 8186, 
        8187, 8188, 8189, 8190, 8191, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 
        8199, 8200, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 
        8211, 8212, 8213, 8214, 8215, 8216, 8217, 8218, 8219, 8220, 8221, 8222, 
        8223, 8224, 8225, 8226, 8227, 8228, 8229, 8230, 8231, 8232, 8233, 8234, 
        8235, 8236, 8237, 8238, 8239, 8240, 8241, 8242, 8243, 8244, 8245, 8246, 
        8247, 8248, 8249, 8250, 8251, 8252, 8253, 8254, 8255, 8256, 8257, 8258, 
        8259, 8260, 8261, 8262, 8263, 8264, 8265, 8266, 8267, 8268, 8269, 8270, 
        8271, 8272, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8280, 8281, 8282, 
        8283, 8284, 8285, 8286, 8287, 8288, 8289, 8290, 8291, 8292, 8293, 8294, 
        8295, 8296, 8297, 8298, 8299, 8300, 8301, 8302, 8303, 8304, 8305, 8306, 
        8307, 8308, 8309, 8310, 8311, 8312, 8313, 8314, 8315, 8316, 8317, 8318, 
        8319, 8320, 8321, 8322, 8323, 8324, 8325, 8326, 8327, 8328, 8329, 8330, 
        8331, 8332, 8333, 8334, 8335, 8336, 8337, 8338, 8339, 8340, 8341, 8342, 
        8343, 8344, 8345, 8346, 8347, 8348, 8349, 8350, 8351, 8352, 8353, 8354, 
        8355, 8356, 8357, 8358, 8359, 8360, 8361, 8362, 8363, 8364, 8365, 8366, 
        8367, 8368, 8369, 8370, 8371, 8372, 8373, 8374, 8375, 8376, 8377, 8378, 
        8379, 8380, 8381, 8382, 8383, 8384, 8385, 8386, 8387, 8388, 8389, 8390, 
        8391, 8392, 8393, 8394, 8395, 8396, 8397, 8398, 8399, 8400, 8401, 8402, 
        8403, 8404, 8405, 8406, 8407, 8408, 8409, 8410, 8411, 8412, 8413, 8414, 
        8415, 8416, 8417, 8418, 8419, 8420, 8421, 8422, 8423, 8424, 8425, 8426, 
        8427, 8428, 8429, 8430, 8431, 8432, 8433, 8434, 8435, 8436, 8437, 8438, 
        8439, 8440, 8441, 8442, 8443, 8444, 8445, 8446, 8447, 8448, 8449, 8450, 
        8451, 8452, 8453, 8454, 8455, 8456, 8457, 8458, 8459, 8460, 8461, 8462, 
        8463, 8464, 8465, 8466, 8467, 8468, 8469, 8470, 8471, 8472, 8473, 8474, 
        8475, 8476, 8477, 8478, 8479, 8480, 8481, 8482, 8483, 8484, 8485, 8486, 
        8487, 8488, 8489, 8490, 8491, 8492, 8493, 8494, 8495, 8496, 8497, 8498, 
        8499, 8500, 8501, 8502, 8503, 8504, 8505, 8506, 8507, 8508, 8509, 8510, 
        8511, 8512, 8513, 8514, 8515, 8516, 8517, 8518, 8519, 8520, 8521, 8522, 
        8523, 8524, 8525, 8526, 8527, 8528, 8529, 8530, 8531, 8532, 8533, 8534, 
        8535, 8536, 8537, 8538, 8539, 8540, 8541, 8542, 8543, 8544, 8545, 8546, 
        8547, 8548, 8549, 8550, 8551, 8552, 8553, 8554, 8555, 8556, 8557, 8558, 
        8559, 8560, 8561, 8562, 8563, 8564, 8565, 8566, 8567, 8568, 8569, 8570, 
        8571, 8572, 8573, 8574, 8575, 8576, 8577, 8578, 8579, 8580, 8581, 8582, 
        8583, 8584, 8585, 8586, 8587, 8588, 8589, 8590, 8591, 8592, 8593, 8594, 
        8595, 8596, 8597, 8598, 8599, 8600, 8601, 8602, 8603, 8604, 8605, 8606, 
        8607, 8608, 8609, 8610, 8611, 8612, 8613, 8614, 8615, 8616, 8617, 8618, 
        8619, 8620, 8621, 8622, 8623, 8624, 8625, 8626, 8627, 8628, 8629, 8630, 
        8631, 8632, 8633, 8634, 8635, 8636, 8637, 8638, 8639, 8640, 8641, 8642, 
        8643, 8644, 8645, 8646, 8647, 8648, 8649, 8650, 8651, 8652, 8653, 8654, 
        8655, 8656, 8657, 8658, 8659, 8660, 8661, 8662, 8663, 8664, 8665, 8666, 
        8667, 8668, 8669, 8670, 8671, 8672, 8673, 8674, 8675, 8676, 8677, 8678, 
        8679, 8680, 8681, 8682, 8683, 8684, 8685, 8686, 8687, 8688, 8689, 8690, 
        8691, 8692, 8693, 8694, 8695, 8696, 8697, 8698, 8699, 8700, 8701, 8702, 
        8703, 8704, 8705, 8706, 8707, 8708, 8709, 8710, 8711, 8712, 8713, 8714, 
        8715, 8716, 8717, 8718, 8719, 8720, 8721, 8722, 8723, 8724, 8725, 8726, 
        8727, 8728, 8729, 8730, 8731, 8732, 8733, 8734, 8735, 8736, 8737, 8738, 
        8739, 8740, 8741, 8742, 8743, 8744, 8745, 8746, 8747, 8748, 8749, 8750, 
        8751, 8752, 8753, 8754, 8755, 8756, 8757, 8758, 8759, 8760, 8761, 8762, 
        8763, 8764, 8765, 8766, 8767, 8768, 8769, 8770, 8771, 8772, 8773, 8774, 
        8775, 8776, 8777, 8778, 8779, 8780, 8781, 8782, 8783, 8784, 8785, 8786, 
        8787, 8788, 8789, 8790, 8791, 8792, 8793, 8794, 8795, 8796, 8797, 8798, 
        8799, 8800, 8801, 8802, 8803, 8804, 8805, 8806, 8807, 8808, 8809, 8810, 
        8811, 8812, 8813, 8814, 8815, 8816, 8817, 8818, 8819, 8820, 8821, 8822, 
        8823, 8824, 8825, 8826, 8827, 8828, 8829, 8830, 8831, 8832, 8833, 8834, 
        8835, 8836, 8837, 8838, 8839, 8840, 8841, 8842, 8843, 8844, 8845, 8846, 
        8847, 8848, 8849, 8850, 8851, 8852, 8853, 8854, 8855, 8856, 8857, 8858, 
        8859, 8860, 8861, 8862, 8863, 8864, 8865, 8866, 8867, 8868, 8869, 8870, 
        8871, 8872, 8873, 8874, 8875, 8876, 8877, 8878, 8879, 8880, 8881, 8882, 
        8883, 8884, 8885, 8886, 8887, 8888, 8889, 8890, 8891, 8892, 8893, 8894, 
        8895, 8896, 8897, 8898, 8899, 8900, 8901, 8902, 8903, 8904, 8905, 8906, 
        8907, 8908, 8909, 8910, 8911, 8912, 8913, 8914, 8915, 8916, 8917, 8918, 
        8919, 8920, 8921, 8922, 8923, 8924, 8925, 8926, 8927, 8928, 8929, 8930, 
        8931, 8932, 8933, 8934, 8935, 8936, 8937, 8938, 8939, 8940, 8941, 8942, 
        8943, 8944, 8945, 8946, 8947, 8948, 8949, 8950, 8951, 8952, 8953, 8954, 
        8955, 8956, 8957, 8958, 8959, 8960, 8961, 8962, 8963, 8964, 8965, 8966, 
        8967, 8968, 8969, 8970, 8971, 8972, 8973, 8974, 8975, 8976, 8977, 8978, 
        8979, 8980, 8981, 8982, 8983, 8984, 8985, 8986, 8987, 8988, 8989, 8990, 
        8991, 8992, 8993, 8994, 8995, 8996, 8997, 8998, 8999, 9000, 9001, 9002, 
        9003, 9004, 9005, 9006, 9007, 9008, 9009, 9010, 9011, 9012, 9013, 9014, 
        9015, 9016, 9017, 9018, 9019, 9020, 9021, 9022, 9023, 9024, 9025, 9026, 
        9027, 9028, 9029, 9030, 9031, 9032, 9033, 9034, 9035, 9036, 9037, 9038, 
        9039, 9040, 9041, 9042, 9043, 9044, 9045, 9046, 9047, 9048, 9049, 9050, 
        9051, 9052, 9053, 9054, 9055, 9056, 9057, 9058, 9059, 9060, 9061, 9062, 
        9063, 9064, 9065, 9066, 9067, 9068, 9069, 9070, 9071, 9072, 9073, 9074, 
        9075, 9076, 9077, 9078, 9079, 9080, 9081, 9082, 9083, 9084, 9085, 9086, 
        9087, 9088, 9089, 9090, 9091, 9092, 9093, 9094, 9095, 9096, 9097, 9098, 
        9099, 9100, 9101, 9102, 9103, 9104, 9105, 9106, 9107, 9108, 9109, 9110, 
        9111, 9112, 9113, 9114, 9115, 9116, 9117, 9118, 9119, 9120, 9121, 9122, 
        9123, 9124, 9125, 9126, 9127, 9128, 9129, 9130, 9131, 9132, 9133, 9134, 
        9135, 9136, 9137, 9138, 9139, 9140, 9141, 9142, 9143, 9144, 9145, 9146, 
        9147, 9148, 9149, 9150, 9151, 9152, 9153, 9154, 9155, 9156, 9157, 9158, 
        9159, 9160, 9161, 9162, 9163, 9164, 9165, 9166, 9167, 9168, 9169, 9170, 
        9171, 9172, 9173, 9174, 9175, 9176, 9177, 9178, 9179, 9180, 9181, 9182, 
        9183, 9184, 9185, 9186, 9187, 9188, 9189, 9190, 9191, 9192, 9193, 9194, 
        9195, 9196, 9197, 9198, 9199, 9200, 9201, 9202, 9203, 9204, 9205, 9206, 
        9207, 9208, 9209, 9210, 9211, 9212, 9213, 9214, 9215, 9216, 9217, 9218, 
        9219, 9220, 9221, 9222, 9223, 9224, 9225, 9226, 9227, 9228, 9229, 9230, 
        9231, 9232, 9233, 9234, 9235, 9236, 9237, 9238, 9239, 9240, 19763, 
        19764, 19765, 19766, 19767, 19768, 19769, 19770, 19771, 19772, 19773, 
        )), ))
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1)
    xyp = session.XYPlot('XYPlot-1')
    chartName = xyp.charts.keys()[0]
    chart = xyp.charts[chartName]
    pth = session.paths['Path-1']
    xy1 = xyPlot.XYDataFromPath(path=pth, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE, 
        removeDuplicateXYPairs=True, includeAllElements=False)
    c1 = session.Curve(xyData=xy1)
    chart.setValues(curvesToPlot=(c1, ), )
    session.charts[chartName].autoColor(lines=True, symbols=True)
    session.viewports['Viewport: 1'].setValues(displayedObject=xyp)



def CloseODBforMemory():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.odbs['F:/Rhys/AbaqusModels/BoundaryHold_Automate/Mesh17-Final_ID_0p05_IT_0p5_IR_0p5_IE_10x_LT_10_LE_25kPa.odb'].close(
        )


def Macro2():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Mesh17-Final'].parts['Layer']
    s = p.features['Partition face-1'].sketch
    mdb.models['Mesh17-Final'].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models['Mesh17-Final'].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    s2.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s2, 
        upToFeature=p.features['Partition face-1'], filter=COPLANAR_EDGES)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.312324, 
        farPlane=0.318976, width=0.0115603, height=0.00707018, cameraPosition=(
        0.0115965, 0.00438046, 0.31565), cameraTarget=(0.0115965, 0.00438046, 
        0))
    d[5].setValues(value=0.0005, )
    s2.unsetPrimaryObject()
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models['Mesh17-Final'].sketches['__edit__']
    p = mdb.models['Mesh17-Final'].parts['Layer']
    p.regenerate()


def Macro3():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.Job(name='Job-2', model='Mesh17-Final', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=5, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=False, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=2, 
        numDomains=2, numGPUs=0)




def Amplitude():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.models['Change_LE'].TabularAmplitude(name='Amp-2', timeSpan=STEP, 
        smooth=SOLVER_DEFAULT, data=((0.0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.3, 
        0.8), (1.0, 1.0)))


def ampchange():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.models['Change_LE'].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, 
        smooth=SOLVER_DEFAULT, data=((0.0, 0.1), (0.2, 0.3), (0.6, 0.8), (1.0, 
        1.0)))


def ChangeHeldPoints():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p1 = mdb.models['Change_LE'].parts['Layer']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Change_LE'].parts['Layer']
    p.features['Datum plane-1'].setValues(offset=0.105)
    p = mdb.models['Change_LE'].parts['Layer']
    p.regenerate()
    p1 = mdb.models['Change_LE'].parts['Sample']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Change_LE'].parts['Sample']
    p.features['Datum plane-1'].setValues(offset=0.065)
    p = mdb.models['Change_LE'].parts['Sample']
    p.regenerate()


def ChangeStartIncrement():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF, adaptiveMeshConstraints=ON)
    mdb.models['Change_LE'].steps['Compress'].setValues(initialInc=0.001)


def Change_Realistic():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy
    import itertools
    import os, time
    import threading
    import signal
    from abaqus import getWarningReply, YES, NO
    # ENTER ALL VALUES AS DOUBLES, NO INTEGERS ALLOWED
    Preload = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3] # Strain in the layer
    PreloadN = ['2p5%','5%','7p5%','10%','12p5%','15%','17p5%','20%','22p5%','25%','27p5%','30%']
    #LayerElasticity = 25000.0*1.0
    #LayerElasticityN = '25kpa'
    #BackgroundElasticity = 25000.0
    InclusionElasticity = 25000.0*10.0 # Doesn't converge after 7.5x background
    InclusionElasticityN = '10x'
    LayerThickness = 0.01
    LayerThicknessN = '10'
    SampleThickness = 0.005
    InclusionDepth = 0.001 # Minimum 60e-06
    InclusionDepthN = '1'
    InclusionThickness = 0.003 # Minimum 10e-06, NOT allowed to have inclusion touching bottom of sample, need at least 10e-06 clearance
    InclusionThicknessN = '3'
    InclusionRadius = 0.003 # Minimum 10e-06
    InclusionRadiusN = '3'
    OutputNames = ['LBU2','P1U2','P2U2','LTU2','LBS2']
    ModelName = 'Realistic'
    Directory = 'F:/Rhys/AbaqusModels/BoundaryHold_Investigations/'
    
    # Change Sample
    p = mdb.models[ModelName].parts['Sample']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    if SampleThickness-InclusionDepth-InclusionThickness == 0: # If Inclusion touches the bottom of the sample
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['BackgroundL']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['InclusionL']) # Inclusion section node set
        d[2].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[0].setValues(value=50e-06, ) # Small multiple of mesh size
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    else:
        mdb.models[ModelName].parts['Sample'].sectionAssignments[0].setValues(region=p.sets['Background']) # background section node set
        mdb.models[ModelName].parts['Sample'].sectionAssignments[1].setValues(region=p.sets['Inclusion']) # Inclusion section node set
        d[2].setValues(value=SampleThickness-InclusionDepth-InclusionThickness, ) # Distance from bottom of sample to bottom of inclusion
        d[0].setValues(value=InclusionThickness, ) # Inclusion Thickness
        d[1].setValues(value=InclusionRadius, ) # Inclusion Radius
    p.features['Partition face-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    p.features['Datum plane-1'].setValues(offset=0.05) # Change sample hold point
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Inclusion'].elastic.setValues(table=((InclusionElasticity, 0.45), )) # Set Inclusion Material Properties
    mdb.models[ModelName].materials['Background'].elastic.setValues(table=((BackgroundElasticity, 0.45), )) # Set Background Material Properties

    # Change Layer
    p = mdb.models[ModelName].parts['Layer']
    s = p.features['Shell planar-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s1 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    d[6].setValues(value=LayerThickness, ) # Layer Thickness
    p.features['Shell planar-1'].setValues(sketch=s1)
    del mdb.models[ModelName].sketches['__edit__']
    s = p.features['Partition face-1'].sketch
    mdb.models[ModelName].ConstrainedSketch(name='__edit__', objectToCopy=s)
    s2 = mdb.models[ModelName].sketches['__edit__']
    g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
    if LayerThickness > 0.00125:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.0005, ) # first to second particle layer
    elif LayerThickness > 0.00076:
        d[0].setValues(value=0.0005, ) # bottom to first particle layer
        d[1].setValues(value=0.00001, ) # first to second particle layer
    else:
        d[0].setValues(value=0.0001, ) # bottom to first particle layer
        d[1].setValues(value=0.0001, ) # first to second particle layer
    if LayerThickness == 0.004:
        d[5].setValues(value=0.0006, )
    else:
        d[5].setValues(value=0.0005, )
    p.features['Partition face-1'].setValues(sketch=s2)
    del mdb.models[ModelName].sketches['__edit__']
    if LayerElasticity >= BackgroundElasticity:
        p.features['Datum plane-1'].setValues(offset=0.145) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.005)
    else:
        p.features['Datum plane-1'].setValues(offset=0.055) # Change layer hold point
        mdb.models[ModelName].steps['Compress'].setValues(initialInc=0.0001) # If layer is soft smaller increments needed
    p.regenerate()
    p.generateMesh()
    mdb.models[ModelName].materials['Layer'].elastic.setValues(table=((LayerElasticity, 0.45), )) # Set Layer Material Properties
    
    # Regenerate Assembly
    a = mdb.models[ModelName].rootAssembly
    a.regenerate()
    # Boundary Conditions, gradient of position with time equation derived in notebook
    # B_t = -Preload[-1]*SampleThickness*LayerElasticity/BackgroundElasticity
    # T_t = -Preload[-1]*(LayerThickness+SampleThickness*LayerElasticity/BackgroundElasticity)
    # Boundary Conditions, gradient of position with time equation derived in notebook, for sample strain controlled
    B_t = -Preload[-1]*SampleThickness
    #T_t = -Preload[-1]*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-Preload[-1])
    #mdb.models[ModelName].boundaryConditions['LayerTop'].setValuesInStep(stepName='Compress', u2=T_t)
    mdb.models[ModelName].boundaryConditions['LayerBottomFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    mdb.models[ModelName].boundaryConditions['SampleTopFarEdge'].setValuesInStep(stepName='Compress', u2=B_t)
    # Output Times Calc and Set
    OutputTimes = (numpy.array(Preload)/Preload[-1]).tolist()
    # Amplitude Tupple
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*BackgroundElasticity/LayerElasticity*numpy.log(1-numpy.array(Preload))
    #Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*(1-numpy.array(Preload))**(BackgroundElasticity/LayerElasticity) - LayerThickness
    Delta_T = -numpy.array(Preload)*SampleThickness + LayerThickness*numpy.power((1.0-numpy.array(Preload)),(BackgroundElasticity/LayerElasticity)) - LayerThickness
    AmplitudeTupple = ()
    AmplitudeTupple = AmplitudeTupple + ((0,0,),)
    for i in range(0,len(OutputTimes)):
        AmplitudeTupple = AmplitudeTupple + ((OutputTimes[i],Delta_T[i],),)
    mdb.models[ModelName].amplitudes['Amp_LayerTop'].setValues(timeSpan=STEP, smooth=SOLVER_DEFAULT, data=AmplitudeTupple)
    #Old #OutputTimes = (numpy.array(Preload)*(Layers[n2][0]*BackgroundElasticity+SampleThickness*Layers[n2][1])/(Preload[-1]*(Layers[n2][0]+SampleThickness)*BackgroundElasticity)).tolist()
    OutputTimesTupple = ()
    for i in range(0,len(OutputTimes)):
        OutputTimesTupple = OutputTimesTupple + ((OutputTimes[i],),)
    mdb.models[ModelName].TimePoint(name='TimePoints-1', points=OutputTimesTupple)
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-1')
        
    # Make Job
#    JobName = 'Submit_One_NewMesh_2ndNewLE_'+ModelName+'_ID_'+InclusionDepthN+'_IT_'+InclusionThicknessN+'_IR_'+InclusionRadiusN+'_IE_'+InclusionElasticityN+'_LT_'+LayerThicknessN+'_LE_'+LayerElasticityN
#    j = mdb.Job(name=JobName, 
#        model=ModelName, 
#        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
#        queue=None, memory=50, memoryUnits=PERCENTAGE, 
#        getMemoryFromAnalysis=False, explicitPrecision=SINGLE, 
#        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
#        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
#        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
#        numDomains=4, numGPUs=0)
#    j.submit(consistencyChecking=OFF)

def Individual_node_boundary_condition():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.300225, 
        farPlane=0.300255, width=0.000148209, height=0.0001, 
        viewOffsetX=-0.0749312, viewOffsetY=0.00200283)
    a = mdb.models['Submit_One_Mesh17-Final_ID_0p5_IT_3_IR'].rootAssembly
    n1 = a.instances['SAMPLE-1'].nodes
    nodes1 = n1.getSequenceFromMask(mask=('[#0:195 #400 ]', ), )
    region = regionToolset.Region(nodes=nodes1)
    mdb.models['Submit_One_Mesh17-Final_ID_0p5_IT_3_IR'].DisplacementBC(
        name='Ind_Node', createStepName='Compress', region=region, u1=UNSET, 
        u2=0.0005, ur3=UNSET, amplitude='LINEAR', fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)
    a = mdb.models['Submit_One_Mesh17-Final_ID_0p5_IT_3_IR'].rootAssembly
    n1 = a.instances['LAYER-1'].nodes
    nodes1 = n1.getSequenceFromMask(mask=('[#0:55 #200000 ]', ), )
    region = regionToolset.Region(nodes=nodes1)
    mdb.models['Submit_One_Mesh17-Final_ID_0p5_IT_3_IR'].DisplacementBC(
        name='Nex_Node', createStepName='Compress', region=region, u1=UNSET, 
        u2=0.0004, ur3=UNSET, amplitude='LINEAR', fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)



