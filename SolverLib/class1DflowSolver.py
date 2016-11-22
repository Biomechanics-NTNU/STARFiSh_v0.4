import sys,os
import numpy as np

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
topFolder = cur + "/../"
#topFolder should point to where LOGcurrentWaveSpeed.txt is
# located, and if it changes location, renaming it is appropriate.

sys.path.append(os.path.join(cur,'..'))
import logging
logger = logging.getLogger('starfish')

import UtilityLib.classStarfishBaseObject as cSBO

import UtilityLib.progressBar as cPB

from NetworkLib.classBoundaryConditions import VaryingElastance, Valve
from NetworkLib.classVascularNetwork import VascularNetwork

from classBoundarys import Boundary

import classSystemEquations
import classConnections
import classFields
from classCommunicators import *
from UtilityLib import classRuntimeMemoryManager
from UtilityLib.moduleHelperFunctions import memoryUsagePsutil

import gc

class FlowSolver(cSBO.StarfishBaseObject):
    def __init__(self,vascularNetwork, quiet=False):
        """
        Constructor
        """

        if vascularNetwork == None: raise ValueError("ERROR: No vascularNetwork given!")
        assert isinstance(vascularNetwork, VascularNetwork) 
        # the vascular network to solve                    
        self.vascularNetwork = vascularNetwork
        self.vascularNetwork.quiet = quiet
        self.quiet = quiet

        self.vessels = self.vascularNetwork.vessels
        self.fields = {}
        self.dataHandler = None
        # the boundarys of the network { vesselID : [<instance>::classBoundary_02(Characteristics.py), .. ]}
        # 1 boundary for each start/end-vessel except if only 1 vessel in the network
        self.boundarys = {}
        # the system Equations of each vessel { vesselID : <instance>::classSystemEquations(SystemEquations) }
        self.systemEquations = {}
        # the connections of the Network { motherVesselID : <instance>::classConnections, ...}
        self.connections  = {}
        # the communicator objects of the vascularNetwork {communicatorID : <instance>::classCommunicator}
        self.communicators = {}

        self.venousPool = self.vascularNetwork.venousPool
        #

        # list of numerical objects (field,connection,boundary objects as in the traversing list)
        self.numericalObjects = []
        # time step
        self.dt = None
        # number of Timesteps
        self.nTSteps = None

        # total Simulation time
        self.totalTime = None

        # timestep Counter
        self.currentTimeStep = [0]

        self.cycleMode = False

        # Initialize indices to track where in memory the current solution is stored
        self.memoryOffset = [0]
        self.currentMemoryIndex = [0]

        # Set div output
        self.output = {}

        self.simplifyEigenvalues = self.vascularNetwork.simplifyEigenvalues

        self.riemannInvariantUnitBase = self.vascularNetwork.riemannInvariantUnitBase

        self.solvingSchemeConnections = self.vascularNetwork.solvingSchemeConnections

        # bool for cfl meshing with automatic grid adaption
        self.automaticGridAdaptation = self.vascularNetwork.automaticGridAdaptation

        # rigidAreas True: A = A0 False: A = A(P)
        self.rigidAreas = self.vascularNetwork.rigidAreas

        #define solve function
        self.solve = self.MacCormack_Field

        # initialize system
        self.vascularNetwork.currentMemoryIndex = self.currentMemoryIndex
        self.vascularNetwork.initialize(initializeForSimulation = True)

        self.initializeTimeVariables(quiet)

        self.vascularNetwork.initializeNetworkForSimulation()

        self.initializeSystemEquations()
        self.initializeBoundarys()
        self.initializeConnections()
        self.initializeFields()
        
        if hasattr(self.venousPool, "initializeWithFlowSolver"):
            self.venousPool.initializeWithFlowSolver(self)
        
        self.initializeCommunicators()
        self.initializeNumericalObjectList()

        if self.quiet==False:
            self.initOutput() # feedback

    """
    ########################################################################################
    # initialisation Methods
    ########################################################################################
    """
    def calcTimeStep(self,dz,c,CFL):
        return (CFL*dz)/c

    def initializeTimeVariables(self, quiet):
        """
        initialize time variable dt and Tstep
        """
        self.totalTime = self.vascularNetwork.totalTime

        initialValues = self.vascularNetwork.initialValues

        dt_min,dz_min,c_max,gridNodens = [],[],[],[]
        #create waveSpeed Log file
        logfileData = {}

        dt = self.totalTime
        for vessel in self.vessels.itervalues():
        # Calculate time variables
            #estimate initial pressure
            p0,p1 = initialValues[vessel.Id]['Pressure']

            initialPressure = np.linspace(p0,p1,vessel.N)

            A0_max = max(vessel.A(initialPressure))
            #c_high = vessel.c(A0_max,vessel.initialPressure)

            #c_high1 = vessel.c(A0_max,0)
            Compliance = vessel.C(initialPressure)
            c_high = vessel.c(A0_max, Compliance)

            dz_low = min(vessel.dz)
            dt = self.calcTimeStep(dz_low,c_high,self.vascularNetwork.CFL)
            c_max = np.append(c_max,c_high)
            dt_min = np.append(dt_min,dt)
            dz_min = np.append(dz_min,dz_low)
            gridNodens = np.append(gridNodens,vessel.N)

            logfileData[vessel.Id] = [max(c_high),min(c_high),min(dt),vessel.dz,vessel.N]

        # Set time variables
        self.dt = min(dt_min)
        # calculate time steps
        self.nTSteps = int(np.ceil(self.totalTime/self.dt))
        # calculate time steps for initialisation phase
        nTstepsInitPhase = 0
        if self.vascularNetwork.initialisationPhaseExist:
            initPhaseTimeSpan = self.vascularNetwork.initPhaseTimeSpan
            nTstepsInitPhase = int(np.ceil(initPhaseTimeSpan/self.dt))
        # correct time steps
        self.nTSteps += nTstepsInitPhase

        # update vascular network variables
        self.vascularNetwork.update({'dt':self.dt,
                                     'nTSteps': self.nTSteps,
                                     'nTstepsInitPhase': nTstepsInitPhase})

        self.output['dz_min']     = min(dz_min)
        self.output['c_min']      = min(c_max)
        self.output['c_max']      = max(c_max)
        self.output['gridNodens'] = sum(gridNodens)

        self.output['CFLcorrect'] = []

        ###########
        #### grid correction methods

        automaticGridCorrection = {}

        logfile = open(str(topFolder+'LOGcurrentWaveSpeed.txt'),'wb')
        logfile2 = open(str(topFolder+'LOGproposedGrid.txt'),'wb')
        CFL = self.vascularNetwork.CFL
        for vesselT,data in logfileData.iteritems():
            #number of deltaX
            Nnew = int((sum(data[3])*CFL/(self.dt*data[0])))

            L = sum(data[3])
            #calculate number of dx: M = L/dx
            Mnew = int(L*CFL/self.dt/data[0])

            dz_new = L/Mnew
            dt_new = self.calcTimeStep(dz_new,data[0],CFL)

            while dt_new < self.dt:
                Mnew = Mnew-1
                dz_new = L/Mnew
                dt_new = self.calcTimeStep(dz_new,data[0],CFL)

            # calculate gridpoints N = N+1
            Nnew = Mnew+1
            if Nnew == data[4]:
                logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || already best N', '\n']))
                #self.output['CFLcorrect'].append(' '.join([str(vesselT).rjust(2).ljust(15),' no correction']))
            else:
                logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || prop: N%2.0f'%Nnew,'    dN %.0f'%(Nnew-data[4]),'    dtNew %.4f'%(dt_new*1.e3),'   CFL: %2.3f'%(data[0]*self.dt/dz_new), '\n']))
                self.output['CFLcorrect'].append(' '.join([str(vesselT).ljust(3),'|',str(int(data[4])).ljust(3),'->',str(Nnew).ljust(3),'|', '%2.3f'%(data[0]*self.dt/min(data[3])),'->', '%2.3f'%(data[0]*self.dt/dz_new),'|']))
                automaticGridCorrection[vesselT] = int(Nnew)
            logfile2.write(''.join([str(int(Nnew)),'\n']))
        logfile.close()
        logfile2.close()

        if self.output['CFLcorrect'] != []:
            if quiet == False:
                logger.info('=====================================')
                logger.info('___CFL-correction: grid adaptation___')
                logger.info('Id  | gridPoints |       CFL      |')
                logger.info('    | now -> new |   now -> new   |')
                for CFLcorr in self.output['CFLcorrect']:
                    logger.info('%s' % (CFLcorr))

            if automaticGridCorrection != {}:
                gridCorrection = 'ohYesDoItPlease'
                
                # TODO: comment grid correction in again !!
                #if self.automaticGridAdaptation == True: gridCorrection = 'y'
                #while  gridCorrection not in (' ','','y','Y','n','N'): 
                #    gridCorrection = raw_input('Do you whish to adapt grid? (yes [<ENTER>,<y>,<Y>]/no [<n>,<N>])')

                if gridCorrection in (' ','','y','Y'):
                    #if quiet == False: logger.info(' proceed with: grid aptation for vessels {} \n'.format(automaticGridCorrection.keys()))
                    arrayN = 0
                    for vesselId,Nnew in automaticGridCorrection.iteritems():
                        self.vessels[vesselId].update({'N':Nnew})
                        self.vessels[vesselId].initialize({})

                        newOutput = self.output['CFLcorrect'][arrayN].split(' no')[0]
                        self.output['CFLcorrect'][arrayN] = ' '.join([newOutput,'yes'])
                        arrayN = arrayN+1
        gridNodes = 0
        for vessel in self.vascularNetwork.vessels.itervalues():
            gridNodes += vessel.N

        self.output['gridNodens'] = int(gridNodes)

    def initializeSystemEquations(self):
        """
        initialize system Equations
        """
        for vesselId,vessel in self.vessels.iteritems():
            self.systemEquations[vesselId] = classSystemEquations.System(vessel,
                                                    self.simplifyEigenvalues,
                                                    self.riemannInvariantUnitBase,
                                                    self.currentTimeStep,
                                                    self.dt)
            # initialize system equations
            self.systemEquations[vesselId].updateSystem(self.vessels[vesselId].Psol[0],
                                                        self.vessels[vesselId].Qsol[0],
                                                        self.vessels[vesselId].Asol[0])

    def initializeBoundarys(self):
        """
        initialize boundarys
        """
        if len(self.vessels) == 1:
            rootId = self.vascularNetwork.root
            bcList0 = []
            bcList1 = []
            for bc in self.vascularNetwork.boundaryConditions[rootId]:
                if bc.position == 0:     bcList0.append(bc)
                elif bc.position == -1:  bcList1.append(bc)
            self.boundarys[rootId] = [Boundary( self.vessels[rootId],
                                                bcList0,
                                                self.rigidAreas,
                                                self.dt,
                                                self.currentMemoryIndex,
                                                self.currentTimeStep,
                                                self.nTSteps,
                                                self.systemEquations[rootId]),
                                      Boundary( self.vessels[rootId],
                                                bcList1,
                                                self.rigidAreas,
                                                self.dt,
                                                self.currentMemoryIndex,
                                                self.currentTimeStep,
                                                self.nTSteps,
                                                self.systemEquations[rootId])]
            self.output['BndrNR'] = 2
        else:
            for vesselId,boundaryConditions in self.vascularNetwork.boundaryConditions.iteritems():
                self.boundarys[vesselId] = [  Boundary( self.vessels[vesselId],
                                                        boundaryConditions,
                                                        self.rigidAreas,
                                                        self.dt,
                                                        self.currentMemoryIndex,
                                                        self.currentTimeStep,
                                                        self.nTSteps,
                                                        self.systemEquations[vesselId])]

            self.output['BndrNR'] = len(self.boundarys)

    def initializeConnections(self):
        """
        initialize Connections of the network
        by traversing the network tree
        """
        treeList = self.vascularNetwork.treeTraverseList

        for leftMother,rightMother,leftDaughter,rightDaughter  in self.vascularNetwork.treeTraverseConnections:
            ## link
            if rightMother == None and rightDaughter == None:
                self.connections[leftMother] = classConnections.Link(  self.vessels[leftMother],
                                                      self.systemEquations[leftMother],
                                                      self.vessels[leftDaughter],
                                                      self.systemEquations[leftDaughter],
                                                      self.currentMemoryIndex,
                                                      self.dt,
                                                      self.rigidAreas,
                                                      self.solvingSchemeConnections)
            ## bifurcation
            elif rightMother == None:
                self.connections[leftMother] = classConnections.Bifurcation(  self.vessels[leftMother],
                                                             self.systemEquations[leftMother],
                                                             self.vessels[leftDaughter],
                                                             self.systemEquations[leftDaughter],
                                                             self.vessels[rightDaughter],
                                                             self.systemEquations[rightDaughter],
                                                             self.currentMemoryIndex,
                                                             self.dt,
                                                             self.rigidAreas,
                                                             self.solvingSchemeConnections)
            ## anastomosis
            elif rightDaughter == None:
                anastomosisId = leftMother
                if treeList.index(leftMother) > treeList.index(rightMother):
                    anastomosisId = rightMother
                self.connections[anastomosisId] = classConnections.Anastomosis(self.vessels[leftMother],
                                                             self.systemEquations[leftMother],
                                                             self.vessels[rightMother],
                                                             self.systemEquations[rightMother],
                                                             self.vessels[leftDaughter],
                                                             self.systemEquations[leftDaughter],
                                                             self.currentMemoryIndex,
                                                             self.dt,
                                                             self.rigidAreas,
                                                             self.solvingSchemeConnections)

    def initializeFields(self):
        """
        creates field numerical objects for each vessel in the network
        """
        for vesselId,vessel in self.vessels.iteritems():
            self.fields[vesselId] = classFields.Field(  vessel,
                                            self.currentMemoryIndex,
                                            self.dt,
                                            self.systemEquations[vesselId],
                                            self.rigidAreas,
                                            self.vascularNetwork.solvingSchemeField)

    def initializeCommunicators(self):
        """ TODO: Document what this should do. Currently a bit hard wired"""
        for comId, comData in self.vascularNetwork.communicators.iteritems():
            data = {'Pressure': self.vessels[comData['vesselId']].Psol,
                    'Flow'    : self.vessels[comData['vesselId']].Qsol,
                    'Area'    : self.vessels[comData['vesselId']].Asol,
                    }

            comData['data']           = data

            comData['currentMemoryIndex'] = self.currentMemoryIndex
            comData['currentTimeStep']    = self.currentTimeStep
            comData['dt']                 = self.dt

            self.communicators[comId] = eval(comData['comType'])(comData) # call the constructor

    def initializeNumericalObjectList(self):
        """
        ## fill numObjectList (self.currentTimeStepumericalObjects) traversing the treeList
        # 1. add root boundary
        # 2  add vessels
        # 3  add connection or distal boundary condition
        # 4. repeat 2,3 for the whole tree
        # 5. add communicators
        # 6. add blocking Wait if multiprocessing
        """

        # get tree traversing list
        treeList = self.vascularNetwork.treeTraverseList

        singleVesselNetwork = False
        if len(self.vessels) == 1:
            singleVesselNetwork = True

        for vesselId in treeList:
            int(vesselId)
            ## check if root add BC
            try:
                if vesselId == self.vascularNetwork.root:
                    self.numericalObjects.append(self.boundarys[vesselId][0])
            except Exception: self.warning("old except: pass #1 clause in c1dFlowSolv.initializeNumObjList", oldExceptPass= True)

            ## add field
            self.numericalObjects.append(self.fields[vesselId])

            ## try add Connection
            try: self.numericalObjects.append(self.connections[vesselId])
            except Exception: self.warning("old except: pass #2 clause in c1dFlowSolv.initializeNumObjList", oldExceptPass= True)

            ## try add distal BC
            try:
                if vesselId in self.vascularNetwork.boundaryVessels:
                    if singleVesselNetwork == False:
                        self.numericalObjects.append(self.boundarys[vesselId][0])
                    else:
                        self.numericalObjects.append(self.boundarys[vesselId][1])
            except Exception: self.warning("old except: pass #3 clause in c1dFlowSolv.initializeNumObjList", oldExceptPass= True)

        for communicator in self.communicators.itervalues():
            self.numericalObjects.append(communicator)
            try:    communicator.startRealTimeVisualisation()
            except Exception: self.warning("old except: pass #4 clause in c1dFlowSolv.initializeNumObjList", oldExceptPass= True)

        if self.venousPool:
            self.numericalObjects.append(self.venousPool)

        self.vascularNetwork.runtimeMemoryManager.registerGlobalTimeSteps(self.currentMemoryIndex,self.currentTimeStep)
        
        
        self.numericalObjects.append(self.vascularNetwork)

        # TODO: This depends on sequential iteration through the numerical objects list!!!!
        self.numericalObjects.append(self.vascularNetwork.runtimeMemoryManager)
        
        self.memoryOffset = self.vascularNetwork.runtimeMemoryManager.memoryOffset

    def initOutput(self):
        """
        initialize solution matrices
        """
        logger.info('=====================================')
        logger.info('___________Time variables ___________')
        logger.info('%-20s %2.3f' % ('totaltime (sec)',self.totalTime))
        logger.info('%-20s %2.3f' % ('dt (ms)',self.dt*1.0E3))
        logger.info('%-20s %4d' % ('nTSteps',self.nTSteps))
        logger.info('___________Div variables ____________')
        logger.info('%-20s %2.1f' % ('Q init (ml s-1)',self.vascularNetwork.initialValues[self.vascularNetwork.root]['Flow']*1.e6))
        logger.info('%-20s %2.1f' % ('P init (mmHg)',self.vascularNetwork.initialValues[self.vascularNetwork.root]['Pressure'][0]/133.32))
        try: logger.info('%-20s %2.1f' % ('R_cum (mmHg s ml-1)',self.vascularNetwork.Rcum[self.vascularNetwork.root]/133.32*1.e-6))
        except Exception: self.warning("old except: pass clause in c1dFlowSolv.initOutput", oldExceptPass= True)
        logger.info('%-20s %2.1f' % ('CFL init max',self.vascularNetwork.CFL))
        logger.info('%-20s %2.1f' % ('dz min (mm)',self.output['dz_min']*1.0E3))
        logger.info('%-20s %2.1f' % ('c min (m/s)',self.output['c_min']))
        logger.info('%-20s %2.1f' % ('c max (m/s)',self.output['c_max']))
        logger.info('%-20s %4d' % ('Grid nodens',self.output['gridNodens']))
        logger.info('___________Num variables ____________')
        logger.info('%-20s %4d' % ('NumObj',len(self.fields)+len(self.boundarys)+len(self.connections)))
        logger.info('%-20s %4d' % ('NumFields',len(self.fields)))
        logger.info('%-20s %4d' % ('NumConnections',len(self.connections)))
        logger.info('%-20s %4d' % ('NumBoundarys',len(self.boundarys)))
        logger.info('%-20s %4d' % ('NumCommunicators',len(self.communicators)))
        logger.info('%-20s %4d' % ('NumObj calls',len(self.numericalObjects)*self.nTSteps))
        logger.info('%-20s %4d' % ('used Memory (Mb)',memoryUsagePsutil()))
        logger.info('===================================== \n')


    """
    ########################################################################################
    # Solver Methods:
    #
    #    MacCormack_Field
    #
    ########################################################################################
    """

    def MacCormack_Field(self):
        """
        MacCormack solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.

        This method is solving the system by looping through the defined network
        imposing the boundary conditions based on Riemann Invariants and then solving the vessels,
        conncetions, bifucations with a predictor-corrector step method
        """
        if self.quiet == False:
            logger.info("Solving system ...")
            progressBar = cPB.ProgressBar(35,self.nTSteps, subpressPrint = self.quiet)

        reflectionCoefficientCount = 0
        maxRef = 0
    
        if self.cycleMode == False:
            # original
            for n in xrange(self.nTSteps):
                self.currentTimeStep[0] = n
                self.currentMemoryIndex[0] = n - self.memoryOffset[0]
                
                if self.quiet == False:
                    progressBar.progress(n)
                
                for numericalObject in self.numericalObjects:
                    try:
                        numericalObject()
                    except Exception:
                        # Save the Solution data for debugging
                        logger.critical("Exception caught in  {} by MacCormack_Field attempting to save solution data file...".format(numericalObject))
                        self.vascularNetwork.runtimeMemoryManager.flushSolutionMemory()
                        self.vascularNetwork.saveSolutionData()
                        logger.critical("Success in saving solution data file. Reraising Exception")
                        raise # TODO: why does self.exception() not force the program to quit?
                        # self.exception()
                
                if self.quiet == False:
                    progressBar.progress(n)
                
        ## to be concentrated with original cycle mode !!
        else:
            # steady state variables
            P_lastCycle  = {}
            Q_lastCycle  = {}
            A_lastCycle  = {}

            for vesselId,vessel in self.vessels.iteritems():

                initialValues = self.vascularNetwork.initialValues
                p0,p1 = initialValues[vesselId]['Pressure']
                Qm    = initialValues[vesselId]['Flow']

                P_lastCycle[vesselId]  = np.ones((self.nTSteps,vessel.N))
                Q_lastCycle[vesselId]  = np.ones((self.nTSteps,vessel.N))
                A_lastCycle[vesselId]  = np.ones((self.nTSteps,vessel.N))


            for cycle in xrange(self.numberCycles-1):
                logger.info(' solving cycle {}'.format(cycle+1))
                # 1. solve cycle
                for n in xrange(self.nTSteps-1):

                    self.currentTimeStep[0] = n
                    for numericalObject in self.numericalObjects:
                        numericalObject()
                # 2. check for steady state
                if self.quiet == False:
                    for vesselId in self.vessels.keys():
                        #Perror =  np.sum(np.sqrt((np.divide((P_lastCycle[vesselId]-self.P[vesselId]),P_lastCycle[vesselId]))**2.0))/self.P[vesselId].size
                        #Qerror =  np.sum(np.sqrt((np.divide((Q_lastCycle[vesselId]-self.Q[vesselId]),Q_lastCycle[vesselId]))**2.0))/self.Q[vesselId].size
                        #Aerror =  np.sum(np.sqrt((np.divide((A_lastCycle[vesselId]-self.A[vesselId]),A_lastCycle[vesselId]))**2.0))/self.A[vesselId].size

                        Perror =  np.max(np.sqrt((np.divide((P_lastCycle[vesselId]-self.P[vesselId]),P_lastCycle[vesselId]))**2.0))
                        #Qerror =  np.max(np.sqrt((np.divide((Q_lastCycle[vesselId]-self.Q[vesselId]),Q_lastCycle[vesselId]))**2.0))
                        Aerror =  np.max(np.sqrt((np.divide((A_lastCycle[vesselId]-self.A[vesselId]),A_lastCycle[vesselId]))**2.0))

                        P_lastCycle[vesselId] = self.P[vesselId].copy()
                        #Q_lastCycle[vesselId] = self.Q[vesselId].copy()
                        A_lastCycle[vesselId] = self.A[vesselId].copy()
                        logger.debug('{} {}'.format(Perror, Aerror))

                # 3. rehash solution arrays if not last cycle
                if cycle is not self.numberCycles-1:
                    for vesselId in self.vessels.keys():
                        self.P[vesselId][0]    = self.P[vesselId][-1]
                        self.Q[vesselId][0]    = self.Q[vesselId][-1]
                        self.A[vesselId][0]    = self.A[vesselId][-1]


        if self.quiet == False: logger.info("\nSystem solved!")

        BVcheck = False
        if BVcheck == True:

            logger.info('\n=====================================')
            logger.info('___Blood volume - consistency check___')

            logger.info('Vessels  init (ml)  sol (ml)  diff')
            vesselsInit = 0
            vesselsSol  = 0
            for vesselId,vessel in self.vessels.iteritems():
                A1 = self.vessels[vesselId].Asol[0][0:-1]
                A2 = self.vessels[vesselId].Asol[0][1:]
                volumeInit = np.sum(vessel.dz*(A1+A2+np.sqrt(A1*A2))/3.0)*1.e6

                A1 = self.vessels[vesselId].Asol[-1][0:-1]
                A2 = self.vessels[vesselId].Asol[-1][1:]
                volumeSol  = np.sum(vessel.dz*(A1+A2+np.sqrt(A1*A2))/3.0)*1.e6
                diff = volumeInit-volumeSol
                vesselsInit += volumeInit
                vesselsSol  += volumeSol
                logger.info('{:<4}    {:7.2f}    {:7.2f}    {:4.2f}'.format(str(vesselId), volumeInit, volumeSol, diff))

            logger.info(' ----------------------------------- ')
            logger.info('Boundarys        in (ml)   out (ml)  ')
            totalIn = 0
            totalOut = 0
            for boundaryList in self.boundarys.itervalues():
                for boundary in boundaryList:
                    logger.info('{:<12}     {:6.2f}    {:4.2f}'.format(boundary.name,
                                                                 abs(boundary.BloodVolumen[0])*1.e6,
                                                                 abs(boundary.BloodVolumen[1])*1.e6 ))
                    totalIn  += abs(boundary.BloodVolumen[0])*1.e6
                    totalOut += abs(boundary.BloodVolumen[1])*1.e6
            logger.info(' ----------------------------------- ')
            logger.info('total')
            logger.info('  vessels initial  (ml)  {:4.2f}'.format(vesselsInit))
            logger.info('  vessels solution (ml)  {:4.2f}'.format(vesselsSol))
            logger.info('  boundarys in     (ml)  {:4.2f}'.format(totalIn))
            logger.info('  boundarys out    (ml)  {:4.2f}'.format(totalOut))
            logger.info('                    ------------')
            logger.info('  volume diff      (ml)  {:.2f}'.format(vesselsInit - vesselsSol + totalIn - totalOut))

        ## stop realtime visualisation
        for communicator in self.communicators.itervalues():
            try: communicator.stopRealtimeViz()
            except Exception: self.warning("old except: pass #1 clause in c1dFlowSolv.MacCormack_Field", oldExceptPass= True)

        ### garbage collection
        gc.collect()

        del self.numericalObjects
        del self.fields
        del self.connections
        del self.boundarys
        del self.communicators

