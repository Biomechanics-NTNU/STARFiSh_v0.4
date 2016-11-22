import sys
import os
import numpy as np
cur = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(cur+'/../')
import UtilityLib.classStarfishBaseObject as cSBO
import UtilityLib.classConfigurableObjectBase as cCOB


class StaticVenousPressure(cSBO.StarfishBaseObject):
    """
    """
    # defined external data
    externVariables      = {'P0':cCOB.ConfigurableObjectBase.ExtValue(float, unit='Pa'),
                             'pressureGain':cCOB.ConfigurableObjectBase.ExtValue(float, unit=None )}
    externXmlAttributes  = []
    externXmlElements    = ['P0', 'pressureGain']
            
    def __init__(self):
        self.veinId  = 0

        self.pressureGain = 3. # pressure gain between CVP and LAP - Bell paper
        self.P0 = 2.0*133.32
        self.Pmin = 0.0
        self.P = [self.P0]
        self.P_LA =[self.pressureGain*self.P[0]]
        
    def __call__(self):
        # TODO: Make this more efficient
        self.P[0] = self.P0
        self.P_LA[0] = self.pressureGain*self.P[0]
    
    def update(self,dataDict):
        """
        updates the data
        Dict = {'variableName': value}
        """
        for key,value in dataDict.iteritems():
            if hasattr(self,key):
                setattr(self,key,value)
            else:
                logger.debug("StaticVenousPool.update(): wrong key: %s, could not set up venousPool" %key)

class venousPool(cSBO.StarfishBaseObject):
    """
    Very simple model of the venous side, considering the veins as one big compliant reservoir,
    and assuming a pure pressure gain between CVP and LAP
    The Baroreflex regulates the unstretched volume of the venous side, through which the CVP and ultimately
    the LAP are changed

    self.V is the blood volume in the veins
    self.Vusv: unstretched Volume of Veins with zero external pressure
    self.P0: constant
    self.k: constant
    self.pressureGain: pressure gain from CVP (right atrial) to LAP (left atrial) --> a pure gain is used, value according to Bell
    self.P: Central Venouse Pressure i.e. right atrial pressure
    self.P_LA: pressure in left atrium (LAP)
    self.Qin: inflow
    self.Qout: outflow
    """
    # defined external data
    externVariables      = {'P0' : cCOB.ConfigurableObjectBase.ExtValue(float,  unit = 'Pa'),
                            'V0' :  cCOB.ConfigurableObjectBase.ExtValue(float,  unit = 'm3'),
                            'Vusv0' :  cCOB.ConfigurableObjectBase.ExtValue(float,  unit = 'm3'),
                            'k' :  cCOB.ConfigurableObjectBase.ExtValue(float,  unit = None),
                             'pressureGain': cCOB.ConfigurableObjectBase.ExtValue(float,  unit = None )}
    externXmlAttributes  = []
    externXmlElements    = externVariables.keys()
    
    
    solutionMemoryFields    = ["Vusv", "V", "P", "Qin", "Qout", "P_LA"]
    solutionMemoryFieldsToSave = ["Vusv", "V", "P", "Qin", "Qout", "P_LA"]

    def __init__(self):

        self.dt = 0 #will be updated with update method
        self.currentTimeStep = 0 # current time step
        self.currentMemoryIndex = 0
        self.nTSteps = 0

        self.boundarys = {} # make it a dictionary/needs to be initialized in FlowSolver
        self.veinId  = 0


        """
        ### FINDING INITIAL VALUES FOR VOLUME
        V = 3892. # 5600.0e-6 *0.61# * 0.61 # estimated blood volume on venous side under normal conditions
        Vusv0 = 2378. #  ? 3213e-6 # unstretched volume at reference state
        P0 = 2.0 # pressure constant for calculation of P venous
        k = 0.1124

        Vusv = Vusv0
        def pFct(V):
            P = P0*(np.exp(k*(V-Vusv)**1.5/(V)))
            return P

        from scipy import optimize as opt
        fct = lambda V: pFct(V) - 3.0
        opt.brentq(fct, 2400,4000)
        # 2850.912397321067
        fct = lambda V: pFct(V) - 4.0
        opt.brentq(fct, 2400,4000)
        # 3091.6779241832583
        fct = lambda V: pFct(V) - 5.0
        opt.brentq(fct, 2400,4000)
        # 3270.4477485970647
        fct = lambda V: pFct(V) - 6.0
        opt.brentq(fct, 2400,4000)
        # 3414.6023352947177
        fct = lambda V: pFct(V) - 7.0
        opt.brentq(fct, 2400,4000)
        # 3536.118289840244
        fct = lambda V: pFct(V) - 8.0
        opt.brentq(fct, 2400,4000)
        # 3641.5166289669837
        """

        # TODO figure out how to have BRX not need these as maximal values
        self.V0 = 3770.4477485970647e-6 # 3892e-6 # 5600.0e-6 *0.61# * 0.61 # estimated blood volume on venous side under normal conditions
        self.Vusv0 =  3400.e-6 # 2378e-6 #  ? 3213e-6 # unstretched volume at reference state


        self.P0 = 2.0 * 133.322368 # pressure constant for calculation of P venous
        self.Pmin = 0.0
        self.k = 0.1124 #0.1124e-9 # constant

        self.pressureGain = 1.0/0.228 # pressure gain between CVP and LAP - Bell paper

        self.Qin = np.zeros(0) # in and outflow to the venous pool
        self.Qout = np.zeros(0)

        ### vectors for export
        self.V = np.zeros(0)
        self.P = np.array([self.pressureFromVolume(self.V0, self.Vusv0)])
        self.P_LA = np.zeros(0)
        self.Vusv = np.zeros(0)

    def initializeWithFlowSolver(self, flowSolver):
        self.currentTimeStep         = flowSolver.currentTimeStep
        self.currentMemoryIndex      = flowSolver.currentMemoryIndex
        self.boundarys               = flowSolver.boundarys

    def initializeForSimulation(self, vascularNetwork):
        self.dt = vascularNetwork.dt
        self.nTSteps = vascularNetwork.nTSteps
        self.dsetGroup = vascularNetwork.solutionDataFile.create_group('Venous')
        self.allocate(vascularNetwork.runtimeMemoryManager)
        
        self.boundaryCondtions = vascularNetwork.boundaryConditions
        self.Vusv[:] = self.Vusv0
        self.V[0] = self.V0
        self.P[0] = self.pressureFromVolume(self.V[0],self.Vusv[0])
        self.P_LA[0] = self.pressureGain*self.P[0]
        self.Qin[0] = 0.0
        self.Qout[0] = 0.0

    def estimateInflow(self):
        """
        calculate the inflow to the venous side, from terminal boundaries
        """
        nmem = self.currentMemoryIndex[0]
        Qin = 0

        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                if self.boundarys[key][x].position == -1: #distal boundaries
                    bcCondition = self.boundarys[key][x].bcType2[0]
                    # TODO: Add other types
                    if bcCondition.name == 'Windkessel-3Elements':
                        deltaP = self.boundarys[key][x].P[nmem,-1] - bcCondition.venousPressure[nmem]
                        Qin = Qin + deltaP/bcCondition.Rtotal
        self.Qin[nmem] = Qin

    def estimateOutflow(self):
        """
        calculate outflow
        """
        nmem = self.currentMemoryIndex[0]
        Qout = 0
        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                if self.boundarys[key][x].position == 0: #proximal boundaries

                    if self.boundarys[key][x].bcType2[0].name == 'VaryingElastanceHeart':
                        Qout = Qout + self.boundarys[key][x].bcType2[0].mitralQ[nmem]
                    else:
                        Qout = Qout + self.boundarys[key][x].Q[nmem,0]
        self.Qout[nmem] = Qout

    def pressureFromVolume(self,V, Vusv):
        return self.P0*(np.exp(self.k*((V-Vusv)*1e6)**1.5/(V*1e6)))

    def updateVenousPool(self):
        """
        update the state of the venous pool (volume and pressure)
        the model is from Ursino_1999
        """
        
        nmem = self.currentMemoryIndex[0]
        
        self.V[nmem+1] = self.V[nmem] + self.dt*(self.Qin[nmem] - self.Qout[nmem])
        try:
            # TODO: Concurrency issue as Vusv may be modified by BRX on same time step?
            self.P[nmem+1] = self.pressureFromVolume(self.V[nmem+1], self.Vusv[nmem])
        except ValueError:
            if self.V[nmem+1] - self.Vusv[nmem+1] < 0:
                self.exception(
                   "Venous volume {} is lower the unstressed venous volume {} at time {} and time step {}.".format(
                        self.V[nmem+1], self.Vusv[nmem+1], self.dt*self.currentTimeStep[0], self.currentTimeStep[0]))
            else:
                raise

        self.P_LA[nmem+1] = self.pressureGain*self.P[nmem+1]


    def updateBoundaryConditions(self):
        """
        update all the boundary conditions
        """
        n = self.currentMemoryIndex[0]

        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):

                if self.boundarys[key][x].position == 0:
                    if self.boundarys[key][x].type == 'VaryingElastanceHeart':
                        self.boundarys[key][x].atriumPressure[n+1] = self.P_LA[n+1]

                    else: pass

                if self.boundarys[key][x].position == -1:
                        #
                        # Replaces precalculated Pv(t) = Pcv_init + pgh(t) + -Pext by
                        # Pv[n+1] = Pv[n+1] - Pv[0] + CVP[n+1],
                        # where Pv[n+1] - Pv[0] = Pcv_init + pgh[n+1] + -Pext - (Pcv_init + pgh[0] + -Pext)
                        # which leaves Pv[n+1] - Pv[0] = pgh[n+1] - pgh[0], and thus
                        # Pv[n+1] = CVP[n+1] + rho gh[n+1] - rho gh[0], which is correct so long as rho gh[0] = 0 and Pext = 0
                        #

                    if self.boundarys[key][x].type == 'Windkessel-2Elements':
                        self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P[n+1]

                    if self.boundarys[key][x].type == 'Windkessel-3Elements':
                        if n < (np.size(self.boundarys[key][x].bcType2[0].venousPressure)-1):
                            self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P[n+1]

                    if self.boundarys[key][x].type == 'Resistance':
                        self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P[n+1]

    def __call__(self):
        """
        call function for Venous pool, updates the venous pool volume and pressure
        and updates the boundary conditions with new pressure values
        """
        self.estimateInflow()
        self.estimateOutflow()
        self.updateVenousPool()
        self.updateBoundaryConditions()

    def update(self,dataDict):
        """
        updates the data
        Dict = {'variableName': value}
        """
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except Exception:
                self.warning("venousPool.update(): wrong key: %s, could not set up venousPool" %key)

class VenousPoolXMLWrapper(cSBO.StarfishBaseObject):
    # TODO: This is a hack until the top level structure is resolved fully
    externVariables      = {'venousPoolContent':cCOB.ConfigurableObjectBase.ExtObject({'StaticVenousPressure':StaticVenousPressure,
                                         'venousPool':venousPool})}
    externXmlAttributes  = []
    externXmlElements    = ['venousPoolContent']
    
    def __init__(self):
        self.venousPoolContent = None
