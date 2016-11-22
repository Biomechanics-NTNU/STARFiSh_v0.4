import sys
import os
import numpy as np
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')

import UtilityLib.classStarfishBaseObject as cSBO
import NetworkLib.classBoundaryConditions as ccBC


class Boundary(cSBO.StarfishBaseObject):
    def __init__(self, vessel, boundaryConditions, rigidArea, dt, currentMemoryIndex, currentTimeStep, nTsteps, systemEquation):
        """
        Constructor of Boundary
        s
        Initializes one Boundary of a vessel
        input: - boundaryConditions:
                  a list of all boundaryConditions of the vessel at one position i.g 0 or -1

        Note: if network consits of a single vessel make sure the function is called twice
              once for the start and once for the end boundary at 0 and -1!

        variables (most important):
        self.position        :    position of the boundary
        self.ID              :    Vessel ID to which the boundary belongs
        self.duFunction      :    function which gives du-vector back
        self.bcType1         :    list of all type1 functions found in the given boundaryConditions
        self.omegaInput      :    function which gives the _omega of omega_ back depending on self.position
        self.omegaFunctio    :    function with boundaryCondtion of type2 calculating the omega-vector
        """
        self.position = None

        self.name = ' '.join(['Boundary',str(vessel.Id)])
        self.type = ''

        #Function which gives du-vector back
        self.duFunction = None
        self.duVector = np.zeros((nTsteps,2))
        self._du_ = np.empty(2)
        # list of all type1 functions found in the given boundaryConditions, which determine
        # the du-Function
        self.bcType1 = []
        #Function which gives the _omega of omega_ back depending on self.position
        self.omegaInput = None
        #Function with boundaryCondtion of type2 calculating the omega-vector
        self.omegaFunction = None
        #list of all type2 functions found in the given boundaryConditions
        self.bcType2 = []

        #System and Vessel Variables
        self.z  = vessel.z
        self.dt = dt
        self.currentTimeStep = currentTimeStep
        self.currentMemoryIndex = currentMemoryIndex
        self.systemEquation = systemEquation

        #SolutionVariables
        self.P = vessel.Psol
        self.Q = vessel.Qsol
        self.A = vessel.Asol

        self.BloodFlowSep  = np.zeros((nTsteps,2)) ## [ dbloodVolume in, dbloodVolume out] blood volume from one timestep to the next
        self.BloodVolumen  = np.zeros(2)

        posTemp = []
        for bC in boundaryConditions:
            if bC.type == 1:
                self.bcType1.append(bC)
            elif bC.type == 2:
                self.bcType2.append(bC)
            posTemp.append(bC.position)

        #find out position of the Boundary + check if all conditions are defined at the same side
        if sum(posTemp) == 0: self.position = 0
        elif sum(np.array(posTemp)+1) == 0: self.position = -1
        else: raise ValueError("One position of one boundaryCondition is not correct")

        #initialize solution variable
        self.BloodFlowSep[-1][self.position] = vessel.Qsol[0][self.position]

        # 1. Set duFunction
        if len(self.bcType1) == 0:
            # in for callMacCormack normal
            self.duFunction = self.duFunctionZero
            self.duRuntimeEvaluation = False
        elif len(self.bcType1) == 1:
            self.duFunction = self.duFunctionSingle
            self.bcType1ConditionQuantity = self.bcType1[0].conditionQuantity
            self.duRuntimeEvaluation      = self.bcType1[0].runtimeEvaluation
            precribeTotalValues           = self.bcType1[0].prescribeTotalValues
        elif len(self.bcType1) > 1:
            self.duFunction = self.duFunctionMulti
            self.bcType1ConditionQuantity = self.bcType1[0].conditionQuantity
            self.duRuntimeEvaluation      = self.bcType1[0].runtimeEvaluation
            precribeTotalValues           = self.bcType1[0].prescribeTotalValues
            print "WARNING: Multiple Type1 Boundary Conditions are prescribed as influx condition!"
        # prepare duVector from boundaryConditions:
        self.duEvaluateVector()
        # set bcType1 to type instead of list
        try: self.bcType1 = self.bcType1[0]
        except:pass

        # 3. Set Condition there should only one! if None apply the one condition depending on the side
        if len(self.bcType2) == 1:
            self.omegaFunction = self.bcType2[0]
            self.type = self.bcType2[0].name

        elif len(self.bcType2) == 0:
            # set calculation type for Standard Boundary
            self.type = self.bcType1.name
            try:
                if precribeTotalValues == False:
                    self.omegaFunction = ccBC.PrescribedInflux()
                elif precribeTotalValues == True and self.bcType1ConditionQuantity == 'Flow':
                    self.omegaFunction = ccBC.PrescribedTotalFlow()
                elif precribeTotalValues == True and self.bcType1ConditionQuantity == 'Pressure':
                    self.omegaFunction = ccBC.PrescribedTotalPressure()
            except:
                self.omegaFunction = ccBC.PrescribedInflux()
            self.omegaFunction.setPosition(self.position)
        else:
            self.warning("classBoundary: Too many type2-boundary Conditions defined!", noException= True)

        # 4. Define the output of A, dependend if rigidArea
        self.rigidArea = rigidArea
        self.A_nID = vessel.A_nID

        ## 5. Define the call function depending on the solving Scheme
        #if solvingScheme == "MacCormack_Field":
        #    self.__call__ = self.callMacCormackField

    ### Function which calculated du
    def duFunctionZero(self,currentTimeStep,dt):
        """
        Determine the du-vector = [0,0] if no type1 boundaryConditions are given
        """
        return np.zeros(2)

    def duFunctionSingle(self,currentTimeStep,dt):
        """
        Determine the du-vector with the values of the given
        type1 boundaryCondition
        """
        return self.bcType1.calculateDu(currentTimeStep,dt)

    def duFunctionMulti(self,currentTimeStep,dt):
        """
        Determine the summized du-vector with the values of all given
        type1 boundaryConditions
        """
        du = np.zeros(2)
        for bc in self.bcType1:
            du = du + bc.calculateDu(currentTimeStep,dt)
        return du

    def duEvaluateVector(self):
        """
        Pre-Calculate the BoundaryConditions duVector of Type1
        """
        self.duVector = self.duVector*0.
        for bc in self.bcType1:
            self.duVector = self.duVector+bc.calculateDuVector(len(self.duVector),self.dt)

    def duFunctionRuntime(self,currentTimeStep,dt,Z1,Z2,position):

        if self.duRuntimeEvaluation == True:
            duPrescribed = self.duFunction(currentTimeStep,dt)
        else:
            duPrescribed = self.duVector[currentTimeStep]
        # check if influx Values should be prescribed
        if self.bcType1 != []:
            # check if Flow is prescribed
            if self.bcType1ConditionQuantity == 'Flow':
                #position 0: pf = Zf*Qf
                if position == 0:
                    duPrescribed[0]  = Z1*duPrescribed[1]
                #position -1: pb = -Zb*Qb
                elif position == -1:
                    duPrescribed[0] = -Z2*duPrescribed[1]
            # check if Pressure is prescribed
            elif self.bcType1ConditionQuantity == 'Pressure':
                #position 0: Qf = pf/Zf
                if position == 0:
                    duPrescribed[1]  = duPrescribed[0]/Z1
                #position -1: Qb = - pf/Zb
                elif position == -1:
                    duPrescribed[1] = -duPrescribed[0]/Z2

        return duPrescribed

    def __call__(self): #callMacCormackField(self):
        """
        new Boundary method calculates the values at the boundary
        and applies it to the new values
        """

        # create local variables for this timestep
        dt = self.dt
        currentMemoryIndex = self.currentMemoryIndex[0]
        currentTimeStep = self.currentTimeStep[0]

        P = self.P[currentMemoryIndex]
        Q = self.Q[currentMemoryIndex]
        A = self.A[currentMemoryIndex]

        z = self.z
        position = self.position

        # calculate need values
        L,R,LMBD,Z1,Z2,_omega_field = self.systemEquation.updateLARL(P,Q,A,position)

        #calculate the du_vector using given boundaryConditions of type 1
        duPrescribed = self.duFunctionRuntime(currentTimeStep, dt, Z1, Z2, position)

        #calculate the dp and dq using given boundaryConditions of type 2
        dPQ_calc,dBloodVolumen = self.omegaFunction(_omega_field, duPrescribed, R, L, currentMemoryIndex, currentTimeStep, dt, P[position], Q[position], A[position],Z1, Z2)

        # calculate new values for p and q
        P_calc = dPQ_calc[0]+ P[position]
        Q_calc = dPQ_calc[1]+ Q[position]

        # check new p value
        if P_calc < 0:
            raise ValueError("{} calculated negative pressure P_calc = {} at time {} (n {},dt {})".format(self.name, P_calc, currentTimeStep*dt,currentTimeStep,dt))
            #exit()

        # calculate new value for the area
        A_calc = A[position] # assign old value
        if self.rigidArea == False:
            A_calc =  self.A_nID([P_calc],position)

        # apply values to solution array
        self.P[currentMemoryIndex+1][position] = P_calc
        self.Q[currentMemoryIndex+1][position] = Q_calc
        self.A[currentMemoryIndex+1][position] = A_calc

        try: self.BloodFlowSep[currentTimeStep] = self.BloodFlowSep[currentTimeStep-1]+dBloodVolumen
        except Exception: print "passed bloodflow integration"
        try:    self.BloodVolumen = self.BloodVolumen + 0.5*(self.BloodFlowSep[currentTimeStep-1]+self.BloodFlowSep[currentTimeStep])*dt
        except Exception: self.BloodVolumen = self.BloodVolumen + self.BloodFlowSep[currentTimeStep]*dt

