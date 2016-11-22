import numpy as np
import csv
from math import exp
import scipy.optimize as so
import ODESolver as OD

import os, sys
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))

sys.path.append(cur+'/../')

import UtilityLib.classStarfishBaseObject as cSBO
import UtilityLib.moduleFilePathHandler as mFPH

class BoundaryCondition(cSBO.StarfishBaseObject):
    """
    Base-class for all boundary conditions
    """
    position = None
    name = None

    # #
    # if False prescribe influx of Flow and Pressure
    # if True  prescribe total Values of Pessure and Flow
    prescribeTotalValues = False
    # #
    # conditionQuantity <string> should be either 'Flow' or 'Pressure'
    conditionQuantity = None


    def update(self, bcDict):
        """
        updates the updateBoundaryDict data using a dictionary in form of
        bcDict = {'variableName': value}
        """
        for key, value in bcDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key, value)
            except Exception:
                self.warning("boundaryConditions.update(): wrong key: %s, could not set up boundaryCondition" % key)

    def setPosition(self, position):
        """Set the position of the boundaryCondition """
        self.position = position

    def getVariableValue(self, variableName):
        """
        Returns value of variable with name : variableName
        States Error if not such variable
        """
        try:
            return self.__getattribute__(variableName)
        except Exception:
            self.warning("BoundaryCondition.getVariable() : BoundaryCondition has no variable {}".format(variableName))


class BoundaryConditionType2(BoundaryCondition):
    def setPosition(self, position):
        """
        Set the position of the boundaryCondition
        and determines the return function
        based on the position of the bC (0,-1)
        """
        self.position = position
        if self.position == 0:
            self.returnFunction = self.funcPos0
        elif self.position == -1:
            self.returnFunction = self.funcPos1
        else:
            self.returnFunction = None


class generalPQ_BC(BoundaryConditionType2):
    """
    Defines a wrapper for boundary conditions which may be defined as
        f(u,du/dt,t) = 0,
    where u = [P,Q]

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """
    NEWTON_TOL = 1e-6 # This is the accuracy of the estimate of \omega, thus it's units are  order of [Pa] or [m^3/s]
    def funcPos0(self, _domegaField, R, dt, P, Q,A, nmem, n):
        """return function for position 0 at the start
        of the vessel
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        omegaOld_ = 0.0
        # TODO: Generalize residual wrappers?
        # Compared Brentq, fsolve and newton, and newton seems to be the most efficient.
        domega_ = so.newton(self.residualW1Newton, x0 = omegaOld_, tol=self.NEWTON_TOL, args = (_domegaField,P, Q,A, dt, nmem,n, r11, r12, r21, r22))

        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut

    def funcPos1(self, domegaField_, R, dt, P, Q,A, nmem, n):
        """
        return function for position -1 at the end
        of the vessel
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]
        _omegaOld = 0.0

        # Compared Brentq, fsolve and newton, and newton seems to be the most efficient.
        _domega = so.newton(self.residualW2Newton, x0 = _omegaOld,  tol=self.NEWTON_TOL, args = (domegaField_, P, Q,A, dt, nmem, n, r11, r12, r21, r22))

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domega

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()
        return np.dot(R, self.omegaNew), self.dQInOut

    def residualW1Newton(self,domega_, _domegaField, P, Q,A, dt, nmem, n, r11, r12, r21, r22):
        """
        Maps the input state of omega_1 and omega_2 to dP and dQ for the PQ residual equaition.

        Args:
            args = [_domegaField,P, Q, dt,n, r11, r12, r21, r22]

        """
        dP = r11*domega_+ r12*_domegaField
        dQ = r21*domega_+ r22*_domegaField

        return self.residualPQ(dP,dQ,P,Q,A,dt,nmem,n)

    def residualW2Newton(self,_domega, domegaField_,P, Q,A, dt, nmem, n, r11, r12, r21, r22):
        """
        Maps the input state of omega_1 and omega_2 to dP and dQ for the PQ residual equaition.

        Args:
            args = [domegaField_,P, Q, dt,n, r11, r12, r21, r22]

        """
        dP = r11*domegaField_+ r12*_domega
        dQ = r21*domegaField_+ r22*_domega
        return self.residualPQ(dP,dQ,P,Q,A,dt,nmem,n)

    def residualW1(self,domega_, args):
        """
        Maps the input state of omega_1 and omega_2 to dP and dQ for the PQ residual equaition.

        Args:
            args = [_domegaField,P, Q, dt,nmem,n, r11, r12, r21, r22]

        """
        _domegaField,P, Q, A, dt, nmem, n, r11, r12, r21, r22 = args
        if _domegaField>0.0:
            dP = r11*domega_+ r12*_domegaField
            dQ = r21*domega_+ r22*_domegaField
        else:
            dP = r11*domega_
            dQ = r21*domega_

        return self.residualPQ(dP, dQ, P, Q, A, dt, nmem, n)

    def residualW2(self,_domega, args):
        """
        Maps the input state of omega_1 and omega_2 to dP and dQ for the PQ residual equaition.

        Args:
            args = [domegaField_,P, Q, dt,nmem,n, r11, r12, r21, r22]

        """
        domegaField_,P, Q, A, dt, nmem, n, r11, r12, r21, r22 = args
        dP = r11*domegaField_+ r12*_domega
        dQ = r21*domegaField_+ r22*_domega

        return self.residualPQ(dP, dQ, P, Q, A, dt, nmem, n)

    def residualPQ(self,dP,dQ,P,Q,A,dt,nmem,n):
        raise NotImplementedError('Residual equation residualPQ for omega must be implemented.')

########################################################################################
# Type 1 boundary Conditions
########################################################################################

class BoundaryConditionType1(BoundaryCondition):
    """
    Boundary profile - type 1

        creates frame for a periodic signal

        signal type can be quantified with function1 / function2
        which are a function of (t)

        call function input (n,dt)
        gives amplitude value for pressure/flow back
        for the function delta g
    """
    def __init__(self):

        BoundaryCondition.__init__(self)

        self.type = 1
        self.prescribe = 'influx'

        # # variables from xml
        self.amp     = 0
        self.ampConst = 0
        self.Npulse = 1.
        self.freq     = 1.
        self.Tpulse = 0.
        self.Tspace = 0.
        self.runtimeEvaluation = False

        # # evaluated values
        self.lastU                = 0.0
        self.TmeanFlow            = 0.0
        self.initPhaseTimeSpan      = 0.0
        self.MeanFlow             = 0.0
        self.Tperiod              = 0.0
        self.pulseTime            = []
        self.duMatrix             = np.ones(2)
        self.duVector             = np.empty(2)

        self.initialisationPhaseExist = False
        self.nTstepsInitPhase          = 0

        # values to manipulate period during runtime
        self.updateTime = -10.0
        self.freqNew = self.freq
        self.TspaceNew = self.Tspace

    def initialize(self, bcDict):
        """
        updates - the updateBoundaryDict data using a dictionary in from of
                  bcDict = {'variableName': value}
                - time period of one pulse
                - pulseTimes
        """

        self.update(bcDict)

        self.Tperiod = self.Tspace + 1.0 / self.freq

        self.pulseTime = []
        for puls in np.arange(self.Npulse):
            self.pulseTime.append([self.Tpulse + (self.Tperiod) * puls,
                                   self.Tpulse + 1.0 / self.freq + (self.Tperiod) * puls,
                                   self.Tpulse + self.Tperiod + (self.Tperiod) * puls])

        self.pulseTimeInit = [[0,
                               1.0 / self.freq,
                               self.Tperiod]]

        self.lastU = self.calculateOneStep(0, 0)

        if self.prescribe == 'total':
            self.prescribeTotalValues = True
        elif self.prescribe == 'influx':
            self.prescribeTotalValues = False

        if 'Flow' in self.name:
            self.duMatrix = np.array([0, 1])
            self.conditionQuantity = 'Flow'
        elif 'Pressure' in self.name:
            self.duMatrix = np.array([1, 0])
            self.conditionQuantity = 'Pressure'

        # # load file for Flow-From File
        if self.name == 'Flow-FromFile':
            if self.loadedFile == False:
                self.loadFile()

    def updatePeriodRuntime(self, TperiodNew, updateTime):
        """
        This function updates the freq, Tspace and pulsTime for new given Tperiod
        """

        freq = self.freq
        Tspace = self.Tspace
        Tperiod = self.Tperiod
        # # calculate new values
        self.freqNew = Tperiod * freq / TperiodNew
        self.TspaceNew = TperiodNew * Tspace / Tperiod

        manipulateOtherPulses = False
        findFirstPulse = True

        for pulse in self.pulseTime:

            if  updateTime >= pulse[0] and updateTime <= pulse[2] and findFirstPulse == True:
                pStart = pulse[2]
                manipulateOtherPulses = True

                self.updateTime = pulse[2]

            if manipulateOtherPulses == True:
                if findFirstPulse == False:
                    pulse[0] = pStart
                    pulse[1] = pStart + 1.0 / self.freqNew  # - self.TmeanFlow
                    pulse[2] = pStart + TperiodNew  # - self.TmeanFlow
                    pStart = pulse[2]
                else: findFirstPulse = False

    def calculateDuVector(self, nTsteps, dt):
        """
        Function calculates the duVector for a given number of time steps Tsteps and dt
        return: self.duVector
        """
        # create du Vector
        self.duVector = np.zeros((nTsteps, 2))

        initPhase = False
        # check if initphase needs to be added at start
        if self.initialisationPhaseExist:
            initPhase = True

        lastStep = self.calculateOneStep(0, dt, initPhase=initPhase)

        # add rest of simulation time
        for nt in xrange(nTsteps):
            if nt == self.nTstepsInitPhase: initPhase = False
            nextStep = self.calculateOneStep(nt + 1, dt, initPhase=initPhase)
            self.duVector[nt] = nextStep - lastStep
            lastStep = nextStep

        return self.duVector

    def calculateDu(self, n, dt):
        """
        Function calculates the du for a given time step and dt
        return: du
        """
        nextStep = self.calculateOneStep(n + 1, dt)
        du = nextStep - self.lastU
        self.lastU = nextStep

        if n * dt > self.updateTime - dt and n * dt < self.updateTime + dt :

            self.freq = self.freqNew
            self.Tspace = self.TspaceNew
            self.Tperiod = self.TspaceNew + 1.0 / self.freqNew

            self.updateTime = -10

        return du

    def calculateOneStep(self, n, dt, initPhase=False):
        """
        calculates the amplitude value for one time step n and dt

        return: array([ampP,ampQ])
        """
        t = n * dt
        ampT = self.ampConst
        timeInPulseTime = False

        pulseTime = self.pulseTime
        if initPhase == True: pulseTime = self.pulseTimeInit
        # get time slots
        for pTA in pulseTime:
            if pTA[0] <= t and t <= pTA[2]:
                pulsNum = pulseTime.index(pTA)
                pTA1 = pTA[0]
                pTA2 = pTA[1]
                pTA3 = pTA[2]
                timeInPulseTime = True
                break
        if timeInPulseTime == False:
            return ampT * self.duMatrix

        if initPhase == True: t = t + self.TmeanFlow

        try:
            if pTA1 <= t and t <= pTA2:
                t0 = pTA1
                # # function1 according to the signal type
                ampT = ampT + self.function1(t, t0, pulsNum)

            elif pTA2 <= t and t < pTA3:
                t0 = pTA1
                # # function2 according to the signal type in most cases return self.ampConst
                ampT = ampT + self.function2(t, t0, pulsNum)

            t = n * dt
            if self.Tpulse <= t and t < self.Tpulse + 1.5 * dt and self.initialStep == False:
                self.initialStep = True

            return ampT * self.duMatrix

        except Exception:          #TODO This should get some comment to explain what is happening.
            return ampT * self.duMatrix

    def findMeanFlowAndMeanTime(self, givenMeanFlow = None, quiet = False):
        """
        This function calculates the mean flow of the signal self.MeanFlow
        and the first occurence evaluatedTime of the mean flow self.TmeanFlow
        """
        #find meanFlow
        self.initialize({})

        period = self.Tperiod
        totalTime = period+self.Tpulse

        MeanFlow = 123221.00
        dt = 1e-5*2 #0.1ms
        tolerance = 5.e-10

        notConverged = True
        while notConverged == True:
            dt = dt/2.
            nTsteps = int(np.round(totalTime/dt, 0))
            nTstepsStart = int(np.round(self.Tpulse/dt, 0))
            evaluatedFlow = np.zeros(nTsteps-nTstepsStart-1)
            evaluatedTime = np.zeros(nTsteps-nTstepsStart-1)

            for n in xrange(nTstepsStart,nTsteps-1):
                flow = self.calculateOneStep(n,dt)
                evaluatedFlow[n-nTstepsStart] = flow[1]
                evaluatedTime[n-nTstepsStart]= (n-nTstepsStart)*dt

            evaluatedMeanFlow = np.mean(evaluatedFlow)
            #print dt, evaluatedMeanFlow ,    MeanFlow-evaluatedMeanFlow, nTsteps-nTstepsStart
            if evaluatedMeanFlow < min(evaluatedFlow) or evaluatedMeanFlow > max(evaluatedFlow):
                print "Meanflow is wrong"
                exit()
            if abs(MeanFlow-evaluatedMeanFlow) < tolerance:
                notConverged = False
            
            MeanFlow = evaluatedMeanFlow

        if givenMeanFlow == None:
            self.MeanFlow = evaluatedMeanFlow
        else:
            self.MeanFlow = givenMeanFlow
            difference = abs(givenMeanFlow - evaluatedMeanFlow)
            tmpstring = ("given meanFlow given {} differs \n from meanFlow of boundaryCondition {}".format(givenMeanFlow*1.e6,evaluatedMeanFlow*1.e6)
                         + "\n(evaluated with integral over one period)."
                         + "\nThe difference is {} ml s-1".format(difference*1.e6))
            self.warning(tmpstring,  noException= True, quiet= quiet)
#            if quiet == False:
#                print """\n  WARNING:  given meanFlow given {} differs from meanFlow of boundaryCondition {}
#                (evaluated with integral over one period).
#                The difference is {} ml s-1 \n""".format(givenMeanFlow*1.e6,evaluatedMeanFlow*1.e6,difference*1.e6)
# TODO wrote this before convert to spaces


        
        self.TmeanFlow = 0
        self.initPhaseTimeSpan = 0
        # find mean flow time
        if givenMeanFlow != 0 or evaluatedFlow[0] != 0:
            for n in xrange(len(evaluatedFlow)):
                tn = evaluatedTime[n]
                qn = evaluatedFlow[n]
                if n < len(evaluatedFlow)-1:
                    ## someting for the end
                    qn1 = evaluatedFlow[n+1]
                    tn1 = evaluatedTime[n+1]
                    if (qn-self.MeanFlow)*(qn1-self.MeanFlow) < 0:
                        qm = self.MeanFlow
                        linearInt = tn+(qm-qn)/(qn1-qn)*(tn1-tn)
                        self.TmeanFlow = linearInt
                        break
                else:
                    self.warning("no corresponding time found to set mean flow {}, simulation is initialized with no time shift".format(self.MeanFlow*1.e6))

        if self.TmeanFlow != 0:
            self.initPhaseTimeSpan = self.Tperiod - self.TmeanFlow
        self.initialize({})
        

        if quiet == False:
            print '====================================='
            print '___BC _Type1: mean flow evaluation___'
            print 'meanFlow evaluated (ml/s)  {:.6}'.format(str(evaluatedMeanFlow*1.e6).ljust(5))
            print 'meanFlowTime (s)           {:.6}'.format(str(self.TmeanFlow).ljust(5))
            print 'initPhaseTimeSpan          {:.6}'.format(str(self.initPhaseTimeSpan).ljust(5))
            print 'total volume/period (ml)   {:.6}'.format(str(np.sum(evaluatedFlow*1.e6)*dt).ljust(5))

        return  self.MeanFlow, self.initPhaseTimeSpan

    def function1(self, t, t0, pulsNum):
        """
        amplitude function caracterized by signal type
        """
        pass
    def function2(self, t, t0, pulsNum):
        """
        amplitude function caracterized by signal type
        """
        return 0

class RampMean(BoundaryConditionType1):
    """
    Boundary profile - type 1

        ramps to a mean amplitude self.amp starting from self.ampConst
    """

    def function1(self, t, t0, pulsNum):
        return self.amp * pow(np.sin(np.pi * (t - self.Tpulse) / (1. / self.freq * 2.0)), 2.0)

    def function2(self, t, t0, pulsNum):
        return self.amp

class Sinus(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates in a periodic sinus signal
    """
    def function1(self, t, t0, pulsNum):
        return self.amp * np.sin(2 * np.pi * self.freq * (t - t0))


class Sinus2(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates in a periodic sinus-squared signal
    """
    def function1(self, t, t0, pulsNum):
        return self.amp * pow(np.sin(np.pi * self.freq * (t - t0)), 2.0)


class ExpFunc(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates an exponentialfunction according to section 3.1 Comparison1Dscheme
    """
    def function1(self, t, t0, pulsNum):

        return 0.000001*np.exp(-10000*pow((t-0.05),2))


class CCAInflow(BoundaryConditionType1):
    """
    Creates the inflow used in test Of Common Carotid artery section 3.2 in Comparison1Dscheme
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.T = 1.1

    def function1(self, t, t0, pulsNum):

        return 1e-6*(6.5+3.294*np.sin(2*np.pi*t/self.T-0.023974)+1.9262*np.sin(4*np.pi*t/self.T-1.1801)-1.4219*np.sin(6*np.pi*t/self.T+0.92701)-0.66627*np.sin(8*np.pi*t/self.T-0.24118)-0.33933*np.sin(10*np.pi*t/self.T-0.27471)-0.37914*np.sin(12*np.pi*t/self.T-1.0557)+0.22396*np.sin(14*np.pi*t/self.T+1.22)+0.1507*np.sin(16*np.pi*t/self.T+1.0984)+0.18735*np.sin(18*np.pi*t/self.T+0.067483)+0.038625*np.sin(20*np.pi*t/self.T+0.22262)+0.012643*np.sin(22*np.pi*t/self.T-0.10093)-0.0042453*np.sin(24*np.pi*t/self.T-1.1044)-0.012781*np.sin(26*np.pi*t/self.T-1.3739)+0.014805*np.sin(28*np.pi*t/self.T+1.2797)+0.012249*np.sin(30*np.pi*t/self.T+0.80827)+0.0076502*np.sin(32*np.pi*t/self.T+0.40757)+0.0030692*np.sin(34*np.pi*t/self.T+0.195)-0.0012271*np.sin(36*np.pi*t/self.T-1.1371)-0.0042581*np.sin(38*np.pi*t/self.T-0.92102)-0.0069785*np.sin(40*np.pi*t/self.T-1.2364)+0.0085652*np.sin(42*np.pi*t/self.T+1.4539)+0.0081881*np.sin(44*np.pi*t/self.T+0.89599)+0.0056549*np.sin(46*np.pi*t/self.T+0.17623)+0.0026358*np.sin(48*np.pi*t/self.T-1.3003)-0.0050868*np.sin(50*np.pi*t/self.T-0.011056)-0.0085829*np.sin(52*np.pi*t/self.T-0.86463))


class AortaInflow(BoundaryConditionType1):
    """
    Creates the inflow used in test Of Single Aorta used in section 3.3 in Comparison1Dscheme
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.T = 0.9550

    def function1(self, t, t0, pulsNum):
#         print "self.T is: ", self.T
#         print "t is: ", t
#         print "return is: ", 4.0*np.arctan(1.0)*0.0126*0.0126*(0.20617+0.37759*np.sin(2*np.pi*t/self.T+0.59605)+0.2804*np.sin(4*np.pi*t/self.T-0.35859)+0.15337*np.sin(6*np.pi*t/self.T-1.2509)-0.049889*np.sin(8*np.pi*t/self.T+1.3921)+0.038107*np.sin(10*np.pi*t/self.T-1.1068)-0.041699*np.sin(12*np.pi*t/self.T+1.3985)-0.020754*np.sin(14*np.pi*t/self.T+0.72921)+0.013367*np.sin(16*np.pi*t/self.T-1.5394)-0.021983*np.sin(18*np.pi*t/self.T+0.95617)-0.013072*np.sin(20*np.pi*t/self.T-0.022417)+0.0037028*np.sin(22*np.pi*t/self.T-1.4146)-0.013973*np.sin(24*np.pi*t/self.T+0.77416)-0.012423*np.sin(26*np.pi*t/self.T-0.46511)+0.0040098*np.sin(28*np.pi*t/self.T+0.95145)-0.0059704*np.sin(30*np.pi*t/self.T+0.86369)-0.0073439*np.sin(32*np.pi*t/self.T-0.64769)+0.0037006*np.sin(34*np.pi*t/self.T+0.74663)-0.0032069*np.sin(36*np.pi*t/self.T+0.85926)-0.0048171*np.sin(38*np.pi*t/self.T-1.0306)+0.0040403*np.sin(40*np.pi*t/self.T+0.28009)-0.0032409*np.sin(42*np.pi*t/self.T+1.202)-0.0032517*np.sin(44*np.pi*t/self.T-0.93316)+0.0029112*np.sin(46*np.pi*t/self.T+0.21405)-0.0022708*np.sin(48*np.pi*t/self.T+1.1869)-0.0021566*np.sin(50*np.pi*t/self.T-1.1574)+0.0025511*np.sin(52*np.pi*t/self.T-0.12915)-0.0024448*np.sin(54*np.pi*t/self.T+1.1185)-0.0019032*np.sin(56*np.pi*t/self.T-0.99244)+0.0019476*np.sin(58*np.pi*t/self.T-0.059885)-0.0019477*np.sin(60*np.pi*t/self.T+1.1655)-0.0014545*np.sin(62*np.pi*t/self.T-0.85829)+0.0013979*np.sin(64*np.pi*t/self.T+0.042912)-0.0014305*np.sin(66*np.pi*t/self.T+1.2439)-0.0010775*np.sin(68*np.pi*t/self.T-0.79464)+0.0010368*np.sin(70*np.pi*t/self.T-0.0043058)-0.0012162*np.sin(72*np.pi*t/self.T+1.211)-0.00095707*np.sin(74*np.pi*t/self.T-0.66203)+0.00077733*np.sin(76*np.pi*t/self.T+0.25642)-0.00092407*np.sin(78*np.pi*t/self.T+1.3954)-0.00079585*np.sin(80*np.pi*t/self.T-0.49973))
        return 4.0*np.arctan(1.0)*0.0126*0.0126*(0.20617+0.37759*np.sin(2*np.pi*t/self.T+0.59605)+0.2804*np.sin(4*np.pi*t/self.T-0.35859)+0.15337*np.sin(6*np.pi*t/self.T-1.2509)-0.049889*np.sin(8*np.pi*t/self.T+1.3921)+0.038107*np.sin(10*np.pi*t/self.T-1.1068)-0.041699*np.sin(12*np.pi*t/self.T+1.3985)-0.020754*np.sin(14*np.pi*t/self.T+0.72921)+0.013367*np.sin(16*np.pi*t/self.T-1.5394)-0.021983*np.sin(18*np.pi*t/self.T+0.95617)-0.013072*np.sin(20*np.pi*t/self.T-0.022417)+0.0037028*np.sin(22*np.pi*t/self.T-1.4146)-0.013973*np.sin(24*np.pi*t/self.T+0.77416)-0.012423*np.sin(26*np.pi*t/self.T-0.46511)+0.0040098*np.sin(28*np.pi*t/self.T+0.95145)-0.0059704*np.sin(30*np.pi*t/self.T+0.86369)-0.0073439*np.sin(32*np.pi*t/self.T-0.64769)+0.0037006*np.sin(34*np.pi*t/self.T+0.74663)-0.0032069*np.sin(36*np.pi*t/self.T+0.85926)-0.0048171*np.sin(38*np.pi*t/self.T-1.0306)+0.0040403*np.sin(40*np.pi*t/self.T+0.28009)-0.0032409*np.sin(42*np.pi*t/self.T+1.202)-0.0032517*np.sin(44*np.pi*t/self.T-0.93316)+0.0029112*np.sin(46*np.pi*t/self.T+0.21405)-0.0022708*np.sin(48*np.pi*t/self.T+1.1869)-0.0021566*np.sin(50*np.pi*t/self.T-1.1574)+0.0025511*np.sin(52*np.pi*t/self.T-0.12915)-0.0024448*np.sin(54*np.pi*t/self.T+1.1185)-0.0019032*np.sin(56*np.pi*t/self.T-0.99244)+0.0019476*np.sin(58*np.pi*t/self.T-0.059885)-0.0019477*np.sin(60*np.pi*t/self.T+1.1655)-0.0014545*np.sin(62*np.pi*t/self.T-0.85829)+0.0013979*np.sin(64*np.pi*t/self.T+0.042912)-0.0014305*np.sin(66*np.pi*t/self.T+1.2439)-0.0010775*np.sin(68*np.pi*t/self.T-0.79464)+0.0010368*np.sin(70*np.pi*t/self.T-0.0043058)-0.0012162*np.sin(72*np.pi*t/self.T+1.211)-0.00095707*np.sin(74*np.pi*t/self.T-0.66203)+0.00077733*np.sin(76*np.pi*t/self.T+0.25642)-0.00092407*np.sin(78*np.pi*t/self.T+1.3954)-0.00079585*np.sin(80*np.pi*t/self.T-0.49973))

class AoBifInflow(BoundaryConditionType1):
    """
    Creates the inflow used in test Of Single Aorta used in section 3.4 in Comparison1Dscheme
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.T = 1.1

    def function1(self, t, t0, pulsNum):

        return (7.9853e-06+2.6617e-05*np.sin(2*np.pi*t/self.T+0.29498)+2.3616e-05*np.sin(4*np.pi*t/self.T-1.1403)-1.9016e-05*np.sin(6*np.pi*t/self.T+0.40435)-8.5899e-06*np.sin(8*np.pi*t/self.T-1.1892)-2.436e-06*np.sin(10*np.pi*t/self.T-1.4918)+1.4905e-06*np.sin(12*np.pi*t/self.T+1.0536)+1.3581e-06*np.sin(14*np.pi*t/self.T-0.47666)-6.3031e-07*np.sin(16*np.pi*t/self.T+0.93768)-4.5335e-07*np.sin(18*np.pi*t/self.T-0.79472)-4.5184e-07*np.sin(20*np.pi*t/self.T-1.4095)-5.6583e-07*np.sin(22*np.pi*t/self.T-1.3629)+4.9522e-07*np.sin(24*np.pi*t/self.T+0.52495)+1.3049e-07*np.sin(26*np.pi*t/self.T-0.97261)-4.1072e-08*np.sin(28*np.pi*t/self.T-0.15685)-2.4182e-07*np.sin(30*np.pi*t/self.T-1.4052)-6.6217e-08*np.sin(32*np.pi*t/self.T-1.3785)-1.5511e-07*np.sin(34*np.pi*t/self.T-1.2927)+2.2149e-07*np.sin(36*np.pi*t/self.T+0.68178)+6.7621e-08*np.sin(38*np.pi*t/self.T-0.98825)+1.0973e-07*np.sin(40*np.pi*t/self.T+1.4327)-2.5559e-08*np.sin(42*np.pi*t/self.T-1.2372)-3.5079e-08*np.sin(44*np.pi*t/self.T+0.2328))


class ExperimentalInflow(BoundaryConditionType1):
    """
    Creates the inflow used in Benchmark experimental test used in section 3.5 in Comparison1Dscheme
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.T = 0.827

    def function1(self, t, t0, pulsNum):

        return (3.1199+7.7982*np.sin(2*np.pi*t/self.T+0.5769)+4.1228*np.sin(4*np.pi*t/self.T-0.8738)-1.0611*np.sin(6*np.pi*t/self.T+0.7240)+0.7605*np.sin(8*np.pi*t/self.T-0.6387)-0.9148*np.sin(10*np.pi*t/self.T+1.1598)+0.4924*np.sin(12*np.pi*t/self.T-1.0905)-0.5580*np.sin(14*np.pi*t/self.T+1.042)+0.3280*np.sin(16*np.pi*t/self.T-0.5570)-0.3941*np.sin(18*np.pi*t/self.T+1.2685)+0.2833*np.sin(20*np.pi*t/self.T+0.6702)+0.2272*np.sin(22*np.pi*t/self.T-1.4983)+0.2249*np.sin(24*np.pi*t/self.T+0.9924)+0.2589*np.sin(26*np.pi*t/self.T-1.5616)-0.1460*np.sin(28*np.pi*t/self.T-1.3106)+0.2141*np.sin(30*np.pi*t/self.T-1.1306)-0.1253*np.sin(32*np.pi*t/self.T+0.1552)+0.1321*np.sin(34*np.pi*t/self.T-1.5595)-0.1399*np.sin(36*np.pi*t/self.T+0.4223)-0.0324*np.sin(38*np.pi*t/self.T+0.7811)-0.1211*np.sin(40*np.pi*t/self.T+1.0729))/1000./60.

class Adan55InflowFromfile(BoundaryConditionType1):
    """
    Creates the inflow used in Benchmark experimental test used in section 3.5 in Comparison1Dscheme
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)
        mypath= '/home/Fredrik/fredrik/Master/Comparios1D/AoBif/Matlab/'
        Flowarray=np.array([])
        Timearray=np.array([])
        # # additional variables for this function
        with open('/home/Fredrik/Code/NetworkFiles/Adan55/Adan55Inflow.csv') as csvfile:
            reader = csv.reader(csvfile)
            for rows in reader:
                Timearray=np.append(Timearray,[float(rows[0])])
                Flowarray=np.append(Flowarray,[float(rows[1])])
        self.Timearray=Timearray
        self.Flowarray=Flowarray*1e-6
        self.T=1

    def function1(self, t, t0, pulsNum):
        t=t-self.T*pulsNum
        Q=np.interp(t,self.Timearray,self.Flowarray)
        #print "t, Q: ",t,", ",Q
        #print "pulsNum: ", pulsNum
        return Q

class PhysiologicalFunction(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates a similar heart-outflow signal as found in
        Stergiopulos et al. 1992

        set together from 4 continous functions sin2,sin2,cos,sin2
        to lead to a continous function
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.lowPoint = 0.1739
        self.fracSin2 = 0.36
        self.fracCos = 0.43
        self.fracRes = 1.0 - (self.fracSin2 + self.fracCos)

    def function1(self, t, t0, pulsNum):
        if t < self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq) + self.fracSin2 * 1 / self.freq:
            ampT = self.amp * pow(np.sin(np.pi * self.freq / (2 * self.fracSin2) * (t - (self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq)))), 2.0)
        elif t < self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq) + (self.fracSin2 + self.fracCos) * 1 / self.freq:
            ampT = self.amp * (np.cos((np.pi * self.freq / (2.0 * self.fracCos) * (t - self.fracSin2 / self.freq - (self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq))))))
        elif t < self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq) + (self.fracSin2 + self.fracCos + self.fracRes * 3.0 / 8.0) * 1 / self.freq:
            ampT = self.amp * (pow(np.sin((t - (self.fracSin2 + self.fracCos) / self.freq - (self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq)) - (1.0 / (self.freq) * self.fracRes) * 3. / 8.) * np.pi * self.freq / (1.5 * self.fracRes)), 2.0) * 2 * self.lowPoint - self.lowPoint)
        else:
            ampT = self.amp * (pow(np.sin((t - (self.fracSin2 + self.fracCos) / self.freq - (self.Tpulse + pulsNum * (self.Tspace + 1.0 / self.freq)) - (1.0 / (self.freq) * self.fracRes) * 3. / 8.0) * np.pi * self.freq / (self.fracRes * 1.25)), 2.0) * self.lowPoint - self.lowPoint)
        return ampT

class expVelocity(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates a single gaussian peak signal as found
        in D.Xiu and S.P.Sherwin 2007
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.gaussC = 5000
        self.area = 1

        # # set to duMatrix
        self.duMatrix[1] = 1


    def function1(self, t, t0, pulsNum):
        return self.amp * exp(-self.gaussC * (t - t0 + 0.5 / self.freq) ** 2) * self.area

class Fourier(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates fourier signal similar to the flow of the heart
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        self.scale = 1
        self.harm = [[0.2465, 0.1975, 0.2624, 0.2132, 0.0424, 0.0722, 0.0411, 0.0093, 0.0141, 0.0044],
                     [0.0, 2.2861, 4.5723, 6.8584, 9.1445, 11.4307, 13.7168, 16.0029, 18.2891, 20.5752],
                     [0.0, 2.5010, -2.9986, -0.2689, 1.4904, -2.9854, -0.0317, 1.5292, -2.5394, 0.5771]]

        # # set variables
        self.freq = self.harm[1][1]

    def function1(self, t, t0, pulsNum):
        amp = 0
        for i in range(len(self.harm[0])):
            amp += self.harm[0][i] * np.sin(2.0 * np.pi * self.harm[1][i] * (t - t0) + self.harm[2][i] + np.pi / 2.0)
        # amp = amp*self.scale
        return amp


class PhysiologicalData(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates signal based on measured physiological values
        (source unknown: values in NetworkLib/physiologicalData.py)
    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables for this function
        from physiologicalData import aorticFlowPressure
        self.data = aorticFlowPressure()

        tDataMax = max(self.data.t)
        tDataMin = min(self.data.t)
        self.timeData = np.linspace(tDataMin, tDataMax, 2000)


    def function1(self, t, t0, pulsNum):
        self.timeSim = np.linspace(0.0, 1. / self.freq, 2000)
        tReset = t - t0
        tInter = np.interp(tReset, self.timeSim, self.timeData)
        Q = np.interp(tInter, self.data.t, self.data.Q)
        # P = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.P)
        return Q * self.duMatrix

class FlowFromFile(BoundaryConditionType1):
    """
    Boundary profile - type 1

        creates signal based on data values stored in *.csv file
        saved in network directory

        csv delimiter ;
        first line of the colums must be t and Q respectively:
        t;Q
        ..; ..

        expected units: t [s]; Q[m3 / s]

    """
    def __init__(self):

        BoundaryConditionType1.__init__(self)

        # # additional variables fill in with data in the readXML file
        self.networkDirectory = ''
        self.filePathName = ''
        self.networkDirectory = ''
        self.dataTime = []
        self.dataFlow = []
        self.loadedFile = False

    def loadFile(self, inflowFilePath):

        try:
            # set the path relative to THIS file not the executing file!
            
            #pathAndFilename = mFPH.getFilePath('inflowFile', networkName, dataNumber, mode, exception)
            reader = csv.DictReader(open(self.filePathName, 'rb'), delimiter=';')
        except Exception:
            self.exception("boundaryConditions.FlowFromFile could not open file <<{}>> with boundary values, system exit".format(self.filePathName.split('/')[-1]))
        try:
            dataTime = []
            dataFlow = []
            for row in reader:
                dataTime.append(float(row['t']))
                dataFlow.append(float(row['Q']))
            self.dataTime = np.asarray(dataTime)
            self.dataFlow = np.asarray(dataFlow)
            self.loadedFile = True
        except Exception: self.warning("Old except: pass clause in ccBC.FlowFromFile.loadFile")


    def function1(self, t, t0, pulsNum):
        # set up simulation time (needed here if freq changes after init known to slow down)
        self.timeData = np.linspace(min(self.dataTime), max(self.dataTime), 2000)
        self.timeSim = np.linspace(0.0, 1. / self.freq, 2000)
        # reset time to fit in array [0 , 1/freq]
        tReset = t - t0
        # interpolate time to dataTime
        tInter = np.interp(tReset, self.timeSim, self.timeData)
        # interpolate the Q value from data and dataTime
        Q = np.interp(tInter, self.dataTime, self.dataFlow)
        # P = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.P)
        return Q * self.duMatrix


########################################################################################
# Type 2 boundary Conditions
########################################################################################
class PrescribedInflux(BoundaryConditionType2):
    """
    Boundary profile - type 2

    calculates omega if influx of pressure or flow is prescribed with a type1 condition

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2
        self.name = "standardType2"

        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, L, R)

    def funcPos0(self, _domegaField, duPrescribed, L, R):
        """
        return function for position 0 at the start
        of the vessel
        """
        domegaPrescribed_ = np.dot(L[0], duPrescribed)

        self.omegaNew[0] = domegaPrescribed_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut


    def funcPos1(self, domegaField_, duPrescribed, L, R):
        """
        return function for position -1 at the end
        of the vessel
        """
        _domegaPrescribed = np.dot(L[1], duPrescribed)

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domegaPrescribed

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew) , self.dQInOut


class PrescribedTotalFlow(BoundaryConditionType2):
    """
    Boundary profile - type 2

    calculates omega if total flow is prescribed with a type1 condition

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2
        self.name = "standardType2"

        self.R = np.empty((2, 2))
        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.duNew = np.empty((2))
        self.dQInOut = np.zeros((2))

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, L)


    def funcPos0(self, _domegaField, duPrescribed, L):
        """
        return function for position 0 at the start
        of the vessel
        """
        """
            # 1. of L set values l11 = 0 l12 = 1.
            # 2. and make inverse
            # det = l11*l22 - l12*l21 = -l21
            self.R[0][0] =  L[1][1] / -L[1][0]
            self.R[0][1] =  1./L[1][0]     #-1 / -L[1][0]
            self.R[1][0] =  1.            #-L[1][0] / -L[1][0]
            self.R[1][1] =  0
            # 3. set omega vector
            self.omegaNew[0] = duPrescribed[1]
            self.omegaNew[1] = _domegaField
            # 4. du = self.R * self.omegaNew
            du = np.dot(self.R,self.omegaNew)
        How it is done:
            inserted R in dot product
        """
        self.duNew[0] = -L[1][1] / L[1][0] * duPrescribed[1] + 1. / L[1][0] * _domegaField
        self.duNew[1] = duPrescribed[1]
        self.dQInOut[0] = duPrescribed[1]
        return self.duNew, self.dQInOut

    def funcPos1(self, domegaField_, duPrescribed, L):
        """
        return function for position -1 at the end
        of the vessel
        """
        """
        How to do it:
            # 1. of L set values l21 = 0 l22 = 1.
            # 2. and make inverse
            # det = l11*l22 - l12*l21 = l11
            self.R[0][0] =  1. / L[0][0] # 1/l11
            self.R[0][1] =  -L[0][1]/L[0][0] # -l12/l11
            self.R[1][0] =  0
            self.R[1][1] =  1.
            # 3. set omega vector
            self.omegaNew[0] = domegaField_
            self.omegaNew[1] = duPrescribed[1]
            # 4. du = self.R * self.omegaNew
            du = np.dot(self.R,self.omegaNew)
        How it is done:
            inserted R in dot product
        """
        self.duNew[0] = domegaField_ / L[0][0] - L[0][1] / L[0][0] * duPrescribed[1]
        self.duNew[1] = duPrescribed[1]

        self.dQInOut[1] = duPrescribed[1]
        return self.duNew, self.dQInOut

class PrescribedTotalPressure(BoundaryConditionType2):
    """
    Boundary profile - type 2

    calculates omega if total pressure is prescribed with a type1 condition

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2
        self.name = "standardType2"

        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.R = np.empty((2, 2))
        self.duNew = np.empty((2))
        self.dQInOut = np.zeros(2)

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, L)


    def funcPos0(self, _domegaField, duPrescribed, L):
        """
        return function for position 0 at the start
        of the vessel
        """
        """
        How to do it:
            # 1. of L set values l11 = 0 l12 = 1.
            # 2. and make inverse
            # det = l11*l22 - l12*l21 = l22
            self.R[0][0] =  1.
            self.R[0][1] =  0
            self.R[1][0] =  -L[1][0]/L[1][1]
            self.R[1][1] =  1./L[1][1]
            # 3. set omega vector
            self.omegaNew[0] = duPrescribed[0]
            self.omegaNew[1] = _domegaField
            # 4. du = self.R * self.omegaNew
            du = np.dot(self.R,self.omegaNew)
        How it is done:
            inserted R in dot product
        """
        self.duNew[0] = duPrescribed[0]
        self.duNew[1] = -L[1][0] / L[1][1] * duPrescribed[0] + _domegaField / L[1][1]
        self.dQInOut[0] = self.duNew[1]
        return self.duNew, self.dQInOut

    def funcPos1(self, domegaField_, duPrescribed, L):
        """
        return function for position -1 at the end
        of the vessel
        """
        """
        How to do it:
            # 1. of L set values l11 = 0 l12 = 1.
            # 2. and make inverse
            # det = l11*l22 - l12*l21 = -l12
            self.R[0][0] =  0.
            self.R[0][1] =  1.
            self.R[1][0] =  1./L[0][1]
            self.R[1][1] =  L[0][0] / -L[0][1]

            # 3. set omega vector
            self.omegaNew[0] = domegaField_
            self.omegaNew[1] = duPrescribed[0]
            # 4. du = self.R * self.omegaNew
            du = np.dot(self.R,self.omegaNew)
        How it is done:
            inserted R in dot product
        """
        self.duNew[0] = duPrescribed[0]
        self.duNew[1] = domegaField_ / L[0][1] - L[0][0] / L[0][1] * duPrescribed[0]
        self.dQInOut[1] = self.duNew[1]
        return self.duNew, self.dQInOut

class ReflectionCoefficient(BoundaryConditionType2):
    """
    Boundary profile - type 2

    Terminal reflection
    only in combination with prescribed influx condition (Type1) or alone

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """
    def __init__(self):
        self.type = 2
        self.Rt = 0

        self.returnFunction = None

        self.omegaNew = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem, n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, L, R)

    def funcPos0(self, _domegaField, duPrescribed, L, R):
        """
        return function for position 0 at the start
        of the vessel
        """
        domegaPrescribed_ = np.dot(L[0], duPrescribed)
        domegaReflected_ = _domegaField * self.Rt + domegaPrescribed_

        self.omegaNew[0] = domegaReflected_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut


    def funcPos1(self, domegaField_, duPrescribed, L, R):
        """
        return function for position -1 at the end
        of the vessel
        """
        _domegaPrescribed = np.dot(L[1], duPrescribed)
        _domegaReflected = domegaField_ * self.Rt + _domegaPrescribed

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domegaReflected

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew), self.dQInOut


class ReflectionCoefficientTimeVarying(BoundaryConditionType2):
    """
    Boundary profile - type 2

    Terminal reflection
    only in combination with prescribed influx condition (Type1) or alone

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """
    def __init__(self):
        self.type = 2

        # # reflection Values
        self.RtOpen = 0.2
        self.RtClosed = 0.8

        # # Timing
        self.Topen1 = 0.0  # verschiebung nach rechts von 0.5/freq
        self.Topen2 = 0.0  # verschiebung nach links von 0.5/freq

        self.Tclosed1 = -0.0  # verschiebung nach rechts von Tspace anfang
        self.Tclosed2 = 0.0  # verschiebung nach links von Tspace ende

        ### need this from boundaryConditionType1 ####
        # # is set with update method in vascularNetwork in calculateInitialValues
        self.Tpulse = 0.0
        self.Tspace = 0.0
        self.freq = 1.
        self.Npulse     = 1.
        self.TmeanFlow = 0.0
        # ##

        self.pulseTimeR = []

        self.returnFunction = None

        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, L, R, n, dt)

    def funcPos0(self, _domegaField, duPrescribed, L, R, n, dt):
        """
        return function for position 0 at the start
        of the vessel
        """
        Rt = self.calcRt(n, dt)
        domegaPrescribed_ = np.dot(L[0], duPrescribed)
        domegaReflected_ = _domegaField * Rt + domegaPrescribed_

        self.omegaNew[0] = domegaReflected_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut


    def funcPos1(self, domegaField_, duPrescribed, L, R, n, dt):
        """
        return function for position -1 at the end
        of the vessel
        """
        Rt = self.calcRt(n, dt)
        _domegaPrescribed = np.dot(L[1], duPrescribed)
        _domegaReflected = domegaField_ * Rt + _domegaPrescribed

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domegaReflected

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew), self.dQInOut



    def update(self, bcDict):
        """
        updates the updateBoundaryDict data using a dictionary in from of
        bcDict = {'variableName': value}
        """
        # TODO: this should refer to BoundaryCondition's update function instead of doing the exact same thing.
        for key, value in bcDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key, value)
            except: pass
                # print 'ValueError: wrong key: %s, could not set up boundaryCondition' %key

        self.type = 2
        self.name = 'ReflectionCoefficientTimeVarying'

        self.pulseTimeR = []

        # let the valve be open 20 + of the pulse
        # self.Topen1 = - (0.5/self.freq)*0.1
        # self.Topen2 =   (0.5/self.freq)*0.1

        self.Tperiod = self.Tspace + 1.0 / self.freq

        for puls in np.arange(self.Npulse):
            self.pulseTimeR.append([self.Tpulse + (self.Tperiod) * puls + self.Tclosed2 - self.TmeanFlow,
                                    self.Tpulse + 0.5 / self.freq + (self.Tperiod) * puls - self.TmeanFlow + self.Topen1,
                                    self.Tpulse + 0.5 / self.freq + (self.Tperiod) * puls - self.TmeanFlow + self.Topen2,
                                    self.Tpulse + 1.0 / self.freq + (self.Tperiod) * puls - self.TmeanFlow + self.Tclosed1,
                                    self.Tpulse + self.Tperiod + (self.Tperiod) * puls - self.TmeanFlow + self.Tclosed2,
                                    self.Tpulse + self.Tperiod + (self.Tperiod) * puls + self.Tclosed2 ])


    def calcRt(self, n, dt):
        t = n * dt
        RC = self.RtClosed

        # get time slots
        for pTA in self.pulseTimeR:
            if pTA[0] < t and t < pTA[5]:
                # pulsNum = self.pulseTimeR.index(pTA)
                pTA1 = pTA[0]
                pTA2 = pTA[1]
                pTA3 = pTA[2]
                pTA4 = pTA[3]
                pTA5 = pTA[4]
                pTA6 = pTA[5]
                break
        try:
            if (pTA1 < t and t < pTA2):
                RC = (self.RtOpen - self.RtClosed) / (pTA2 - pTA1) * (t - pTA1) + self.RtClosed
            elif pTA2 <= t and t <= pTA3:
                RC = self.RtOpen
            elif pTA3 < t and t < pTA4:
                RC = -(self.RtOpen - self.RtClosed) / (pTA4 - pTA3) * (t - pTA3) + self.RtOpen
            elif (pTA5 < t and t < pTA6):
                # RC = self.RtOpen
                RC = (self.RtOpen - self.RtClosed) / (pTA2 - pTA1) * (t - pTA5) + self.RtClosed
            else: pass
            return RC
        except Exception: #TODO: should be explained.
            return RC


class Resistance(BoundaryConditionType2):
    """
    Boundary profile - type 2

    signle Resistance element

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2
        self.Rc = 1

        self.venousPressure = 20.*133  # not needed

        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, R, Z1)

    def funcPos0(self, _domegaField, R, Z1):
        """
        return function for position 0 at the start
        of the vessel
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        Rc = self.Rc
        if Rc == 'VesselImpedance':
            Rc = Z1

        domega_ = -(r22 * Rc + r12) / (r21 * Rc + r11) * _domegaField

        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut


    def funcPos1(self, domegaField_, R, Z1):
        """
        return function for position -1 at the end
        of the vessel
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        Rc = self.Rc
        if Rc == 'VesselImpedance':
            Rc = Z1

        _domega = (r11 - r21 * Rc) / (r22 * Rc - r12) * domegaField_

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domega

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew), self.dQInOut


class Windkessel2DAE(generalPQ_BC):
    """
    Boundary profile - type 2

    2 Element Windkessel solved using as a DAE using Backwards Euler

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    solutionMemoryFields    = ["pressure", "volume", "Qin", "Qout"]
    solutionMemoryFieldsToSave = ["pressure", "volume", "Qin", "Qout"]

    def __init__(self):
        self.type = 2
        self.Rc = 1
        self.C = 0

        self.venousPressure = 7.*133.32
        self.pressure = np.zeros(0)
        self.volume = np.zeros(0)
        self.Qout = np.zeros(0)
        self.Qin = np.zeros(0)
        self.vesselId = None
        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def initializeSolutionVectors(self, runtimeMemoryManager, solutionDataFile):
        """
        Initializes the solutionMemoryFields
        """
        self.dsetGroup = solutionDataFile.create_group('WK2-'+str(self.vesselId))
        self.allocate(runtimeMemoryManager)

        self.pressure[0] = 100.*133.32
        self.volume[0] = 0.0

        if isinstance(self.venousPressure, float):
            Pv = self.venousPressure
        else:
            Pv = self.venousPressure[0]

        self.Qout[0] = (self.pressure[0]-Pv)/self.Rc

        self.Qin[0] = 0.0

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):

        self.pressure[nmem] = P
        self.Qin[nmem] = Q

        # TODO: Polymorphic scalar/array variables?
        if isinstance(self.venousPressure, float):
            Pv = self.venousPressure
        else:
            Pv = self.venousPressure[nmem]

        self.Qout[nmem] = (P - Pv)/self.Rc
        self.volume[nmem] = self.volume[nmem-1] + (Q - self.Qout[nmem])*dt

        return self.returnFunction(_domegaField_, R, dt, P, Q, nmem, n)

    def residualPQ(self, dP, dQ, P, Q, dt,nmem, n):
        return self.C*dP/dt + (P+dP-self.venousPressure[n])/self.Rc - (Q+dQ)

class Windkessel3DAE(generalPQ_BC):
    """
    Boundary profile - type 2

    3 Element Windkessel solved using as a DAE using Backwards Euler

    """

    solutionMemoryFields    = ["pressure", "volume", "Qin", "Qout"]
    solutionMemoryFieldsToSave = ["pressure", "volume", "Qin", "Qout"]

    def __init__(self):
        self.type = 2
        self.Z0 = None
        self.Z = 'VesselImpedance'
        self.Zrt = self.Z0
        self.firstRun = False
        self.Rc = None
        self.Rcrt = 1.0
        self.Rtotal = None
        self.firstRun = False
        self.C = 0

        self.venousPressure = 7.*133.32
        self.Pv = self.venousPressure
        self.pressure = np.zeros(0)
        self.volume = np.zeros(0)
        self.Qout = np.zeros(0)
        self.Qin = np.zeros(0)
        self.vesselId = None
        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def initializeSolutionVectors(self, runtimeMemoryManager, solutionDataFile):
        """
        Initializes the solutionMemoryFields
        """
        self.dsetGroup = solutionDataFile.create_group('WK3-'+str(self.vesselId))
        self.allocate(runtimeMemoryManager)

        self.pressure[0] = 100.*133.32
        self.volume[0] = 0.0

        if isinstance(self.venousPressure, float):
            self.Pv = self.venousPressure
        else:
            self.Pv = self.venousPressure[0]

        self.Qout[0] = 0.0

        self.Qin[0] = 0.0

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):


        if self.Z == 'VesselImpedance' and self.firstRun == False:
            self.Zrt = Z1
            # # Z not time-varying activate this lines:
            # self.Z0 = Z1
            # self.firstRun = True
        elif self.firstRun == True:
            self.Zrt = self.Z0
        else:
            self.Zrt= self.Z

        self.Rcrt = self.Rc
        if self.Rc == None:
            self.Rcrt = self.Rtotal - self.Zrt

        self.pressure[nmem] = P - Q*self.Zrt
        self.Qin[nmem] = Q

        # TODO: Polymorphic scalar/array variables?
        if isinstance(self.venousPressure, float):
            self.Pv = self.venousPressure
        else:
            self.Pv = self.venousPressure[nmem]

        self.Qout[nmem] = (P-Q*self.Zrt - self.Pv)/self.Rcrt
        self.volume[nmem] = self.volume[nmem-1] + (Q - self.Qout[nmem])*dt

        return self.returnFunction(_domegaField_, R, dt, P, Q, nmem, n)

    def residualPQ(self, dP, dQ, P, Q, dt,nmem, n):
        dPW = P + dP - self.pressure[nmem] - (Q+dQ)*self.Zrt
        return self.C*dPW/dt + (self.pressure[nmem]+dPW-self.Pv)/self.Rcrt - (Q+dQ)


class Windkessel2(BoundaryConditionType2):
    """
    Boundary profile - type 2

    2 Element Windkessel

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2

        self.Rc = 1
        self.C = 0

        self.venousPressure = 7.*133

        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L,nmem,  n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, R, dt, P, Q, n)

    def funcPos0(self, _domegaField, R, dt, P, Q, n):
        """return function for position 0 at the start
        of the vessel
        """

        """Knut Petter: The old version of the WK2 scheme is shown in the comments so you can compare it to the new version,
         the new scheme is based on the "half step central difference"-scheme shown in my master thesis, which was tested and gave good results.

        Old version:
        P0 = du[0] #/2.take this out is old ##??
        r21,r22 = R[1][0],R[1][1]
        dw_,_dw = dO[self.position][0],dO[self.position][1]

        if self.Rc == None: self.Rc = -1./R[1][1]
        tau = 2.*self.Rc*self.C
        taudt = tau/dt
        denom = 1.+taudt+r21*self.Rc
        domega_ = ((-1.-taudt-r22*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
        return np.array([domega_ , _domega_])

        """
        """ Version does not take hole R matrix into account
        r21,r22 = R[1][0],R[1][1]
        taudt = 2*self.Rc*self.C/dt
        a = 1 + taudt + self.Rc*r21
        b = 1 + taudt + self.Rc*r22
        c = 2*(self.Rc*Q + P - self.venousPressure)

        domega_ = (-b*_domegaField_ - c)/a
        #return np.array([domega_ , _domega])
        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domegaField_
        return self.omegaNew
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        taudt = self.Rc * self.C / dt
        a = -self.Rc * r21 - (1. + 2.*taudt) * r11
        b = (2.*taudt + 1.) * r12 + self.Rc * r22

        domega_ = (2 * (self.Rc * Q + (P - self.venousPressure[n])) + b * _domegaField) / a

        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut


    def funcPos1(self, domegaField_, R, dt, P, Q, n):
        """
        return function for position -1 at the end
        of the vessel
        """

        """Old version:
        P0 = du[0] #/2.take this out is old ##??
        r21,r22 = R[1][0],R[1][1]
        dw_,_dw = dO[self.position][0],dO[self.position][1]

        if self.Rc == None: self.Rc = 1./R[1][0]
        tau = 2.*self.Rc*self.C
        taudt = tau/dt
        denom = 1.+taudt-r22*self.Rc
        _domega = ((-1.-taudt+r21*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
        return np.array([_domega_ , _domega])
        """
        """ Version does not take hole R matrix into account
        r21,r22 = R[1][0],R[1][1]
        taudt = 2*self.Rc*self.C/dt
        a = taudt - self.Rc*r21 - 1
        b = taudt - self.Rc*r22 - 1
        c = 2*(self.Rc*Q - (P - self.venousPressure))

        _domega = (-a*domega_ + c)/b
        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega
        return self.omegaNew
        #return np.array([domega_ , _domega])
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        taudt = self.Rc * self.C / dt
        a = self.Rc * r21 - (1. + 2.*taudt) * r11
        b = (2.*taudt + 1.) * r12 - self.Rc * r22
        _domega = (2 * (self.Rc * Q - (P - self.venousPressure[n])) + a * domegaField_) / b


        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domega

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew), self.dQInOut

    def solveW2(self,_domega, args):

        domegaField_,P, Q, dt,n, r11, r12, r21, r22 = args

        dP = r11*domegaField_+ r12*_domega
        dQ = r21*domegaField_+ r22*_domega

        return self.f(dP,dQ,P,Q,dt,n)

    def f(self,dP,dQ,P,Q,dt,n):
        return self.C*dP/dt + (P+dP-self.venousPressure[n])/self.Rc - (Q+dQ)

class Windkessel3(BoundaryConditionType2):
    """
    Boundary profile - type 2

    3 Element Windkessel

    call function input:
     _domega_,dO,du,R,L,n,dt
    returns the domega-vector with (domega_ , _domega) based on the input values
    and its returnFunction
    """

    def __init__(self):
        self.type = 2

        # parameters in xml
        self.Rc = None
        self.Z = 'VesselImpedance'
        self.C = 1
        self.Rtotal = None
        # parameters for calculation
        self.Z0 = None  # vesselimpedance at first step
        self.RcNum = None

        self.venousPressure = [0]  # 7.*133.
        self.returnFunction = None
        self.omegaNew = np.empty((2))
        self.dQInOut = np.empty((2))

        self.firstRun = False

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem, n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, R, dt, P, Q, Z1, Z2, n)

    def funcPos0(self, _domegaField, R, dt, P, Q, Z1, Z2, n):
        """
        return function for position 0 at the start
        of the vessel

        Old version
        venousPressure = P

        r21,r22 = R[1][0],R[1][1]
        dw_,_dw = dO[self.position][0],dO[self.position][1]

        if self.R1 == None: self.R1 = -1./r22
        if self.Rc == None: self.Rc = self.RT - self.R1
        tau = 2.*self.Rc*self.C
        taudt = tau/dt
        denom = taudt + 1. + r21*(self.R1*taudt + self.Rc + self.R1)
        domega_ = ((1.+r21*self.R1)*taudt*dw_ + (1.+r22*self.R1)*taudt*_dw + (-1.-taudt-r22*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom

        Old Version: does not take hole R matrix into account
        r21,r22 = R[1][0],R[1][1]

        Z = self.Z
        if self.Z == 'VesselImpedance':# and self.firstRun == False:
            Z = Z1 #+Z2
            #self.Z = Z
            #self.firstRun = True
        Rc = self.Rc
        if self.Rc == None:
            Rc = self.Rtotal - Z
            ## Rc not time-varying activate this line:
            #self.Rc = Rc
        C = self.C

        taudt = 2*Rc*C/dt

        a = taudt + 1
        b = (taudt +1)*Z + Rc
        c = 2*Q*(Rc + Z) - 2*(P-self.venousPressure)

        domega_ = -((a + b*r21)*_domega_ - c)/(a + b*r22)


        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega_
        return self.omegaNew
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]




        if self.Z == 'VesselImpedance' and self.firstRun == False:
            Z = Z1
            # # Z not time-varying activate this lines:
            self.Z0 = Z1
            self.firstRun = True
        elif self.firstRun == True:
            Z = self.Z0
        else:
            Z = self.Z

        Rc = self.Rc
        if self.Rc == None:
            Rc = self.Rtotal - Z
            # # Rc not time-varying activate this line:
#            self.Rc = Rc
            # #

        C = self.C

        a = Z * Rc * C + 0.5 * dt * (Z + Rc)
        b = 0.5 * dt + Rc * C

        domega_ = (dt * (P - self.venousPressure[n]) + dt * (Z + Rc) * Q + _domegaField * (a * r22 + b * r12)) / (-a * r21 - b * r11)

        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domegaField

        self.dQInOut = R[:][1] * self.omegaNew

        return np.dot(R, self.omegaNew), self.dQInOut



    def funcPos1(self, domegaField_, R, dt, P, Q, Z1, Z2, n):
        """
        return function for position -1 at the end
        of the vessel  (newer version from Knut Petter ??? )

        Knut Petter: Yes this is a newer version based on the "half-step central difference"-scheme (shown for WK2 in my master thesis, but not WK3)
        I am also putting the old version back into the comments so the two can be compared.


        Old Version:
        P0 = du[0] #/2.take this out is old ##??

        r21,r22 = R[1][0],R[1][1]
        dw_,_dw = dO[self.position][0],dO[self.position][1]

        if self.R1 == None: self.R1 = 1./r21
        if self.Rc == None: self.Rc = self.RT - self.R1
        tau = 2*self.Rc*self.C
        taudt = tau/dt
        denom =  taudt + 1. - r22*(self.R1*taudt + self.Rc + self.R1)
        _domega = ((1.-r21*self.R1)*taudt*dw_ + (1.-r22*self.R1)*taudt*_dw + (-1.-taudt+r21*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom
        return np.array([_domega_ , _domega])

        Old Version Not taking total R-matrix into account:
        taudt = 2*Rc*C/dt

        a = taudt + 1
        b = (taudt +1)*Z + Rc
        c = 2*Q*(Rc + Z) - 2*(P-self.venousPressure)

        _domega = -((a + b*r21)*domega_ - c)/(a + b*r22)

        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega
        return self.omegaNew
        """
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]


        if self.Z == 'VesselImpedance' and self.firstRun == False:
            Z = Z1
            # # Z not time-varying activate this lines:
            self.Z0 = Z1
            self.firstRun = True
        elif self.firstRun == True:
            Z = self.Z0
        else:
            Z = self.Z

        Rc = self.Rc
        if self.Rc == None:
            Rc = self.Rtotal - Z
            # # Rc not time-varying activate this line:
            # self.Rc = Rc

        C = self.C

        a = Z * Rc * C + 0.5 * dt * (Z + Rc)
        b = 0.5 * dt + Rc * C

        _domega = (dt * (P - self.venousPressure[n]) - dt * (Z + Rc) * Q + domegaField_ * (b * r11 - a * r21)) / (a * r22 - b * r12)

        self.omegaNew[0] = domegaField_
        self.omegaNew[1] = _domega

        self.dQInOut = (R[:][1] * self.omegaNew)[::-1].copy()

        return np.dot(R, self.omegaNew), self.dQInOut



class L_network(BoundaryConditionType2):
    """
    Boundary profile - type 2

    L-network

    call function input:
     _domega_,dO,du,R,L,n,dt
    return:
    the domega-vector with (domega_ , _domega) based on the input values and its returnFunction
    """

    def __init__(self):
        self.type = 2

        self.C = 0
        self.R1 = 1

        self.returnFunction = None
        self.omegaNew = np.empty((2))

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem, n, dt, P, Q, A, Z1, Z2):
        return self.returnFunction(_domegaField_, duPrescribed, R, L, n, dt)

    def funcPos0(self, _domega_, du, R, dt):
        """
        return function for position 0 at the start
        of the vessel
        """
        raise NotImplementedError("boundaryCondition Lnet is not implemented correct!")
        #print "ERROR: boundaryCondition Lnet is not implemented correct!"
        #exit()

        dQ0 = du[1] / 2.0

        r21, r22 = R[1][0], R[1][1]
        dw_, _dw = 23, 23  # # dO[self.position][0],dO[self.position][1]
        if self.R1 == None: self.R1 = -1. / R[1][1]
        tau = 2.*self.R1 * self.C
        taudt = tau / dt
        denom = taudt + taudt * r21 * self.R1 + r21 * self.R1
        domega_ = ((1. + r21 * self.R1) * taudt * dw_ + (1. + r22 * self.R1) * taudt * _dw + (-taudt - taudt * r22 * self.R1 - r22 * self.R1) * _domega_ + self.R1 * dQ0) / denom
        # return np.array([domega_ , _domega_])
        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega_

        return np.dot(R, self.omegaNew)

    def funcPos1(self, _domega_, du, R, dt):
        """
        return function for position -1 at the end
        of the vessel
        """
        raise NotImplementedError("boundaryCondition Lnet is not implemented correct!")

        #print "ERROR: boundaryCondition Lnet is not implemented correct!"
        #exit()

        dQ0 = du[1] / 2.0
        r21, r22 = R[1][0], R[1][1]
        dw_, _dw = 23, 23  # dO[self.position][0],dO[self.position][1]
        if self.R1 == None: self.R1 = 1. / R[1][0]
        tau = 2.*self.R1 * self.C
        taudt = tau / dt
        denom = -taudt + taudt * r22 * self.R1 + r22 * self.R1
        _domega = ((-1. + r21 * self.R1) * taudt * dw_ + (-1. + r22 * self.R1) * taudt * _dw + (taudt - taudt * r21 * self.R1 - r21 * self.R1) * _domega_ + self.R1 * dQ0) / denom

        self.omegaNew[0] = _domega_
        self.omegaNew[1] = _domega

        return np.dot(R, self.omegaNew)


class VaryingElastance(BoundaryConditionType2):

    """
    An implementation of a time-varying elastance model of the left ventricle
     based on the modfied varying elastance equation including a source
     resistance K (as proposed by Shroff), and a parametrized time varying
     elastance function (as given by Stergiopulos).

    The general shape of the elastance function is given by three shape
    parameters. Various conditions of heart contractility are then created
    by scaling this function using the parameters:
    T - Heart period
    Emax - Maximum elastance
    Emin - Minimum elastance
    Tpeak - Time to peak elastance

    The equation also requires the two constants:
    V0 - Volume axis intercept
    K - Source resistance

    NB! The source resistance (K) was introduced to the implementation as
    an experiment. It modifies the elastance curve depending on ventricular
    outflow, so that it becomes dependent on the afterload of the heart.
    Introducing the source resistance did produce a load dependence (shown
    by curved isochrones in p-v loops), however the results are in no way to
    be trusted since the modified varying elastance equation was not intended
    to be used together with the specific varying elastance curve shape used
    here (Stergiopulos). A proper implementation of the source resistance
    requires a different curve shape. The parameter K is therefore set to zero
    by default, but it should perhaps be removed from the code altogether?

    Currently only the return method "def funcPos0" has been implemented so
    that the boundary condition can only be put at the proximal end of a blood
    vessel. It is fairly straightforward to implement funcPos1 if necessary,
    this does however require a lot of duplicated code.
    """
    def __init__(self):
        self.type = 2

        self.subiterations = 0

        self.omegaNew = np.empty((2))

        # Default parameters
        self.T = 1
        self.Emax = 2.31 * 133.3e6
        self.Emin = 0.06 * 133.3e6
        self.Tpeak = 0.4

        self.V0 = 20e-6

        self.K = 0.0

        """Shape parameters"""
        self.alpha = 1.672
        self.n1 = 1.32
        self.n2 = 21.9

        # n-1 values
        self.aorticFlowPreviousTimestep = None

        self.system = {'both open':np.array([0, 1, 2]), 'mitral open': np.array([0, 1]), 'aortic open':np.array([1, 2])}

        self.newCycle = False
        self.cycleNumber = 0
        self.num = 0
        self.atriumPressure = 7.5 * 133.32  # Pressure in the atrium ## venouse pressure?!

        self.x0 = np.array([0.0, 0.0, 0.0])  # Initial values for the iterative solver

        self.mitral = None  # mitral valve
        self.aortic = None  # aortic valve
        self.initializeValves()  # intialize valves

        self.dQInOut = np.empty((2))


    def initializeSolutionVectors(self, Tsteps):
        """Initializes some solution vectors storing pressure, flow and volume of the ventricle, as well as opening and closing state

        NB! This method is not called from the class constructor, but is called externally by the initializeSolutionMatrices method in the solver,
        this is a bit messy, but was the easiest way to do it since the BC is initiated before the number of time steps is known.
        """

        """ Initialize Solution Vectors """

        print """ Initialize Solution Vectors """

        self.mitral.initializeSolutions(Tsteps+1)
        self.aortic.initializeSolutions(Tsteps+1)

        self.pressure = np.zeros(Tsteps+1)
        self.volume = np.zeros(Tsteps+1)
        self.mitralQ = np.zeros(Tsteps+1)
        self.Elastance = np.zeros(Tsteps+1)  # New
        self.Flow = np.zeros(Tsteps+1)
        self.DtFlow = np.zeros(Tsteps+1)
        self.Turb = np.zeros(Tsteps+1)
        self.Inert = np.zeros(Tsteps+1)
        self.InbyTurb = np.zeros(Tsteps+1)
        self.deltaP = np.zeros(Tsteps+1)
        self.aortaP = np.zeros(Tsteps+1)

        """ Initial conditions in the ventricle"""
        self.pressure[0] = self.atriumPressure
        self.volume[0] = self.atriumPressure / self.E(0) + self.V0


    def initializeValves(self):
        """Mitral valve parameters"""

        self.mitral_annulus_area = 0.0006

        mitral_M_st = 1
        mitral_M_rg = 0.0
        mitral_delta_p_open = 0
        mitral_delta_p_close = 0
        mitral_K_v_open = 0.3
        mitral_K_v_close = 0.4


        """Aortic valve parameters"""
        aortic_M_st = 1
        aortic_M_rg = 0
        aortic_delta_p_open = 0 * 133.32
        aortic_delta_p_close = 0 * 133.32  # 2mmHg
        aortic_K_v_open = 0.12
        aortic_K_v_close = 0.12


        # Create valves
        self.mitral = Valve(mitral_M_st, mitral_M_rg, mitral_delta_p_open, \
                        mitral_delta_p_close, mitral_K_v_open, mitral_K_v_close)
        self.aortic = Valve(aortic_M_st, aortic_M_rg, aortic_delta_p_open, \
                        aortic_delta_p_close, aortic_K_v_open, aortic_K_v_close)

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem,  n, dt, P, Q, A, Z1, Z2):

        self.newCycle = False
        self.updateValves(P, n, dt)  # Update the state of the mitral and aortic valve at timestep n + 1
        self.startNewCycleIfCriteriaIsMet(n, dt)
        self.funcPos0(_domegaField_, R, n, dt, P, Q, A)  # Compute the riemann variant going into the vessel save in omegaNew

        self.dQInOut = R[:][1] * self.omegaNew
        # calculate du and return this!
        return np.dot(R, self.omegaNew), self.dQInOut


    def updateValves(self, P, n, dt):
        mitralPressureDifference = self.atriumPressure - self.pressure[n]
        self.mitral.updateValveState(mitralPressureDifference, n, dt)

        aorticPressureDifference = self.pressure[n] - P
        self.aortic.updateValveState(aorticPressureDifference, n, dt)


    def getCycleTime(self, n, dt):
        return self.num * dt

    def startNewCycleIfCriteriaIsMet(self, n, dt):
        if self.getCycleTime(n + 1, dt) > self.T:
            self.cycleNumber += 1
            self.num = 0
            self.newCycle = True

    def funcPos0(self, _domega, R, n, dt, Pn, Qn, A):

        # Qn1 == value at old time step
        # change to self.aorticPressurePreviousTimestep ...

        Qn1 = self.aorticFlowPreviousTimestep

        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]

        LdivB = self.aortic.LdivideB(A, n + 1)
        B = self.aortic.computeB(A, n + 1)
        L = self.aortic.computeL(A, n + 1, B, self.aortic.state[n - 1])
        mitrL = self.mitral.computeL(self.mitral_annulus_area, n + 1, B, self.aortic.state[n - 1])
        mitrLdivB = self.mitral.LdivideB(A, n + 1)
        mitrB = self.mitral.computeB(self.mitral_annulus_area, n + 1)  #
        mitrQn = self.mitralQ[n]
        mitrQn1 = self.mitralQ[n - 1]
        venoP = self.atriumPressure
        t = self.getCycleTime(n + 1, dt)
#         ttemp = t-dt
#         print "n is",n
#         print "self num is", self.num
#         print "time is",t
#         print "numtime is", self.num*dt
        E = self.E(t)
        Vn = self.volume[n]
        self.Elastance[n + 1] = E / 133.3e6
        self.Flow[n] = Qn * 1e6
        self.aortaP[n] = Pn
#         self.DtFlow[n]=(Qn-Qnold)/dt
        ventrPn = self.pressure[n]
        if self.cycleNumber == 3:
            self.T = 0.7
            self.Tpeak = 0.43 * self.T

        if self.cycleNumber == 6:
            self.T = 1
            self.Tpeak = 0.43 * self.T

        if self.cycleNumber == 10:
            self.T = 0.8
            self.Tpeak = 0.43 * self.T

        if self.cycleNumber == 11:
            self.T = 0.7
            self.Tpeak = 0.43 * self.T

        if self.cycleNumber == 12:
            self.T = 0.6
            self.Tpeak = 0.43 * self.T

        if self.cycleNumber == 14:
            self.T = 0.5
            self.Tpeak = 0.43 * self.T

        if B:

            self.Turb[n] = abs(Qn) * Qn * B / 133
            self.Inert[n] = L * (Qn - Qn1) / (dt * 133)
            self.deltaP[n] = abs(Qn) * Qn * B / 133 + L * (Qn - Qn1) / (dt * 133)



        """ Because of the large differences in magnitude between pressure and flow in the currently used dimensions some attemts were made to scale the
        variables using for example these scaling factors for B,p, q,and omega"""
        B_ref = 1060 / (2 * A ** 2)
        n_p = 1.  # self.Emax*self.V0
        n_q = 1.  # (n_p/B_ref)**0.5
        n_o = 1.  # n_q/r21
#         if Qn<-20e-6:
#         B=None

        args = dt, mitrLdivB, mitrB, LdivB, L, mitrL, B, mitrQn1, mitrQn, ventrPn, venoP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p, B_ref


        """The following section computes the increment domega_ which goes into the vessel from the ventricle, """

        if not B and not mitrB:
            """both valves are closed: """


            self.mitralQ[n + 1] = 0
            domega_ = _domega
            self.volume[n + 1] = self.volume[n] - 0.5 * (Qn - mitrQn) * dt
            self.pressure[n + 1] = E * (self.volume[n + 1] - self.V0)
            if Qn == 0:
                domega_ = _domega
            else:
                domega_ = (-0.5 * Qn - r22 * _domega) / r21  # Correction for non-zero flow at full aortic valve closure

            self.x0 = np.array([0, 0, domega_])
        else:
            if not mitrB:
                """only the aortic valve is open"""
                x = self.newtonSolver(self.x0, args, partialSystem='aortic open')
                self.mitralQ[n + 1] = 0
                self.pressure[n + 1] = self.pressure[n] + x[0] * n_p
                domega_ = x[1] * n_o
                self.x0 = np.concatenate((np.array([0]), x))

            elif not B:
                """only  the mitral valve is open"""
                x = self.newtonSolver(self.x0, args, partialSystem='mitral open')
                self.mitralQ[n + 1] = self.mitralQ[n] + x[0] * n_q
                self.pressure[n + 1] = self.pressure[n] + x[1] * n_p


                if Qn == 0:
                    domega_ = _domega
                else:
                    domega_ = (-0.5 * Qn - r22 * _domega) / r21  # Correction for non-zero flow at full aortic valve closure

                self.x0 = np.concatenate((x, np.array([0])))
            else:
                """both valves are open"""
                x = self.newtonSolver(self.x0, args)
                self.mitralQ[n + 1] = self.mitralQ[n] + x[0] * n_q
                self.pressure[n + 1] = self.pressure[n] + x[1] * n_p
                domega_ = x[2] * n_o
                self.x0 = x

            dQ = r21 * domega_ + r22 * _domega
            self.volume[n + 1] = Vn - (Qn + 0.5 * (-self.mitralQ[n] - self.mitralQ[n + 1] + dQ)) * dt

        self.aorticFlowPreviousTimestep = Qn
#         Qnold=Qn
        self.num = self.num + 1
        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega




    def funcPos1(self, _domega, R, L, n, dt, P, Q, A):
        pass

    def E(self, t):
        """Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
           according to Tpeak, T, Emax and Emin. """
        a1 = 0.708 * self.Tpeak
        a2 = 1.677 * a1

        n1, n2 = self.n1, self.n2
        shapeFunction1 = (t / (a1)) ** n1 / (1 + (t / (a1)) ** n1)
        shapeFunction2 = (1 + (t / (a2)) ** n2) ** (-1)
        return (self.Emax - self.Emin) * self.alpha * shapeFunction1 * shapeFunction2 + self.Emin

    def newtonSolver(self, x0, args, partialSystem='both open'):
        """Solves the partial or full equation system of the varying elastance model"""
        maxIterations = 20
        iterations = 0

        xn = x0[self.system[partialSystem]]
        res = self.solverResiduals(xn, *args, partialSystem=partialSystem)
        # error = np.linalg.norm(res, 2)

        while True:
            # print x0
            iterations += 1
            J_inv = self.solverInverseJacobian(xn, *args, partialSystem=partialSystem)



            x = xn - np.dot(J_inv, res).T


            error = np.linalg.norm(x - xn, 2) / np.linalg.norm(xn, 2)
            if error < 0.0001:
                break
            xn = x
            res = self.solverResiduals(x, *args, partialSystem=partialSystem)

            if iterations > maxIterations:
                x *= 0
                break

        return x

    def solverResiduals(self, x_partial, dt, mitrLdivB, mitrB, LdivB, L, mitrL, B, mitrQn1, mitrQn, ventrPn, atrP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p, B_ref, partialSystem=np.array([0, 1, 2])):  # dt, mitrL, mitrB,L,B, mitrQn1, mitrQn, venoP, E, Vn, Qn, Qn1, r21, r22, Pn, _domega):
        """Computes  are the resisduals of the functions f1,f2 and f3, they are defined as functions that are only called when they are needed. The
        argument partialSystem determines which of the residuals are computed and returned."""

        x = np.array([0.0, 0.0, 0.0])
        x[self.system[partialSystem]] += x_partial
        dQm, dPv, domega_ = x

#         Knut Petters version:
        """
        def f1():
            a = mitrQn/n_q + dQm
            return a*abs(a) + mitrLdivB/(2*n_q*dt)*(3*dQm + (mitrQn1 - mitrQn)/n_q) + (n_p*dPv + ventrPn - atrP)/(mitrB*n_q**2)
        def f2():
            return E/n_p*(Vn - (Qn - mitrQn + 0.5*(n_q*domega_ + r22*_domega - n_q*dQm))*dt - self.V0)*(1-self.K*(Qn + n_q*domega_ + r22*_domega)) - ventrPn/n_p - dPv
        def f3():
            a = (Qn + r22*_domega)/n_q + domega_
            return a*abs(a) + LdivB/(2*n_q*dt)*(3*domega_ + (3*r22 - Qn + Qn1)/n_q) + (n_q/r21*domega_ + _domega + Pn - ventrPn - n_p*dPv)/(B*n_q**2)
        """
        # Simple, K is set to 0
        def f1():
            a = mitrQn / n_q + dQm
            return mitrB * a * abs(a) + (0.5 / dt) * mitrL * (3 * dQm + (mitrQn1 - mitrQn) / n_q) + (n_p * dPv + ventrPn - atrP)
        def f2():
            return E / n_p * (Vn - (Qn - mitrQn + 0.5 * (n_q * domega_ * r21 + r22 * _domega - n_q * dQm)) * dt - self.V0) * (1 - self.K * (Qn + n_q * domega_ * r21 + r22 * _domega)) - ventrPn / n_p - dPv

        def f3():
            a = (Qn + r22 * _domega) / n_q + domega_ * r21
            return B * a * abs(a) + (0.5 / dt) * L * (3 * domega_ * r21 + (3 * r22 * _domega - Qn + Qn1) / n_q) + (r11 * domega_ + _domega * r12 + Pn - ventrPn - n_p * dPv)

        def f1simple():
            a = mitrQn + dQm
            return mitrB * a * abs(a) + 0.5 * mitrL * (3 * dQm + mitrQn1 - mitrQn) / dt - atrP + ventrPn + dPv

        def f2simple():
            return E * (Vn - 0.5 * dt * (r21 * domega_ + r22 * _domega + 2 * Qn - 2 * mitrQn - dQm) - self.V0) - ventrPn - dPv

        def f3simple():
            a = r21 * domega_ + r22 * _domega + Qn
            return B * a * abs(a) + (L / (2 * dt)) * (3 * r21 * domega_ + 3 * r22 * _domega - Qn + Qn1) - ventrPn - dPv + Pn + r11 * domega_ + _domega * r12

        functions = np.array([f1simple, f2simple, f3simple])
        return np.array([f() for f in functions[self.system[partialSystem]]])

    def solverInverseJacobian(self, x_partial, dt, mitrLdivB, mitrB, LdivB, L, mitrL, B, mitrQn1, mitrQn, ventrPn, venoP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p, B_ref, partialSystem=np.array([0, 1, 2])):
        """Computes the inverse Jacobian of the system. The components a1, a2, a3, a4, a5 and a6 are declared using strings, and evaluated using eval()
         only when needed. (Not sure how smart this is). Reduces lines of code, but is probably slower."""


        x = np.array([0, 0, 0])
        x[self.system[partialSystem]] += x_partial
        dmQ, dvP, domega_ = x

        if partialSystem == 'mitral open':
            # simple -> K is not part of the jacobi
            a1 = mitrB * 2 * (mitrQn / n_q + dmQ) * np.sign(mitrQn / n_q + dmQ) + 1.5 * mitrL / dt
            a2 = 0.5 * E * dt * (1 - self.K * (r21 * domega_ + r22 * _domega + Qn))

            a1simple = mitrB * 2 * abs(mitrQn + dmQ) + 1.5 * mitrL / dt
            a2simple = 0.5 * E * dt

            J_inv = np.array([[1, 1],
                              [a2simple, -a1simple]]) / (a1simple + a2simple)
            return J_inv
        elif partialSystem == 'aortic open':
            a3 = -0.5 * E * dt * (1 - self.K * (r21 * domega_ + r22 * _domega + Qn)) - r21 * E * self.K * (Vn - (Qn - mitrQn + 0.5 * (n_q * domega_ * r21 + r22 * _domega - n_q * dmQ)) * dt - self.V0)
            a4 = B * r21 * 2 * (Qn + r21 * domega_ + r22 * _domega) * np.sign(Qn + r21 * domega_ + r22 * _domega) + 1.5 * L * r21 / dt + r11

            a3simple = -0.5 * E * dt * r21
            a4simple = 2 * r21 * B * abs(r21 * domega_ + r22 * _domega + Qn) + 1.5 * L * r21 / (dt) + r11

            J_inv = np.array([[a4simple, -a3simple],
                             [1, -1]]) / (-a4simple + a3simple)
            return J_inv
        else:
            a1 = mitrB * 2 * (mitrQn / n_q + dmQ) * np.sign(mitrQn / n_q + dmQ) + 1.5 * mitrL / dt
            a2 = 0.5 * E * dt * (1 - self.K * (r21 * domega_ + r22 * _domega + Qn))
            a3 = -0.5 * E * dt * (1 - self.K * (r21 * domega_ + r22 * _domega + Qn)) - r21 * E * self.K * (Vn - (Qn - mitrQn + 0.5 * (n_q * domega_ * r21 + r22 * _domega - n_q * dmQ)) * dt - self.V0)
            a4 = B * r21 * 2 * (Qn + r21 * domega_ + r22 * _domega) * np.sign(Qn + r21 * domega_ + r22 * _domega) + 1.5 * L * r21 / dt + r11

            J_inv = np.array([[ -a4 + a3, -a4, a3  ],
                              [ -a2 * a4, a1 * a4, -a1 * a3 ],
                              [  -a2, a1, a1 + a2 ]]) / (a1 * a3 - a1 * a4 - a2 * a4)


            return J_inv


        # Version from Master Thesis Knut Petter:
        """Version from Master Thesis Knut Petter
        expressions = np.array([
        '2*(mitrQn/n_q + dmQ)*np.sign(mitrQn/n_q + dmQ) + 1.5*mitrLdivB/(n_q*dt)',
        'B_ref/mitrB',
        '0.5*n_q/n_p*E*dt*(1-self.K*(Qn + n_q*domega_ + r22*_domega))',
        '-0.5*n_q/n_p*E*dt*(1-self.K*(Qn + n_q*domega_ + r22*_domega)) - E*self.K*n_q/n_p*(Vn -(Qn - mitrQn + 0.5*(n_q*domega_ + r22*_domega - n_q*dmQ))*dt - self.V0)',
        '-B_ref/B',
        '2*((Qn + r22*_domega)/n_q + domega_)*np.sign((Qn + r22*_domega)/n_q + domega_) + 1.5*LdivB/(n_q*dt) + 1/(r21*B*n_q)'])






        if partialSystem == 'mitral open':
            a1, a2, a3 = [eval(e) for e in expressions[[0,1,2]]]

            J_inv = np.array([[-1,  -a2],
                              [-a3, a1]])/(-a1-a2*a3)
            return J_inv
        elif partialSystem == 'aortic open':
            a4, a5,a6 = [eval(e) for e in expressions[[3,4,5]]]
            J_inv =np.array([[a6, -a4],
                             [-a5, -1]])/(-a6 - a4*a5)
            return J_inv
        else:
            a1, a2, a3, a4, a5, a6 = [eval(e) for e in expressions]

            J_inv = np.array([[ a6+a4*a5,  a2*a6,  -a2*a4   ],
                              [    a3*a6, -a1*a6,   a1*a4   ],
                              [   -a3*a5,  a1*a5,   a1+a2*a3]])/(a3*a2*a6 + a1*a6 + a1*a4*a5)


            return J_inv

        """


class Valve:
    """ A valve model meant to be used with the varying elastance boundary condition. Two instances of the class are then made representing the mitral and
    aortic valves. The valve instance has a variable state, which is updated every timestep using the updateValveState method and is then able to return
    updated values of the coefficients B and L, used by the varying elastance solver.
    """

    def __init__(self, M_st, M_rg, delta_p_open, delta_p_close, K_v_open, K_v_close):

        """
        M_st - Determines the effective orifice area (EOA) at maximum opening. M_st ~= 1 => Healthy valve
                                                                               M_st  < 1 => Stenosed valve
                                                                               (If set to 1 there is no pressure difference, meaning that the valve can not close
                                                                               once it has fully opened, 0.99 has been used)

        M_rg - Determines the EOA at minimum opening. M_rg = 0 => Healthy valve
                                                      M_rg > 0 => Regurgitant valve

        """
        self.M_st = M_st
        self.M_rg = M_rg
        self.delta_p_open = delta_p_open
        self.delta_p_close = delta_p_close
        self.K_v_open = K_v_open
        self.K_v_close = K_v_close

        self.rho = 1060  # # default value is update while start up in classVascularNetwork.initialize()

    def initializeSolutions(self, Tsteps):
        """This method is called by the initializeSolutionVectors method in the VaryingElastance-instance to initialize the vector containing the state variable,
        the state variable is stored so that it can be used in visualization to show the opening and closing of the valve """
        self.state = np.zeros(Tsteps)

    def computeB(self, A, n):
        """ Returns the turbulent resistance coefficient B, used in computing the pressure difference across the valve"""
        A_eff = self.effectiveOrificeArea(A, n)



        if A_eff == 0:
            B = None
        elif A / A_eff > 1e4:
            B = None
        else:
#             B = 5*0.5*self.rho*(1/A_eff - 1/A)**2
#             B = 0.5*self.rho*(1/A_eff - 1/A)  # expression: Masterthesis Knut Petter
            B = self.rho / (2 * A_eff ** 2)  # paper Mynard et al. 2012
        return B

    def computeL(self, A, n, B, state):
        """ Returns the inertance coefficient L,  used in computing the pressure difference across the valve"""
        A_s = self.effectiveOrificeArea(A, n)
        if B:
            # leff=0.008+0.01*(1-state) # aortic
            leff = 0.018  # 2cm
        else:
            leff = 0.01  # mitral

        if A_s == 0:
            L = None
        elif A / A_s > 1e4:
            L = None
        else:
#             L = 4*np.pi*self.rho*(1/A_s - 1/A)**0.5 # expression: Masterthesis Knut Petter
            L = self.rho * leff / A_s  # paper Mynard et al. 2012
        return L


    def LdivideB(self, A, n):
        A_s = self.effectiveOrificeArea(A, n)
        if A_s == 0:
            return None
        elif A / A_s > 1e4:
            return None
        else:
            return 4 * np.pi * (1 / A_s - 1 / A) ** (-1.5)

    def effectiveOrificeArea(self, A, n):
        """ Computes the effective orifice area (A_eff) by interpolating between maximum and minimum area."""
        A_max = self.M_st * A
        A_min = self.M_rg * A

        return (A_max - A_min) * self.state[n] + A_min

    def A_max(self, A):
        return self.M_st * A

    def A_min(self, A):
        return self.M_rg * A

    def updateValveState(self, delta_p, n, dt):
        """The opening state at timestep n+1 is computed by explicit numerical solution of the opening and closing rate equations (Mynard) """

        if delta_p > self.delta_p_open:
            """The valve starts to open if the pressure difference is positive and above the opening threshold pressure"""
            if self.state[n] == 1.0:
                self.state[n + 1] = 1.0
            else:
                self.state[n + 1] = self.state[n] + (1 - self.state[n]) * self.K_v_open * (delta_p - self.delta_p_open) * dt
                if self.state[n + 1] > 1.0:
                    self.state[n + 1] = 1.0

        elif delta_p < -self.delta_p_close:
            """ The valve starts to close """
            if self.state[n] == 0.0:
                self.state[n + 1] = 0.0
            else:
                self.state[n + 1] = self.state[n] + self.state[n] * self.K_v_close * (delta_p - self.delta_p_close) * dt
                if self.state[n + 1] < 0:
                    self.state[n + 1] = 0.0
        else:
            """ The case where -delta_p_close < delta_p < delta_p_open, the valve state stays unchanged """
            self.state[n + 1] = self.state[n]

class VaryingElastanceSimple(BoundaryConditionType2):
    """
    An implementation of a time-varying elastance model of the left ventricle,ing based on the modfied varying elastance equation including a source resistance K
    (as proposed by Shroff), and a parametrized time varying elastance function (as given by Stergiopulos).

    The general shape of the elastance function is given by three shape parameters. Various conditions of heart contractility are then created by scaling this
    function using the parameters:
    T - Heart period
    Emax - Maximum elastance
    Emin - Minimum elastance
    Tpeak - Time to peak elastance

    The equation also requires the two constants:
    V0 - Volume axis intercept
    K - Source resistance

    NB! The source resistance (K) was introduced to the implementation as an experiment. It modifies the elastance curve depending on ventricular outflow,
    so that it bself.mitral = None # mitral valve
        self.aortic = None # aortic valve
        self.initializeValves() # intialize valves becomes dependent on the afterload of the heart. Introducing the source resistance did produce a load dependence
    (shown by curved isochrones in p-v loops), however the results are in no way to be trusted since the modified varying elastance equation was not intended
    to be used together with the specific varying elastance curve shape used here (Stergiopulos). A proper implementation of the source resistance
    requires a different curve shape. The parameter K is therefore set to zero by default, but it should perhaps be removed from the code altogether??

    Currently only the return method "def funcPos0" has been implemented so that the boundary condition can only be put at the proximal end of a blood vessel.
    It is fairly straightforward to implement funcPos1 if necessary, this does however require a lot of duplicated code.   """

    solutionMemoryFields    = ["pressure", "atriumPressure", "volume", "mitralQ", "Elastance", "Flow", "Flow2", "deltaP", "aortaP"]
    solutionMemoryFieldsToSave = ["pressure", "atriumPressure", "volume", "mitralQ", "Elastance", "Flow", "Flow2", "deltaP", "aortaP"]

    def __init__(self):
        self.type = 2

        self.subiterations = 0

        self.omegaNew = np.empty((2))

        # Default parameters (NOT Time VARYING)
        self.T = 1
        self.Emax = 2.31 * 133.3e6
        self.Emin = 0.06 * 133.3e6
        self.Tpeak = 0.4

        # BRX update variables
        self.T_BRX = self.T
        self.Emax_BRX = self.Emax
        self.Emin_BRX = self.Emin

        # Runtime Variables
        self.rt_T = self.T
        self.rt_Emax = self.Emax
        self.rt_Emin = self.Emin
        self.rt_Tpeak = self.Tpeak

        self.V0 = 20e-6

        self.K = 0.0
        # self.Rv = 0.00007500616 * 133 / (10 ** -6) # From Lau and Figueroa paper
        self.Rv = 0.003*133.32*1e6 #0.005 * 133 / (10 ** -6)

        ## Shape parameters
        self.alpha = 1.672
        self.n1 = 1.32
        self.n2 = 21.9
        self.R11 = None
        self.R12 = None
        self.R21 = None
        self.R22 = None
        self.DtW2 = None

        # Cycle management for restarting isovolumetric phase
        self.newCycle = False
        self.cycleNumber = 0
        self.num = 0

        self.atriumPressure0 = 7.5 * 133.32  # TODO: Fix this: Pressure in the atrium ## venouse pressure?!

        self.dQInOut = np.empty((2))

        self.dsetGroup = None

        self.pressure = np.zeros(0)
        self.atriumPressure =  np.zeros(0)
        self.volume = np.zeros(0)
        self.mitralQ = np.zeros(0)
        self.Elastance = np.zeros(0)
        self.Flow = np.zeros(0)
        self.Flow2 = np.zeros(0)
        self.DtFlow = np.zeros(0)
        self.deltaP = np.zeros(0)
        self.aortaP = np.zeros(0)

    def update(self, bcDict):
        super(VaryingElastanceSimple,self).update(bcDict)

        # BRX update variables
        self.T_BRX = self.T
        self.Emax_BRX = self.Emax
        self.Emin_BRX = self.Emin
        # Runtime Variables
        self.rt_T = self.T
        self.rt_Emax = self.Emax
        self.rt_Emin = self.Emin
        self.rt_Tpeak = self.Tpeak


    def initializeSolutionVectors(self, runtimeMemoryManager, solutionDataFile):
        """Initializes some solution vectors storing pressure, flow and volume
        of the ventricle, as well as opening and closing state

        NB! This method is not called from the class constructor, but is
        called externally by the initializeSolutionMatrices method in the
        solver, this is a bit messy, but was the easiest way to do it since the
        BC is initiated before the number of time steps is known.
        """

        """ Initialize Solution Vectors """

        print """ Initialize Solution Vectors """

        self.dsetGroup = solutionDataFile.create_group('Heart')
        self.allocate(runtimeMemoryManager)

        """ Initial conditions in the ventricle"""
        self.pressure[0] = self.atriumPressure0
        self.volume[0]   = self.atriumPressure0 / self.E(0) + self.V0
        self.atriumPressure[0] = self.atriumPressure0

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem,  n, dt, P, Q, A, Z1, Z2):

#         self.updateValves(P, n, dt)
        self.newCycle = False                     # Update the state of the mitral and aortic valve at timestep n + 1
        self.startNewCycleIfCriteriaIsMet(n, dt)
        self.funcPos0(_domegaField_, R, nmem, n, dt, P, Q, A)  # Compute the riemann variant going into the vessel save in omegaNew

        self.dQInOut = R[:][1] * self.omegaNew
        # calculate du and return this!
        return np.dot(R, self.omegaNew), self.dQInOut

    def getCycleTime(self, n, dt):
        return self.num * dt

    def startNewCycleIfCriteriaIsMet(self, n, dt):
        if self.getCycleTime(n + 1, dt) > self.T:
            self.cycleNumber += 1
            self.num = 0
            self.newCycle = True
            self.rt_T = self.T_BRX
            self.rt_Tpeak = 0.4*self.rt_T
            self.rt_Emax = self.Emax_BRX
            self.rt_Emin = self.Emin_BRX

    def funcPos0(self, _domega, R, nmem, n, dt, Pn, Qn, A):

        # Qn1 == value at old time step
        # change to self.aorticPressurePreviousTimestep ...

#         Qn1 = self.aorticFlowPreviousTimestep

        L = np.linalg.inv(R)
        L11, L12, L21, L22 = L[0][0], L[0][1], L[1][0], L[1][1]
        omegaprevious_ = L11 * Pn + L12 * Qn
        r11, r12, r21, r22 = R[0][0], R[0][1], R[1][0], R[1][1]
    #     deltatdiff = 0.00001
#         mitrQn = self.mitralQ[n]
#         mitrQn1 = self.mitralQ[n-1]
        venoP = self.atriumPressure[nmem]
        t = self.getCycleTime(n + 1, dt)
        t2 = self.getCycleTime(n, dt)
#         ttemp = t-dt
        E = self.E(t)
        e2 = self.E(t)


    #     dE= (self.E(t+deltatdiff) -E)/deltatdiff
        Vn = self.volume[nmem]
        self.R11 = r11
        self.R12 = r12
        self.R21 = r21
        self.R22 = r22
        if _domega<0.0:
            self.DtW2 = _domega / dt
        else:
            self.DtW2 = 0.0

        self.Elastance[nmem + 1] = E
        self.Flow[nmem] = Qn
        self.aortaP[nmem] = Pn
#         self.DtFlow[n]=(Qn-Qnold)/dt
        ventrPn = self.pressure[nmem]

        def diastole(u, t):
            """Differential equations during diastole u[0]=V, u[1]=Pv, return dV/dt and dP/dt"""
            u = u.reshape(2)
            deltatdiff = 0.00001
            E = self.E(t)
            dE = (self.E(t + deltatdiff) - E) / deltatdiff
            DV = (venoP - u[1]) / self.Rv
            DP = dE * (u[0] - self.V0) + E * (venoP - u[1]) / self.Rv
            return [DV, DP]

        def systole(u, t):
            """Differential equations during systole u[0]=V, u[1]=Pv=Pa, u[2]=Q return dV/dt and dP/dt, dQ/dt"""
            u = u.reshape(3)
            deltatdiff = 0.00001
            Etemp = self.E(t)
            dE = (self.E(t + deltatdiff) - Etemp) / deltatdiff

            DV = -u[2]
            DP = (dE * (u[0] - self.V0) - Etemp * u[2])
    #         #DW1 = ((dE*(u[0]-self.V0) - Etemp*u[3] - self.R12*self.DtW2)/self.R11)
            DQ = (self.R21 / self.R11) * (dE * (u[0] - self.V0) - Etemp * u[2]) + self.DtW2 * (self.R22 - (self.R12 * self.R21) / self.R11)
            return [DV, DP, DQ]

        solverSys = OD.ForwardEuler(systole)  # Runga Kutta solver for ejection phase, Systole

        solverdias = OD.ForwardEuler(diastole)  # Runga Kutta solver for diastole

        if (venoP > ventrPn):  # Diastolecondition
            solverdias.set_initial_condition([Vn, ventrPn])
            if (t2) < 0:
                t_pointsd = np.linspace(0, dt, 2)

            else:
                t_pointsd = np.linspace(t2, t2 + dt, 2)


            udiastole, td = solverdias.solve(t_pointsd)
            udiastole = udiastole[1]
            V = udiastole[0]
            Pv = udiastole[1]
            self.volume[nmem + 1] = V
            self.pressure[nmem + 1] = Pv
            if Qn == 0:
                domega_ = _domega
            else:
                domega_ = (-0.5 * Qn - r22 * _domega) / r21



        elif (Qn >= -1e-15 and ventrPn - Pn > (-0.0000001)):  # sysstolecondition
            solverSys.set_initial_condition([Vn, ventrPn, self.Flow2[nmem]])
            t_pointss = np.linspace(t2, t2 + dt, 2)
            uSystole, ts = solverSys.solve(t_pointss)
            uSystole = uSystole[1]
            V = uSystole[0]
            P = uSystole[1]
            Q = uSystole[2]
            self.Flow2[nmem + 1] = Q
            omeganew_ = L11 * P + L12 * Q



            self.volume[nmem + 1] = V
            self.pressure[nmem + 1] = P
            domega_ = omeganew_ - omegaprevious_

        else:  # isocondition
            self.volume[nmem + 1] = Vn
            self.pressure[nmem + 1] = E * (Vn - self.V0)

            domega_ = _domega


        self.num = self.num + 1
        self.omegaNew[0] = domega_
        self.omegaNew[1] = _domega


    def funcPos1(self, _domega, R, L, n, dt, P, Q, A):
        pass

    def flushSolutionData(self, saving, nDB, nDE, nSB, nSE, nSkip):

        if saving:
            self.dsetGroup['pressure'][nDB:nDE] = self.pressure[nSB:nSE:nSkip]
            self.dsetGroup['volume'][nDB:nDE] = self.volume[nSB:nSE:nSkip]
            self.dsetGroup['mitralQ'][nDB:nDE] = self.mitralQ[nSB:nSE:nSkip]
            self.dsetGroup['Elastance'][nDB:nDE] = self.Elastance[nSB:nSE:nSkip]
            self.dsetGroup['Flow'][nDB:nDE] = self.Flow[nSB:nSE:nSkip]
            self.dsetGroup['Flow2'][nDB:nDE] = self.Flow2[nSB:nSE:nSkip]
            self.dsetGroup['deltaP'][nDB:nDE] = self.deltaP[nSB:nSE:nSkip]
            self.dsetGroup['aortaP'][nDB:nDE] = self.aortaP[nSB:nSE:nSkip]

    def E(self, t):
        """
        Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
        according to Tpeak, T, Emax and Emin.
        """
        a1 = 0.708 * self.rt_Tpeak
        a2 = 1.677 * a1

        n1, n2 = self.n1, self.n2
        shapeFunction1 = (t / (a1)) ** n1 / (1 + (t / (a1)) ** n1)
        shapeFunction2 = (1 + (t / (a2)) ** n2) ** (-1)
        return (self.rt_Emax - self.rt_Emin) * self.alpha * shapeFunction1 * shapeFunction2 + self.rt_Emin

class VaryingElastanceSimpleDAE(generalPQ_BC):
    """
    TODO: Convert this to Napoleon compliant docstring
    An implementation of a time-varying elastance model of the left ventricle,ing based on the modfied varying elastance equation including a source resistance K
    (as proposed by Shroff), and a parametrized time varying elastance function (as given by Stergiopulos).

    The general shape of the elastance function is given by three shape parameters. Various conditions of heart contractility are then created by scaling this
    function using the parameters:
    T - Heart period
    Emax - Maximum elastance
    Emin - Minimum elastance
    Tpeak - Time to peak elastance

    The equation also requires:
    V0 - Volume axis intercept

    NB! The source resistance (K) was introduced to the implementation as an experiment. It modifies the elastance curve depending on ventricular outflow,
    so that it bself.mitral = None # mitral valve
        self.aortic = None # aortic valve
        self.initializeValves() # intialize valvesecomes dependent on the afterload of the heart. Introducing the source resistance did produce a load dependence
    (shown by curved isochrones in p-v loops), however the results are in no way to be trusted since the modified varying elastance equation was not intended
    to be used together with the specific varying elastance curve shape used here (Stergiopulos). A proper implementation of the source resistance
    requires a different curve shape. The parameter K is therefore set to zero by default, but it should perhaps be removed from the code altogether??

    Currently only the return method "def funcPos0" has been implemented so that the boundary condition can only be put at the proximal end of a blood vessel.
    It is fairly straightforward to implement funcPos1 if necessary, this does however require a lot of duplicated code.   """
    solutionMemoryFields    = ["pressure", "volume", "mitralQ", "aorticQ","Elastance","atriumPressure"]
    solutionMemoryFieldsToSave = ["pressure", "volume", "mitralQ", "aorticQ","Elastance","atriumPressure"]
    def __init__(self):
        self.type = 2

        self.subiterations = 0
        self.residualName = "FlowResidual"

        self.omegaNew = np.empty((2))

        # Default parameters (NOT Time VARYING)
        self.T = 1
        self.Emax = 2.31 * 133.3e6
        self.Emin = 0.06 * 133.3e6
        self.Tpeak = 0.4

        # BRX update variables
        self.T_BRX = self.T
        self.Emax_BRX = self.Emax
        self.Emin_BRX = self.Emin

        # Runtime Variables
        self.rt_T = self.T
        self.rt_Emax = self.Emax
        self.rt_Emin = self.Emin
        self.rt_Tpeak = self.Tpeak

        self.V0 = 20e-6

        self.Rv = 0.003*133.32*1e6 #0.005 * 133 / (10 ** -6)

        ## Shape parameters
        self.alpha = 1.672
        self.n1 = 1.32
        self.n2 = 21.9
        self.a1 = 0.708 * self.rt_Tpeak
        self.a2 = 1.677 * self.a1


        # Cycle management for restarting isovolumetric phase
        self.newCycle = False
        self.cycleNumber = 0
        self.num = 0
        self.atriumPressure = np.array([7.5 * 133.32])  # TODO: Fix this: Pressure in the atrium ## venouse pressure?!

        self.dQInOut = np.empty((2))

        self.dsetGroup = None

        self.pressure = np.zeros(0)
        self.volume = np.zeros(0)
        self.mitralQ = np.zeros(0)
        self.aorticQ = np.zeros(0)
        self.Elastance = np.zeros(0)

    def update(self, bcDict):
        super(VaryingElastanceSimpleDAE,self).update(bcDict)

        # BRX update variables
        self.T_BRX = self.T
        self.Emax_BRX = self.Emax
        self.Emin_BRX = self.Emin
        # Runtime Variables
        self.rt_T = self.T
        self.rt_Emax = self.Emax
        self.rt_Emin = self.Emin
        self.rt_Tpeak = self.Tpeak

        self.a1 = 0.708 * self.rt_Tpeak
        self.a2 = 1.677 * self.a1

    def initializeSolutionVectors(self, runtimeMemoryManager, solutionDataFile):
        """Initializes some solution vectors storing pressure, flow and volume
        of the ventricle, as well as opening and closing state

        NB! This method is not called from the class constructor, but is
        called externally by the initializeSolutionMatrices method in the
        solver, this is a bit messy, but was the easiest way to do it since the
        BC is initiated before the number of time steps is known.
        """

        self.dsetGroup = solutionDataFile.create_group('Heart')
        self.allocate(runtimeMemoryManager)

        """ Initial conditions in the ventricle"""
        self.Elastance[0] = self.E(0)
        self.atriumPressure[::] = 7.5*133.32
        self.pressure[0] = self.atriumPressure[0]
        self.volume[0]   = 145.*1e-6 # self.atriumPressure[0] / self.E(0) + self.V0

    def __call__(self, _domegaField_, duPrescribed, R, L, nmem, n, dt, P, Q, A, Z1, Z2):
        self.newCycle = False                     # Update the state of the mitral and aortic valve at timestep n + 1
        self.startNewCycleIfCriteriaIsMet(n, dt)
        self.Z1 = Z1
        self.Z2 = Z2
        self.funcPos0(_domegaField_, R, dt, P, Q, A, nmem, n)  # Compute the riemann variant going into the vessel save in omegaNew

        self.dQInOut = R[:][1] * self.omegaNew

        # Increment period time step
        self.num = self.num + 1

        # calculate du and return this!
        return np.dot(R, self.omegaNew), self.dQInOut

    def getCycleTime(self, n, dt):
        return self.num * dt

    def startNewCycleIfCriteriaIsMet(self, n, dt):
        if self.getCycleTime(n + 1, dt) > self.T:
            self.cycleNumber += 1
            self.num = 0
            self.newCycle = True
            self.rt_T = self.T_BRX
            self.rt_Tpeak = 0.39 # 0.4*self.rt_T
            self.rt_Emax = self.Emax_BRX
            self.rt_Emin = self.Emin_BRX

    def FlowImpedanceResidual(self, dP, dQ, P, Q, A, dt, nmem, n):
        # Assume pressure and derive flow as impedance
        dVImposed = -dt*(2*Q + dP/self.Z1)/2
        self.volume[nmem+1] = self.volume[nmem] + dVImposed
        self.pressure[nmem+1] = (self.Elastance[nmem+1])*(self.volume[nmem+1]-self.V0)
        residual = self.pressure[nmem+1] - (P+dP)
        return residual

    def PressureResidual(self, dP, dQ, P, Q, A, dt, nmem, n):
        # Assume Flow and derive pressure
        dVImposed = -dt*(2*Q + dQ)/2
        self.volume[nmem+1] = self.volume[nmem] + dVImposed
        self.pressure[nmem+1] = (self.Elastance[nmem+1])*(self.volume[nmem+1]-self.V0)
        residual = self.pressure[nmem+1] - (self.pressure[nmem]+dP)
        return residual

    def FlowResidual(self, dP, dQ, P, Q, A, dt, nmem, n):
        # Assume Pressure and derive Q
        dVImposed = -dt*(2*Q + dQ)/2
        self.pressure[nmem+1] = P + dP
        self.volume[nmem+1] = (P+dP)/self.Elastance[nmem+1] + self.V0
        residual = self.volume[nmem+1] - self.volume[nmem] - dVImposed
        return residual

    def EllipsoidalPressureResidual(self, dP, dQ, P, Q, A, dt, nmem, n):
        ## Ellipsoidal and Bernoulli's law
        # P,dP are P2
        dVImposed = -dt*(2*Q + dQ)/2
        Vh = self.volume[nmem] + dVImposed
        self.volume[nmem+1] =Vh
        Ph = self.Elastance[nmem+1]*(self.volume[nmem+1]-self.V0)
        self.pressure[nmem+1] = Ph
        a = 4.5e-2 # length of left ventricle in cm
        h = 0.15*a
        Ah = 3*Vh*(2*a-h)*h/(4*a**3)
        rho = 1050.0
        q2 =  (Q+dQ)**2
        residual = (P+dP - Ph)  +rho*q2/2*(1./Ah -1./A)
        return residual

    def ForwardFlow(self, _domegaField, R, dt, P, Q, A, nmem, n):
        # Explicit calculation of hte forward flow and pressure
        self.omegaNew[1] = _domegaField
        dVImposed = -dt*Q
        self.volume[nmem+1] = self.volume[nmem] + dVImposed
        self.pressure[nmem+1] = (self.Elastance[nmem+1])*(self.volume[nmem+1]-self.V0)
        dP = self.pressure[nmem+1]-self.pressure[nmem]
        dQ = dP/self.Z1
        self.omegaNew[0] = dP
        self.dQInOut = dQ # R[:][1] * self.omegaNew
        u =  np.dot(R[0][:], self.omegaNew[0])
        dQInOut = self.dQInOut
        return u, dQInOut

    def diastoleStep(self,dt,nmem):
        self.volume[nmem+1] = self.volume[nmem] + dt*(self.atriumPressure[nmem]-self.pressure[nmem])/self.Rv
        self.pressure[nmem+1] = self.Elastance[nmem+1]*(self.volume[nmem+1]-self.V0)

    def isovolumetricStep(self,dt,nmem):
        self.volume[nmem+1] = self.volume[nmem]
        self.pressure[nmem+1] = self.Elastance[nmem+1]*(self.volume[nmem+1]-self.V0)

    def residualPQ(self, dP, dQ, P, Q, A, dt, nmem, n):
        # Only Valid during systole
        fct = getattr(self,self.residualName)
        residual = fct(dP, dQ, P, Q, A, dt, nmem, n)
        return residual

    def reflection(self,_domegaField, R, dt, P, Q, n):
        domegaReflected_ = _domegaField
        self.omegaNew[0] = domegaReflected_
        self.omegaNew[1] = _domegaField
        self.dQInOut = R[:][1] * self.omegaNew

    def dampedReflection(self,_domegaField, R, dt, P, Q, nmem, n):
        """
        TODO: What is the physical interpretation of this? Is it nonconservative in volume? momentum?
        """
        domegaReflected_ = (-0.5 * Q - R[1][1] * _domegaField) / R[1][0]
        self.omegaNew[0] = domegaReflected_
        self.omegaNew[1] = _domegaField
        self.dQInOut = R[:][1] * self.omegaNew
        u =  np.dot(R, self.omegaNew)
        return u, self.dQInOut

    def funcPos0(self, _domegaField, R, dt, P, Q, A, nmem, n):
        """return function for position 0 at the start
        of the vessel
        """
        self.aorticQ[nmem] = Q
        ventrPn = self.pressure[nmem]
        venoP = self.atriumPressure[nmem]
        t = self.getCycleTime(nmem + 1, dt)
        self.Elastance[nmem + 1] = self.E(t)

        # TODO: fix how this inherits from generalizedPQ_BC
        if (Q >= -1e-15 and ventrPn - P > (-0.0001)):  # systolecondition
            # TODO: Ignore return values as they are stored as class members
            if self.residualName in ["ForwardFlow"]: #Explicit Forms
                fct = getattr(self,self.residualName)
                u, dQInOut = fct(_domegaField, R, dt, P, Q, A, nmem, n)
            else: #Implicit Form
                u, dQInOut = super(VaryingElastanceSimpleDAE,self).funcPos0(_domegaField, R, dt, P, Q, A, nmem, n)

        else:
            u, dQInOut = self.dampedReflection(_domegaField,R,dt,P,Q,nmem,n)
            # TODO: Anyway to factor these functions so they are only active when the heart is "coupled" to the arteries?
            if (venoP > ventrPn):
                self.diastoleStep(dt, nmem)
            else:
                self.isovolumetricStep(dt,nmem)

        return u, dQInOut

    def E(self, t):
        """Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
           according to Tpeak, T, Emax and Emin. """
        a1 = self.a1
        a2 = self.a2
        n1 = self.n1
        n2 = self.n2
        shapeFunction1 = (t / (a1)) ** n1 / (1 + (t / (a1)) ** n1)
        shapeFunction2 = (1 + (t / (a2)) ** n2) ** (-1)
        return (self.rt_Emax - self.rt_Emin) * self.alpha * shapeFunction1 * shapeFunction2 + self.rt_Emin

    def dEdt(self, t):
        """
        Derivative of Elastance Function
        """
        Emax = self.Emax
        Emin = self.Emin
        a1 = self.a1
        a2 = self.a2
        n1 = self.n1
        n2 = self.n2
        alpha = self.alpha


        cg = -(Emax - Emin) * alpha * ((t / a1) ** (2 * n1) * (t / a2) ** n2 \
                * n2 - (t / a2) ** n2 * (t / a1) ** n1 * n1 + (t / a1) ** n1 \
                * (t / a2) ** n2 * n2 - (t / a1) ** n1 * n1) \
                / (1 + (t / a1) ** n1) ** 2 / (1 + (t / a2) ** n2) ** 2 \
                / t
        return cg
