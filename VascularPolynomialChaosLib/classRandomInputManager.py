import numpy as np

import itertools

from testBaseClass import TestBaseClass 

import classRandomInput 
import classCorrelationMatrices

class RandomInputManager(TestBaseClass):
    
    externVariables      = {'randomInputs' : TestBaseClass.ExtDict('randomInput',TestBaseClass.ExtObject({'ParametricRandomInput':classRandomInput.ParametricRandomInput,
                                                                                                          'GeneralRandomInput':classRandomInput.GeneralRandomInput})),
                            'correlation'  : TestBaseClass.ExtObject({'CorrelationMatrix':classCorrelationMatrices.CorrelationMatrix}, optional = True)
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['randomInputs',
                            'correlation']
    
    def __init__(self):        
        '''
        The RandomInputManager class is as in the network xml file
        the container of correlation matrix and random inputs.
        In addition it has methods to sorts out the connection between random inputs
        and connects the random inputs via the update functions to the vascular network variables
        or other random inputs. (see method: linkRandomInputUpdateFunctions)
        '''        
        self.randomInputs = {} # randomInput as they stand in xml
        
        self.correlation       = None
        self.correlationMatrix = None
        
        self.randomInputsExtDist = [] # randomInput which have a external distribution assosiated
        self.randomInputDimension = 0
                                                
    def initialize(self, vascularNetwork):
        '''
        Method to initialize the random inputs and
        defined them into the different types
        
        Args:
            vascularNetwork (instance of VascularNetwork): 
                An instance of the vascular network class to link the 
                random inputs to the corresponding update functions
        '''
        
        randomInputParameters = []
        
        for name,randomInput in self.randomInputs.iteritems():
            # set name            
            randomInput.name = name
            if randomInput.randomInputType == 'parametricRandomInput':
                if randomInput.parameter not in randomInputParameters:
                    randomInputParameters.append(randomInput.parameter)
                else:
                    raise ValueError("assoziated parameter <{}> of randomInput {} is doubled defined!".format(randomInput.location,randomInput.name))
                
        self.linkRandomInputUpdateFunctions(vascularNetwork)
            
        self.checkCorrelationMatrix()
            
    def checkCorrelationMatrix(self):
        '''
        This functions initiates the assembly of the correlation matrix and 
        ensures that the defined correlation matrix is consistent and matches with
        the defined random inputs.
        '''
        # get extern random input name list
        definedBases = []
        for randomInput in self.randomInputsExtDist:
            definedBases.append(randomInput.name) 
        
        self.correlationMatrix = self.correlation.assembleCorrelationMatrix(definedBases)
            
    def linkRandomInputUpdateFunctions(self, vascularNetwork):
        '''
        This function establisehes the links of the random inputs.
        
        ParametricRandomInputs can be in the form of:
        
            a+b*Uniform(0,1) / a+b*Normal(0,1)
            or
            a+b*Z1 where Z1 must be a separate defined general random input 
        
        In addition this function creates a randomInputvector
        which holds all random inputs, with a propability denistiy function 
        which is modeled in stochastic sense.
        
        eg. 4 defined random inputs:
        Z1 (parametric input): target(vessel_1_betaHayashi), 2 + 0.5*Z4
        Z2 (parametric input): target(vessel_2_betaHayashi), 6 + 0.5*Z4
        Z3 (parametric input): target(vessel_1_radiusA), 2 + 0.5* Uniform(0,1)
        Z4 (general random input): 4 + 5.5* Normal(0,1)
        
        The random variables have the connection:
        
        vessel_1_betaHayashi    <-- Z1 <-- Z4
        vessel_2_betaHayashi    <-- Z2 <-- Z4
        vessel_1_radiusA        <-- Z3
        
        the randomInputvector would now be [Z3, Z4], as Z1 and Z2 are dependent on
        Z4 and have no independet distibution which is stochstically modelled.
        In this way many parameters can be controlled by 1 random variable, e.g. age.
        
        Args:
            vascularNetwork (instance of VascularNetwork): 
                An instance of the vascular network class to link the 
                random inputs to the corresponding update functions
        '''
        
        randomInputMap = {}
        for randomInputName,randomInput in self.randomInputs.iteritems():
            dist = randomInput.distributionType 
            
            if randomInput.randomInputType == 'parametricRandomInput':
                # check distribution
                loc = randomInput.parameter.split('_')
                objType = loc[0]
                
                if randomInput.variableName == []:
                    randomInput.variableName.append(loc[-1])
                
                if objType == "boundaryCondition":
                    for bc in vascularNetwork.boundaryConditions[int(loc[2])]:
                        if bc.getVariableValue('name') == loc[1]:
                            randomInput.updateMethods = {randomInput.variableName[0]:
                                                         bc.update}
            
                elif objType == "vessel":
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 vascularNetwork.vessels[int(loc[1])].update}
                
                elif objType == "example_location": # Example for how to add additional types of locations
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 getattr(vascularNetwork,objType)[int(loc[1])].update}
                
                elif objType == "measurementRoutine":
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 vascularNetwork.measurementRoutine.setVariablesDict}
                
                else: raise ValueError("RandomInputManager: Parameter {} of randomInput {} is linkable.".format(randomInput.location,randomInputName))
                if randomInput.updateMethods == {}: break
                
            # find links between random inputs
            if dist in ['Normal','Uniform']:
                # append random input to random input vector
                self.randomInputsExtDist.append(randomInput)
            else:
                if dist not in randomInputMap:
                    randomInputMap[dist] = {randomInputName : randomInput.passRealisationToAssosiatedObj}
                else:
                    randomInputMap[dist][randomInputName] = randomInput.passRealisationToAssosiatedObj
        
        # link general random inputs 
        for dist,updateMethods in randomInputMap.iteritems():
            # find random input with dist:
            for randomInputName in self.randomInputs.iterkeys():
                if dist == randomInputName:
                    generalRandomInput = self.randomInputs[randomInputName]
                    generalRandomInput.updateMethods = updateMethods
                    generalRandomInput.variableName = updateMethods.keys()
                    #self.randomInputsExtDist.append(generalRandomInput)
        
        self.randomInputDimension = len(self.randomInputsExtDist)
                    
    def saveRealisationLog(self, evaluationLogFile, networkName, dataNumber, caseName): 
        '''
        This Function saves a log file with the realisations passed for each sample iteration
        
        Args:
            evaluationLogFile (file): the file in which the log should be saved
            networkName (str) : name of the current network
            dataNumber (str) : data number of the current simulation case
            caseName (str) : case name of the current stochastic case
        '''
        
        logData = np.array([randomInput.updateLog for randomInput in self.randomInputs.values()]).transpose()
        
        logfile = open(evaluationLogFile, "wb")
        logfile.write(''.join(['Stochastic simulation ',str(networkName),' DataNumber ',dataNumber,'\n','\n']))
        logfile.write(''.join(['uqsaCase :', caseName,'  number of evaluations: ',str(len(logData)),'\n','\n']))
        
        for info in self.generateInfo():
            logfile.write(''.join([info,'\n']))
        randomInputName = [randomInput.name for randomInput in self.randomInputs.values()]
        logfile.write(''.join(['\n{:<10}'.format(''), '{:20}'.format('RandomVariableId'),'\n','\n']))
        logfile.write(''.join(['{:<10}'.format('EvalNr'), ''.join(['{:20}'.format(j) for j in ['{:<5}'.format(k) for k in randomInputName]]),'\n',]))     
        for i,logLine in enumerate(logData):
            logfile.write(''.join(['{:<10}'.format(i), ''.join(['{:20}'.format(j) for j in ['{:<.5}'.format(k) for k in logLine]]),'\n']))
        
        
    def printOutInfo(self):
        '''
        Function to print out random variable informations
        '''
        for info in self.generateInfo():
            print info
            
    def generateInfo(self):
        '''
        Function which generates info of the random inputs
        
        Returns:
            randomInputManagerInfo (list): combined info of all random inputs as str
        '''
        randomInputManagerInfo = []
        randomInputManagerInfo.append("\n Defined Random Inputs\n")
        randomInputManagerInfo.append( '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution"))
        randomInputManagerInfo.append( "-------------------------------------------------------------------- \n")
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputs.values()]))
        for info in randomInputInfos:
            randomInputManagerInfo.append( info)
        randomInputManagerInfo.append( "\n Defined Random Variables\n")
        randomInputManagerInfo.append( '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution"))
        randomInputManagerInfo.append( "-------------------------------------------------------------------- \n")
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputsExtDist]))
        for info in randomInputInfos:
            randomInputManagerInfo.append(info)
        return randomInputManagerInfo
    
    def deleteAllRandomInputs(self):
        '''
        Function removes all existing random inputs
        '''
        for rI in self.randomInputsList:
            del rI
        self.randomInputsList = []
        self.randomInputsExtDist = []
        self.map = {}
        self.randomInputDimension = 0
                
