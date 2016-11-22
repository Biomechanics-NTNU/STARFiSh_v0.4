
from testBaseClass import TestBaseClass

import moduleFilePathHandlerVPC as mFPH_VPC

import h5py
import numpy as np
import chaospy as cp

class SampleManager(TestBaseClass):
    
    externVariables      = { 'dependentCase' : TestBaseClass.ExtValue(bool),
                             'samplingMethod'  : TestBaseClass.ExtValue(str, strCases = ['K','R','L','S','H','M','C','NC','G','RG']),
                             } 
    
    externXmlAttributes  = []
    
    externXmlElements    = ['samplingMethod',
                            'dependentCase']
    
    
    def __init__(self):
                
        self.samplingMethod = None
        self.dependentCase = False
        
        
        self.randomVariableNames = None #randomVariableNames
        self.abcSample = False
        
        self.samples = None
        self.samplesDependent = None
        self.currentSampleSize = None #currentMaxSample
        
        
    def sampleSize(self):
        '''
        returns the max number of samples used for this case
        '''
        return self.currentMaxSample
        
    def getSample(self, sampleIndex):
        '''
        Returns sample for a certain sampleIndex
        
        dependenSample or independenSample is chooses dependent on the case definitions
        '''
        # TODO: if random sample, ensure that sub set samples of sample space are randomly choosen
        if sampleIndex in xrange(self.currentSampleSize):
            if self.dependentCase == False:
                sample = self.samples[sampleIndex]
            else:
                sample = self.samplesDependent[sampleIndex]
        else:
            raise ValueError('sampleIndex {} out of range 0:{}'.format(sampleIndex,self.currentSampleSize))
                             
        return sample
    
    def getSampleMatrices(self, sampleSize, abcSample, offset=0):
        '''
        Function which returns the samples, and abcSamples as defined by the input
        '''
        #TODO: check if samplesSize and offsett+sampleSize is within range of samples
        if abcSample == False:
            dependentSamples = False
            if self.dependentCase == True:
                dependentSamples = self.samplesDependent[offset:sampleSize+offset]
            return self.samples[offset:sampleSize+offset], dependentSamples
        else:
            #TODO: implement abc sampling hash if monte carlo is used
            raise NotImplementedError("MC sensitivity, abc-sample hash case not implemented yet, exit")
            
    
    def createSamples(self, distributionManager, randomVariableNames, maxSampleSize, abcSample, sampleSaveFile):
        '''
        create samples for the defined distribution for given samplesSize and sampleMethod
        using the chaospy toolbox
        Input
        
            sampleMethod : str
                (from chaospy)
                Alternative sampling techniques
            
                Normal sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "K"     Korobov             no
                "R"     (Pseudo-)Random     no
                "L"     Latin hypercube     no
                "S"     Sobol               yes
                "H"     Halton              yes
                "M"     Hammersley          yes / (nested only in the
            
                Grided sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "C"     Chebyshev nodes     maybe
                "G"     Gaussian quadrature no
                "E"     Gauss-Legendre      no
            
            
        '''
        if abcSample == False:
            # creating samples
            samples = distributionManager.jointDistribution.sample(maxSampleSize,self.samplingMethod).transpose()   
            # reshape samples
            if distributionManager.distributionDimension == 1:
                samples = samples.reshape(maxSampleSize,1)
            # create dependent samples if correlation exists
            samplesDependent = None
            if self.dependentCase == True:
                samplesDependent = self.createDependentSamples(distributionManager, samples)
        else:
            #TODO: implement abc sampling hash if monte carlo is used
            raise NotImplementedError("MC sensitivity, abc-sample hash case not implemented yet, exit")
        
        self.samples = samples
        self.samplesDependent = samplesDependent
        self.randomVariableNames = randomVariableNames
        self.currentSampleSize = maxSampleSize
        self.abcSample = abcSample
        
        self.saveSamples(sampleSaveFile)
  
    def createDependentSamples(self,distributionManager, samples):
        '''
        Create dependen samples if dependen case = True and dependent random variables true
        '''
        if distributionManager.jointDistributionDependent != None:
            samplesDependent = distributionManager.jointDistributionDependent.inv(distributionManager.jointDistribution.fwd(samples.T)).transpose()
        else:
            raise ValueError("uqsaMethodPolynomialChaos.createDependentSamples(), cannot create dependent samples as distributionManager.jointDistributionDependent is not defined!")
        return samplesDependent
    
    def createSamplesMonteCarlo(self, distributionManager):
        '''
        create samples
        
        samplesA, samplesB and samplesC
        '''
        
        distDim = distributionManager.distributionDimension
        
        if distDim == 1:
            print "WARNING: random distribution dimension == 1, no sensitivity analysis possible!"
            self.sensitivityAnalysis = False
        
        ## if sensitivityAnalysis is false, only UQ is done with one matrix
        if self.sensitivityAnalysis == False:
            samples = distributionManager.jointDistribution.sample(self.sampleSize,self.samplingMethod).transpose()   
            # reshape samples
            if distDim == 1:
                samples = samples.reshape(self.sampleSize,1)
                
            return samples,None
        else:
            # createHash table
            self.createSampleMatixHashTable(distDim)
            
            samples = distributionManager.jointDistribution.sample(self.sampleSize*2,self.samplingMethod).transpose()
            
            samplesA = samples[self.matrixHash['A'][0]:self.matrixHash['A'][1]]
            samplesB = samples[self.matrixHash['B'][0]:self.matrixHash['B'][1]]
                
            samplesTest = np.sum((samplesA-samplesB).ravel())
            if samplesTest == 0:
                print "WARNING: samplesA and samplesB are the same!"
        
            for i in xrange(distDim):
                samplesC = samplesB.copy()
                samplesC[:,i] = samplesA[:,i].copy()
                
                samples = np.vstack([samples,samplesC])
            
            return samples,None
    
    def createSampleMatixHashTable(self, distDim):
        '''
        Creates hash table for the sample matrix
        '''
        self.matrixHash = {}
        self.matrixHash['A'] = [0,self.sampleSize]
        self.matrixHash['B'] = [self.sampleSize,2*self.sampleSize]
        for d in xrange(distDim):
            self.matrixHash[''.join(['C',str(d)])] = [self.sampleSize*(d),self.sampleSize*(d+1),self.sampleSize*(d+2),self.sampleSize*(d+3)]
    
    def evaluateLoadedSamples(self, sampleFile, maxSampleSize, randomVariableNames):
        '''
        load data and evaluate match with the defined case
        '''
        # load data
        loadedData = self.loadSamples(sampleFile)
        # TODO: type conversion should be moved to load function
        loadedData['randomVariableNames'] = loadedData['randomVariableNames'].tolist()
        loadedData['abcSample'] = bool(loadedData['abcSample'])
        
        # check if data is correct
        if loadedData['randomVariableNames'] != randomVariableNames:
            raise ValueError("Loaded randomVariableNames {} for samples does not match with defined random variable vector {}".format(loadedData['randomVariableNames'], randomVariableNames))
                     
        if len(loadedData['randomVariableNames']) != len(loadedData['samples'][0]):
            raise ValueError("Created sample matrix does not match with defined random variable vector {}".format(self.randomVariableNames))
                         
        if loadedData['samplingMethod'] != self.samplingMethod:
            raise ValueError("Loaded {}: {}  does not match with defined variable {}".format('samplingMethod', loadedData['samplingMethod'],self.samplingMethod))
          
        if loadedData['abcSample'] != self.abcSample:
            raise ValueError("Loaded {}: {}  does not match with defined variable {}".format('abcSample', loadedData['abcSample'],self.abcSample))
           
           
        # TODO: if random sample, ensure that sub set samples of sample space are randomly choosen          
        # crop sampling solution data
        if maxSampleSize < len(loadedData['samples']):
            loadedData['samples'] = loadedData['samples'][:maxSampleSize]
            # TODO: dependent Samples
        elif maxSampleSize > len(loadedData['samples']):
            raise ValueError("Loaded sample size {} is to small for max sample size {} in defined uqsaCases, rerun sample creation! and simulate extra cases!".format(len(loadedData['samples']), maxSampleSize) )
        
        # update loaded data
        self.setVariablesDict(loadedData)
        self.currentSampleSize = maxSampleSize
        
        
    ##---- read write data
    def loadSamples(self, sampleFile):
        '''
        load the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''
        ## TODO: rewrite in more base-class manner
        # TODO: datatypeconversion!!
        # TODO: optional bool for object types
        
        hdf5SaveGroup = h5py.File(sampleFile,'r')
        
        loadedData = {}
        
        ## pure variables
        # in class definition
        variablesToLoad = ["samples",
                           'samplesDependent']
        # functionality
        for variableName in variablesToLoad:
            if variableName in hdf5SaveGroup.keys(): 
                loadedData[variableName] = hdf5SaveGroup[variableName][:]
            #else:
            #    raise ValueError("Could not load variable {} from file {}".format(variableName, hdf5SaveGroup))
        ## attributes
        # in class definition
        attributesToLoad = ['randomVariableNames',
                            'dependentCase',
                            'samplingMethod',
                            'abcSample']
        # functionality
        for attributeName in attributesToLoad:
            if attributeName in hdf5SaveGroup.attrs.keys(): 
                loadedData[attributeName] = hdf5SaveGroup.attrs.get(attributeName)
            else:
                raise ValueError("Could not load variable {} from file {}".format(attributeName, hdf5SaveGroup))
            
        hdf5SaveGroup.close()
        
        return loadedData
      
    def saveSamples(self, sampleFileName):
        '''
        save the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''        
        hdf5SaveGroup = h5py.File(sampleFileName,'w')
        
        ## pure variables
        # in class definition
        variablesToSave = ["samples",
                           'samplesDependent']
        # functionality
        for variableName in variablesToSave:
            variableValue = self.getVariable(variableName)
            if variableValue != None: 
                hdf5SaveGroup.create_dataset(variableName, data=variableValue)
        
        ## attributes
        # in class definition
        attributesToSave = ['randomVariableNames',
                            'dependentCase',
                            'samplingMethod',
                            'abcSample']
        # functionality
        for attributeName in attributesToSave:
            attributeValue = self.getVariable(attributeName)
            if attributeValue != None: 
                hdf5SaveGroup.attrs.create(attributeName, data = attributeValue)
        
        hdf5SaveGroup.flush()
        hdf5SaveGroup.close()  