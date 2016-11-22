import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 

import moduleFilePathHandlerVPC as mFPH_VPC

import UtilityLib.progressBar as cPB

import classLocationOfInterestManager
import classUqsaMethods
import classSampleManager

import shutil

import time
import numpy as np
import multiprocessing

class UqsaCase(TestBaseClass):
    
    externVariables      = { 'createSample'             : TestBaseClass.ExtValue(bool),
                             'createEvaluationXmlFiles' : TestBaseClass.ExtValue(bool),  
                             'simulateEvaluations'      : TestBaseClass.ExtValue(bool),  
                             'preProcessData'           : TestBaseClass.ExtValue(bool),  
                             'postProcessing'           : TestBaseClass.ExtValue(bool),  
                             'localEvaluation'          : TestBaseClass.ExtValue(bool),             
                             'multiprocessing'          : TestBaseClass.ExtValue(bool),  
                             'numberOfProcessors'       : TestBaseClass.ExtValue(int),
                             'simulateEvaluationNumbers': TestBaseClass.ExtValue(int, multiVar=True),
                             'sampleManager'            : TestBaseClass.ExtObject({'sampleManager':classSampleManager.SampleManager}),
                             'uqsaMethods'              : TestBaseClass.ExtDict('uqsaMethod', TestBaseClass.ExtObject({'uqsaMethodPolynomialChaos':classUqsaMethods.UqsaMethodPolynomialChaos,
                                                                                                                       'uqsaMethodMonteCarlo'     :classUqsaMethods.UqsaMethodMonteCarlo,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirLR': classUqsaMethods.UqsaMethodPolynomialChaosDepDirLR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirQR': classUqsaMethods.UqsaMethodPolynomialChaosDepDirQR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirQL': classUqsaMethods.UqsaMethodPolynomialChaosDepDirQL,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirLRorder': classUqsaMethods.UqsaMethodPolynomialChaosDepDirLRorder,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFLR' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFLR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFR' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFL' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFL,
                                                                                                                       'uqsaMethodMonteCarloParametrizedBootstrapping' : classUqsaMethods.UqsaMethodMonteCarloParametrizedBootstrapping},
                                                                                                                        )),
                             'locationOfInterestManager' : TestBaseClass.ExtObject({'LocationOfInterestManager':classLocationOfInterestManager.LocationOfInterestManager}),
                           } 
    
    externXmlAttributes  = []
    
    externXmlElements    = ['createSample',
                            'createEvaluationXmlFiles',
                            'simulateEvaluations',
                            'preProcessData' ,
                            'postProcessing',
                            'localEvaluation',                 
                            'multiprocessing',     
                            'numberOfProcessors',   
                            'simulateEvaluationNumbers',
                            'sampleManager',
                            'uqsaMethods',
                            'locationOfInterestManager']
    
    
    def __init__(self):
        
        self.networkName = None
        self.dataNumber  = None
        
        ### data read in from file
        ##control variables
        # create samples ( TRUE == create and save, FALSE == load existing)
        self.createSample     = True
        # ceate simulation case files       
        self.createEvaluationXmlFiles = True
        #  run simulations 
        self.simulateEvaluations    = True
        self.localEvaluation        = True #TODO: add functions for server
        self.multiprocessing        = True
        self.numberOfProcessors     = 12
        self.simulateEvaluationNumbers = []
        # pre process data for all quantitiy of interest
        self.preProcessData   = True
        self.preProcessData   = True
        # post processing - uncertainty quantification and sensitivity analysis
        self.postProcessing  = True
        
        # sample manager for the case
        self.sampleManager = None
        
        ## uqsa method instance of this case to use
        self.uqsaMethods = None # {}
        
        ## location of interests i.e. datastructure including all evaluated data y = f(z)
        self.locationOfInterestManager = None
        
        ### data assoziated during run time
        ## samples of Z
        
        
    def initialize(self,networkName, dataNumber):
        '''
        Initialize case class
        '''
        self.networkName = networkName
        self.dataNumber  = dataNumber
        
        self.locationOfInterestManager.initialize()
        
        for uqsaMethodName in self.uqsaMethods.iterkeys():
            print "Info classUQSACase 96: running ", uqsaMethodName
        
    def aquireSamples(self, distributionManager, randomInputsExtDist):
        '''
        Function that envokes either sample creation of loading depending on the defined control variable        
        '''
        randomVariableNames = []
        for randomInput in randomInputsExtDist:
            randomVariableNames.append(randomInput.name) 
        
        maxSampleSize = 0
        abcSample = False
        #find out maximum numbers of samples needed and if ABC sample is needed
        for uqsaMethod in self.uqsaMethods.itervalues():
            sampleSizeCurrent,abcSampleCurrent = uqsaMethod.evaluateSamplesSize(distributionManager.distributionDimension)   
            
            if sampleSizeCurrent > maxSampleSize:
                maxSampleSize = sampleSizeCurrent
            if abcSampleCurrent == True:
                abcSample = True
                
        if self.createSample == True:
            sampleSaveFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, 
                                                  mode = "write", caseName = self.sampleManager.samplingMethod)
            self.sampleManager.createSamples(distributionManager,randomVariableNames, 
                                             maxSampleSize, abcSample, sampleSaveFile)
            
        else:
            sampleLoadFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, 
                                                  mode = "read", caseName = self.sampleManager.samplingMethod)
            self.sampleManager.evaluateLoadedSamples(sampleLoadFile, maxSampleSize, randomVariableNames)
                            
        self.locationOfInterestManager.sampleSize = self.sampleManager.currentSampleSize
    
        
    def createEvaluationCaseFiles(self): 
        '''
        
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
        '''
        self.evaluationCaseFiles = [] # list of dict:  [ caseFileDict1,caseFileDict2 ..]  for each evaluation
        
        #TODO: replace create evaluation files or save them to disc!!!
        if self.simulateEvaluations == True or self.preProcessData == True or self.createEvaluationXmlFiles == True:
        
            sampleSize = self.sampleManager.currentSampleSize
            
            print "Create evaluation case file list"
            progressBar = cPB.ProgressBar(35, sampleSize)
            
            for simulationIndex in xrange(sampleSize):
                                        
                networkXmlFileLoad = mFPH_VPC.getFilePath('uqsaEvaluationNetworkXmlFile', self.networkName, 
                                                           self.dataNumber, 'write',
                                                           caseName = self.sampleManager.samplingMethod, 
                                                           evaluationNumber=simulationIndex)
                networkXmlFileSave = networkXmlFileLoad
                pathSolutionDataFilename = mFPH_VPC.getFilePath('uqsaEvaluationSolutionDataFile', 
                                                                self.networkName, self.dataNumber, 'write',
                                                                 caseName = self.sampleManager.samplingMethod, 
                                                                 evaluationNumber=simulationIndex)
                
                caseFileDict1= {'simulationIndex': simulationIndex,
                                'networkName': self.networkName,
                                'dataNumber': self.dataNumber,
                                'networkXmlFileLoad': networkXmlFileLoad,
                                'networkXmlFileSave': networkXmlFileSave,
                                'pathSolutionDataFilename': pathSolutionDataFilename}
                
                self.evaluationCaseFiles.append(caseFileDict1)
            
                progressBar.progress(simulationIndex) 
            
        
    def getSimulatioNBatchFileList(self):
        '''
        Returns simulation batch file list, as defined in configs
        '''
        
        startIndex = 0
        endIndex   = int(self.sampleManager.currentSampleSize)
        
        newRange = self.simulateEvaluationNumbers
        if len(newRange) == 2:
            # check if the indices are avaliable
            if all([i in xrange(self.sampleManager.currentSampleSize) for i in newRange]):
                if newRange[0] < newRange[1]:
                    startIndex = newRange[0]
                    endIndex   = newRange[1]
                    
        batchFileList = self.evaluationCaseFiles[startIndex:endIndex+1]
    
        return batchFileList
        
    def preprocessSolutionData(self):
        '''
        envoke preprocessing of solution data
        '''
        if self.preProcessData == True:
            #TODO rename caseName!!
            caseName = '_'.join([self.sampleManager.samplingMethod])
            preprocessedSolutionData = mFPH_VPC.getFilePath('preprocessedDataFile', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName = caseName )
            
            simulationTimeFileSave  = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName =  caseName)
            simulationTimeFileLoad = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "read", caseName = caseName, exception = 'No')
            
            self.locationOfInterestManager.preprocessSolutionData(self.evaluationCaseFiles,
                                                                  preprocessedSolutionData,
                                                                  simulationTimeFileSave,
                                                                  simulationTimeFileLoad)
                    
    def quantifyUncertaintyAndAnalyseSensitivtiy(self, distributionManager):
        '''
        evnoke uq sa process
        '''
        if self.postProcessing == True:
            
            # copy of preprocessed data file 
            caseName = '_'.join([self.sampleManager.samplingMethod])
            preprocessedSolutionData = mFPH_VPC.getFilePath('preprocessedDataFile', self.networkName, self.dataNumber, 
                                                     mode = "read", caseName = caseName)
            uqsaSolutionDataFile = mFPH_VPC.getFilePath('uqsaSolutionDataFile', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName = caseName)
            shutil.copy(preprocessedSolutionData,uqsaSolutionDataFile)
            
            # open solution file
            self.locationOfInterestManager.openQuantityOfInterestFile(uqsaSolutionDataFile, mode = 'r+')
            
            # loop through data objects
            for qoi in self.locationOfInterestManager.getQoiIterator():
                
                timeStartBatch = time.time()
                
                ## ranges 
                
                ### Normal distribution
                ## last wave start to first discontinuity
                ## 0.2572545014772901 - 0.2709600872555526
                #basis = np.linspace(0.258,0.269,100)
                
                ## first discontinuity to first wave return
                ## 0.2709600872555526 - 0.29017871810043244
                #basis = np.linspace(0.271,0.29,100)
                
                ### Uniform distribution
                
                ## last wave start to first discontinuity
                #basis = np.linspace(0.258,0.277,100)
                
                ## first discontinuity to last discontinuity
                #basis = np.linspace(0.27890,0.28824,100)
                
                #basis = np.array([0.280])
                
                # all
                #basis = np.linspace(0.025,0.045,500)
                
                # just discontinuity
                
                #basis = np.linspace(0.031,0.039,50)
                
                numerOfSamples = 20
                
                # called now
                #lowerEnd = 0.031 
                
                # called now2
                #lowerEnd = 0.0318
                
                #upperEnd = 0.035
                
                #called full
                lowerEnd = 0.032
                upperEnd = 0.038
                
                import chaospy as cp
                u = cp.Uniform()
                basis = lowerEnd+(upperEnd-lowerEnd)*np.sort(u.sample(numerOfSamples,'H'))
                                
                #basis = np.linspace(0.03170,0.035,22)
                
                #basis = np.linspace(0.034,0.035,25)
                
                
                # half
                
                #basis = np.linspace(0.025,0.035,100)
                
                # no sicontinuity
                #basis = np.linspace(0.025,0.03,25)
                
                
                print "hashDataForGivenBases {}".format(basis)
                #qoi.hashDataForGivenBases(basis, self.sampleManager.currentSampleSize)
                
                timeBatchJob= time.time()-timeStartBatch
                minutesBatch = int(timeBatchJob/60.)
                secsBatch = timeBatchJob-minutesBatch*60.
                print '====================================='
                print 'total runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
                print '====================================='
                print
                
                multiprocessingUQSA = False
                
                if multiprocessingUQSA == True:
                    self.multiprocessingUQSA(qoi, distributionManager)
                
                else:
                    
                    timeStartTotal = time.time()
                    
                    for uqsaMethodName,uqsaMethod in self.uqsaMethods.iteritems():
                        
                        timeStartBatch = time.time()
                        print "calculate uqsa measure for {}".format(uqsaMethodName)
                        
                        uqsaMeasures = uqsaMethod.calculateStatistics(distributionManager, self.sampleManager, qoi)
                        qoi.addUqsaMeasures(uqsaMethodName, uqsaMeasures)
                        
                        timeBatchJob= time.time()-timeStartBatch
                        minutesBatch = int(timeBatchJob/60.)
                        secsBatch = timeBatchJob-minutesBatch*60.
                        print '====================================='
                        print 'runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
                        print '====================================='
                        print
            
                    timeTotal= time.time()-timeStartTotal
                    minutesTotal = int(timeTotal/60.)
                    secsTotal = timeTotal-minutesTotal*60.
                    print '====================================='
                    print 'total runtime:  {} min {} sec'.format(minutesTotal,secsTotal)
                    print '====================================='
                    print
            
            self.locationOfInterestManager.closeAndSaveQuantityOfInterestFile()
    
    
    
    def multiprocessingUQSA(self,qoi, distributionManager):
        '''
        Run all uqsa methods in as a local multiprocess
        '''
        print "running Multiprocessing UQSA"
        timeStartBatch = time.time()
        
        # create batch list for all jobs    
        batchList = []
        for uqsaMethodName,uqsaMethod in self.uqsaMethods.iteritems():
            batchList.append([uqsaMethodName,uqsaMethod,distributionManager,self.sampleManager,qoi])         
        # run jobs
        pool = multiprocessing.Pool( multiprocessing.cpu_count())
        pool.imap(self.batchJobUQSA,batchList)
        pool.close() 
        pool.join()
        
        timeBatchJob= time.time()-timeStartBatch
        minutesBatch = int(timeBatchJob/60.)
        secsBatch = timeBatchJob-minutesBatch*60.
        print '====================================='
        print 'total runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
        print '====================================='
        print
        
        
        
    
    def batchJobUQSA(self, args):
        '''
        batch job for local uqsa multiprocessing
        '''
        
        uqsaMethodName,uqsaMethod,distributionManager,sampleManager,qoi = args
        
        uqsaMeasures = uqsaMethod.calculateStatistics(distributionManager, sampleManager, qoi)
        qoi.addUqsaMeasures(uqsaMethodName, uqsaMeasures)
            
        
        
    
    