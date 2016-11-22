


import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))


from classLocationOfInterest import LocationOfInterest

#sys.path.append(''.join([cur,'/../UtilityLib']))
#import moduleXML as mXML

from UtilityLib import moduleXML
import UtilityLib.progressBar as cPB

import numpy as np
import h5py

from testBaseClass import TestBaseClass 


class LocationOfInterestManager(TestBaseClass):
    '''
    
    '''
    externVariables      = {'locationsOfInterest'   : TestBaseClass.ExtDict('locationOfInterest',TestBaseClass.ExtObject({'LocationOfInterest':LocationOfInterest})),
                            'evaluateSimulationTime': TestBaseClass.ExtValue(bool)
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['evaluateSimulationTime',
                            'locationsOfInterest']
    
    # in class definition
    ## pure variables
    variablesHdf5Memory = ["simulationTime"]
    ## dictionary with objects to load
    objectDictsHdf5Memory = ["locationsOfInterest"]
    
    def __init__(self):
        
        self.locationsOfInterest = {}
        self.sampleSize = None
            
        self.evaluateSimulationTime = False
        self.simulationTime = None
            
    def initialize(self):
        '''
        Function to initialize all locations of interest which creates quantity of intertes objects
        '''
        for locationOfInterest in self.locationsOfInterest.itervalues():
            locationOfInterest.initialize()
                
    def addLocationOfInterest(self,locationId, locationName, quantitiesOfInterestToProcess, xVal, confidenceAlpha):
        '''
        
        '''
        self.locationsOfInterest[locationId] = LocationOfInterest(locationName,quantitiesOfInterestToProcess, xVal, confidenceAlpha)
                                
    def preprocessSolutionData(self,
                               evaluationCaseFiles,
                               preprocessedSolutionData,
                               simulationTimeFileSave,
                               simulationTimeFileLoad):
        '''
        load all samples and pass data to locations of interest
        
        find simulation time array with largest dt
        
        invoke time cropping and post processing for the data        
        '''
        # loop through all solution data files 
        # find time array
        #solutionTime min
        nPoints = 1e32
        timeS = 0
        timeE = 1e32
        
        if self.evaluateSimulationTime == True:
            
            print "estimate simulation time of all simulations"
            progressBar = cPB.ProgressBar(35, len(evaluationCaseFiles))
            
            for batchData in evaluationCaseFiles:
                networkName              = batchData['networkName']
                dataNumber               = batchData['dataNumber']
                networkXmlFileLoad       = batchData['networkXmlFileLoad']
                pathSolutionDataFilename = batchData['pathSolutionDataFilename']
                
                vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
                vascularNetwork.linkSolutionData()
                
                cPoints = len(vascularNetwork.simulationTime)
                ctimeS  = min(vascularNetwork.simulationTime)
                ctimeE  = max(vascularNetwork.simulationTime)
                
                if cPoints < nPoints:
                    nPoints = cPoints
                if ctimeS > timeS:
                    timeS = ctimeS
                if ctimeE < timeE:
                    timeE = ctimeE
                                                    
                vascularNetwork.solutionDataFile.close()
                del vascularNetwork
            
                progressBar.progress(evaluationCaseFiles.index(batchData))  
            
            self.simulationTime = np.linspace(timeS, timeE, nPoints)
            # save the simulationTime
            self.saveQuantitiyOfInterestData(simulationTimeFileSave) 
            
        else:
            print "load estimated simulation time of all simulations"
            if simulationTimeFileLoad != None:
                self.loadQuantitiyOfInterestData(simulationTimeFileLoad)
            else:
                raise ValueError('simulationTime hdf5 file for case does not exist! {}'.format(simulationTimeFileLoad))
        
        
        self.openQuantityOfInterestFile(preprocessedSolutionData)
        print "estimate solution data for quantities of interest"
        progressBar = cPB.ProgressBar(35, len(evaluationCaseFiles))
        
        # pass the data to the locationsOfInterests which will load the information needed
        for batchData in evaluationCaseFiles:
            simulationIndex          = batchData['simulationIndex']
            networkName              = batchData['networkName']
            dataNumber               = batchData['dataNumber']
            networkXmlFileLoad       = batchData['networkXmlFileLoad']
            pathSolutionDataFilename = batchData['pathSolutionDataFilename']
            
            vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
            vascularNetwork.linkSolutionData()
            for locationOfInterest in self.locationsOfInterest.values():
                locationOfInterest.preprocessSolutionData(vascularNetwork,self.simulationTime, self.sampleSize, simulationIndex)
        
            progressBar.progress(evaluationCaseFiles.index(batchData))    
            
        # save hdf5 file
        
        # second postprocessing find extrema if needed also for variables defined over space
        ## TODO: fix data saving methods as it will not work now
        for locationOfInterest in self.locationsOfInterest.values():
            locationOfInterest.preprocessSolutionDataExtremaAndInflectionPoints(self.simulationTime, self.sampleSize)
            
        #for locationOfInterest in self.locationsOfInterest.values():
        #    locationOfInterest.preprocessSolutionDataTrajectory(self.simulationTime, self.sampleSize)
        
        self.closeAndSaveQuantityOfInterestFile()
        
    def getQoiIterator(self):
        '''
        Function that creates and iterator element which iterates through all stored qoi data elements
        '''
        qoiList = []
        for locationOfInterest in self.locationsOfInterest.values():
            for qoi in locationOfInterest.quantitiesOfInterest.itervalues():
                qoiList.append(qoi)
        return qoiList
    
    ### --- data saving and writing procedure
    
    ## wrappers for load/save of memory objects at a certain point in time
    def loadQuantitiyOfInterestData(self, hdf5File):
        '''
        function taylored to open hdf5 file, save all data into it and close it afterwards
        '''
        self.openHdf5File(hdf5File, mode='r+')
        baseGroupName = 'LocationOfInterestManager' 
        self.updateHdf5Groups(baseGroupName)
        self.loadDataHdf5()
        self.closeHdf5File()
    
    def saveQuantitiyOfInterestData(self, hdf5File):
        '''
        function taylored to open hdf5 file, load all data into it and close it afterwards
        '''
        self.openHdf5File(hdf5File, mode='w')
        baseGroupName = 'LocationOfInterestManager' 
        self.updateHdf5Groups(baseGroupName)
        self.saveDataHdf5()
        self.closeHdf5File()
    
    def openQuantityOfInterestFile(self, uqsaSolutionDataFileSave, mode = 'w'):
        
        self.openHdf5File(uqsaSolutionDataFileSave, mode)
        baseGroupName = 'LocationOfInterestManager' 
        self.updateHdf5Groups(baseGroupName)
        
        
    def closeAndSaveQuantityOfInterestFile(self):
        
        baseGroupName = 'LocationOfInterestManager' 
        self.updateHdf5Groups(baseGroupName)
        self.saveDataHdf5()
        self.closeHdf5File()