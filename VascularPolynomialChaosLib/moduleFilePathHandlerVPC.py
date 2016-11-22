#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
#
# moduleFilePathHandler VPC
#
# provide functions to save, load, maipulate, parse etc. pickle files of simulation cases
#
# loadSolutionDataFile(networkName,dataNumbers)
# parseDirectoryForSimulationCases(networkName)
# updateSimulationDescriptions(networkName)
#
#
# created by Vinzenz Eck // vinzenz.g.eck@ntnu.no
##

import cPickle
import os,sys,shutil
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')

from copy import copy as copy 
from pprint import pprint as pp

import UtilityLib.moduleFilePathHandler as mFPH


def getFilePath(fileType, networkName, dataNumber, mode, caseName = "None", evaluationNumber = "None", exception = 'Error'):
    '''
    Function return a requested file path, if this file exists
    
    Args:
        fileType (str): 'vpcConfigXmlFile',
                        'vpcConfigTemplateFile',
                        'vpcNetworkXmlFile',
                        'uqsaSampleFile',
                        'uqsaEvaluationNetworkXmlFile',
                        'uqsaEvaluationSolutionDataFile',
                        'vpcProcessedSolutionDataFile',
                        'evaluationLogFile',
                        'vpcSolutionDataFile'
                        'simulationTime'
                        'preprocessedDataFile'
        
        networkName (str):name of the network file
        dataNumber (str):data number of solution file or xml file
        mode (str): read or write
        exception (str): (for read mode),
                    Error (default): raise error and exit if file is not exiting,
                    Warning: just raise Warning and return with error string,  
                    
    Returns:
        requestedFilePath (str): path to the requested file 
    '''
    existingFileTypes = ['uqsaCaseXmlFile',
                         'uqsaCaseTemplateFile',
                         'vpcNetworkXmlFile',
                         'uqsaSampleFile',
                         'uqsaEvaluationNetworkXmlFile',
                         'uqsaEvaluationSolutionDataFile',
                         'vpcProcessedSolutionDataFile',
                         'evaluationLogFile',
                         'uqsaSolutionDataFile',
                         'simulationTime',
                         'preprocessedDataFile']
    
    if fileType not in existingFileTypes:
        raise ValueError("ERROR: getFilePath, requested file type {}\
                          is not in existingFileTypess {}".format(fileType, existingFileTypes))
        
    # file names
    filenames = {
                 'uqsaCaseTemplateFile'          : 'uqsaCase_xxx.xml',
                 'uqsaCaseXmlFile'               : ''.join([networkName,'_uqsaCase_',dataNumber,'.xml']),
                 'vpcNetworkXmlFile'             : ''.join([networkName,'_vpc_',dataNumber,'.xml']),
                 'uqsaSampleFile'                : ''.join(['samples_',caseName,'.hdf5']),
                 'simulationTime'                : ''.join(['simulationTime_',caseName,'.hdf5']),
                 'uqsaEvaluationNetworkXmlFile'  : ''.join([networkName,'_evaluation_',str(evaluationNumber).zfill(7),'.xml']),
                 'uqsaEvaluationSolutionDataFile': ''.join([networkName,'_evaluation_',str(evaluationNumber).zfill(7),'.hdf5']),
                 'evaluationLogFile'             : ''.join(['evaluationLogFile.txt']),
                 'uqsaSolutionDataFile'          : ''.join([networkName,'_uqsa-SolutionData_',dataNumber,'.hdf5']),
                 'preprocessedDataFile'          : ''.join(['preprocessedData_',caseName,'.hdf5']),
                 }    
        
    ## find requested file name
    requestedFilename  = filenames[''.join([fileType])]
    ## find directory    
    requestedDirectory = getDirectory(''.join([fileType,'Directory']), 
                                      networkName, 
                                      dataNumber,
                                      mode,
                                      caseName = caseName,
                                      exception = exception)
    if requestedDirectory == None:
        if exception == "Warning":
            print "WARNING: moduleFilePathHandler.getFileAndPaths() directory of file '{}' does not exits. Exit()".format(requestedFilename)
            return None
        elif exception == "No":
            pass
        else:
            raise ValueError("ERROR: moduleFilePathHandler.getFileAndPaths() directory of file '{}' does not exits. Exit()".format(requestedFilename))
              
    ## requested file path 
    requestedFilePath = ''.join([requestedDirectory,'/',requestedFilename])
                
    if mode == 'read':    
        # if mode read
        ## ensure that the file exists
        if not os.path.isfile(requestedFilePath):
            if exception == "Warning":
                print "WARNING: moduleFilePathHandler.getFileAndPaths() file '{}' does not exits. Exit()".format(requestedFilePath)
                return None
            elif exception == "No":
                print "raise no exception"
                return None
            else:
                raise ValueError("ERROR: moduleFilePathHandler.getFileAndPaths() file '{}' does not exits. Exit()".format(requestedFilePath))
              
    return requestedFilePath
    
def getDirectory(directoryType, networkName, dataNumber, mode, exception = 'Error', caseName = 'None'):
    '''
    Function returns a requested directory path, if this directory does not exists
    it is created.
    
    Args:
        directoryType (str):  'workingDirectory',
                              'vpcConfigXmlFileDirectory',
                              'vpcNetworkXmlFileDirectory',
                              'uqsaEvaluationNetworkXmlFileDirectory',
                              'uqsaEvaluationSolutionDataFileDirectory',
                              'vpcSampleFileDirectory',
                              'evaluationLogFileDirectory',
                              'vpcSolutionDataFileDirectory',
                              'uqsaCaseTemplateFile',
                              'simulationTimeDirectory'
        networkName (str): name of the network file
        dataNumber (str): data number of solution file or xml file
        mode (str): read or write
        exception (str): (for read mode) 
                    Error (default): raise error and exit if file is not exiting-
                    Warning: just raise Warning and return with error string      
    Returns:
        requestedDirectory (str): path to the requested directory
    '''
    
    existingDirectoryTypes = {'workingDirectory',
                              'uqsaCaseXmlFileDirectory',
                              'vpcNetworkXmlFileDirectory',
                              'uqsaEvaluationNetworkXmlFileDirectory',
                              'uqsaEvaluationSolutionDataFileDirectory',
                              'uqsaSampleFileDirectory',
                              'evaluationLogFileDirectory',
                              'uqsaSolutionDataFileDirectory',
                              'uqsaCaseTemplateFileDirectory',
                              'simulationTimeDirectory',
                              'preprocessedDataFileDirectory'} 
    
    if directoryType not in existingDirectoryTypes:
        raise ValueError("ERROR: getDirectory, requested directoryType {}\
                          is not in existingDirectoryTypes{}".format(directoryType, existingDirectoryTypes))
    ##definitions
    starfishHomeDirectory       = ''.join([cur,'/..'])
    workingDirectory            = mFPH.readConfigFile(['WorkingDirectory'])['WorkingDirectory']
    networkXmlFileDirectory     = ''.join([workingDirectory,'/',networkName])
    vpcCaseDirectory            = ''.join([networkXmlFileDirectory,'/vascularPolynomialChaos_',str(dataNumber)]) 
    vpcOrderMethodDirectory     = ''.join([vpcCaseDirectory,'/','samplingMethod',caseName])
    vpcEvaluationNetDirectory   = ''.join([vpcOrderMethodDirectory,'/evaluationNetworkFiles'])
    vpcEvaluationSolDirectory   = ''.join([vpcOrderMethodDirectory,'/evaluationSolutionData'])
    vpcConficTemplateDicrectory = ''.join([starfishHomeDirectory,'/TemplateNetworks','/vascularPolynomialChaos'])
    
    ## look up tables
    # directories
    directories = {
                   'workingDirectory'                       : workingDirectory,
                   # vascular polynomial chaos
                   'uqsaCaseXmlFileDirectory'               : vpcCaseDirectory,
                   'vpcNetworkXmlFileDirectory'             : vpcCaseDirectory,
                   'uqsaSampleFileDirectory'                : vpcOrderMethodDirectory,
                   'uqsaEvaluationNetworkXmlFileDirectory'  : vpcEvaluationNetDirectory,
                   'uqsaEvaluationSolutionDataFileDirectory': vpcEvaluationSolDirectory,
                   'simulationTimeDirectory'                : vpcOrderMethodDirectory,
                   'evaluationLogFileDirectory'             : vpcOrderMethodDirectory,
                   'preprocessedDataFileDirectory'          : vpcOrderMethodDirectory,
                   'uqsaSolutionDataFileDirectory'          : vpcOrderMethodDirectory,
                   'uqsaCaseTemplateFileDirectory'          : vpcConficTemplateDicrectory
                   }
    
    requestedDirectory = os.path.normpath(directories[directoryType])
    
    # if mode write
    if mode == 'write':
        ## ensure that the directory exists
        if not os.path.exists(requestedDirectory):
            os.makedirs(requestedDirectory)  
    if mode == 'read':
        if not os.path.exists(requestedDirectory):
            requestedDirectory = None
    
    return requestedDirectory
