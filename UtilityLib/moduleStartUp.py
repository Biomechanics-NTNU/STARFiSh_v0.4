#!/usr/bin/env python
# -*- coding: utf-8 -*-


#############################################################################
#
# moduleStartUp
#
# includes option parser function and 
# provides functions to ask user for need data such as NetworkName, dataNumber Visualisation cases
# if not given with options while staring Main or Visualistion or VascularPolynomialChaos
#
# created by Vinzenz Eck // vinzenz.g.eck@ntnu.no
##

import os,shutil,sys
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(os.path.join(*[cur,'..']))

import inspect
from pprint import pprint as pp

import UtilityLib.moduleFilePathHandler as mFPH
import VascularPolynomialChaosLib.moduleFilePathHandlerVPC as mFPH_VPC

from optparse import OptionParser

def parseOptions(activeOptions, visualisationOnly = False, vascularPolynomialChaos = False):
    """
    parse options for the code

    enable options to parse with str in activeOptions
    input:     activeOptions := list with str

            'f' : parse filename
            'n' : dataNumber
            'd' : simulation description
            'r' : resimulated case
            's' : save
            'v' : visualisation type
            'c' : connect visualisations
            'w' : set working directory
            'p' : open workind directory settings

            visualisationOnly := bool if True proposal of visualisation cases are made if non is given

    return: dict with options and arguments

    Usage e.g., :

        optionsDict = parseOptions(['f','n','d','s','v','r','w','p'])

        networkName           = optionsDict['networkName']
        save                  = optionsDict['save']
        dataNumber            = optionsDict['dataNumber']
        simulationDescription = optionsDict['simulationDescription']
        vizOutput             = optionsDict['vizOutput']
        resimulate            = optionsDict['resimulate']
    """
    
    parser = OptionParser()
    
    for activeOption in activeOptions:
        if activeOption == 'f':
            parser.add_option("-f", "--file", dest='networkName',
                              help = "open file with given network name")
        elif activeOption == 'n':
                parser.add_option("-n", "--dataNumber", dest='dataNumber',
                                  help = "number of the solution data (last number in filename), default = 999, max 999;")
        elif activeOption == 's':
            parser.add_option("-s", "--save", action="store_true", dest="save", 
                              help = "if set solution data is saved")
        elif activeOption == 'd':
                parser.add_option('-d', '--description', dest='description',
                                  help = "simulation case description; NB: no space subported")
        elif activeOption == 'v':
                parser.add_option('-v', '--vizBool', dest='vizBool',
                                  help = "choose visualisation mode, 0: no visualisation, 1: 2d and 3d, 2: 2d plots, 3: 3d visualisation")
        elif activeOption == 'c':       
                parser.add_option("-c", "--connect", action="store_true",  dest='connect',
                                  help="connect to 3dViz (True) or not (False); currently not working")
        elif activeOption == 'r':
            parser.add_option("-r", "--resimulate", action="store_true", dest="resimulate", 
                              help = "resimulate case with same network saved in datanumber file, 0 = False, 1 = True")
        elif activeOption == 'w':
            parser.add_option("-w", "--workingDirectory", dest="workingDirectory", 
                              help = "set the absolute path of your working Directory where you the networkfiles are stored (If fresh installed use this option to set your first directory)")
        elif activeOption == 'p':
            parser.add_option("-p", "--workingDirectorySettings", action="store_true", dest="workingDirectorySettings", 
                              help = "open working directory settings menu")

    (options, args) = parser.parse_args()
    optionsDict = options.__dict__
        
    networkName = None
    save        = False
    dataNumber  = False
    dataSetNumber = False
    vizOutput   = 'non'
    simulationDescription = None
    connect     = False
    resimulate  = False
    
    
    for option,optionArgument in optionsDict.iteritems():
        # -f network name
        if option == 'networkName':
            if optionArgument != None:
                networkName = optionArgument
        # -n dataNumber
        elif option == 'dataNumber':
            dataNumber,dataSetNumber = evaluateDataNumber(optionArgument)  
            if dataNumber != None:
                save = True     
        # -s save
        elif option == 'save':
            if optionArgument != None:
                save = optionArgument
            #save solution data and the vascularNetwork in c pickle, if no save take temporary slot 999
            #if save == False:
            #    dataNumber = '999'
        # -d simulation Description
        elif option == 'description':
            if optionArgument != None:
                simulationDescription = optionArgument
        # -v visialisation type        
        elif option == 'vizBool':
            
            if optionArgument == None:
                optionArgument = defineVisualisation()
            
            if optionArgument == 1:
                vizOutput = "2D+3D"
            elif optionArgument == 2:
                vizOutput = "2D"
            elif optionArgument == 3:
                vizOutput = "3D"
        # -c connect
        elif option == 'connect':
            if optionArgument != None:
                connect = optionArgument
        # -r resimulate
        elif option == 'resimulate':
            if optionArgument != None:
                    resimulate = optionArgument
        elif option == 'workingDirectory':
            if optionArgument != None:
                insertWorkingDirectory(optionArgument)
                exit()
        elif option == 'workingDirectorySettings':
            if optionArgument != None:
                # Settings moved there to prefent circular imports
                mFPH.workingDirectorySettings()
                exit()
                    
    # catch up non given but necessary options
    ## simulation and visualisation
    if visualisationOnly == False and vascularPolynomialChaos == False:
        if networkName == None:
            networkName = chooseNetwork()
        if simulationDescription == None:
            simulationDescription = defineSimulationDescription()
        if dataNumber == None:
            dataNumber = defineDataNumber(networkName)
    ## visualisation only
    if visualisationOnly == True and (dataSetNumber == None or networkName == None) and vascularPolynomialChaos == False: 
            print "\n  No networkName passed, choose between all available networks:"
            print "  (NB: use -f networkName to define a specific file you want to open) \n"
            networkName,dataNumber = chooseSolutionDataCase()
    ## polynomial chaos
    if vascularPolynomialChaos == True:
        #TODO Critical Need to check for network template here?
        if networkName == None:
            networkName = chooseNetwork(showTemplates = False)
        if dataSetNumber == None:
            dataNumber = chooseUQSACaseFile(networkName)
        
    del parser
    
    return {'networkName'           : networkName,
            'save'                  : save,
            'dataNumber'            : dataNumber,
            'dataSetNumber'         : dataSetNumber,
            'vizOutput'             : vizOutput,
            'simulationDescription' : simulationDescription,
            'connect'               : connect,
            'resimulate'            : resimulate}

def prettyPrintList(title, listToPrint, indexOffSet = 0):
    """
    Function to pretty print a list to STDOUT with numbers to choose from
    """
    print title
    for index,listElement in enumerate(listToPrint):
        print "   [ {:3} ] - {}".format(index+indexOffSet,listElement)

def userInputEvaluationInt(maxBound, minBound=0, question = "    insert your choice, (q)-quit: "):
    '''
    Question user to isert an integer number between minBound and maxBound
    '''
    appropriateInputList = [str(int(i+minBound)) for i in xrange(maxBound-minBound)]
    userInput = "NONE"
    appropriateInputList.append('q')
    print ""
    while userInput not in appropriateInputList:
        userInput = raw_input(question)
    print ""
    if userInput == 'q': exit()
    else: return int(userInput)


def chooseNetwork(showTemplates = True):
    """
    console Interface to choose a VascularNetwork for simulation / vascularPolynomial Chaos
    """
    dirNamesTemplate = []
    prettyPringOffset = 0
    if showTemplates:
        # network templates
        templatePath = mFPH.getDirectory('networkXmlFileTemplateDirectory','','','read')
        dirNamesTemplate = [d for d in os.listdir(templatePath) if '.' not in  d]
        prettyPrintList("\n Template Networks: \n",dirNamesTemplate)
        prettyPringOffset = len(dirNamesTemplate)
    # working directory
    workingDirectoryPath = mFPH.getDirectory('workingDirectory','','','read')
    dirWorkingDirectory = [d for d in os.listdir(workingDirectoryPath) if '.' not in  d]
    prettyPrintList("\n WorkingDirectory Networks: \n", dirWorkingDirectory, indexOffSet = prettyPringOffset)
    
    dirNames = dirNamesTemplate+dirWorkingDirectory    
    userInput = userInputEvaluationInt(len(dirNames), 0)
    print ""
    print '====================================='
    return dirNames[userInput]

# TODO: (einar) fix exception variable
def evaluateDataNumber(dataNumberString, exception = "Error"):
    """
    Function to evaluate DataNumbers given as a string
    Max lenght of dataNumber = 3
    Input:
        dataNumber String sparate datanumbers with ,
        (e.g. '12,23' or '4')
    output:
        dataNumber = String (of the first one)
        dataSetNumber = [ String, ... ,String] 
        (e.g. dataNumber = '012' dataSetNumber = ['012','023']
              dataNumber = '004' dataSetNumber = ['004'] 
    """
    dataSetNumber = None
    dataNumber = None
    if dataNumberString != None:
        dataNumber = dataNumberString
        if ',' in dataNumber:
            dataNumbers = dataNumber.split(',')
            dataSetNumber = []
            for dataNum in dataNumbers:
                if len(dataNum) < 4:
                    dataSetNumber.append(dataNum.zfill(3))
                else:
                    if exception == "Error":
                        raise ValueError('moduleStartUp.evaluateDataNumber. Datanumber {} to high! system exit'.format(dataSetNumber))
                        
                    elif exception == 'Warning':
                        print 'moduleStartUp.evaluateDataNumber. Datanumber {} to high'.format(dataSetNumber)
                        return False,False
                    else:
                        raise Exception
                    
            dataNumber = dataSetNumber[0]
        else:
            dataSetNumber = [dataNumber.zfill(3)]
            dataNumber = dataNumber.zfill(3)
            if len(dataNumber) > 3:
                if exception == "Error":
                    raise ValueError('moduleStartUp.evaluateDataNumber. Datanumber {} to high! system exit'.format(dataSetNumber))
                    
                elif exception == 'Warning':
                    print 'moduleStartUp.evaluateDataNumber. Datanumber {} to high'.format(dataSetNumber)
                    return False,False
                else:
                    raise Exception   
                
    return dataNumber,dataSetNumber
  
  
def defineVisualisation():
    """
    Function to ask for visualisation type
    
    returns: visualisationBool (int)
    """
    listToPrint = ['No visualisation', '2d and 3d', '2d visualisation', '3d visualisation']
    prettyPrintList("\n        Choose visualisation mode:",listToPrint)
    userInput = userInputEvaluationInt(4, 0, "     What to do? (q)-quit: ")
    return userInput
    
    
  
def defineSimulationDescription():
    simulationDescription = str(raw_input("\n  Type in description of the simulation case: "))
    try:
        simulationDescription = str(simulationDescription)
    except:
        print 'ERROR: no String Convertable input given, system exit'
        exit()
    if simulationDescription in ['',' ']: simulationDescription = '-'
    return simulationDescription
    
def defineDataNumber(networkName):
    
    print "No solution-datanumber defined (3 characters)! It needs to be defined"

    #TODO CRITICAL this fails for template networks
    existingDataNumbers = findExistingDataNumbers(networkName)
    correctDataNumber = False
    
    prettyPrintList("\n        Existing dataNumbers for this network: {}".format(networkName),existingDataNumbers)
       
    while correctDataNumber == False:
        dataNumber = str(raw_input("\n  Please enter datanumber (3 characters): "))
        testedFailed = False
        try:
            dataNumber,xx = evaluateDataNumber(dataNumber)
        except ValueError as e:
            print e
            testedFailed = True
            
        if dataNumber in existingDataNumbers:
            #userInput = str(raw_input
            listToPrint = [" keep this datanumber and overwrite simulation case",
                           " enter new data number"]
            prettyPrintList("\n        Simulation case with datanumber exits already".format(networkName),listToPrint)
            userInput = userInputEvaluationInt(2, 0, "     What to do? ")
            if userInput == 0:
                correctDataNumber = True
            
            elif userInput == 0:
                dataNumber = False
                correctDataNumber = False
                
        elif testedFailed == False:
            correctDataNumber = True
            
    return dataNumber     

            
def findExistingDataNumbers(networkName):
    """
    evaluates exiting data numbers of a network
    
    Args: networkName (str) := networkName
    
    Returns: existingDataNumbers (list)
    """
    # TODO: better way to handle this case?
    destinationNetworkName = networkName
    if "_template" in networkName:
        destinationNetworkName = ''.join(networkName.split("_template"))
        try:
            _ = mFPH.getFilePath('simulationDescriptionFile', destinationNetworkName, "xxx", 'read', exception = 'No')
        except:
            _ = mFPH.getFilePath('simulationDescriptionFile', destinationNetworkName , "xxx", 'write')#, exception = 'No')

    simulationCaseDict = mFPH.getSimulationCaseDescriptions(destinationNetworkName)#, exception = 'No')
    print(simulationCaseDict)
    networkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory',destinationNetworkName,'xxx','read')
    
    existingDataNumbers = []
    
    for root, dirs, files in os.walk(networkDirectory):
        for file in files:
            if ".hdf5" in file and "_SolutionData_" in file:
                solutionDataFile = file.split('.')[0]
                dataNumber = solutionDataFile.split('_SolutionData_')[-1]
                if len(dataNumber) == 3:
                    existingDataNumbers.append(dataNumber)
    return existingDataNumbers
    
    
def chooseSolutionDataCase():
    """
    console Interface to choose a vascular1DFlow simulation case for e.g. Visualisation
    Output:
        networkName <string>
        dataNumber  <string>
    """
    workingDirectory = mFPH.getDirectory('workingDirectory','','','read')
    networkCases = [d for d in os.listdir(workingDirectory) if '.' not in  d]
        
    fileNameDataNumber = []
    
    indexOffSet = 0
    for networkName in networkCases:
        simulationCaseDict = mFPH.getSimulationCaseDescriptions(networkName )#, exception = 'No')
        networkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory',networkName,'xxx','read')
        
        listToPrint = []
        if simulationCaseDict != None:
            first = True
            for root, dirs, files in os.walk(networkDirectory):
                for file in files:
                    if ".hdf5" in file and "uqsa" not in file:
                        solutionDataFile = file.split('.')[0]
                        dataNumber = solutionDataFile.split('_SolutionData_')[-1]
                        if len(dataNumber) == 3:
                            if dataNumber not in simulationCaseDict:
                                description =  "'{}' not listed in simulation descriptions of network '{}'.".format(dataNumber,networkName)
                            else:
                                description = simulationCaseDict[dataNumber]
                                
                            listToPrint.append("{} : {}".format(dataNumber, description))
                            fileNameDataNumber.append([networkName,dataNumber])
    
        prettyPrintList("\n        {}".format(networkName),listToPrint, indexOffSet = indexOffSet)
        indexOffSet = indexOffSet+len(listToPrint)
        
    if len(fileNameDataNumber) == 0:
        print "No solutionCases available, system exit"
        exit()
    
    question  = "  Choose simulation case you want to open according to its number, (q)-quit: "
    userInput = userInputEvaluationInt(len(fileNameDataNumber), 0, question)
    networkName = fileNameDataNumber[userInput][0]
    dataNumber  = fileNameDataNumber[userInput][1]
    
    return networkName,dataNumber


def chooseUQSACaseFile(networkName):
    """
    console Interface to choose a vascularPolynomialChaos Config File
     including the possility to create a template Config File
    Input:
        networkName of VascularNetwork
    Output:
        networkName,dataNumber of the Config File (networkName should be the same)
    """
    from VascularPolynomialChaosLib import classUqsaCase
    
    workingDirectory = mFPH.getDirectory('workingDirectory','','','read')
    networkDirectory = os.path.join(*[workingDirectory,networkName])
    
    filesNetworkDir = os.listdir(networkDirectory)
    filenames = []
    
    for fileNetworkDir in filesNetworkDir:
        # check if polychaos directory:
        if 'vascularPolynomialChaos' in fileNetworkDir:
            allFilenames = os.listdir(os.path.join(*[workingDirectory,networkName,fileNetworkDir]))
            for filename in allFilenames:
                if ".xml" in filename and "uqsaCase" in filename:
                    filenames.append(filename)
            
    print "\n  No dataNumber for UQSAcase file passed, choose between all available UQSAcase files:"
    print "  (NB: use -n dataNumber to define a specific UQSAcase you want to open)\n"
    print "   [   0 ] - Create new template uqsa case file and exit"
    if filenames != []:
        prettyPrintList('',filenames, indexOffSet = 1)
        
    question  = "  Choose Option or Config-File you want to open according to its number:, (q)-quit: "
    userInput = userInputEvaluationInt(1+len(filenames), 0, question)
        
    if userInput in [0] :
        
        userInputDataNumber = 'xxxx'
        dataNumber = False
        while dataNumber == False:
            userInputDataNumber = str(raw_input("\n  Insert dataNumber for polynomial Chaos case (overwrites if existing): "))
            dataNumber = evaluateDataNumber(userInputDataNumber, exception = "Warning")[0]
        
        # create template configuration
        configurationFilePathTemplate = mFPH_VPC.getFilePath('uqsaCaseTemplateFile', networkName, dataNumber, 'read')
        
        uqsaCase = classUqsaCase.UqsaCase()
        uqsaCase.loadXMLFile(configurationFilePathTemplate)
        
        # TODO: check if file exits and ask if it should be overwritten        
        configurationFilePath = mFPH_VPC.getFilePath('uqsaCaseXmlFile', networkName, dataNumber, 'write')
        uqsaCase.writeXMLFile(configurationFilePath)
        # copy network file
        toCopyFile = mFPH.getFilePath('networkXmlFile', networkName, 'xxx','write')
        destinationFile = mFPH_VPC.getFilePath('vpcNetworkXmlFile', networkName, dataNumber,'write')
        shutil.copy(toCopyFile, destinationFile)
        print "files created!, exit()"
        exit()
    
    else:
        networkName = filenames[userInput-1]
        dataNumber = networkName.split('.')[0].split('_')[-1]
        
    return dataNumber   
    
    
