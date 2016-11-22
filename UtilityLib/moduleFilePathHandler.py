#####
# all functions to save, load, maipulate, parse etc. pickle files of simulation cases
# .v1d
####

#############################################################################
#
# moduleFilePathHandler
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
sys.path.append(''.join([cur,'/../']))

import ConfigParser

try:
    from lxml import etree
except:
    # TODO: This produces an error!! with the pretty_print argument
    from xml.etree import ElementTree as etree


# TODO: (einar) Rename input variable exception and corresponding strings
def getFilePath(fileType, networkName, dataNumber, mode, exception = 'Error'):
    """
    Function return a requested file path, if this file exists

    Args:
        
        fileType (str):
            'randomVariableCSVFile'
            'vesselCSVFile',
            'boundaryCSVFile',
            'networkXmlFileTemplate',
            'networkXmlFile',
            'solutionFile',
            'configFile',
            'simulationDescriptionFile',
            'vncRescentNetworksFile',
            'vncNetworkGraphFile

        networkName (str): name of the network file

        dataNumber (int): data number of solution file or xml file

        mode (str): 'read' or 'write'


    Returns:
        string: file path

    Raises:
        Error (default): (in read mode) raise error and exit if file is not exiting
        Warning: (in read mode) just raise Warning and return with error string
    """    
    
    existingFileTypes = ['randomVariableCSVFile',
                         'vesselCSVFile',
                         'boundaryCSVFile',
                         'networkXmlFileTemplate',
                         'networkXmlFile',
                         'solutionFile',
                         'configFile',
                         'simulationDescriptionFile',
                         'vncRescentNetworksFile',
                         'vncNetworkGraphFile']
    
    if fileType not in existingFileTypes:
        raise ValueError("ERROR: getFilePath, requested file type {}\
                          is not in existingFileTypess {}".format(fileType, existingFileTypes))
    
    if '_template' in networkName:
        if fileType == 'networkXmlFile':
            fileType = 'networkXmlFileTemplate'
    
    # file names
    filenames = {
                 'configFile'                : 'STARFiSh.config',
                 'vesselCSVFile'             : ''.join([networkName,'.csv']),
                 'boundaryCSVFile'           : ''.join([networkName,'BC.csv']),
                 'networkXmlFileTemplate'    : ''.join([networkName,'.xml']),
                 'networkXmlFileXXX'         : ''.join([networkName,'.xml']),
                 'networkXmlFileSim'         : ''.join([networkName,'_SolutionData_',dataNumber,'.xml']),
                 'solutionFile'              : ''.join([networkName,'_SolutionData_',dataNumber,'.hdf5']),
                 'simulationDescriptionFile' : ''.join(['simulationCaseDescriptions.txt']),
                 'vncRescentNetworksFile'    : '.recentNetworkNames.pickle',
                 'vncNetworkGraphFile'       : ''.join([networkName,'.png']),
                 'randomVariableCSVFile'     : ''.join([networkName,'RV.csv']),
                 }    
                
    ## if fileType=networkXmlFile check if master (dn = XXX) or simulated is meant
    if fileType == 'networkXmlFile':
        if dataNumber == 'xxx':
            fileType = ''.join([fileType,'XXX'])
        else:
            fileType = ''.join([fileType,'Sim'])
    
    ## find requested file name
    requestedFilename  = filenames[''.join([fileType])]
    ## find directory
    requestedDirectory = getDirectory(''.join([fileType,'Directory']), networkName, dataNumber, mode, exception)
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
                #print "raise no exception"
                return None
            else:
                raise ValueError("ERROR: moduleFilePathHandler.getFileAndPaths() file '{}' does not exits. Exit()".format(requestedFilePath))
              
    return requestedFilePath

# TODO: (einar) fix exception variable in function
def getDirectory(directoryType, networkName, dataNumber, mode, exception = 'Error'):
    """
    Function returns a requested directory path, if this directory does not exists
    it is created.

    directoryType:
        'randomVariableCSVFileDirectory'
        'workingDirectory',
        'configFileDirectory',
        'vesselCSVFileDirectory',
        'boundaryCSVFileDirectory',
        'networkXmlFileTemplateDirectory',
        'networkXmlFileXXXDirectory',
        'networkXmlFileSimDirectory',
        'solutionFileDirectory',
        'screenshotDirectory',
        'movieDirectory',
        'simulationDescriptionFileDirectory',
        'vncRescentNetworksFileDirectory',
        'vncNetworkGraphFileDirectory'


    networkName:

        name of the network file

    dataNumber:

        data number of solution file or xml file

    mode:

        read or write

    exception: (for read mode)

        Error (default): raise error and exit if file is not exiting
        Warning: just raise Warning and return with error string
    """
    
    existingDirectoryTypes = {
                              'workingDirectory',
                              'configFileDirectory',
                              'vesselCSVFileDirectory',
                              'boundaryCSVFileDirectory',
                              'networkXmlFileTemplateDirectory',
                              'networkXmlFileXXXDirectory',
                              'networkXmlFileSimDirectory',
                              'solutionFileDirectory',
                              'screenshotDirectory',
                              'movieDirectory',
                              'simulationDescriptionFileDirectory',
                              'vncRescentNetworksFileDirectory',
                              'randomVariableCSVFileDirectory',
                              'vncNetworkGraphFileDirectory'} 
    
    if directoryType not in existingDirectoryTypes:
        raise ValueError("ERROR: getDirectory, requested directoryType {}\
                          is not in existingDirectoryTypes{}".format(directoryType, existingDirectoryTypes))
    ##definitions
    starfishHomeDirectory = ''.join([cur,'/..'])
    
    if directoryType != 'configFileDirectory': 
        # load working directory from config file
        workingDirectory = readConfigFile(['WorkingDirectory'])['WorkingDirectory']
    else:
        workingDirectory = starfishHomeDirectory
           
    networkXmlFileTemplateDirectory = ''.join([starfishHomeDirectory,'/TemplateNetworks/',networkName])
    networkXmlFileDirectory         = ''.join([workingDirectory,'/',networkName])
    solutionFileDirectory           = ''.join([networkXmlFileDirectory,'/SolutionData_',str(dataNumber)])
    movieDirectory              = ''.join([solutionFileDirectory,'/Movies'])
    screenshotDirectory         = ''.join([solutionFileDirectory,'/Screenshots'])
    ## look up tables
    # directories
    directories = {
                   'workingDirectory'                   : workingDirectory,
                   'configFileDirectory'                : starfishHomeDirectory,
                   'vesselCSVFileDirectory'             : networkXmlFileDirectory,
                   'boundaryCSVFileDirectory'           : networkXmlFileDirectory,
                   'randomVariableCSVFileDirectory'     : networkXmlFileDirectory,
                   'networkXmlFileTemplateDirectory'    : networkXmlFileTemplateDirectory,
                   'networkXmlFileXXXDirectory'         : networkXmlFileDirectory,
                   'networkXmlFileSimDirectory'         : solutionFileDirectory,
                   'solutionFileDirectory'              : solutionFileDirectory,
                   'simulationDescriptionFileDirectory' : networkXmlFileDirectory,
                   # 3d viz
                   'screenshotDirectory'                : screenshotDirectory,
                   'movieDirectory'                     : movieDirectory,
                   # vnc
                   'vncRescentNetworksFileDirectory'    : workingDirectory,
                   'vncNetworkGraphFileDirectory'       : networkXmlFileDirectory
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

def createWorkingCopyOfTemplateNetwork(templateNetworkName, destinationNetworkName = None):
    """
    Function which copys all data from a template network wtih template network name into
    the working directory.
    It uses the same name of the template network if destinationNetworkName is not defined
    """
    if destinationNetworkName == None:
        destinationNetworkName = templateNetworkName.split('_template')[0]
    
    pathTemplateNetwork     = getDirectory('networkXmlFileTemplateDirectory', templateNetworkName, 'xxx', 'read')
    pathDestinationNetwork  = getDirectory('networkXmlFileXXXDirectory', destinationNetworkName, 'xxx', 'write')
    
    #loop through files
    for file in os.listdir(pathTemplateNetwork):
        # remove _template from name
        renamedFile = ''.join(file.split('_template'))
        # check if new name needs to be applied
        oldName = templateNetworkName.split('_template')[0]
        if oldName in renamedFile:
            renamedFile = ''.join([destinationNetworkName,renamedFile.split(oldName)[-1]])

        shutil.copy(os.path.join(*[pathTemplateNetwork,file]), os.path.join(*[pathDestinationNetwork,renamedFile]))
        
        newFilePath = os.path.join(*[pathDestinationNetwork,renamedFile])
        if ".xml" in newFilePath:
            setFlowFromFilePathToAbsolute(newFilePath, pathDestinationNetwork)

        
    return destinationNetworkName


def readConfigFile(options):
    """
    Function to read from options from STARFiSh.ini config file
    input:
        options = list with options
                existing options: 'WorkingDirectory'
    output:
        configurations = dict with {option: configuration from file}
    """ 
    
    config = ConfigParser.ConfigParser()
    
    filePath = getFilePath('configFile', "", '', 'read',exception = 'No')
    if filePath == None:
        saveConfigFile({'WorkingDirectory': '', 'knownWorkingDirectories':''})
        filePath = getFilePath('configFile', "", '', 'read')
    
    config.read(filePath)
        
    
    configurations = {}
    
    for option in options:    
        if option == 'WorkingDirectory':
            #try:
            workingDirectory = config.get('Directory Paths', option)
            if os.path.isdir(workingDirectory) == False:
                print Warning("\n ERROR WorkingDirectory {} does not exist \n".format(workingDirectory))
                
                knownWorkingDirectories = config.get('Directory Paths', 'knownWorkingDirectories').split(',')
                if workingDirectory in knownWorkingDirectories:            
                    knownWorkingDirectories.remove(workingDirectory) 
                    
                saveConfigFile({'knownWorkingDirectories':','.join(knownWorkingDirectories)})
                workingDirectory = ''
           # except:
            #    workingDirectory = None
           #     raise ValueError("ERROR pathAndFilenameHandler.readConfigFile reading WorkingDirectory failed ini file corrupted, exit()")
            if workingDirectory == '':
                print Warning("ERROR pathAndFilenameHandler.readConfigFile reading WorkingDirectory failed: no path defined")
                
               
                
                workingDirectorySettings(searchKnowWorkingDirectories = False)
            
            configurations['WorkingDirectory'] = workingDirectory
            
        
        if option == 'knownWorkingDirectories':
            try:
                knownWorkingDirectories = config.get('Directory Paths', option)
            except:
                knownWorkingDirectories = None
                raise ValueError("ERROR pathAndFilenameHandler.readConfigFile reading <knownWorkingDirectories> failed ini file corrupted, exit()")
            if knownWorkingDirectories == '':
                configurations['knownWorkingDirectories'] = []
            else:
                knownWorkingDirectories = knownWorkingDirectories.split(',')
                configurations['knownWorkingDirectories'] = knownWorkingDirectories
            
    return configurations    

def saveConfigFile(configurations):
    """
    Function to save configurations to options in the STARFiSh.ini config file
    The file will update the existing config file, it is not neccessary to pass all 
    configurations which exist.
    input:
        configurations = dict with {option: configuration from file}
        
    """ 
    # open config to get current states
    existingOptions = ['WorkingDirectory','knownWorkingDirectories']
    
    Config = ConfigParser.ConfigParser()
    
    configFilePath = getFilePath('configFile', '','',  'read',exception = 'No')
    
    if configFilePath is not None:  #  file exists
        Config.read(configFilePath)
    else: #  file does not exist
        Config.add_section('Directory Paths')
        
    for option,config in configurations.iteritems(): 
            if option in existingOptions:
                Config.set('Directory Paths', option, config)
                    
    with open(getFilePath('configFile', '','',  'write'), 'wb') as configfile:
        Config.write(configfile)
    
def updateKnownWorkingDirectories():
    """
    Function which updates the known working directories by adding the current working directory 
    list of known working directories
    """
    # 1. get known working directories
    try:
        knownWorkingDirectories = readConfigFile(['knownWorkingDirectories'])['knownWorkingDirectories']
    except ValueError:
        knownWorkingDirectories = []
    # 2. get current working directory
    currentWorkingDirectory = readConfigFile(['WorkingDirectory'])['WorkingDirectory']
    
    # 3. check if current working directory is not already in known working directories
    if currentWorkingDirectory not in knownWorkingDirectories:
        knownWorkingDirectories.append(currentWorkingDirectory)
        
        saveConfigFile({'knownWorkingDirectories':','.join(knownWorkingDirectories)})
    
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


def workingDirectorySettings(searchKnowWorkingDirectories = True):
    '''
    working directory settings
    '''
    
    prettyPrintList(' Working directory settings menu',['add working directory','switch to another known working directory'])
    if searchKnowWorkingDirectories == True: 
        updateKnownWorkingDirectories()
        print "\n current working directory: {} ".format(readConfigFile(['WorkingDirectory'])['WorkingDirectory'])
    userInput = userInputEvaluationInt(2)
    if userInput == 0:
        insertWorkingDirectory(None)
    elif userInput ==1:
        knownWorkingDirectories = readConfigFile(['knownWorkingDirectories'])['knownWorkingDirectories']
        prettyPrintList(' List of all known working directories:',knownWorkingDirectories)
        userInput2 = userInputEvaluationInt(len(knownWorkingDirectories))
        saveConfigFile({'WorkingDirectory': knownWorkingDirectories[userInput2]})

def insertWorkingDirectory(optionArgument):
    
    print "Setting new working directory"
    
    if optionArgument == None:
        optionArgument = ""
        loopQuestion = True
        defaultWD = os.path.expanduser(os.path.join('~','starfish_working_directory')) 
        print("Enter an absolute or relative path to set the working directory, or (q) quit")
        input_prompt="[Press enter to use the default path {}]: ".format(defaultWD)
        optionArgument = raw_input(input_prompt) or defaultWD
        
    if optionArgument != "q":
        optionArgument = os.path.expanduser(optionArgument)
    
        if os.path.isdir(optionArgument):
            saveConfigFile({'WorkingDirectory':optionArgument})
            updateKnownWorkingDirectories()
            print "   working directory set!"
        else:
            print "  working directory does not exist! try to create folder"
            try:
                os.mkdir(optionArgument)
                saveConfigFile({'WorkingDirectory':optionArgument})
                updateKnownWorkingDirectories()
                print "   created working directory folder successfully"
                print "   working directory set!"
            except:
                print "  WARNING: moduleStartUp.insertWorkingDirectory() could not set WorkingDirectory {} directory does not exists!".format(optionArgument)
        
    
def updateSimulationDescriptions(networkName, currentDataNumber, currentDescription):
    """
    Function to update the text-file with the simulation description for the given network:
    Input:
        networkName <String> (= name of the network to process)
    
    workflow:
        1. open all pickle files and write out the simulation descriptions and datanumbers
        2. write information into file
    """
    # open File
    #simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'read')#, exception = 'No')
    try:
        simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'read', exception = 'No')
        simCaseDescFile = open(simCaseDescFilePath, 'r')
    except:
        simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'write')#, exception = 'No')
    
        simCaseDescFile = file(simCaseDescFilePath, 'w+')
        simCaseDescFile.write("DataNumber   Description \n")
        simCaseDescFile.close()
        simCaseDescFile = open(simCaseDescFilePath, 'r')
            
    alreadyWritten = False
    # check if datanumber in file:
    simCaseDescriptionFileLines = simCaseDescFile.readlines()
    simCaseDescFile.close() 
    
    for Line in simCaseDescriptionFileLines:
        if currentDataNumber in Line: 
            if currentDescription not in Line or currentDescription == '':
                index = simCaseDescriptionFileLines.index(Line)
                simCaseDescriptionFileLines.remove(Line)
                simCaseDescriptionFileLines.insert(index,"  {} {} \n".format(currentDataNumber.ljust(10),currentDescription))
            alreadyWritten = True
            
    if alreadyWritten == False:
        simCaseDescriptionFileLines.append("  {} {} \n".format(currentDataNumber.ljust(10),currentDescription))
    
    #print simCaseDescriptionFileLines
    
    simCaseDescFile = open(simCaseDescFilePath, 'w')
    simCaseDescFile.writelines(simCaseDescriptionFileLines)
    simCaseDescFile.close()

def getSimulationCaseDescriptions(networkName, exception = 'Warning'):
    """

    """
    simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, 'xxx', 'read', exception = 'No')
    try:
        simCaseDescFile = open(simCaseDescFilePath, 'r')
    except:
        if exception == 'Warning':
            print "WARNING getSimulationCaseDescriptions() simulation description file of network {} does not exist!".format(networkName)
        elif exception == 'No':
            pass
        else: raise ValueError(exception)
        return None
    simCaseDescriptionFileLines = simCaseDescFile.readlines()
    simCaseDescFile.close() 
    
    dataNumberDescriptionDict = {}
    for Line in simCaseDescriptionFileLines:
        sol = Line.split('{:8}'.format(' '))
        if len(sol)>1:
            #dataNumber
            dataNumberDescriptionDict[sol[0].split(' ')[-1]] = sol[1].split(' \n')[0] 
    
    return dataNumberDescriptionDict
    
   
def regenerateSimulationCaseDescription(networkName):
    """
    Function to regenerate the simulation case descriptions of a network
    """
    
        
def loadExternalDataSet(fileName):
    """
    Function to open external (preprocessed) DataSets with the file ending *.v1dfExD
    Input:
        fileName <string> (with the total path)
    Output:
        extData <dict> 
    """
    
    externalData = {'PressureTime':None,'Pressure':None,'PressureUnit':'',
                    'FlowTime':None,'Flow':None,'FlowUnit':'',
                    'AreaTime':None,'Area':None,
                    'Description': ''}
    try:
        externalDataFile = open(fileName,'rb')
        externalData = cPickle.load(externalDataFile)
    except:
        print "Error: no or corrupted external data-file found"
    
    return externalData

def setFlowFromFilePathToAbsolute(fileName, pathDestinationNetwork):
    """
    Function to change filePathName in Flow-FromFile xml element to absolute path after copying template network
    Input:
        fileName <string> (abs path of the copied xml fie)
        pathDestinationNetwork <string> (abs path of the directory of the xml file)
    Output:
        extData <dict> 
    """
    
    parser = etree.XMLParser(encoding='iso-8859-1')
    tree = etree.parse(fileName, parser)
    root = tree.getroot()
    
    for XmlElement in root:
        
        if XmlElement.tag == 'boundaryConditions':
            for vessel in XmlElement:
                for bc in vessel:
                    if bc.tag == 'Flow-FromFile':
                        for bcTag in bc:
                            if bcTag.tag == 'filePathName':
                                oldFilePath = bcTag.text
                                bcTag.text = pathDestinationNetwork + oldFilePath
    
    tree.write(fileName)
    
    

                                
                                



## Old unneeded functions replaced by functions in classVascularNetwork due to new data handling

# def parseDirectoryForSimulationCases(networkName):
#     """
#     Function to search for all simulationCases of a given network
#     Input:
#         networkName <String> (= name of the network to process) 
#     Output:
#         simulationCases <Dict> = { dataNumber1<String:SimCase1<String>, ...
#                                    dataNumberN<String:SimCaseN<String>}
#     """
#     path = ''.join([cur,'/../','NetworkFiles/',networkName,'/','SolutionData','/'])
#     simulationCases = {}
#     
#     for dirName, dirNames, fileNames in os.walk(path):
#         for fileName in fileNames:
#             if ".hdf5" in fileName and "polyChaos" not in fileName and "pure" not in fileName:
#                 dataNumber = fileName.split('.')[0].split('_SolutionData_')[-1]
#                 simulationCases[dataNumber] = ''.join([path,fileName])
#                 
#     return simulationCases
