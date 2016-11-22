import csv
from numpy import sqrt,pi

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))

#sys.path.append(''.join([cur,'/../NetworkLib/']))
import NetworkLib.classBoundaryConditions as ccBC

from constants import newestNetworkXml as nxml
from constants import variablesDict

import moduleXML as mXML
import moduleFilePathHandler as mFPH


def writeVesselDataToCSV(networkName, vessels, delimiter=';'):
    """
    Functions writes vessel data to *.csv file 
    
    input:
        networkName <string>
        vessels     <dict>    := vessels dict of class vascularNetwork {vesselId : vesselInstance} 
        delimiter   <string>  (default = ';')
    
    """
    tags = []
    for tag in nxml.vesselAttributes: tags.append(tag)
    for vesselElement in nxml.vesselElements:
        if vesselElement == 'compliance':
            for specificCompElements in nxml.vesselElementReference[vesselElement].values():
                for tag in specificCompElements:
                    if tag not in tags:
                        tags.append(tag)
        else: 
            for tag in nxml.vesselElementReference[vesselElement]:
                tags.append(tag)
                
    ## openFile and create writer
    vesselCSVFile = mFPH.getFilePath('vesselCSVFile', networkName, 'xxx', 'write')
    writer = ccBC.csv.DictWriter(open(vesselCSVFile,'wb'),tags,delimiter=delimiter)
    
    # write first row == tags
    firstRow = {}
    for item in tags: firstRow[item] = item
    writer.writerow(firstRow)
    
    # write unit row
    unitRow = {}
    for tag in tags:
        try:
            unitRow[tag] = ''.join(['#',variablesDict[tag]['unitSI']])
        except: unitRow[tag] = ''
    unitRow['Id'] = 'unit'
    
    # write all data
    writer.writerow(unitRow)
    data = [] 
    for vessel in vessels.itervalues():
        vesselDict = {}
        for tag in tags:
            vesselDict[tag] = vessel.getVariableValue(tag)
        data.append(vesselDict)
    writer.writerows(data)

def readVesselDataFromCSV(networkName, delimiter=';'):
    """
    Functions loads vessel data from *.csv file inclusive polynomial chaos definitions
    
    input:
        networkName 
        delimiter   (default = ';')
    
    return:
            dict := {'vesselData': vesselData} which is used by vascular network to
                    update its vessel data with the function vascularNetwork.updateNetwork(dataDict)
    """
        
    vesselCSVFile = mFPH.getFilePath('vesselCSVFile', networkName, 'xxx', 'read', exception = 'Warning')
    
    if vesselCSVFile == None:
        return None
    
    
    with open(vesselCSVFile,'rb') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=";,")
        csvfile.seek(0)
        
        if dialect.delimiter == ',':
            print """\n WARNING: mCSV92 detected delimiter ',': This delimiter might lead to wrong data as
             it is used as well as decimal place separator in some languages! Check the loaded values carefully!"""
            
        #reader = ccBC.csv.DictReader(csvfile,dialect)
        reader = csv.DictReader(csvfile,delimiter = dialect.delimiter)
        # hash data with in dictionary and separate units
        columUnits = {}
        vesselData = {}
        
        for row in reader:
            Id = row.pop('Id')
            if Id == 'unit': columUnits = row
            else:
                Id = int(Id)
                vesselData[Id] = row
            
        variablesToDiscard = []
        for Id,data in vesselData.iteritems():
    #         polyChaos = {}
            for variable,variableValueStr in data.iteritems():
                # check if value is defined
                if variableValueStr not in ['', None]:
                    #find units 
                    if '#' in columUnits[variable]: #  '#' in variable or 
                        nothing,variableUnit = columUnits[variable].split('#',1)
                    # convert variables to corret unit and type
                    data[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit)
                else: variablesToDiscard.append([Id,variable]) # find out variables which have no values
        # remove variables which have no values 
        for Id,variableToDiscard in variablesToDiscard:
            del vesselData[Id][variableToDiscard]
        
        return {'vesselData':vesselData}    

def writeRandomInputstoCSV(networkName, randomInputManager, delimiter = ';'):
    '''
    Function writes random variable data to *.csv  file
    '''
    # TODO: add units
    # evaluate tags
    variablesToSaveTags = []
    variablesToSaveTags.extend(nxml.generalRandomInputsAttributes)
    for listValue in nxml.randomInputsReference.itervalues():
        variablesToSaveTags.extend(listValue)
    variablesToSaveTags.extend(['randomInputType','location'])
    # create writer 
    randomVariableCSVFile = mFPH.getFilePath('randomVariableCSVFile', networkName,'xxx', 'write', exception = 'Warning')
    writer = csv.DictWriter(open(randomVariableCSVFile,'wb'),variablesToSaveTags, delimiter=delimiter)
    
    firstRow = {key: key for key in variablesToSaveTags}
    writer.writerow(firstRow)
    
    for randomInput in randomInputManager.randomInputsList:
        rowDict = {}
        for variablesToSaveTag in variablesToSaveTags:
            rowDict[variablesToSaveTag] = randomInput.getVariableValue(variablesToSaveTag)
        writer.writerow(rowDict)
        
def readRandomInputsfromCSV(networkName, randomInputManager, delimiter = ';'):
    '''
    Function reads the random variable data from a *.csv file
    all existing rv definitions in the randomInputManager will be erased
    '''
    
    # prepare randomInputManager
    randomInputManager.deleteAllRandomInputs()
    
    # create reader
    randomVariableCSVFile = mFPH.getFilePath('randomVariableCSVFile', networkName,'xxx', 'read', exception = 'Warning')
    reader = csv.DictReader(open(randomVariableCSVFile,'rb'), delimiter=delimiter)
        
    # evaluate tags
    variablesToLoadTags = []
    variablesToLoadTags.extend(nxml.generalRandomInputsAttributes)
    for listValue in nxml.randomInputsReference.itervalues():
        variablesToLoadTags.extend(listValue)
    variablesToLoadTags.extend(['randomInputType','location'])
    # add random variables
    for row in reader:
        dataDict = {}
        for variable in variablesToLoadTags: 
            if variable in row.keys():
                # TODO: add units
                variableUnit = None 
                # save converted CSV-value
                dataDict[variable] = mXML.loadVariablesConversion(variable, row[variable], variableUnit)
            else: print "WARNING: mCSV.readRandomInputsfromCSV(), no variable {} defined in the csv file but needed to create proper working random input".format(variable)
        randomInputManager.addRandomInput(dataDict)

def writeBCToCSV(networkName, boundaryConditionDict, boundaryConditionPolyChaos, delimiter=';'):
    """
    Functions writes boundaryCondition data to *.csv file 
    
    input:
        networkName <string>
        boundaryConditionDict      <dict>  := boundaryConditionDict dict of class VascularNetwork {vesselId : boundaryConditionInstance} 
        boundaryConditionPolyChaos <dict>  := boundaryConditionPolyChaos dict of class VascularNetwork 
        delimiter   <string>  (default = ';')
    
    """
    # find all polychaos tags
    # TODO: read write polynomial chaos variables
#     polyChaosTags = {}
#     for id,bcPolyChaosList in boundaryConditionPolyChaos.iteritems():
#         for bcPolyChaosDict in bcPolyChaosList:
#             for variable,interval in bcPolyChaosDict.iteritems():
#                 if variable != 'name': polyChaosTags[variable] = len(interval)     
                
    # find all tags which are known for all boundary conditions
    tagsBCType1 = ['Id','boundaryType']
    tagsBCType2 = [] 
    for boundaryCondition,elementTags in nxml.boundaryConditionElements.iteritems():
        # None - bc for tag evaluation
        if 'None' not in boundaryCondition and '_' not in boundaryCondition:
            #check out type of BoundaryCondition:
            bcType = eval(nxml.bcTagsClassReferences[boundaryCondition])().getVariableValue('type')
            for elementTag in elementTags:
                if bcType == 1:
                    if elementTag not in tagsBCType1:
                        tagsBCType1.append(elementTag)
#                         if elementTag in polyChaosTags.keys():
#                             for count in range(polyChaosTags[elementTag]):
#                                 tagsBCType1.append(''.join([elementTag,'-pC',str(int(count)+1)]))
                elif bcType == 2:
                    if elementTag not in tagsBCType2:
                        tagsBCType2.append(elementTag)
#                         if elementTag in polyChaosTags.keys():
#                             for count in range(polyChaosTags[elementTag]):
#                                 tagsBCType2.append(''.join([elementTag,'-pC',str(int(count)+1)]))
        
    tagsBCType1.extend(tagsBCType2)
    tags = tagsBCType1
    
    boundaryCSVFile = mFPH.getFilePath('boundaryCSVFile', networkName, 'xxx', 'write')
    writer = ccBC.csv.DictWriter(open(boundaryCSVFile,'wb'),tags,delimiter=delimiter)
    
    # write Tag row
    firstRow = {}
    for item in tags:
        firstRow[item] = item
    writer.writerow(firstRow)
    
    # write unit row
    unitRow = {}
    for tag in tags:
        try:
            if '-pC' in tag: 
                tagUnit = tag.split('-pC')[0]
                unitRow[tag] = ''.join(['#',variablesDict[tagUnit]['unitSI']])
            else:
                unitRow[tag] = ''.join(['#',variablesDict[tag]['unitSI']])
        except: unitRow[tag] = ''
    unitRow['Id'] = 'unit'
    writer.writerow(unitRow)
    
    ## fill out data of defined boundary conditions
    for Id,boundaryConditions in boundaryConditionDict.iteritems():
        for boundaryCondition in boundaryConditions:
            boundaryType  = boundaryCondition.getVariableValue('name')
            dataRow       = {}
            dataRow['Id'] = Id
            dataRow['boundaryType'] = boundaryType
            for variable in nxml.boundaryConditionElements[boundaryType]:
                dataRow[variable] = boundaryCondition.getVariableValue(variable)
#             try:
#                 for bcPolyChaosDict in boundaryConditionPolyChaos[Id]:
#                     bcPolyChaosDict.pop('name')
#                     for variable,interval in bcPolyChaosDict.iteritems():
#                         for count,value in enumerate(interval): 
#                             dataRow[''.join([variable,'-pC',str(int(count)+1)])] = value
#             except: pass               
                    
            writer.writerow(dataRow) 
    
    
def readBCFromCSV(networkName, delimiter=';'):
    """
    Functions loads boundaryCondition data from \*.csv file inclusive polynomial chaos definitions

    Args:
        networkName (str): The name of the network
        delimiter (str): Delimiter (default = ';')

    Returns:
        VascularNetwork.boundaryConditionDict
            A description of the VascNw.BCD instance returned

        VascularNetwork.boundaryConditionPolyChaos
            A description of the VascNw.BCPC instance returned

    """
    
    boundaryCSVFile = mFPH.getFilePath('boundaryCSVFile', networkName, 'xxx', 'read', exception = 'Warning')
    reader = ccBC.csv.DictReader(open(boundaryCSVFile,'rb'),delimiter=delimiter)
    
    # hash data with in dictionary and separate units
    columUnits = {}
    boundaryData = []
    for row in reader:
        if row['Id'] == 'unit': columUnits = row
        else:  boundaryData.append(row)
    
    boundaryConditionPolyChaos = {}       
    BCconditionData = {}
    for bcData in boundaryData:
        Id = int(bcData['Id'])
        # TODO: read write polynomial chaos variables
        # create class instance
        boundaryType = bcData['boundaryType']
        try: boundaryInstance = eval(nxml.bcTagsClassReferences[boundaryType])()
        except: 'ERROR moduleCSV.readBCFromCSV: boundaryType <<{}>> does not exist'.format(boundaryType)
        boundaryDataDict = {'name':boundaryType}
        polyChaos = {}
        for variable,variableValueStr in bcData.iteritems():
            if variable in nxml.boundaryConditionElements[boundaryType]: 
                try: variableUnit = columUnits[variable]
                except: variableUnit = None 
                # save converted XML-value
                if variable == 'filePathName':
                    path = ''.join([boundaryCSVFile])
                    if path not in variableValueStr:  variableValueStr = variableValueStr.join([path,''])
                print variableValueStr
                if variableValueStr != '':
                    try:
                        boundaryDataDict[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit)
                    except:
                        pass
            if '-pC' in variable and variableValueStr != '':
                polyChaos['name'] = boundaryType
                variable,number = variable.split('-pC')
                if variable in polyChaos.keys():
                    polyChaos[variable] = ' '.join([polyChaos[variable],variableValueStr])
                else: polyChaos[variable] = variableValueStr
        # convert polynomial chaos variables to corret unit and type
        if polyChaos != {}:
            for variable,variableValueStr in polyChaos.iteritems():
                try:
                    variableUnit = columUnits[variable].split('#',1)
                    polyChaos[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                except: pass          
            if Id not in boundaryConditionPolyChaos.keys(): 
                boundaryConditionPolyChaos[Id] =[polyChaos]
            else: boundaryConditionPolyChaos[Id].append(polyChaos)           
            
        boundaryInstance.update(boundaryDataDict)
        
        if Id not in BCconditionData.keys(): BCconditionData[Id] = [boundaryInstance]
        else: BCconditionData[Id].append(boundaryInstance)
    
    return BCconditionData, boundaryConditionPolyChaos
  
