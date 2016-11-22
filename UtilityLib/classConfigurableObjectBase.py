import os,sys

from UtilityLib.constants import unitsDictSI as unitsDict
from UtilityLib.moduleHelperFunctions import getGitHash

try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree
        
class ConfigurableObjectBase(object):
     
    externVariables      = {}
    externXmlAttributes  = []
    externXmlElements    = []
    
    class ExtValue(object):
        '''
        Class to describe and external "value"-variable
        
        Args:
            variableType (type, list of types): type of the external value-variable, accepted types = [bool,int,float,str,None]
            unit (str): str of the SI-unit to which the variable value is automatically converted to 
            strCases (list): necessary if variableType is str, list of strings with all possible cases the variable str can be, if anything is allowed, this need to be set to ['anything'] 
            multiVar (bool): if True: variable is treated as a list of variables of the given type, values can be both separated wiht ' ', and ',' in the XML-file
            
        Raises:
            ValueError: raise error if an unappropriate variableType is passed
            ValueError: raise error if no strCases is passed, when variableType includes str
        '''
        def __init__(self, variableType, unit = None, strCases = None, multiVar = False, optional = False):
            
            self.optional = optional 
            if type(variableType) is not list:
                variableType = [variableType]
                
                #raise ValueError("ERROR: extValue in {}, variableType is not a list:  <<{}>>.".format(self.__class__.__name__,variableType))
                  
            for varType in variableType:
                if varType not in [bool,int,float,str,None]:
                    raise ValueError("ERROR: extValue in {} has non supported type <<{}>>.".format(self.__class__.__name__,varType))
            
            self.variableType = variableType
            self.unit     = unit
            self.multiVar = multiVar
            
            if str in self.variableType:
                if strCases is not None:
                    self.strCases = strCases
                else: 
                    raise ValueError("ERROR: extValue of type <str> in {} has no <strCases> -list defined.".format(self.__class__.__name__,variableType))
                
    class ExtDict(object):
        '''
        Class to describe an external dictionary-variable
        
        Args:
            dictObjName (str): name of the external dictionary value-element in the xml tag
            dictObjType (ExtValue, ExtObject): external variable type of the dictionary value-element
            
        Raises:
            ValueError: (TODO: NOT IMPLEMENTED) raise error if an unappropriate dictObjType is passed
        '''
        def __init__(self, dictObjName, dictObjType, optional = False):
            
            self.optional = optional 
            self.variableType = ['dict']
            self.dictObjName = dictObjName
            
            self.dictObjType = dictObjType
            
            #TODO: approptiate testing for dictObjType  
#             print type(dictObjType)
#             if type(dictObjType) is type(object):
#             else: raise ValueError("ERROR: ExtDict in {}, dictObjType is not a class instance: <<{}>>".format(self.__class__.__name__,dictObjType))
#             
    class ExtObject(object):
        '''
        Class to describe an external object-variable
        
        Args:
            classCases (dict): dictionary containing all possible object types including the variable can be with the following format
            {'xmlTagOfObject': objectClass}
            
        Raises:
            ValueError: raise error if an classCases is not of dype dict
        '''
        def __init__(self, classCases, optional = False):
            
            self.optional = optional 
            self.variableType = ['object']
            
            if type(classCases) is dict:            
                self.classCases = classCases
            else: raise ValueError("ERROR: ExtObject in {}, classCases is not a dictionary: <<{}>>".format(self.__class__.__name__, classCases))
        
                    
    
    def writeXMLFile(self, filePathName):
        '''
        Function to write a xml file as defined within the class structures with
        externVariables, externXmlAttributes, externXmlElements
        
        Args:
             filePathName (str): path and filename to xml file
        '''
        
        fileName = filePathName.split('/')[-1].split('.')[0]
        root = etree.Element(fileName, gitHash = getGitHash())
        xmlFile = etree.ElementTree(root)
        xmlNodeSelf = etree.SubElement(root, self.__class__.__name__, {'class':self.__class__.__name__})  
        self.writeDataToXmlNode(xmlNodeSelf)
        xmlFile.write(filePathName,encoding='iso-8859-1',pretty_print = True)
       
    def writeDataToXmlNode(self,xmlNode):
        '''
        Function which writes the data defined with 'self.externXmlAttributes' and 'self.externXmlElements' to and         xmlNode.
        
        Args:
            xmlNode (etree.element): XML node element from etree parser.
            
        Raises:
            KeyError: If externXmlElement in 'self.externXmlElement' is not appropriate defined in 'self.externalVariables'
        '''        
        # 1. write attributes
        for attribute in self.externXmlAttributes:
            xmlNode.set(attribute, str(self.getVariable(attribute)))
            
        # 2. loop through data to save
        for externXmlElement in self.externXmlElements: 
            # check if lists are proper defined
            if externXmlElement in self.externVariables:
                # get the extern variable definition class
                externVariable  = self.externVariables[externXmlElement]
                # create xml node
                externXmlNode  = etree.SubElement(xmlNode, externXmlElement)  
                
                variableValue = self.getVariable(externXmlElement)
                if variableValue != None:
                    ## find out what type the externXmlElement variable is:
                    # 2.1 if dict variable -> writeExtValueXml
                    if isinstance(externVariable,self.ExtDict):
                        self.writeExtDictXml(externXmlNode,  variableValue, externVariable)
                    # 2.2 if object variable -> writeExtDictXml       
                    elif isinstance(externVariable,self.ExtObject):
                        self.writeExtObjectXml(externXmlNode, variableValue, externVariable)
                    # 2.3 if value variable -> writeExtObjectXml
                    elif isinstance(externVariable,self.ExtValue):
                        self.writeExtValueXml(externXmlNode, variableValue, externVariable)
                else:
                    if externVariable.optional == False:
                        raise ValueError("""ERROR: try to write <<{}>> to xml-node {},
               however the value of <<{}>> is None""".format(externXmlElement, xmlNode, externXmlElement))
                    
            else: raise KeyError("""ERROR: try to write <<{}>> to xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(externXmlElement, xmlNode, externXmlElement))
            
    def writeExtDictXml(self,externXmlNode, dictToWrite, externVariable):
        '''
        Writes dictionary xml node and envokes the dictionary data to be written
        
        Note:
                    
        Args:
            externXmlNode (etree.element)  : node in the xml file
            externXmlElement (dict)        : dict to be written
            externVariable (ExtDict)       : instance of ExtDict defining the variable properties
        '''
        # iterate through dictionary
        for key,value in dictToWrite.iteritems():
            # find out the name of the variable and create a xml node
            externXmlDictNode  = etree.SubElement(externXmlNode, externVariable.dictObjName) 
            # write id
            externXmlDictNode.set('Id', str(key))
            # if object variable -> writeExtDictXml       
            if isinstance(externVariable.dictObjType,self.ExtObject):
                self.writeExtObjectXml(externXmlDictNode, value, externVariable.dictObjType)
            # if value variable -> writeExtObjectXml
            elif isinstance(externVariable.dictObjType,self.ExtValue):
                self.writeExtValueXml(externXmlDictNode, value, externVariable.dictObjType)
        
    def writeExtObjectXml(self,externXmlNode, classToWrite, externVariable):
        '''
        Writes object xml node and envokes the object to write it's own data by calling object.writeDataToXmlNode(externXmlNode)
        
        Note:
                    
        Args:
            externXmlNode (etree.element)  : node in the xml file
            classToWrite (object)          : class instance to be written
            externVariable (ExtObject)     : instance of ExtObject defining the variable properties
        '''
        classNameDefined = False
        for definedClassName, definedClasses in externVariable.classCases.iteritems():
            if isinstance(classToWrite,definedClasses):
                externXmlNode.set('class',definedClassName)
                classNameDefined = True
        if classNameDefined == False:
            className = classToWrite.__class__.__name__
            raise KeyError("""ERROR: try to write class <<{}>>,
       however <<{}>> is not defined in externVariable.classCases: {}""".format(classToWrite, className,externVariable.classCases ))
        classToWrite.writeDataToXmlNode(externXmlNode)
        
    def writeExtValueXml(self,externXmlNode, variableValues, externVariable):
        '''
        Writes value xml node including units if existing object.writeDataToXmlNode(externXmlNode)
        
        Args:
            externXmlNode  : node in the xml file
            variableValues : str of the variable values
            externVariable : instance of ExtObject defining the variable properties
        '''
        # set unit if existing
        if externVariable.unit != None: externXmlNode.set('unit', externVariable.unit)
        # write variable value
        if externVariable.multiVar == True:
            externXmlNode.text = ' '.join(str(i) for i in variableValues)
        else:
            externXmlNode.text = str(variableValues)
    
    def loadXMLFile(self,filePathName):
        '''
        Function to load a xml file as defined within the class structures with externVariables, externXmlAttributes, externXmlElements
        
        Args:
             filePathName (str): path and filename to xml file
        
        Raises: ValueError if etree.parserError occures including line number of the xml error
        '''
        parser = etree.XMLParser(encoding='iso-8859-1')
        tree = etree.parse(filePathName, parser)
        try:
            parser = etree.XMLParser(encoding='iso-8859-1')
            tree = etree.parse(filePathName, parser)
        except (etree.ParseError, ImportError) as e:
            if isinstance(e, etree.ParseError):
                raise ValueError(" Error in XML file on line: {}".e.position[0])     
        root = tree.getroot()
        
        # try to get the xmlNodeSelf corresponding to the node of the file
        for childNode in root.iterchildren():
            if 'class' in childNode.attrib:
                if childNode.attrib['class'] == self.__class__.__name__:
                    self.loadDataFromXmlNode(childNode)
                    return
        
        raise ValueError("No xml root-child tag for class '{}' defined in xml file, could not load xmlFile {}".format(self.__class__.__name__,filePathName))
        
    def loadDataFromXmlNode(self, xmlNode):
        '''
        Function which loads the data defined with 'self.externXmlAttributes' and 'self.externXmlElements' from a passed xmlNode.
        
        Args:
            xmlNode (etree.element): XML node element from etree parser.
            
        Raises:
            KeyError: If xml node attribute is not defined in 'self.externVariables'
            KeyError: If xml node attribute type is not defiened as self.ExtValue in 'self.externVariables'
            ValueError: If xml node does not have a corresponding xml element as defined in 'self.externXmlElements'
            KeyError: If externXmlElement in 'self.externXmlElements' is not appropriate defined in 'self.externalVariables'
        '''
        newData = {}
        # load object attributes
        for attribute in self.externXmlAttributes:            
            # check if lists are proper defined
            if attribute in self.externVariables:
                if isinstance(self.externVariables[attribute], self.ExtValue):
                    newData[attribute] = self.loadVariableConversion(xmlNode.attrib[attribute], '', self.externVariables[attribute])
                else:raise KeyError("""ERROR: try to read attribute <<{}>> of xml-node {},
                 however <<{}>> is not defined as class-instance <self.ExtValue>  in self.externVariables""".format(attribute, xmlNode, attribute))   
            else: raise KeyError("""ERROR: try to read attribute <<{}>> of xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(attribute, xmlNode, attribute))
            
        ## TODO: check if there is not-needed information in the xml file
        for externXmlElement in self.externXmlElements: 
            # check if lists are proper defined
            if externXmlElement in self.externVariables:
                # get the extern variable definition class
                externVariable  = self.externVariables[externXmlElement]
                                
                # try to get the corresponding xml element:
                try: externXmlNode = xmlNode.findall(''.join(['.//',externXmlElement]))[0]
                except IndexError:
                    if externVariable.optional == False:
                        if externXmlElement in self.externVariables: variableType = self.externVariables[externXmlElement].variableType
                        else: variableType = "No entry defined for This element"
                        raise ValueError("""ERROR loadNetworkFromXML():
                                      variable "{}" of {} is not defined.
                                      (Hint:{}) , system exit!""".format(externXmlElement, xmlNode, variableType))
                    else: 
                        # if the element is optional then continue without reading it!
                        externXmlNode = None
                    
                if externXmlNode != None:
                    ## find out what type the externXmlElement variable is:
                    if isinstance(externVariable,self.ExtDict): 
                        newData[externXmlElement] = self.loadExtDictXml(externXmlNode, externXmlElement, externVariable)
                            
                    elif isinstance(externVariable,self.ExtObject):
                        newData[externXmlElement] = self.loadExtObjectXml(externXmlNode, externXmlElement, externVariable)
                        
                    elif isinstance(externVariable,self.ExtValue):
                        newData[externXmlElement] = self.loadExtValueXml(externXmlNode, externXmlElement, externVariable)
                    
            else: raise KeyError("""ERROR: try to read <<{}>> of xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(externXmlElement, xmlNode, externXmlElement))
            
        self.setVariablesDict(newData)
     
    
    def loadExtDictXml(self,externXmlNode, externXmlElement, externVariable):
        '''
        Loads dictionary xml node and envokes the dictionary data to be loaded (either object or value)
        
        Args:
            externXmlNode (etree.element): node of the xml file with dict definition
            externXmlElement (str) : str of the variable name 
            externVariable := instance of ExtDict defining the variable properties
        '''
        elementDictData= {}
                
        for dictXmlNode in externXmlNode.getchildren():
        
            dictXmlElement = dictXmlNode.tag
            
            if dictXmlElement == externVariable.dictObjName:
                dictVariable = externVariable.dictObjType
                # check if Id == dict.key is defined and unique
                if 'Id' in dictXmlNode.attrib:
                    dictXmlNodeId = dictXmlNode.attrib['Id']
                    if dictXmlNodeId  not in elementDictData:
                        
                        # check if ExtObject is expected or if ExtValue is expected
                        if isinstance(dictVariable,self.ExtObject):
                            elementDictData[dictXmlNodeId] = self.loadExtObjectXml(dictXmlNode,dictXmlElement,dictVariable)
                            
                        elif isinstance(dictVariable,self.ExtValue):
                            elementDictData[dictXmlNodeId] = self.loadExtValueXml(dictXmlNode,dictXmlElement,dictVariable)
                            
                        else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
                   however <{}> is not defined as class-instance <self.ExtValue> or <self.ExtObject> in self.externVariables""".format(dictXmlElement, externXmlElement, dictXmlElement))
                        
                    else: raise KeyError("""ERROR try to read <<{}>> of xml-node {}, however dict-key Id=<<{}>> is defined multiple times""".format(dictXmlElement, externXmlElement,dictXmlNodeId))
                    
                else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
                   however attribute <<Id>> is not defined in the XML-tag""".format(dictXmlElement, externXmlElement))
                
            else: print """WARNING: try to read xml-node <<{}>> as dict element for <<{}>>,
        however this sub-type is not defined as variable-type of <<{}>>. Skipping xml-node""".format(dictXmlElement,externXmlElement,externXmlElement)
                
                
        return elementDictData
            
            
    def loadExtObjectXml(self,externXmlNode, externXmlElement, externVariable):
        '''
        Loads class xml node and envokes the class to read its data (either object or value)
        
        Args:
            externXmlNode (etree.element): node of the xml file with class definition
            externXmlElement (str) : str of the variable name 
            externVariable := instance of ExtDict defining the variable properties
        '''
        # check if class in attribute
        if 'class' in externXmlNode.attrib:
            classXmlNode = externXmlNode.attrib['class']
            if classXmlNode in externVariable.classCases:
                # create extObject class with constructor
                extObject = externVariable.classCases[classXmlNode]()
                # update class with external data
                extObject.loadDataFromXmlNode(externXmlNode)
                # return extObject 
                return extObject
            else:
                raise ValueError("""ERROR try to read <<{}>> of xml-node {},
               however class-tpye <<{}>> is not defined ExtObject.classCases""".format(externXmlElement, externXmlNode, classXmlNode))
        
        else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
               however attribute <<class>> is not defined in the XML-tag""".format(externXmlElement, externXmlNode))
            
        
    def loadExtValueXml(self,externXmlNode, externXmlElement, externVariable):
        '''
        Loads value xml node. Checks the externXmlNode.text string and evaluates it corresponding to the properties in externVariable
        
        Args:
            externXmlNode (etree.element): node of the xml file with class definition
            externXmlElement (str) : str of the variable name 
            externVariable := instance of ExtDict defining the variable properties
        '''
        ## retrieve data from xml load
        # get variable value str from xml node
        variableValueStr = externXmlNode.text
        # get unit
        if 'unit' in externXmlNode.attrib: 
            variableUnit = externXmlNode.attrib['unit']
        else: variableUnit = None 
        
        return self.loadVariableConversion(variableValueStr, variableUnit, externVariable)
        
    def loadVariableConversion(self, variableValueStr, variableUnit, externVariable):
        '''
        Converts a value-str of a variable into the type defined by externVariable and the converts into the unit defined by externVariable from given unit (variableUnit) 
        
        Args:
            variableValueStr (str): value-str of the variable to be converted
            variableUnit (str) : str of the variable unit passed in 'variableValueStr'
            externVariable := instance of ExtDict defining the variable properties
        
        Returns:
            converted variable (str,bool,int,float,None,list): returns the converted variable
        '''
        # start conversion process
        multiVariable = False
        variableValue = 'notConvertable'
        convertError = []
        variableTypes = externVariable.variableType
        
        # check if variable is a multiple variable (means list of variables)
        if externVariable.multiVar:
            variableValueStrings = (variableValueStr.replace(',',' ')).split() 
            multiVariable = True
            variableValues = []
        else:
            variableValueStrings = [variableValueStr]
                
        # start conversion loop over variable and types
        for variableValueString in variableValueStrings:
            for variableType in variableTypes:
                
                if variableType in [float,int]:
                    try: variableValue = float(eval(variableValueString))
                    except ValueError: convertError.append('float') 
                    
                    if externVariable.unit != None and variableUnit != None:
                        if ' ' in variableUnit:
                            variableUnits = variableUnit.split(' ')
                            for variableUnit in variableUnits: variableValue = variableValue*unitsDict[variableUnit]
                        else: variableValue = variableValue*unitsDict[variableUnit]
                        
                    if variableType == int: 
                        try: variableValue = int(variableValue)
                        except ValueError: convertError.append('int') 
                        
                elif variableType == bool: 
                    if variableValueString == 'False': variableValue = False
                    elif variableValueString == 'True': variableValue = True
                    else: convertError.append('bool') 
                    
                elif variableType == str:
                    if variableValueString in externVariable.strCases: variableValue = variableValueString
                    elif externVariable.strCases[0] == 'anything': variableValue = variableValueString
                    else: convertError.append(''.join(['str == ',str(externVariable.strCases)]))
                
                elif variableType == None:
                    if variableValueString == 'None' or variableValueString == '' or variableValueString == None:
                        variableValue = None
                    else: convertError.append('None')
                        
            if variableValue == 'notConvertable':
                raise ValueError("""ERROR: {}.loadExtValueXml():
                      Cannot convert given value "{}" of variable "{}"
                      to {}!
                      Check if it is of type {}, system exit!""".format(self.__class__.__name__,
                                                                        variableValueString,
                                                                        variableValueStr,
                                                                        convertError,
                                                                        variableTypes))
            if multiVariable == False: return variableValue
            else: variableValues.append(variableValue)
        
        return variableValues   
    
        
    def setVariablesDict(self, dataDict):
        '''
        updates the class data using a dictionary in from of dataDict = {'variableName': value}
        
        Args:
            dataDict (dict): dictionary with new data to update class in from of dataDict = {'variableName': value}
        '''
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except KeyError:
                print "WARNING {}.updateData - Wrong key: {}, could not update varibale".format(self.__class__.__name__, key)
                
    def getVariable(self,variableName):
        '''
        getter function for any variable in self
                
        Args:
             variableName (str): str of the variable name self.variableName
        Returns:
            value of self.variableName 
        Raises:
            prints warning/Error if not such variable exist
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            # TODO: exchange with appropriate warning exception
            print "ERROR Vessel.getVariable() : vessel has no variable {}".format(variableName)