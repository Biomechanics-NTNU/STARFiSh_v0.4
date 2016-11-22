
# import dependencies

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )

# functions for vascular1DFlow_v0.2

from NetworkLib.classVascularNetwork import VascularNetwork 
from NetworkLib.classBoundaryConditions import *

import UtilityLib.moduleXML as mXML 
import UtilityLib.moduleCSV as mCSV
import UtilityLib.moduleFilePathHandler as mFPH

from UtilityLib.constants import newestNetworkXml as nxml

### import units of all variables in the SI system
from UtilityLib.constants import variableUnitsSI as variableUnits
from UtilityLib.constants import unitsDictSI as unitsDict

#from modulePickle import loadSolutionDataFile
import UtilityLib.moduleStartUp as mStartUp

### import units of all variales in the Medical System
#from constants import variableUnitsMed as variableUnits
#from constants import unitsDictMed as unitsDict

import pydot
from VncLib import xdot
import gtk

from VncLib.classGraph import Graph
from VncLib.classGraph import MyDotWindow

import cPickle
import pprint as pprint
import numpy as np
import thread
from copy import deepcopy


def enterNetworkName(networkName, recentNetworkNames = None):
    """
    function to evaluate networkName input of user
    
    input:
          networkName <string>      := current existing networkName
          recentNetworkNames <list> := list of recently used networkNames
    
    return:
         networkNameUserInput <string> := if user input == None : networkName passed to function
                                          if user input == number correspoding to recentNetworkNames: 
                                                           the correspinding networkName
                                          if user input == newNetworkName: newNetworkName
    """
    print "     current networkName: ",networkName
    networkNameUserInput = str(raw_input("     enter/change networkName (only! ENTER to use current networkName):\n     "))
    
    if networkNameUserInput == "":
        if networkName != None:
            networkNameUserInput = networkName
        else:
            networkNameUserInput = None
            
    if networkNameUserInput != None:       
        if networkNameUserInput in [str(int(i)) for i in range(0,9)] and len(networkNameUserInput)==1 and recentNetworkNames!= None:
            number = int(networkNameUserInput)
            if len(recentNetworkNames)>= number:
                networkNameUserInput = recentNetworkNames[number-1]
                        
    return networkNameUserInput
    

def findAllDaughters(vascularNetwork, motherVesselID):
    """
    evaluates all daughters of the vessel with the passed motherVesselID
    
    input:
        vascularNetwork <classVascularNetwork Instance>
        motherVesselID <int>
        
    return daughters <list> := list with all daughter ids
    """
    daughters = []
    viz = []
    root = motherVesselID
    if vascularNetwork.vessels[root].leftDaughter != None:
        viz.append(root)
    # loop through tree until all daughters are added to the graph
    while len(viz) != 0:
        # get the mother vessel (already added) and add its daughters
        motherVessel = viz.pop(0)
        # find left daughter
        leftDaughter = vascularNetwork.vessels[motherVessel].leftDaughter
        # add left daughter
        daughters.append(leftDaughter)
        # check if leftDaughter has also daughters 
        if vascularNetwork.vessels[leftDaughter].leftDaughter != None:
            viz.append(leftDaughter)
        # find right daughter
        if vascularNetwork.vessels[motherVessel].rightDaughter != None:
            rightDaughter = vascularNetwork.vessels[motherVessel].rightDaughter
            # add right daughter
            daughters.append(rightDaughter)
            # check if rightDaughter has also daughters
            if vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                viz.append(rightDaughter)
    return daughters
    

def main():
    """
    vnc main function
    
    menu and all options in one big function
    
    need refactoring:
        transformed into class 
        splitted in single member functions
        thus give possibility to make gui
        xdot on click should call function in class vncMain
        vascularNetwork should never be deleted ?!
        modules need to be imported as modules !!
    """
    # set graphs directory
    graphPath = str(cur+'/NetworkFiles/')
    
    # create a new window
    window = MyDotWindow()    
    window.connect('destroy',gtk.main_quit)
    window.show()
    
    #create the main graph instance of a graph class
    mainGraph = Graph()
    
    #START THE MAIN LOOP
    menuInput = ""
    subMenuInput = ''
    
    # create vascularNetwork instance
    vascularNetwork = VascularNetwork()
    networkName = None
    k = None
    while menuInput != "q":
        menuInput = ""
        print ""
        print '====================================='
        print '#    VascularNetworkCreator_v2.1    #'
        print '====================================='
        print " [a] - add vessel to network"
        print " [d] - delete vessel in network"
        print " [n] - new network"
        print " [b] - set boundary conditions"
        print " [f] - set global fluid properties"
        print " [l] - load network"
        print " [s] - save network"
        print " [u] - update XML from CSV file(s)"
        print " [g] - print network graph"
        print " [p] - print network informations"
        print " [q] - quit"
        print ""
        print '===========Parameter Units==========='
        print ' NOTE! All parameter units are in the'
        print ' Pa, seconds, and meters system'
        print ' i.e. not mmHg, seconds, cm, or dynes,'
        print ' seconds, grams'
        print '====================================='
        print '  current network: ', networkName
        while  menuInput not in ("l","b","q","a","s","g","f","d","u",'n','p'):
            menuInput = raw_input("what to do? ")
        
        if menuInput == "a": 
            print "Add new vessel"
            
            existing = False
            vesselId = raw_input(" enter the vessel id:  ")
            while True:
                try:
                    vesselId = int(vesselId)
                    if vesselId not in vascularNetwork.vessels:
                        break
                    else:
                        existing = True
                except ValueError:
                    print "TYPE-ERROR: vessel id must be type(int) not type(string)"
                    vesselId = raw_input(" enter non existing id: ")
                if existing == True:
                    print " the vessel id exists already enter a new one"
                    vesselId = raw_input(" enter non existing id: ")
                    existing = False
            
            
            if vascularNetwork.vessels != {}:
                                    
                existing = False
                mother = raw_input(" enter existing mother id:  ")
                while True:
                    try:
                        mother = int(mother)
                        if mother in vascularNetwork.vessels.keys() and vascularNetwork.vessels[mother].rightDaughter == None:
                            break
                        else:
                            existing = True
                    except ValueError:
                        print "TYPE-ERROR: mother id must be type(int) not type(string)"
                        mother = raw_input(" enter existing mother id:  ")
                    if existing == True:
                        if mother not in vascularNetwork.vessels: print " there exists no vessel with this id"
                        else: print "   only bifurcations possible!"
                        mother = raw_input(" enter existing mother id:  ")
                        existing = False
                    
                if vascularNetwork.vessels[mother].leftDaughter == None: vascularNetwork.vessels[mother].leftDaughter = vesselId
                else: vascularNetwork.vessels[mother].rightDaughter = vesselId
                                
            vascularNetwork.addVessel(vesselId)
            print " define vessel compliance!"
            
            
            inputType = '0'
            print "     available compliance types:"
            print ""
            # get all defined boundaryConditions from constants-dict save as bcTypes
            complianceTypes = nxml.vesselComplianceElements.keys()
            compTypes = ['default (Hayashi)']
            for compType in complianceTypes: 
                compTypes.append(compType)
            # show all compliance types in the compTypes
            index = 0
            for key in compTypes:
                print "       [",str(index).rjust(2),"]    ",key
                index = index+1
            # get user input and check if it was correct to define the bcType 
            existing = False
            inputType = raw_input ("      choose type ")
            while True:
                # check if right input
                try:
                    inputType = int(inputType)
                    if inputType in np.linspace(0,len(compTypes)-1,len(compTypes)):
                        break
                    else:
                        existing = True
                # if not int but string
                except ValueError:
                    print "      TYPE-ERROR: vessel id must be type(int) not type(string)"
                    inputType = (raw_input ("      choose type "))
                # if int but to low or high
                if existing == True:
                    print "       the type does not exist"
                    inputType = (raw_input ("      choose type "))
                    existing = False
            
            compType = compTypes[int(inputType)] 
            
            if compType != 'default (Hayashi)':
                vesselData = {'complianceType':compType}
                nxml.vesselComplianceElements[compType]
                                    
                print ""
                print "      set values for the Compliance: ", compType
                question = True
                for arg in nxml.vesselComplianceElements[compType]:
                    if arg != 'complianceType':
                        currValue = raw_input (str("            set value for "+str(arg)+' '))
                        test = True
                        try: float(currValue)
                        except:
                            print '            VALUE or TYPE ERROR, set to None'
                            test = False
                        if test == True: vesselData[arg] = (float(currValue))
                        else: vesselData[arg] = None
                vascularNetwork.updateNetwork({'vesselData':  { vesselId : vesselData}})
            
            mainGraph.update_graph(vascularNetwork, window)
                    
        if menuInput == "d":
            print "Delete a vessel and all its daugthers"
            if vascularNetwork.vessels.keys() != []:
                
                existing = False
                vesselId = raw_input(" enter existing vessel id: ")
                while True:
                    try:
                        vesselId = int(vesselId)
                        if vesselId in vascularNetwork.vessels:
                            break
                        else: 
                            existing = True
                    except ValueError:
                        print "TYPE-ERROR: vessel id must be type(int) not type(string)"
                        vesselId = raw_input(" enter existing vessel id: ")
                    if existing == True:
                        print " the vessel does not exist"
                        vesselId = raw_input(" enter existing vessel id: ")
                        existing = False
                
                #travers the tree starting with the vessel and collect all ids
                toDelete = findAllDaughters(vascularNetwork,vesselId)
                
                toDelete.append(vesselId)
                
                for vesselToDelete in toDelete:
                    vascularNetwork.deleteVessel(vesselToDelete)
                
                # empty the graph to redraw it
                if vascularNetwork.vessels.keys() == []:
                    mainGraph.update_graph(None, window)
                    vascularNetwork = VascularNetwork()
                else:
                    mainGraph.update_graph(vascularNetwork, window)
                
            else:
                print " there are no vessels to delete"
                
        elif menuInput == "n":
            print "new network"
            question = raw_input(" are u sure to delete all current data? [y] - yes: ")
            if question == 'y':
                # delete vascularNetwork
                del vascularNetwork
                # create vascularNetwork instance
                mainGraph.update_graph(None, window)
                vascularNetwork = VascularNetwork()
                    
        elif menuInput == "p":
            vascularNetwork.showVessels()
            vascularNetwork.showNetwork()
            vascularNetwork.randomInputManager.printOutInfo()
                        
        elif menuInput == "g":
            print mainGraph.getGraph()
            
        elif menuInput == "b":
            subMenuInput = ''
            
            while  subMenuInput not in ['1','2','3','b']:
                if vascularNetwork.getVariableValue('vessels') == {}:
                    print " there are no vessels defined and thus no boundarys available";break
                else:
                    # evaluate boundarys in Network
                    boundarys = []
                    notDefinedBoundarys = []
                    
                    boundarys.extend(vascularNetwork.boundaryVessels)
                    if vascularNetwork.root != None and vascularNetwork.root not in boundarys: 
                        boundarys.append(vascularNetwork.root)
                        
                    boundarysSaved = vascularNetwork.boundaryConditions.keys()
                    
                    # update saved boundary conditions
                    for boundarysCurrent in boundarys:
                        if boundarysCurrent not in boundarysSaved:
                            print " boundary added to vascularNetwork"
                            vascularNetwork.boundaryConditions[boundarysCurrent] = []
                            boundarysSaved.append(boundarysCurrent)
                            
                        if vascularNetwork.boundaryConditions[boundarysCurrent] == []:
                            notDefinedBoundarys.append(boundarysCurrent)
                        
                    nonBoundarys = list(set(boundarys).symmetric_difference(set(boundarysSaved)))
                    for nonBoundary in nonBoundarys:
                        print " boundary removed from vacularNetwork"
                        del(vascularNetwork.boundaryConditions[nonBoundary])
                        if nonBoundary in notDefinedBoundarys: notDefinedBoundarys.remove(nonBoundary)
                        
                    vascularNetwork.evaluateConnections()
                    print ""
                    print "    sub menu: boundary conditions"
                    print ""
                    print "     [1] - show  boundary conditions"
                    print "     [2] - add   boundary condition "
                    print "     [3] - del   boundary condition "
#                     print "     [4] - load  boundary conditions from CSV"
#                     print "     [5] - write boundary conditions to CSV"
                    print "     [b] - back to the main menu"
                    print ""     
                    subMenuInput = raw_input("     what to do? ") 
                    
                    if subMenuInput == '1':
                        print "     boundary conditions"
                        pprint.pprint(vascularNetwork.boundaryConditions)
                        subMenuInput = ''
                        
                    elif subMenuInput == '2' and vascularNetwork.root != []:
                        
                        print "     add   boundary condition"
                        print ""
                        
                        definedBoundarys = list(set(notDefinedBoundarys).symmetric_difference(set(vascularNetwork.boundaryConditions.keys())))
                        print "     vessels with defined boundary condition:"
                        print "       ",'  '.join(str(i) for i in definedBoundarys)
                        
                        print "     vessels with undefined boundary condition:"
                        print "       ",'  '.join(str(i) for i in notDefinedBoundarys)
                        print ""
                                            
                        existing = False
                        vesselId = raw_input(" enter existing vessel id: ")
                        while True:
                            try:
                                vesselId = int(vesselId)
                                if vesselId in vascularNetwork.vessels:
                                    break
                                else:
                                    existing = True
                            except ValueError:
                                print " TYPE-ERROR: vessel id must be type(int) not type(string)"
                                vesselId = raw_input(" enter existing vessel id: ")
                            if existing == True:
                                print " the vessel does not exist"
                                vesselId = raw_input(" enter existing vessel id: ")
                                existing = False
                                                
                        inputType = '0'
                        print "     add boundary condition type:"
                        print ""
                        # get all defined boundaryConditions from constants-dict save as bcTypes
                        bcTypesAll = nxml.bcTagsClassReferences.keys()
                        bcTypes = []
                        for bcType in bcTypesAll: 
                            if "_" is not bcType[0]:
                                bcTypes.append(bcType)
                        bcTypes.sort()
                        # show all boundaryConditions in the bcTypes
                        index = 0
                        for key in bcTypes:
                            print "       [",str(index).rjust(2),"]    ",key
                            index = index+1
                        # get user input and check if it was correct to define the bcType 
                        existing = False
                        inputType = raw_input ("      choose type ")
                        while True:
                            # check if right input
                            try:
                                inputType = int(inputType)
                                if inputType in np.linspace(0,len(bcTypes)-1,len(bcTypes)):
                                    break
                                else:
                                    existing = True
                            # if not int but string
                            except ValueError:
                                print "      TYPE-ERROR: vessel id must be type(int) not type(string)"
                                inputType = (raw_input ("      choose type "))
                            # if int but to low or high
                            if existing == True:
                                print "       the type does not exist"
                                inputType = (raw_input ("      choose type "))
                                existing = False
                        
                        bcType = bcTypes[int(inputType)] 
                        
                        boundaryInstance = eval(nxml.bcTagsClassReferences[bcType])()
                        boundaryDataDict = {}
                        boundaryDataDict['name']= bcType
                        
                        print ""
                        print "      set values for the BC condition: ", bcType
                        print "          enter 'b' for the first value to skip this procedure"
                        question = True
                        for arg in nxml.boundaryConditionElements[bcType]:
                            if question == True: 
                                currValue = raw_input (str("            set value for "+str(arg)+' '))
                                if currValue == 'b': question=False
                                test = True
                                try: float(currValue)
                                except:
                                    print '            VALUE or TYPE ERROR, set to None'
                                    test = False
                                if test == True: boundaryDataDict[arg] = (float(currValue))
                                else: boundaryDataDict[arg] = None
                            else: boundaryDataDict[arg] = None
                        if len(vascularNetwork.boundaryConditions.keys()) == 1:
                            print "      set position of the BC condition"
                            position = '2'
                            while position not in ['1','0']:
                                position = raw_input ("          enter '0' for the start or '1' for the end of the vessel ")
                            if position == '1':
                                #bcType = ''.join(['_',bcType])
                                boundaryDataDict['name']= ''.join(['_',bcType])
                        
                        print bcType,boundaryDataDict
                        
                        boundaryInstances = []
                        boundaryInstance.update(boundaryDataDict)
                        boundaryInstances.append(boundaryInstance)
                        
                        if vesselId not in vascularNetwork.boundaryConditions.keys():
                            vascularNetwork.boundaryConditions[vesselId] = boundaryInstances
                        else:
                            vascularNetwork.boundaryConditions[vesselId].extend(boundaryInstances)
                        
                        #if vascularNetwork.getVariableValue('root') != None:
                        mainGraph.update_graph(vascularNetwork, window)
                        subMenuInput = ''
                            
                    
                    elif subMenuInput == '3' and vascularNetwork.root != []:
                        print "     delete boundary condition"
                        print ""
                        pprint.pprint(vascularNetwork.boundaryConditions)
                        
                        vesselId = -1
                        while vesselId not in vascularNetwork.boundaryConditions.keys():
                            vesselId = int(raw_input ("      choose vessel id "))
                        
                        bcs = vascularNetwork.boundaryConditions[vesselId]
                        if bcs != []:                 
                            print ""
                            index = 0
                            for bc in bcs:
                                print "       [",str(index).rjust(2),"]    ",bc.name
                                index = index+1
                            print ""
                            inType = '0'
                            while inType not in np.linspace(0,len(bcs)-1,len(bcs)):
                                inType = int(raw_input ("      choose condition to delete "))
                             
                            print ""
                            print "     boundary condition ",bcs[inType]," removed!"
                            print ""
                            vascularNetwork.boundaryConditions[vesselId].remove(bcs[inType])
                            
                            mainGraph.update_graph(vascularNetwork, window)
                        else:
                            print "     nothing to delete!"
                        subMenuInput = ''
                    
#                     elif subMenuInput == '4' and vascularNetwork.root != []:
#                         print "     load  boundary conditions from CSV"
#                         print ""
#                         networkName = enterNetworkName(networkName)
#                         boundaryConditions,boundaryConditionPolyChaos = mCSV.readBCFromCSV(networkName)
#                         vascularNetwork.update({'boundaryConditions':boundaryConditions,
#                                                 'boundaryConditionPolyChaos':boundaryConditionPolyChaos})
#                         
#                         mainGraph.update_graph(vascularNetwork, window)
#                                                                 
#                     elif subMenuInput == '5' and vascularNetwork.root != []:
#                         print "     write boundary conditions to CSV"
#                         print ""
#                         networkName = enterNetworkName(networkName)
#                         boundaryConditions = vascularNetwork.getVariableValue('boundaryConditions')
#                         boundaryConditionPolyChaos = deepcopy(vascularNetwork.getVariableValue('boundaryConditionPolyChaos'))
#                         mCSV.writeBCToCSV(networkName, boundaryConditions, boundaryConditionPolyChaos)    
                                                
                    elif subMenuInput == 'b':
                        break
        
        elif menuInput == "f":
            
            subMenuInput = ''
            while  subMenuInput not in ["1","2","b"]:
                print ""
                print "    sub menu: set global fluid properties"
                print ""
                print "     [1] - set all"
                print "     [2] - set individual"
                print "     [b] - back to the main menu"
                print ""
                print "    current fluid properties:"
                for key,value in vascularNetwork.globalFluid.iteritems():
                    print "     {0:20}     {1:10}  {2:10}".format(key,value,variableUnits[key])
                print ""
                subMenuInput = raw_input("    what to do? ")
                
                if subMenuInput == '1':
                    print "     set all fluid properties"
                    
                    for key in vascularNetwork.globalFluid.keys():
                        inputType = "1"
                        typeFalse = False
                        while typeFalse == False:
                            try: 
                                inputType = raw_input("      type value for property: {} with unit {} : ".format(key,variableUnits[key]))
                                inputType = float(inputType)
                                typeFalse = True
                                vascularNetwork.globalFluid[key] = inputType
                            except: pass
                    subMenuInput = ''
                    
                elif subMenuInput == '2':
                    print "     set individual fluid property:"
                    i = 0
                    properties = vascularNetwork.globalFluid.keys()
                    for property_i in properties:
                        print "          [",i,'] - ',property_i
                        i = 1+i
                    inputType = 0
                    while inputType not in [str(i) for i in range(0,len(properties))]:
                        inputType = raw_input("      choose property to set ")
                    
                    inputProperty = ""   
                    inputFalse = False
                    while inputFalse == False:
                        try: 
                            inputProperty = raw_input("      type value for property: {} with unit {} : ".format(properties[int(inputType)],variableUnits[properties[int(inputType)]]))
                            inputProperty = float(inputProperty)
                            vascularNetwork.globalFluid[properties[int(inputType)]] = inputProperty
                            inputFalse = True
                        except: pass
                    subMenuInput = ''
                     
                elif subMenuInput == 'b':
                    break
                
            
        elif menuInput == "l":
            try:
                vncRescentNetworksFile = open(mFPH.getFilePath('vncRescentNetworksFile', 'networkName', 'xxx', 'read'),'rb')
                recentNetworkNames = cPickle.load(vncRescentNetworksFile)
                vncRescentNetworksFile.close()
            except:
                recentNetworkNames = []
                
            subMenuInput = ''
            
            while  subMenuInput not in ["1","2","3",'5','6',"b"]: #["1","2","3",'4','5','6',"b"]:
                                
                print ""
                print "    sub menu: load data"
                print ""
                print "     [1] - load network from XML"
                print "     [2] - load template network"
                print "     [3] - load vessel data from CSV"
                #print "     [4] - load vessel data and boundary conditions from CSV"
                print "     [4] - "
                print "     [5] - load network from SolutionData"
                #print "     [6] - load random inputs from CSV"
                print "     [b] - back to the main menu"
                print ""
                
                subMenuInput = raw_input("what to do? ")
                
                
                if subMenuInput in ["1","2","3",'5']:
                    print ""
                    print "         recently used networks"
                    i = 1
                    for name in recentNetworkNames:
                        print "          [",i,'] - ',name
                        i = 1+i
                    print ""
                
                if subMenuInput == '1':
                    print "     load from XML\n"
                    
                    networkName = enterNetworkName(networkName,recentNetworkNames = recentNetworkNames)
                    if networkName == None:break
                    
                    # delete the old network
                    del vascularNetwork
                    
                    #load the new network
                    try:
                        vascularNetwork = mXML.loadNetworkFromXML(networkName)
                    except ValueError as e:
                        mainGraph.update_graph(None, window)
                        vascularNetwork = VascularNetwork()
                        print "\n  could not load network, it does not exist or the file is not up-to-date! \n"
                        print(str(e))
                        #if networkName in recentNetworkNames:
                        #    recentNetworkNames.remove(networkName)
                        networkName = None
                        
                        #vncRescentNetworksFile = open(mFPH.getFilePath('vncRescentNetworksFile', 'networkName', 'xxx', 'write'),'wb')
                        # store pickle
                        #cPickle.dump(recentNetworkNames, vncRescentNetworksFile, protocol=2)
                        #vncRescentNetworksFile.close()
                        
                        break
                    
                    if networkName != None:
                        mainGraph.update_graph(vascularNetwork, window)
                    break
                
                elif subMenuInput == '2':
                    print "\n     load template network:\n"
                    
                    # network templates
                    templatePath = mFPH.getDirectory('networkXmlFileTemplateDirectory','','','read')
                    
                    dirNamesTemplate = [d for d in os.listdir(templatePath) if '.' not in  d]
                    for index,dirName in enumerate(dirNamesTemplate):
                        print "        [ {:3} ] - {}".format(index,dirName)
                    print ""
                    indexChoosen = None
                    while  indexChoosen not in [str(i) for i in xrange(len(dirNamesTemplate))]:
                        indexChoosen = raw_input("     Insert index of network: ")
                                                   
                    templateNetworkName = dirNamesTemplate[int(indexChoosen)]
                    #load the new network
                    
                    del vascularNetwork
                    # delete the old network
                    try:
                        vascularNetwork = mXML.loadNetworkFromXML(templateNetworkName)
                        vascularNetwork.name = None
                        networkName = None
                    except:
                        mainGraph.update_graph(None, window)
                        vascularNetwork = VascularNetwork()
                        print "\n  could not load network, it does not exist! \n"
                        if networkName in recentNetworkNames:
                            recentNetworkNames.remove(networkName)
                        networkName = None
                        vncRescentNetworksFile = open(mFPH.getFilePath('vncRescentNetworksFile', 'networkName', 'xxx', 'write'),'wb')
                        # store pickle
                        cPickle.dump(recentNetworkNames, vncRescentNetworksFile, protocol=2)
                        vncRescentNetworksFile.close()
                        break
                    mainGraph.update_graph(vascularNetwork, window)
                    break
                
                elif subMenuInput == '3':
                    print "     load vessel data from CSV - non existing vessels are added automatically"
                    print ""
                    networkName = enterNetworkName(networkName, recentNetworkNames = recentNetworkNames)
                    if networkName == None:break
                    
                    vesselData = mCSV.readVesselDataFromCSV(networkName)
                    if vesselData == None:
                        print "\n  could not load network csv, it does not exist! \n"
                    else:
                        vascularNetwork.updateNetwork(vesselData)
                        mainGraph.update_graph(vascularNetwork, window)
                    break
                
#                 elif subMenuInput == '4':
#                     print "     load vessel data and boundary conditions from CSV"
#                     networkName = enterNetworkName(networkName, recentNetworkNames = recentNetworkNames)
#                     if networkName == None:break
#                     
#                     vascularNetwork.updateNetwork(mCSV.readVesselDataFromCSV(networkName))
#                     boundaryConditions,boundaryConditionPolyChaos = mCSV.readBCFromCSV(networkName)
#                     vascularNetwork.update({'boundaryConditions':boundaryConditions,
#                                             'boundaryConditionPolyChaos':boundaryConditionPolyChaos})
#                     
#                     mainGraph.update_graph(vascularNetwork, window)
#                     break
                
                elif subMenuInput == '5':
                    print "     load network from SolutionData"
                    try:
                        networkName,dataNumber = mStartUp.chooseSolutionDataCase()
                        
                        vascularNetwork = mXML.loadNetworkFromXML(networkName,dataNumber = dataNumber)
                        vascularNetwork.name = None
                        networkName = None
                        mainGraph.update_graph(vascularNetwork, window)
                                                    
                    except: print "\n ERROR occured could not open requested network, file does not exist or is out dated"
                    break
                
                
#                 elif subMenuInput == '6':
#                     print "     load random input from CSV"
#                     networkName = enterNetworkName(networkName)
#                     mCSV.readRandomInputsfromCSV(networkName, vascularNetwork.randomInputManager)
#                     vascularNetwork.randomInputManager.linkRandomInputUpdateFunctions(vascularNetwork)
                
                elif subMenuInput == 'b':
                    break
            
            if networkName != None:    
                if networkName not in recentNetworkNames: recentNetworkNames.insert(0,networkName)
                else: 
                    recentNetworkNames.remove(networkName)
                    recentNetworkNames.insert(0,networkName)
                if len(recentNetworkNames) > 9: recentNetworkNames.pop(-1)
                            
                vncRescentNetworksFile = open(mFPH.getFilePath('vncRescentNetworksFile', 'networkName', 'xxx', 'write'),'wb')
                # store pickle
                cPickle.dump(recentNetworkNames, vncRescentNetworksFile, protocol=2)
                vncRescentNetworksFile.close()
            
        elif menuInput == "s":
            subMenuInput = ''
            print ""
            print "    sub menu: save data"
            print ""
            print "     [1] - write to XML"
            print "     [2] - write vessel data to CSV"
            #print "     [3] - write vessel data and boundary conditions to CSV"
            print "     [3] - "
            print "     [4] - write graph to .png"
            #print "     [5] - write random input data to CSV"
            print "     [b] - back to the main menu"
            print ""
            while subMenuInput not in ["1","2",'4',"b"]: #["1","2","3",'4','5',"b"]:
                subMenuInput = raw_input("what to do? ")
                     
                if subMenuInput == '1':
                    print "     write to XML"
                    networkName = enterNetworkName(networkName)
                    vascularNetwork.name = networkName
                    if networkName == None:break
                    mXML.writeNetworkToXML(vascularNetwork)
                    break
                    
                elif subMenuInput == '2':
                    print "     write vessel data to CSV"
                    networkName = enterNetworkName(networkName) 
                    if networkName == None:break
                    mCSV.writeVesselDataToCSV(networkName, vascularNetwork.vessels)
                    break
                
#                 elif subMenuInput == '3':
#                     print "     write vessel data and boundary conditions to CSV"
#                     networkName = enterNetworkName(networkName) 
#                     if networkName == None:break
#                     boundaryConditions = vascularNetwork.getVariableValue('boundaryConditions')
#                     boundaryConditionPolyChaos = deepcopy(vascularNetwork.getVariableValue('boundaryConditionPolyChaos'))
#                     # write data
#                     mCSV.writeVesselDataToCSV(networkName, vascularNetwork.vessels)
#                     mCSV.writeBCToCSV(networkName, boundaryConditions,boundaryConditionPolyChaos)  
#                     break
                
                elif subMenuInput == '4':
                    print "     write graph to .png"
                    vncNetworkGraphFile = mFPH.getFilePath('vncNetworkGraphFile', networkName, 'xxx', 'write')
                    #mainGraph.graph.write(graphPath+networkName+'/'+pictureName+'.dot')
                    mainGraph.graph.write_png(vncNetworkGraphFile)
                    break
                
#                 elif subMenuInput == '5':
#                     print "     write random input data to CSV"
#                     networkName = enterNetworkName(networkName) 
#                     mCSV.writeRandomInputstoCSV(networkName, vascularNetwork.randomInputManager)                
#                     break
                
                if subMenuInput == 'b':
                    break
        
        elif menuInput == "u":
            try:
                vncRescentNetworksFile = open(mFPH.getFilePath('vncRescentNetworksFile', 'networkName', 'xxx', 'read'),'rb')
                recentNetworkNames = cPickle.load(vncRescentNetworksFile)
                vncRescentNetworksFile.close()
            except:
                recentNetworkNames = []
                        
            subMenuInput = ''
            print ""
            print "    sub menu: update XML from CSV"
            print ""
            print "     load from XML"
            print ""
            print "         recently used networks"
            i = 1
            for name in recentNetworkNames:
                print "          [",i,'] - ',name
                i = 1+i
            print ""
            networkName = enterNetworkName(networkName, recentNetworkNames = recentNetworkNames)
            if networkName == None:break
            # delete the old network
            del vascularNetwork
            #load the new network
            try:
                vascularNetwork = mXML.loadNetworkFromXML(networkName)
            except ValueError as e:
                mainGraph.update_graph(None, window)
                vascularNetwork = VascularNetwork()
                print "\n  could not load network, it does not exist or the file is not up-to-date! \n"
                print(str(e))
                #if networkName in recentNetworkNames:
                #    recentNetworkNames.remove(networkName)
                networkName = None
            
            if networkName == None:
                mainGraph.update_graph(None, window)
            else:
                mainGraph.update_graph(vascularNetwork, window)
                
            if networkName is not None:
                print "     load vessel data from CSV - non existing vessels are added automatically"
                vesselData = mCSV.readVesselDataFromCSV(networkName)
                if vesselData == None:
                    print "\n  could not load network csv, it does not exist! \n"
                else:
                    vascularNetwork.updateNetwork(vesselData)
                #print "     load boundaryData from csv as well? press [u]" 
                #subMenuInput = raw_input("yes [u]? ")
                #if subMenuInput == 'u':
                #    boundaryConditions,boundaryConditionPolyChaos = mCSV.readBCFromCSV(networkName)
                #    vascularNetwork.update({'boundaryConditions':boundaryConditions,
                 #                           'boundaryConditionPolyChaos':boundaryConditionPolyChaos})          
                
                    mainGraph.update_graph(vascularNetwork, window)
                
                print "     write to XML"
                vascularNetwork.name = networkName
                mXML.writeNetworkToXML(vascularNetwork)
            else:
                print "\n  could not update network \n"
            
    print "bye bye .."

if __name__ == '__main__':
    main()
