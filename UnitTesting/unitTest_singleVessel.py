import sys, os, gc
import numpy as np
import math

cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

import UtilityLib.moduleXML as mXML
import SolverLib.class1DflowSolver as c1dFS

def test_singleVessel():

    # Set which values should be tested, and the threshold for mean square error
    testDict = {"Pressure": 10.0,
                "Area": 10.0,
                "Flow": 10.0}

    # load and run a new simulation of reference
    networkName = "singleVessel"
    dataNumber = "999"
    dataNumberSol = "012"
    networkXmlFileLoad = cur + "/singleVessel/singleVessel.xml"
    networkXmlSolutionLoad = cur + "/singleVessel/singleVessel_SolutionData_012.xml"

    # Temporary files for saving data
    networkXmlFileSave = cur + "/tmp.xml"
    pathSolutionDataFilename = cur + "/tmpSol.hdf5"


    vascularNetworkNew = mXML.loadNetworkFromXML(networkName,
                                                 dataNumber,
                                                 networkXmlFile = networkXmlFileLoad,
                                                 pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.quiet = True
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()
    # link simulation data
    vascularNetworkNew.linkSolutionData()

    # load reference data and link it
    vascularNetworkRef = mXML.loadNetworkFromXML(networkName,
                                                      dataNumber = dataNumberSol,
                                                      networkXmlFile = networkXmlSolutionLoad,
                                                      pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkRef.linkSolutionData()

    # Make a dictionary for the root mean square errors
    RMSEDict = {}

    for vesselId in vascularNetworkNew.vessels:
        RMSEDict[vesselId] = {}
        dataDictNew = {}
        dataDictRef = {}
        keyList = []
        for key in testDict:
            keyList.append(key)

        dataDictNew = vascularNetworkNew.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])
        dataDictRef = vascularNetworkRef.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])

        for key in testDict:
            testVals = np.array(dataDictNew[key])
            refVals = np.array(dataDictRef[key])
            if len(testVals) != len(refVals):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            if len(testVals[0]) != len(refVals[0]):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            SEArrayEntrance = np.zeros(len(testVals))
            SEArrayExit = np.zeros(len(testVals))
            for i in xrange(len(testVals)):
                SEArrayEntrance[i] = (testVals[i][0] - refVals[i][0])*(testVals[i][0] - refVals[i][0])
                SEArrayExit[i] = (testVals[i][1] - refVals[i][1])*(testVals[i][1] - refVals[i][1])
            RMSEDict[vesselId][key] = math.sqrt((np.sum(SEArrayEntrance) + np.sum(SEArrayExit))/(len(testVals)*2))

    ##    Uncomment these if you want to  have a printout of all RMSE values
    #    print "Root Mean Square Error for each part of the network is:"
    #    print RMSEDict

    TooHighError = False
    for vesselId in RMSEDict:
        for key in RMSEDict[vesselId]:
            if testDict[key] < RMSEDict[vesselId][key] :
                TooHighError = True
                print "Error was found to be too high for Vessel {}, key {}, with value {} being above threshold of {}".format(vesselId, key, RMSEDict[vesselId][key], testDict[key])

    assert(not TooHighError)
    if not TooHighError:
        print "\nAll values below error threshold"
        print "Test Successful!"

def test_saveSkipping():

    # Set which values should be tested, and the threshold for mean square error
    testDict = {"Pressure": 10.0,
                "Area": 10.0,
                "Flow": 10.0}

    # load and run a new simulation of reference
    networkName = "singleVessel"
    dataNumber = "999"
    dataNumberSol = "012"
    networkXmlFileLoad = cur + "/singleVessel/singleVessel.xml"
    networkXmlSolutionLoad = cur + "/singleVessel/singleVessel_SolutionData_012.xml"

    # Temporary files for saving data
    networkXmlFileSave = cur + "/tmpNoSkip.xml"
    pathSolutionDataFilename = cur + "/tmpSolNoSkip.hdf5"


    vascularNetworkNew = mXML.loadNetworkFromXML(networkName,
                                                 dataNumber,
                                                 networkXmlFile = networkXmlFileLoad,
                                                 pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.quiet = True
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()

    # modified minSaveDT
    previousDT = vascularNetworkNew.saveDt
    previousCFL = vascularNetworkNew.CFL
    networkXmlFileSave = cur + "/tmpSkip.xml"
    pathSolutionDataFilename = cur + "/tmpSkipSol.hdf5"

    vascularNetworkNew = mXML.loadNetworkFromXML(networkName,
                                                 dataNumber,
                                                 networkXmlFile = networkXmlFileLoad,
                                                 pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.quiet = True
    vascularNetworkNew.minDt = previousDT
    vascularNetworkNew.CFL = previousCFL/2.
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()
    # link simulation data
    vascularNetworkNew.linkSolutionData()

    # load reference data and link it
    vascularNetworkRef = mXML.loadNetworkFromXML(networkName,
                                                      dataNumber = dataNumberSol,
                                                      networkXmlFile = networkXmlSolutionLoad,
                                                      pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkRef.linkSolutionData()

    # Make a dictionary for the root mean square errors
    RMSEDict = {}

    for vesselId in vascularNetworkNew.vessels:
        RMSEDict[vesselId] = {}
        dataDictNew = {}
        dataDictRef = {}
#        print "testing vessel nr. {}".format(vesselId)
        keyList = []
        for key in testDict:
            keyList.append(key)

        dataDictNew = vascularNetworkNew.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])
        dataDictRef = vascularNetworkRef.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])

        for key in testDict:
#            print "testing key {}".format(key)
            testVals = np.array(dataDictNew[key])
            refVals = np.array(dataDictRef[key])
            if len(testVals) != len(refVals):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            if len(testVals[0]) != len(refVals[0]):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            SEArrayEntrance = np.zeros(len(testVals))
            SEArrayExit = np.zeros(len(testVals))
            for i in xrange(len(testVals)):
                SEArrayEntrance[i] = (testVals[i][0] - refVals[i][0])*(testVals[i][0] - refVals[i][0])
                SEArrayExit[i] = (testVals[i][1] - refVals[i][1])*(testVals[i][1] - refVals[i][1])
            RMSEDict[vesselId][key] = math.sqrt((np.sum(SEArrayEntrance) + np.sum(SEArrayExit))/(len(testVals)*2))

##    Uncomment these if you want to  have a printout of all RMSE values
#    print "Root Mean Square Error for each part of the network is:"
#    print RMSEDict

    TooHighError = False
    for vesselId in RMSEDict:
        for key in RMSEDict[vesselId]:
            if testDict[key] < RMSEDict[vesselId][key] :
                TooHighError = True
                print "Error was found to be too high for Vessel {}, key {}, with value {} being above threshold of {}".format(vesselId, key, RMSEDict[vesselId][key], testDict[key])

    assert(not TooHighError)
    if not TooHighError:
        print "\nAll values below error threshold"
        print "Test Successful!"

#        print dataDictNew
#        print dataDictNew
#        print dataDictRef
# root mean square

def test_memoryChunking():
    # Set which values should be tested, and the threshold for mean square error
    testDict = {"Pressure": 10.0,
                "Area": 10.0,
                "Flow": 10.0}

    # load and run a new simulation of reference
    networkName = "singleVessel"
    dataNumber = "999"
    dataNumberSol = "012"
    networkXmlFileLoad = cur + "/singleVessel/singleVessel.xml"
    networkXmlSolutionLoad = cur + "/singleVessel/singleVessel_SolutionData_012.xml"

    # modified maxMemory to force chunking

    networkXmlFileSave = cur + "/tmpChunk.xml"
    pathSolutionDataFilename = cur + "/tmpChunkSol.hdf5"

    vascularNetworkNew = mXML.loadNetworkFromXML(networkName,
                                                 dataNumber,
                                                 networkXmlFile = networkXmlFileLoad,
                                                 pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.maxMemory = 0.00001
    vascularNetworkNew.quiet = True
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()


    # link simulation data
    vascularNetworkNew.linkSolutionData()

    # load reference data and link it
    vascularNetworkRef = mXML.loadNetworkFromXML(networkName,
                                                      dataNumber = dataNumberSol,
                                                      networkXmlFile = networkXmlSolutionLoad,
                                                      pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkRef.linkSolutionData()

    # Make a dictionary for the root mean square errors
    RMSEDict = {}

    for vesselId in vascularNetworkNew.vessels:
        RMSEDict[vesselId] = {}
        dataDictNew = {}
        dataDictRef = {}
#        print "testing vessel nr. {}".format(vesselId)
        keyList = []
        for key in testDict:
            keyList.append(key)

        dataDictNew = vascularNetworkNew.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])
        dataDictRef = vascularNetworkRef.getSolutionData( vesselId,
                                                          keyList,
                                                          vascularNetworkNew.simulationTime, [0.0, 0.2])

        for key in testDict:
#            print "testing key {}".format(key)
            testVals = np.array(dataDictNew[key])
            refVals = np.array(dataDictRef[key])
            if len(testVals) != len(refVals):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            if len(testVals[0]) != len(refVals[0]):
                print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
            SEArrayEntrance = np.zeros(len(testVals))
            SEArrayExit = np.zeros(len(testVals))
            for i in xrange(len(testVals)):
                SEArrayEntrance[i] = (testVals[i][0] - refVals[i][0])*(testVals[i][0] - refVals[i][0])
                SEArrayExit[i] = (testVals[i][1] - refVals[i][1])*(testVals[i][1] - refVals[i][1])
            RMSEDict[vesselId][key] = math.sqrt((np.sum(SEArrayEntrance) + np.sum(SEArrayExit))/(len(testVals)*2))

##    Uncomment these if you want to  have a printout of all RMSE values
#    print "Root Mean Square Error for each part of the network is:"
#    print RMSEDict

    TooHighError = False
    for vesselId in RMSEDict:
        for key in RMSEDict[vesselId]:
            if testDict[key] < RMSEDict[vesselId][key] :
                TooHighError = True
                print "Error was found to be too high for Vessel {}, key {}, with value {} being above threshold of {}".format(vesselId, key, RMSEDict[vesselId][key], testDict[key])
    assert(not TooHighError)

    if not TooHighError:
        print "\nAll values below error threshold"
        print "Test Successful!"

if __name__ == "__main__":
    test_singleVessel()
    test_saveSkipping()
    test_memoryChunking()
