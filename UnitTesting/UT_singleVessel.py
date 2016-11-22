import unittest
import sys, os, gc
import numpy as np
import math

cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

import UtilityLib.moduleXML as mXML
import SolverLib.class1DflowSolver as c1dFS


class setUp_singleVessel(unittest.TestCase):

    def setUp(self):
        # Set which values should be tested, and the threshold for mean square error
        self.testDict = {"Pressure": 10.0,
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
        self.RMSEDict = {}
        self.vesselDict = {}

        for vesselId in vascularNetworkNew.vessels:
            self.RMSEDict[vesselId] = {}
            self.dataDictNew = {}
            self.dataDictRef = {}
#            print "testing vessel nr. {}".format(vesselId)
            self.vesselDict[vesselId] = []
            for key in self.testDict:
                self.vesselDict[vesselId].append(key)

            self.dataDictNew[vesselId] = vascularNetworkNew.getSolutionData( vesselId,
                                                                        self.vesselDict[vesselId],
                                                                        vascularNetworkNew.simulationTime, [0.0, 0.2])
            self.dataDictRef[vesselId] = vascularNetworkRef.getSolutionData( vesselId,
                                                                        self.vesselDict[vesselId],
                                                                        vascularNetworkNew.simulationTime, [0.0, 0.2])

#       If this is needed for other tests
#       uncomment the following:

#       self.vascularNetworkNew = vascularNetworkNew
#       self.vascularNetworkRef = vascularNetworkRef

#       and comment the next three lines
        del vascularNetworkNew
        del vascularNetworkRef
        gc.collect()

    def tearDown(self):
        del self.RMSEDict
        del self.dataDictNew
        del self.dataDictRef
        gc.collect()


class testSingleVessel(setUp_singleVessel):

    def test_errorThreshold(self):
        for vesselId, keylist in self.vesselDict.iteritems():
            for key in keylist:
#                print "testing key {}".format(key)
                testVals = np.array(self.dataDictNew[vesselId][key])
                refVals = np.array(self.dataDictRef[vesselId][key])
                if len(testVals) != len(refVals):
                    print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
                if len(testVals[0]) != len(refVals[0]):
                    print "different number of values for test and reference for vessel {}, key {} !".format(vesselId, key)
                SEArrayEntrance = np.zeros(len(testVals))
                SEArrayExit = np.zeros(len(testVals))
                for i in xrange(len(testVals)):
                    SEArrayEntrance[i] = (testVals[i][0] - refVals[i][0])*(testVals[i][0] - refVals[i][0])
                    SEArrayExit[i] = (testVals[i][1] - refVals[i][1])*(testVals[i][1] - refVals[i][1])
                self.RMSEDict[vesselId][key] = math.sqrt((np.sum(SEArrayEntrance) + np.sum(SEArrayExit))/(len(testVals)*2))

##       Uncomment these if you want to  have a printout of all RMSE values
#        print "Root Mean Square Error for each part of the network is:"
#        print self.RMSEDict

        TooHighError = False
        errorString = ""
        for vesselId in self.RMSEDict:
            for key in self.RMSEDict[vesselId]:
                if self.testDict[key] < self.RMSEDict[vesselId][key] :
                    TooHighError = True
                    errorString = errorString + "\nRMSE was found to be too high for Vessel {}, key {}, with value {} being above threshold of {}".format(vesselId, key, self.RMSEDict[vesselId][key], self.testDict[key])

        self.assertFalse(TooHighError, errorString)
#        if not TooHighError:
#            print "\nAll values below error threshold"
#            print "Test Successful!"


if __name__ == "__main__":
    unittest.main()
