import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')

#sys.path.append(cur+'/../UtilityLib')
import UtilityLib.moduleXML as mXML 
import UtilityLib.moduleStartUp as moduleStartUp

#import UtilityLib.modulePickle as modulePickle
### modulePickle doesn't seem to exist anymore

#sys.path.append(cur+'/../NetworkLib')
from NetworkLib.classVascularNetwork import VascularNetwork 

from optparse import OptionParser
import cPickle
import numpy as np

from pprint import pprint as pp

__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3_dev"


### file to help user post processing

### 1. open solution data file after a simulation

## use -f and -n to define solution file you want to open, or use the inbulid menu

optionsDict = moduleStartUp.parseOptions(['f','n','c'], visualisationOnly = True)
    
networkName           = optionsDict['networkName']
dataSetNumber         = optionsDict['dataNumber']

print dataSetNumber

##  open data file choosed above
try:
    print " Try to open network {} with data number {}".format(networkName, dataSetNumber)
    vascularNetwork = mXML.loadNetworkFromXML(filename = networkName, dataNumber = dataNumber) 
    vascularNetwork.linkSolutionData() 
except:
    print "Error could not open solution data with data number {} of network {}".format(dataSetNumber,networkName)
    exit()
    

### 2. saveData To csv

vesselId = 1
n = 5 #gridnode
filename = "reymondVessel1.csv"

## load data into memory:
##vascularNetwork.laodData 

Psol = vascularNetwork.vessels[vesselId].Psol[:,[n]]
Qsol = vascularNetwork.vessels[vesselId].Qsol[:,[n]]
Asol = vascularNetwork.vessels[vesselId].Asol[:,[n]]

Psol,Qsol,Asol = vascularNetwork.getSolutionData(['Pressure','Qsol','Asol'], vesselId, t1, t2, gridnode)

Psol2,Qsol2,Asol2 = vascularNetwork.getSolutionData(['Pressure','Qsol','Asol'], vesselId, t1, t2, gridnode)

timeValues = vascularNetwork.simulationTime.ravel()

#'P_f [Pa]','P_b [Pa]', 'Q_f [m^3/s]', 'Q_b [m^3/s]'

import csv
with open(filename, 'wb') as f:
    writer = csv.writer(f,delimiter=';')
    writer.writerow(['t [s]', 'P [Pa]','Q [m^3/s]','A [m^2]'])
    for t,p,q,a in zip(timeValues,Psol.ravel(),Qsol.ravel(),Asol.ravel()):
        writer.writerow([t,p,q,a])
