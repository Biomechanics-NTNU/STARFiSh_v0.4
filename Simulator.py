##!/usr/bin/env python
# -*- coding: utf-8 -*- 
import time 
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )
import logging
logger = logging.getLogger('starfish')
logger.setLevel(logging.DEBUG)


import SolverLib.class1DflowSolver as c1DFlowSolv

import UtilityLib.moduleXML as mXML
import UtilityLib.moduleStartUp as mStartUp #import parseOptions
import UtilityLib.moduleFilePathHandler as mFPH

import matplotlib.pyplot as plt

import gc

import subprocess

def main():
    optionsDict = mStartUp.parseOptions(['f','n','d','s','v','r','w','p'])
    
    networkName           = optionsDict['networkName']
    save                  = optionsDict['save']
    dataNumber            = optionsDict['dataNumber']
    simulationDescription = optionsDict['simulationDescription']
    vizOutput             = optionsDict['vizOutput']
    resimulate            = optionsDict['resimulate']
    
    filename = str(networkName+'.xml')
        
    logger.info('____________Simulation_______________')
    logger.info('%-20s %s' % ('Network name',networkName))
    logger.info('%-20s %s' % ('Data number', dataNumber))
    logger.info('%-20s %s' % ('Save simulation', save))
    logger.info('%-20s %s' % ('Case description', simulationDescription))
    logger.info('%-20s %s' % ('Resimulate', resimulate))
    logger.info('%-20s %s' % ('Visualisationmode', vizOutput))
    
    ## check if template
    if '_template' in networkName:
        networkName = mFPH.createWorkingCopyOfTemplateNetwork(networkName)
    
    # load network from the path!
    if resimulate == False:
        vascularNetwork = mXML.loadNetworkFromXML(networkName) # moved to vascularNetowrk constror
    else:
        # resimulate network
        vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber = dataNumber)        
        if simulationDescription == '':
            simulationDescription = vascularNetwork.description
    
    if vascularNetwork == None: exit()
    
    
    vascularNetwork.update({'description':simulationDescription,
                            'dataNumber' :dataNumber})
    
    timeSolverInitStart = time.clock()
    #initialize Solver
    flowSolver = c1DFlowSolv.FlowSolver(vascularNetwork)
    timeSolverInit = time.clock()-timeSolverInitStart
    timeSolverSolveStart = time.clock()
    #solve the system
    flowSolver.solve()
    timeSolverSolve = time.clock()-timeSolverSolveStart
    
    minutesInit = int(timeSolverInit/60.)
    secsInit = timeSolverInit-minutesInit*60.
    minutesSolve = int(timeSolverSolve/60.)
    secsSolve = timeSolverSolve-minutesSolve*60.
    
    
    logger.info('____________ Solver time _____________')
    logger.info('Initialisation: {} min {} sec'.format(minutesInit,secsInit))
    logger.info('Solving:        {} min {} sec'.format(minutesSolve,secsSolve))
    logger.info('=====================================')
    
    vascularNetwork.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetwork, dataNumber = dataNumber) # needs to be moved to vascularNetwork
    
    
    del flowSolver
    gc.collect()
    
    mFPH.updateSimulationDescriptions(networkName, dataNumber, simulationDescription)
    
    gc.collect()
    
    
    if vizOutput == "2D":
        string = ' '.join([sys.executable, cur + '/VisualisationLib/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber)])                
        subprocess.Popen(string, shell=True)
        
    
    if vizOutput == "3D":
        string = ' '.join([sys.executable, cur+'/VisualisationLib/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        subprocess.Popen(string, shell=True)
        
        
    if vizOutput == "2D+3D":
           
        string1 = ' '.join([sys.executable, cur + '/VisualisationLib/class2dVisualisation.py', '-f', vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        string2 = ' '.join([sys.executable, cur + '/VisualisationLib/class3dVisualisation.py', '-f', vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        
        viz2d = subprocess.Popen(string1, shell = True )
        viz3d = subprocess.Popen(string2, shell = True )
        
        while True:
            
            if viz2d.poll() != None:
                viz3d.terminate()
                exit()
                
            if viz3d.poll() != None:
                viz2d.terminate()
                exit()
        
if __name__ == '__main__':
    main()
