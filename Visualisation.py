#!/usr/bin/env python
import subprocess

from optparse import OptionParser
import os,sys
cur = os.path.dirname( os.path.realpath( __file__ ) )

#sys.path.append(cur+'/UtilityLib')
import UtilityLib.moduleStartUp as mStartUp



def main():    
    print ""
    print '====================================='
    print '#    STARFiSh_v0.4 Visualisation    #'
    print '====================================='
    print ""
    
    optionsDict = mStartUp.parseOptions(['f','n','v'], visualisationOnly = True)
    
    networkName    = optionsDict['networkName']
    dataNumber     = optionsDict['dataNumber']
    vizOutput      = optionsDict['vizOutput']
    if vizOutput == 'non':
        vizOutput = "2D+3D"
        

    string1 = ' '.join([sys.executable, cur+'/VisualisationLib/class2dVisualisation.py', '-f', networkName, '-n',dataNumber, '-c']) 
    string2 = ' '.join([sys.executable, cur + '/VisualisationLib/class3dVisualisation.py', '-f', networkName, '-n',dataNumber, '-c True']) 
    
    if vizOutput == "2D":
        
        viz2d = subprocess.Popen(string1, shell = True )
            
    if vizOutput == "3D":
        
        viz3d = subprocess.Popen(string2, shell = True )
        
    if vizOutput == "2D+3D":
        
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
