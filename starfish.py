#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )
import logging
logger = logging.getLogger('starfish')
logger.setLevel(logging.DEBUG)
file_h = logging.FileHandler('starfish.log', mode='w')
file_h.setLevel(logging.DEBUG)
file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_h.setFormatter(file_formatter)
logger.addHandler(file_h)

console_h = logging.StreamHandler(sys.stdout)
console_h.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(message)s')
console_h.setFormatter(console_formatter)
logger.addHandler(console_h)


import Simulator as simulator

try: 
    import Visualisation as visualisationToolBox
    visualization_available = True
except ImportError as e: 
    logger.debug("{}".format(e))
    visualization_available = False

try: 
    import vnc 
    vnc_available = True
except ImportError as e: 
    logger.debug("{}".format(e))
    vnc_available = False

try: 
    import VascularPolynomialChaos as uqsaToolBox
    uqsa_available = True 
except ImportError as e: 
    logger.debug("{}".format(e))
    uqsa_available = False

def main():
        
    print("")
    print('=====================================')
    print('#     STARFiSh_v0.4.19.10.2016      #')
    print('=====================================')
    
    mainMenuInput = ""
    while mainMenuInput not in ['q']:
        print('\n Main menu:\n')
        print " [1] - run simulation | Main.py"
        if vnc_available:
            print " [2] - run vnc (vascular network creator) | vnc.py"
        if visualization_available:
            print " [3] - run visualisation | Visualisation.py"
        if uqsa_available:
            print " [4] - run uncertainty quantification tool box | VascularPolynomialChaos.py"
        print " [q] - quit \n"
        while  mainMenuInput not in ('1','2','3','4','q'):
            mainMenuInput = raw_input("what to do? ")
        
        if mainMenuInput == '1':
            print "\n .. running simulation \n"
            simulator.main()
            mainMenuInput = ""
            
        if mainMenuInput == '2':
            if vnc_available:
                print "\n .. running vnc \n"
                vnc.main()
            else:
                print "vnc not available check log for import error"
            mainMenuInput = ""
            
        if mainMenuInput == '3':
            if visualization_available:
                print "\n .. running visualisation \n"
                visualisationToolBox.main()
            else:
                print "visualization package not available check log for import error"
            mainMenuInput = ""

            
        if mainMenuInput == '4':
            if uqsa_available:
                print "\n .. running uncertainty quantification tool box \n"
                uqsaToolBox.uncertaintyPropagation()
            else:
                print "uqsa toolbox not available check log for import error"
            mainMenuInput = ""
            
    print "bye bye .." 
    
if __name__ == '__main__':
    main()
