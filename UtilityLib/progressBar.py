#!/usr/bin/env python
# encoding: utf-8

import sys
import StringIO

class ProgressBar(object):
    """
    progress bar to show the progress of a for loop iteration process in python
    working also for multiprocessing
    
    with the following look:
    
    [###########         ]
    
    Examples:
        create a progress bar with 20 elements for a loop with 200 iterations 
        and loop over the 200 iterations
    
        >>> progressBar = ProgressBar(20, 200)
        [                    ]
        >>> for i in xrange(200): progressBar.progress(i)
        [###################]
    
        create a progress bar with 2 elements for a loop with 12 iterations 
        and loop over the 12 iterations
    
        >>> progressBar = ProgressBar(2, 12)
        [#]
        >>> for i in xrange(12): progressBar.progress(i)
        [##]
        
        create a progress bar with 4 elements for a loop with 10 iterations 
        check the print out after iteration 5, 7 and 10 respectively
        
        >>> progressBar = ProgressBar(4, 10)
        [   ]
        >>> progressBar.progress(5)
        [#  ]
        >>> progressBar.progress(7)
        [##  ]
        >>> progressBar.progress(10)
        [####]
    """
    def __init__(self, nElements, nIterations, subpressPrint = True):
        """
        Initialization of the the progress bar, prints out the brackets of 
        the progress bar immediately
        
        All outputs to stdout will be captured and displayed after the 
        progress bar is fully loaded.

        Args:
            nElements (int): number of elements of the progress bar
            nIterations (int): number of iterations of the process which 
            progress should be displayed
            subpressPrint (bool): subpress print out statements of function
            in the function which is looped
        """
        
        self.subpressPrint = subpressPrint
        # handle for stdout
        self.stdout = sys.stdout
        # buffer for stdout outputs
        self.outputBuffer = StringIO.StringIO()
        # redirect system stdout to buffer
        sys.stdout = self.outputBuffer
        self.finished = False
        
        self.nElements = nElements
        self.nIterations = float(nIterations)
        self.eCount = 1.

        self.writeProgress(self.eCount)
        

    def clearScreen(self, dCount):


        # write backspaces to the start of the loading bar
        backspacing = ''.join([' \b\b' for i in xrange(dCount)])
        self.stdout.write(backspacing)
        self.stdout.flush()


    def writeProgress(self, eCount):
        '''

        '''
        spaces = int(self.nElements-eCount+1)
        # write empty brackets of the loading bar
        prog = ''.join([''.join(['#' for i in xrange(int(eCount)-1)]),''.join([' ' for i in xrange(spaces)])])
        loadingBar = ''.join(['[',prog,']'])
        self.stdout.write(loadingBar)
        self.stdout.flush()
        

    def progress(self, currentIteration):
        """
        checks the progress of the adjoin process, and advances the progress 
        bar # if necessary
        
        Args:
            currentIteration (int): number of the current iterations of the 
            adjoin process
        """
        if self.finished == False:

            # completed tasks of the process in precent
            loopPercentage = currentIteration/(self.nIterations-1.) 
            # completed shown progress in percent
            barPercentage = self.eCount/self.nElements
            
            while (barPercentage <= loopPercentage):
                
                self.clearScreen(self.nElements+2)

                if self.subpressPrint == False:
                    toPrint = self.outputBuffer.getvalue()
                    self.stdout.write(toPrint)                
                    self.stdout.flush()
                    self.outputBuffer.truncate(0)

                self.eCount += 1.

                self.writeProgress(self.eCount)
                barPercentage = self.eCount/self.nElements
            
            if loopPercentage == 1:

                self.stdout.write("\n")
                self.stdout.flush()
                sys.stdout = self.stdout
                sys.stdout.write(self.outputBuffer.getvalue())
                sys.stdout.write("\n")
                self.outputBuffer.close() 
                self.finished = True
if __name__ == '__main__':
    
    import doctest
    doctest.testmod()
    
        
