import os
import math
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(cur + '/NetworkLib')

class RuntimeMemoryManager(object):
    """
    This class implements a system to limit memory needed to run simulations,
    and to synchronize saving the results of the simulation to data files.

        Initialize RuntimeMemoryManager object for all components.

        Args:
           nSaveBegin := the index in the solution time series where saving begins
           nSaveEnd := the index in the solution time series where saving ends
           nTSteps := the total number of time steps required for the solution
           network := the vascular network object

        Attributes:
           memoryArraySize    := the number of time steps to be stored in the runtime memory
           currentMemoryIndex := the index in vessel memory buffer corresponding to currentTimeStep
           currentTimeStep    := the current position of the solver relative to the total solution


    """

    def __init__(self, nSaveBegin, nSaveEnd, nSaveSkip, nTsteps, maxMemory):
        self.nTSteps = nTsteps
        self.memoryArraySizeTime = None
        self.totalDataPointsPerTimeStep = 0
        self.registeredData = []
        self.chunkCount = 0
        self.nDCurrent = 0
        self.nSaveBegin= nSaveBegin
        self.nSaveEnd = nSaveEnd
        self.nSkipShift = 0
        self.nSaveSkip = nSaveSkip
        self.maxMemory = maxMemory
        self.savedArraySize =  (nSaveEnd-nSaveBegin)//nSaveSkip + 1
        # Updated when called by flow solver
        self.memoryOffset = [0]
        self.currentMemoryIndex = [0]
        self.currentTimeStep = [0]

    def registerGlobalTimeSteps(self, currentMemoryIndex, currentTimeStep):
        # Interface method to link solver time step to memory manager
        self.currentMemoryIndex = currentMemoryIndex
        self.currentTimeStep = currentTimeStep

    def registerDataSize(self,dataSizes):
        """
        Args:
            dataSizes := a tuple with the number of values stored at a single time point
        """
        self.totalDataPointsPerTimeStep += sum(dataSizes)
        estimatedMemorySolutionDataSpace = self.totalDataPointsPerTimeStep*8
        optimalMemorySize = int(math.floor(self.maxMemory * 1024.*1024. / estimatedMemorySolutionDataSpace))
        self.memoryArraySizeTime = max(2, optimalMemorySize)

        # Don't allocate more memory than needed
        if self.memoryArraySizeTime > (self.nTSteps + 1):
            self.memoryArraySizeTime = self.nTSteps + 1

    # Based on the interface implemented in StarfishObjectBase
    def registerSimulationObject(self,simObject):
        solutionMemory, dataBuffers = simObject.getSolutionMemory()
        self.registerSimulationData(solutionMemory, dataBuffers)

    def registerSimulationData(self, solutionMemory, dataBuffers):
        """
        Args:
            solutionMemory := a list of numpy array objects to be registered with the memory manager
            dataBuffers := a list corresponding to the solutinMemoryList with either a dataBuffer or
             the value None if the corresponding vector doesn't need to be saved.
        """
        self.registeredData.extend(zip(solutionMemory,dataBuffers))
        # Each simulation object registers itself here by passing in a tuple of spatial dimensions
        # The

    def registerSaveData(self, solutionMemory,  dataBuffers):
        pass

    def rollSimulationData(self):
        for solutionMemory, dataBuffer in self.registeredData:
            solutionMemory[0] = solutionMemory[-1]


    def flushSolutionMemory(self):
        """
        saving utility function to determine if solution data needs to be sent to the output file,
        and to calculate the correct indices between solution memory and the data output file.

        Explanation of index variables
        nCB,nCE where the beginning and end of the current solution data in memory would lie
         in a full time history of the solution. These are position indices
        nMB,nME what indices of the current memory mark the beginning and end of what should be saved
         nME is a slice index, i.e. position + 1
        nSB,nSE where the beginning and end of the current data to save would lie in the whole of the
         simulation. nSE is a slice index, i.e. position + 1
        nDB,nDE where does the current selection of data belong in the whole of the saved data
        """

        memoryArraySize = self.memoryArraySizeTime;
        offset = (memoryArraySize - 1) * self.chunkCount
        # indices mapping beginning and end of memory to the absolute number time steps in solution
        nCB = offset+1
        nCE = (memoryArraySize - 1) * (self.chunkCount + 1)# nCE == currentTimeStep+1
        nDB = None
        nDE = None
        # check if we need to save
        saving = not(nCE < self.nSaveBegin or nCB > self.nSaveEnd) # not(not saving)

        if saving:
            if self.chunkCount == 0:
                # saving first time point if starting from 0
                nMB = self.nSaveBegin
            elif self.nSaveBegin > nCB:
                nMB = 1 + (self.nSaveBegin-nCB)
            else:
                nMB = 1   + self.nSkipShift

            #determine length to write
            # 1. assume we write out through the end of memory
            # 1.a Accounting for skipping write first value, then as many values as remain divisible by the skip
            nME = memoryArraySize

            # 2. correct assumption 1 if save index is less than the current time step
            if self.nSaveEnd < nCE:
                # set the index to end saving, and include the final time point
                nME -= (nCE - (self.nSaveEnd))


            numWrittenAfterFirst = ((nME - nMB) -1)//self.nSaveSkip
            lengthToWrite = 1 + numWrittenAfterFirst

            # 3. how many entries are skipped after rolling memory
            self.nSkipShift = (self.nSaveSkip - (nME-nMB)%self.nSaveSkip)%self.nSaveSkip

            # Where in the data file to put the data
            nDB = self.nDCurrent
            nDE = self.nDCurrent + lengthToWrite
            self.nDCurrent += lengthToWrite


        for solutionMemory, dataBuffer in self.registeredData:
            if saving and dataBuffer:
                dataBuffer[nDB:nDE] = solutionMemory[nMB:nME:self.nSaveSkip]
            solutionMemory[0] = solutionMemory[-1]
        self.chunkCount += 1

    def runtimeUpdate(self):
        currentMemoryIndex = self.currentMemoryIndex[0]

        if currentMemoryIndex == self.memoryArraySizeTime - 2:
            self.flushSolutionMemory()
            self.memoryOffset[0] = (self.memoryArraySizeTime - 1) * self.chunkCount
        # TODO: This depends on being the last object called for the time step!!!!
        elif self.currentTimeStep[0] == self.nTSteps - 1:
            self.flushSolutionMemory()

    def __call__(self):
        '''
        call function for DataHandler to save the data for each vessel in the network
        '''
        self.runtimeUpdate()




