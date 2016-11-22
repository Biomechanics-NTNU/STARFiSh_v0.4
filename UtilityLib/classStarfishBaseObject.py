# SuperClass.py (name should maybe be changed)
import sys, os
import traceback
import time
import numpy as np
import classConfigurableObjectBase
cur = os.path.dirname( os.path.realpath( __file__ ) )
root = ''.join([cur,'/../'])
logFilePath = root + 'warninglog.txt'

class StarfishBaseObject(classConfigurableObjectBase.ConfigurableObjectBase):
    """Super class every other class in STARFiSh
    will inherit from.

    Contains:
        Global Warning function; warn()

        Global Error Appending function; appendException()

    Planned features:
    Global Update function for dictionaries.
    """
    
    # BEGIN Time Series Memory Management Interface 
    solutionMemoryFields    = []
    solutionMemoryFieldsToSave = []
    
    def allocate(self, runtimeMemoryManager):
        self.createSolutionMemory(runtimeMemoryManager.memoryArraySizeTime)
        self.createFileDataBuffers(runtimeMemoryManager.savedArraySize)
        solMemory, dsets = self.getSolutionMemory()
        runtimeMemoryManager.registerSimulationData(solMemory, dsets)
        

    def getSolutionMemorySizes(self):
        sizes = []
        for key in self.solutionMemoryFields:
            try:
                sizes.append(self.__dict__[key].shape[1])
            except IndexError:
                sizes.append(1)
        return sizes

    def createSolutionMemory(self,memorySize):
        for key in self.solutionMemoryFields:
            shape = (memorySize,) + self.__dict__[key].shape[1::]
            self.__dict__[key] = np.zeros(shape)

    def createFileDataBuffers(self, savedArraySize):
        for key in self.solutionMemoryFields:
            data = self.__dict__[key]
            if key in self.solutionMemoryFieldsToSave:
                try:
                    size = data.shape[1]
                    self.dsetGroup.create_dataset(key,(savedArraySize,size))
                except IndexError:
                    self.dsetGroup.create_dataset(key,(savedArraySize,))
    
    def loadFileDataBuffers(self,dsetGroup,nSelectedBegin,nSelectedEnd,nTStepSpaces):
        for key in self.solutionMemoryFields:
            if key in dsetGroup:
                self.__dict__[key] = dsetGroup[key][nSelectedBegin:nSelectedEnd:nTStepSpaces]
            else:
                self.warning("No data set found for {}".format(key))
                
    def getSolutionMemory(self):
        solutionMemory = []
        dataBuffers = []
        for key in self.solutionMemoryFields:
            solutionMemory.append(self.__dict__[key])
            try:
                dataBuffers.append(self.dsetGroup[key])
            except (KeyError, TypeError)  as e:
                dataBuffers.append(None)

        # Ensure the correct number of entries in both tuples
        assert len(self.solutionMemoryFields) == len(solutionMemory), "Number of fields doesn't match number of arrays"
        assert(len(self.solutionMemoryFields) == len(dataBuffers)), "Number of fields doesn't match number of data sets"
        return solutionMemory, dataBuffers
    # END Time Series Memory Management Interface



    def warning(self, infoString = None, noException = False,
                quiet = False, verbose = False, saveToFile = False,
                oldExceptPass = False):
        """
        Global Warning Function

        Args:
            infoString (String): String describing why the warning was made.
            noException (Optional Bool): Boolean determining if it
                should use the last exception's info or not.
                Defaults to False.
            quiet (Opt. Bool): Whether to print warning at all.
                Defaults to False.
            verbose (Opt. Bool): Whether or not to print full traceback.
                Defaults to False.
            saveToFile (Opt. Bool): Whether or not to save info to file.
                Will save info in STARFiSh/warninglog.txt
                Defaults to False.
        """
        if oldExceptPass:
            quiet = True
            saveToFile = False

        completeString = "Warning: "
        if infoString != None:
            completeString = completeString + infoString
        else:
            completeString = completeString + "no string passed to warn()"

        if not noException:
            try: raise
            except: exceptionInfo = "\n" + traceback.format_exc()
        else: exceptionInfo = " "

        if verbose:
            completeString = completeString + exceptionInfo

        #if not quiet:
            # print completeString

        if saveToFile:
            if not verbose: completeString = completeString + exceptionInfo
            try: f = open(logFilePath, 'a')
            except: f = open(logFilePath, 'w')
            f.write("\n======================")
            f.write("\n{}".format(
                time.strftime("%Y.%b.%d, %H:%M:%S")))
            f.write("\n----------------------\n")
            f.write(completeString)
            f.write("======================\n\n")
            f.close()

    def exception(self, infoString = None):
        """
        The exception appending function.

        Args:
            infoString (string): A string with the information you
                wish to append to the exception's message

        """
        if infoString == None:
            try: raise
            except: return sys.exc_info()[0],sys.exc_info()[1],sys.exc_info()[2]
        else:
            try: raise
            except Exception as e:
                raise type(e), type(e)(e.message + "\nAppended : " + infoString), sys.exc_info()[2]







