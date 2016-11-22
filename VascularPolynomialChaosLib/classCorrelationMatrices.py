from testBaseClass import TestBaseClass 
import numpy as np

class Correlation(TestBaseClass):

    def assembleCorrelationMatrix(self, definedBases):
        pass

class CorrelationMatrix(Correlation):
    
    externVariables      = {'bases': TestBaseClass.ExtValue([str], multiVar=True, strCases = ['anything']),
                            'data' : TestBaseClass.ExtValue([float], multiVar=True)}
    externXmlAttributes  = []
    externXmlElements    = ['bases',
                            'data']
    
    def __init__(self):
        '''
        Init function of the correlation matrix class.
        This class defines a correlation matrix for random variables in the
        xml file.
        After a correlation matrix is read from xml file (str)
        it can be converted to a corresponding numpy array with
        the assembleCorrelationMatrix function.
        
        Attributes:
            bases (str): the name of the random variables which form the bases of the correlation matrix.
            data (str): list of data entries for the correlation matrix (Although the correlation matrix is symmetric
            it must be defined in its full form!

        '''
        self.bases = None
        self.data  = None 
    
    def assembleCorrelationMatrix(self, definedBases):
        '''
        Function assemble a numpy array corresponding to the values given in the xml correlation matrix.
        As the order of the bases in the xml-correlation matrix (in the file) is not neccessary in the same order
        as the bases used in the code (alphabetical) the matrix is reordered accoding to the defined in code ordering.
        
        Args:
            definedBases (list): list of defined bases in code
        
        Returns:
            dataMatrix (np.array): correlation matrix in numpy array form
        '''
        # check if definedBases and bases are the same, if not raise error
        if set(definedBases) == set(self.bases):
            dim = len(self.bases)
            if len(self.data) == dim*dim:
                # create data matrix
                dataMatrix = np.array(self.data).reshape((dim,dim))
                # generate index array
                idx = np.empty(dim, int)
                for i,base in enumerate(self.bases):
                    idx[i] = int(definedBases.index(base))       
                # generate permutation matrix
                p = np.eye(dim)[idx]
                # permutate matrix
                dataMatrix = np.dot(np.dot(p,dataMatrix),p.T)
            else:
                raise ValueError("Defined data {0} for correlation matrix is not well defined, dim {1}x{1} is not met".format(self.data,dim))
        else:
            raise ValueError("Defined bases {} for correlation matrix do not reflect defined random inputs with ass. random variables: {}".format(self.bases, definedBases))
                
        return dataMatrix
        