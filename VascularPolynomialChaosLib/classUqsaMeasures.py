

from testBaseClass import TestBaseClass 
import numpy as np

class UqsaMeasures(TestBaseClass):
    # in class definition
    ## pure variables
    variablesHdf5Memory = ["expectedValue",
                           "variance",
                           "standardDeviation",
                           "conficenceInterval",
                           "conditionalExpectedValue",
                           "conditionalVariance",
                           "firstOrderSensitivities",
                           "totalSensitivities",
                           "numberOfSamples"]
    ## dictionary with objects to load
    objectDictsHdf5Memory = []
    
    def __init__(self):
        '''
        This class is a pure data container class
        for all kind of uncertinaty measures.
        It is connected to the baseclass hdf5 reader.
        '''
        # statistic properties
        self.expectedValue            = None
        self.variance                 = None
        self.standardDeviation        = None 
        
        self.conditionalExpectedValue = None
        self.conditionalVariance      = None
        
        self.conficenceInterval       = None
        self.confidenceAlpha          = None
        
        self.firstOrderSensitivities  = None
        self.totalSensitivities       = None       
        
        self.numberOfSamples          = None