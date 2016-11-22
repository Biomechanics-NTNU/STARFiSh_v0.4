
import chaospy as cp
import h5py


class DistributionManager(object):
    
    def __init__(self, randomInputVector = None):
        '''
        Distribution Manager base class defining a distribution manager.
        Which creates and manages the marginal distributions and
        joint distribution of a defined case.
        
        The distribution manager defined as base class and super classes to
        account for different stochstic packages handeling the distributions.
        
        Nomenclature: 
            as <marginalDistributions> we define the marginal distributions all X in the 
            random inputs which are defined with an actual propability distribution.
        
            the <jointDistribution> is the joint distributions of all marginal distributions.
                                                
        Args:
            randomInputVector (Optional (list)): 
                Defaults to None. The list should contain instances of class RandomInputs
                created from the read in xml-file.
                
        '''
        # distribution
        self.randomInputsExtDist     = randomInputVector
        self.marginalDistributions = []
        self.jointDistribution     = None
        self.jointDistributionDependent = None
        self.distributionDimension = None
        self.dependentCase = False
        
    def passRealisation(self, sample, sampleIndex):
        '''
        Function to pass samples of the random variables to
        the random inputs which are connected to the network parameters
        
        Args:
            sample (np.array): 
                sample drawn from the joint distribution of thr random variables
                
            sampleIndex (int):
                index of the current sample              
        '''
            
        if len(sample) == len(self.randomInputsExtDist):
            print "\nSample number {}".format(sampleIndex)
            print '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","realisation")
            print "--------------------------------------------------------------------"          
            for randomInput,sample_i in zip(self.randomInputsExtDist,sample):
                print "random variable {} with realisation {}".format(randomInput.name,sample_i)
                randomInput.passRealisationToAssosiatedObj(sample_i)
    
    ## methods created by the toolbox-child class implementation
    def createRandomVariables(self):
        '''
        create a random variable vector from
        for the random input variables in the random input variable vector
        and the joint distribution 
        '''
        pass
    
    def createDependentDistribution(self):
        
        '''
        Method creates a dependent distribution 
        '''
        pass
                     
class DistributionManagerChaospy(DistributionManager):
    
    def __init__(self, randomInputVector = None):
        '''
        Distribution Manager base class defining a distribution manager.
        Which creates and manages the marginal distributions and
        joint distribution of a defined case.
        
        This distribution manager is based around the chaospy package
        
        Args:
            randomInputVector (Optional (list)): 
                Defaults to None. The list should contain instances of class RandomInputs
                created from the read in xml-file.
        
        '''
        super(DistributionManagerChaospy, self).__init__(randomInputVector)
        
        
    def createRandomVariables(self):
        '''
        Create marginal and joint distributions for the random inputs 
        '''
        if self.randomInputsExtDist == None:
            print "WARNING: DistributionManager.createRandomVariablesChaospy() no randomInputsExtDist are defined"
            return
        
        evalDistDict = {'Uniform': cp.Uniform, 
                        'Normal' : cp.Normal}
        
        # create marignal distributions
        for randomInput in self.randomInputsExtDist:
            distType = randomInput.distributionType
            if distType in evalDistDict.keys():
                marginalDistribution = evalDistDict[distType]() #eval(distType,{"__builtins__":None},evalDistDict)()
                self.marginalDistributions.append(marginalDistribution)
            #self.jointDistribution = cp.J(self.jointDistribution, marginalDistribution)
        
        # create joint distributions
        self.jointDistribution = cp.J(*self.marginalDistributions)
        self.distributionDimension = len(self.jointDistribution)
    
        
    def createDependentDistribution(self,CorrelationMatrix):
        '''
        Method creates a dependent distribution with Nataf transformation 
        based on a correlation matrix of the dependent random variables.
        
        Args:
            CorrelationMatrix (np.array):
                Correlation Matrix which describes the correlation matrix between the random variables
                defined in the joint distribution.
        
        '''
        self.dependentCase = True
        self.jointDistributionDependent = cp.Nataf(self.jointDistribution, CorrelationMatrix)
    
         
        