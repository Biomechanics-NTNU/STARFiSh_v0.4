

from testBaseClass import TestBaseClass 

class RandomInput(TestBaseClass):
    def __init__(self):
        '''
        class decription of a random input to STARFiSh
        
        nomenclature: 
            a <random input> is the describtion of a random variable Y, composed in the form of
            Y = a + b X, where a,b are parameter of choice, and X can be:
                i.) another random input variable
                ii.) a random variable with a propability distribution Normal(0,1) or Uniform(0,1)
                    on a normalized support.
                
        The random input class is base class for the 2 types of random inputs defined by i.) and ii.)
        named GeneralRandomInput (i) and ParametricRandomInput (ii).
        
        From this defintions, it follows that a GeneralRandomInput can be only connected to other random inputs, where as
        ParametricRandomInput is only connected a variable in the starfish base code.
                
        Every random input is defined in the network xml file and as such these classes are 
        also data objects, created by reading the xml file.
        
        After a random input has calculated its realisation, it can forward this to all connected
        other random inputs, or variables in the starfish base code with the updateMethods defined
        in the updateMethods dictionary, this dictionary is set up by the randomInputManager. 
        
        '''
        self.variableName     = [] # list of variable names (e.g. betaHayashi, Z1) which are associated with the random input
        ##        
        self.name             = None # name of the random input
        ##
        self.distributionType = None # distribution type: Normal(a,b), Uniform(a,b), another random input Z. a+bZ
        self.a                = 0 # 
        self.b                = 0 #
        ##
        self.updateMethods    = {} # method to update connected random inputs or variables in starfish base code
        self.updateLog        = [] # list to save the passed realisations
                
        self.printOutDistributions = {'Uniform': 'U(0,1)',
                                      'Normal' : 'N(0,1)'}
        
    def evaluateRealisationFromSample(self,sample_i):
        '''
        Function to calculate the realisation of the random input defined Y = a + b X,
        the realisation becomes y_i =  a + b x_i
                
        where the x_i (sample_i) comes from the governing distribution: eg. Uniform, Normal or another RandomInput 
        
        Args:
            sample_i (float): realisation x_i of the random variable X
            
        Returns:
            realisation (float): realisation of the random input Y; y_i =  a + b x_i
        '''
        realisation = self.a + self.b * sample_i
        return realisation
        
    def passRealisationToAssosiatedObj(self, input):
        '''
        This function evlauates a given input, a realisation x_i, 
        calculates it's realisation y_i with the evaluateRealisationFromSample function.
        Then it passes the realisation to all connected random inputs or variables in the
        starfish base code with the defined update methods in the updateMethods dictionary.
        
        Args:
            input (float) or (dict):
                if the input is of type float, the 'sample' comes from the distribution manager
                elif the input is of type dict the sample comes from another random input
        '''
        sample_i = None      
        # if input comes from other random Input
        if isinstance(input,dict):
            sample_i = input[self.name]
        elif isinstance(input,float):        
            # if input comes from distribution handler
            sample_i = input        
        if sample_i != None:
            realisation = self.evaluateRealisationFromSample(sample_i)
            if self.updateMethods != {}:
                #for i,loc in enumerate(self.parameter.split('_')):
                #    if i == 0:
                #        print '{:3} | {:20} | {:21} | {:.4} '.format(self.name,self.variableName,loc, realisation)
                #    else:
                #        print '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
                #print  
                for variableIdentifier,updateMethod in self.updateMethods.iteritems():                    
                    updateMethod({variableIdentifier : realisation})
                self.updateLog.append(realisation)
                    
    def generateInfo(self):
        '''
        Function to create a log-string with random input information.
        
        Returns:
            randomInputInfo (str): information about the random input definitions
        '''
        if self.distributionType in self.printOutDistributions:
            dist = self.printOutDistributions[self.distributionType]
        else:
            dist = self.distributionType
        
        randomInputInfo = []
        
        for i,loc in enumerate(self.parameter.split('_')):
            if i == 0:
                info = '{:3} | {:20} | {:21} | {:3} + {:3} {} '.format(self.name,self.variableName,loc, self.a, self.b, dist)
            else:
                info = '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
            randomInputInfo.append(info)
        randomInputInfo.append('\n')
        return randomInputInfo
        
        
class ParametricRandomInput(RandomInput):
    
    externVariables      = {'parameter': TestBaseClass.ExtValue(str,strCases = ['anything']),
                            'distributionType':TestBaseClass.ExtValue(str,strCases = ['anything']),
                            'a':TestBaseClass.ExtValue(float),
                            'b':TestBaseClass.ExtValue(float),
                            } 
    
    externXmlAttributes  = []
    externXmlElements    = ['parameter',
                            'distributionType',
                            'a',
                            'b'
                            ]    
    
    def __init__(self):
        '''
        Class for a parametric random input variable defined as Y = a + b X which is connected to 
        a variable in the starfish base code, e.g. betaHayashi of the hayashi compliance model.
        '''
        super(ParametricRandomInput,self).__init__()
        
        self.randomInputType = 'parametricRandomInput'
        self.parameter  = None
        
        
class GeneralRandomInput(RandomInput):
    
    externVariables   = {'distributionType':TestBaseClass.ExtValue(str,strCases = ['anything']),
                         'a':TestBaseClass.ExtValue(float),
                         'b':TestBaseClass.ExtValue(float),
                        } 
    
    externXmlAttributes  = []
    externXmlElements    = ['distributionType',
                            'a',
                            'b'
                            ]
    
    def __init__(self):
        '''
        Class for a general random input variable Y = a + b X which can be the random variable X of other random inputs.
        '''
        super(GeneralRandomInput,self).__init__()
        
        self.randomInputType = 'generalRandomInput'
        self.parameter = '-'
    