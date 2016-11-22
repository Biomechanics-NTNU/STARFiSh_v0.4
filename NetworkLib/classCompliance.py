import  numpy as np 
import sys, os
# raise errors
np.seterr(invalid='raise')

cur = os.path.dirname(os.path.realpath(__file__))

sys.path.append(cur+'/../')

import UtilityLib.classStarfishBaseObject as cSBO

class Compliance(cSBO.StarfishBaseObject):
    """
    Base-Class for the Compliance-models
    """    
    def __init__(self, rho, As):
        """
        rho  : flunode density
        Ps   : reference pressure of the compliance model
        As   : area at reference pressure of the compliance model 
        beta : Material model parameter beta
        P0   : P0
        """
        self.rho                 = rho
        self.As                 = As                
        self.Ps                 = None
        self.externalPressure     = None
        self.complianceType     = None
        
        self.C0preCalculated     = None 
            
    def A(self, P):
        """
        Function which calculates the area from given pressure
        Input:  Pressure
        Output: Area
        """
        pass
    def A_Node(self, P, node):
        """
        Function which calculates the area from given pressure
        at certain node of the vessel
        Input:  Pressure, node
        Output: Area
        """
        pass
    def C(self, P) : 
        """
        Function which calculates the compliance from given pressure
        Input:  Pressure
        Output: Compliance
        """
        pass
    def C_Node(self, P, node) : 
        """
        Function which calculates the compliance from given pressure
        at certain node of the vessel
        Input:  Pressure, node
        Output: Compliance
        """
        pass
    def C0(self, P):
        """
        Function which calculates the compliance from given reference pressure
        Input:  Pressure ( not used just for consistency needed)
        Output: pre calculated constant Compliance C0vector
        """
        return self.C0preCalculated
    def C0_Node(self, P, node):
        """
        Function which calculates the compliance from given reference pressure
        at certain node of the vessel
        Input:  Pressure ( not used just for consistency needed), node
        Output: pre calculated constant Compliance C0vector[node]
        """
        return self.C0preCalculated[node]    
    
    def update(self, complianceData):
        """
        updates the compliance data using a dictionary in from of 
        complianceData = {'variableName': value}
        """
        for key,value in complianceData.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except Exception:
                self.warning("compliance.updateData Wrong key: {}, could not update varibale".format(key))
#                print "ERROR compliance.updateData Wrong key: {}, could not update varibale".format(key)
    
    def adaptMaterialCoefficient(self, alpha_c, materialCoefficient = None):
        """
        Function which updates the specific stiffness/compliance coefficient of
        the wall model by a coefficient alpha_c
        
        Args alpha_c (float): correction factor alpha_c to update the wall model coefficient
        """
        raise NotImplementedError("adaptMaterialCoefficient is not implemented for this wall model type")

    def getMaterialCoefficient(self):
        raise NotImplementedError("getMaterialCoefficient is not implemented for this wall model type")


class Exponential(Compliance):
    """
    Exponential Compliance Model 
    """
    def __init__(self, rho, As):
        Compliance.__init__(self, rho, As)
        self.betaExponential = None
        self.P0 = 0 #self.Ps*np.exp(self.beta*(self.A0/self.As-1.))
        
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        self.betaExponential = np.ones(len(self.As))*self.betaExponential 
        self.C0preCalculated = self.C(self.Ps)
    
    def A(self, P):
        #return self.As * (np.log((P + self.P0) / self.Ps) / self.betaExponential + 1.)
        P = P-self.externalPressure
        return self.As * (np.log(P/ self.Ps) / self.betaExponential + 1.)
    def C(self, P):
        #return self.As / (self.betaExponential * (P + self.P0))
        P = P-self.externalPressure
        return self.As / (self.betaExponential * P)    
    
    def A_Node(self, P, node):
        #return self.As[node] * (np.log((P[node] + self.P0[node]) / self.Ps) / self.betaExponential[node] + 1.)
        P = P[node]-self.externalPressure
        return self.As[node] * (np.log(P / self.Ps) / self.betaExponential[node] + 1.)
    def C_Node(self, P, node):
        #return self.As[node] / (self.betaExponential[node] * (P[node] + self.P0[node]))
        P = P[node]-self.externalPressure
        return self.As[node] / (self.betaExponential[node] * P)    

   
class Laplace(Compliance):
    """
    Laplace Compliance Model actually Hookean Model
    """
    def __init__(self, rho, As):
        Compliance.__init__(self, rho, As)
        self.betaLaplace   = None
        
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        self.betaLaplace     = np.ones(len(self.As))*self.betaLaplace 
        self.C0preCalculated = self.C(self.Ps)
    
    def A(self, P):
        P = P-self.externalPressure-self.Ps
        return (P/self.betaLaplace + np.sqrt(self.As))** 2.
    def C(self, P):
        P = P-self.externalPressure-self.Ps
        return (2.*(P / self.betaLaplace + np.sqrt(self.As))) / self.betaLaplace
    
    def A_Node(self, P, node):
        P = P[node]-self.externalPressure-self.Ps
        return (P / self.betaLaplace[node] + np.sqrt(self.As[node]))**2.
    def C_Node(self, P, node):
        P = P[node]-self.externalPressure-self.Ps
        return (2.*(P / self.betaLaplace[node] + np.sqrt(self.As[node]))) / self.betaLaplace[node]

    def adaptMaterialCoefficient(self, alpha_c, materialCoefficient = None):
        """
        Function which updates the specific stiffness/compliance coefficient of
        the wall model by a coefficient alpha_c
        
        Args alpha_c (float): correction factor alpha_c to update the wall model coefficient
        """
        if type(alpha_c) in [float,int,np.float64]:
            if materialCoefficient == None:
                self.betaLaplace = self.betaLaplace*alpha_c
            else:
                self.betaLaplace = materialCoefficient*alpha_c
        else:
            raise TypeError("wall model Laplace.adaptMaterialCoefficient() got alpha_c of wrong format, should be float/int")
            
    def getMaterialCoefficient(self):
        return self.betaLaplace

    def estimateMaterialCoefficientEmpiricalAreaRel(self):
        '''
        modified with diameter wavespeed equation found in Reymond et al. 2009
        '''
        # calculate an averaged beta coefficient based on the diameter area relation ship
        Amm = self.As*1e6 
        c = (13.3/(np.sqrt(Amm*4./np.pi)**0.3))
        self.betaLaplace = 2.0*c**2.*self.rho/np.sqrt(self.As)

class Laplace2(Laplace):
    """
    Laplace Compliance Model actually Hookean Model
    same model equations as Laplace, but it calculates the marterial parameter from 
    youngs modulus and wallThickness
    """
    def __init__(self, rho, As):
        Compliance.__init__(self, rho, As)
        self.wallThickness = None
        self.youngModulus  = None
        self.betaLaplace   = None
        
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        self.betaLaplace = (4./3.)*(np.sqrt(np.pi)*self.wallThickness*self.youngModulus)/self.As
        self.C0preCalculated = self.C(self.Ps)
        
        

class Hayashi(Compliance):
    """
    Compliance model found in Hayashi et al. 1993
    """
    def __init__(self, rho, As):
        Compliance.__init__(self, rho, As)
        self.betaHayashi = None
            
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        self.betaHayashi     = np.ones(len(self.As))*self.betaHayashi 
        self.betaHayashi     = np.ones(len(self.As))*self.betaHayashi 
        self.C0preCalculated = self.C(self.Ps)
        
    def A(self, P):
        P = P-self.externalPressure
        return self.As*( 1.0 + np.log( P/self.Ps ) / self.betaHayashi )**2.0    
    def C(self, P):
        P = P-self.externalPressure
        return 2.0* self.As / self.betaHayashi * ( 1.0 + np.log( P/self.Ps ) / self.betaHayashi ) / P
    
    def A_Node(self, P, node):
        P = P[node]-self.externalPressure
        return self.As[node]*( 1.0 + np.log( P/self.Ps ) / self.betaHayashi[node] )**2.0
    def C_Node(self, P, node):
        P = P[node]-self.externalPressure
        return 2.0* self.As[node] / self.betaHayashi[node] * ( 1.0 + np.log( P/self.Ps ) / self.betaHayashi[node] ) / P
    
    def adaptMaterialCoefficient(self, alpha_c, materialCoefficient = None):
        """
        Function which updates the specific stiffness/compliance coefficient of
        the wall model by a coefficient alpha_c
        
        Args alpha_c (float): correction factor alpha_c to update the wall model coefficient
        """
        if type(alpha_c) in [float,int,np.float64]:
            if materialCoefficient == None:
                self.betaHayashi = self.betaHayashi*alpha_c
            else:
                self.betaHayashi = materialCoefficient*alpha_c
            
        else:
            raise TypeError("wall model Hayashi.adaptMaterialCoefficient() got alpha_c of wrong format, should be float/int")
            
    def getMaterialCoefficient(self):
        return self.betaHayashi
    
    def estimateMaterialCoefficientEmpiricalAreaRel(self):
        '''
        modified with diameter wavespeed equation found in Reymond et al. 2009
        '''
        # calculate an averaged beta coefficient based on the diameter area relation ship
        Amm = self.As*1e6 
        c = (13.3/(np.sqrt(Amm*4./np.pi)**0.3))
        self.betaHayashi = c**2.*2.*self.rho/self.Ps
    
class HayashiEmpirical(Compliance):
    """
    Compliance model found in Hayashi et al. 1993
    modified with diameter wavespeed equation found in Reymond et al. 2009
    """
    def __init__(self, rho, As):
        Compliance.__init__(self, rho, As)
        self.betaHayashi = None
        self.rho = rho
             
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        
        #print "DB using area relation for beta", 
        self.C0preCalculated = self.C(self.Ps)
        
#         print "beta hayashi", self.betaHayashi
#         print "beta HAy mean", np.mean((13.3/(np.sqrt(Amm*4./np.pi)**0.3))**2.*2.*self.rho/self.Ps*self.betaHayashi)
#         print "beta HAy mean+0.05", np.mean((13.3/(np.sqrt(Amm*4./np.pi)**0.3))**2.*2.*self.rho/self.Ps*(self.betaHayashi+0.05))
#         print "beta HAy mean+0.4", np.mean((13.3/(np.sqrt(Amm*4./np.pi)**0.3))**2.*2.*self.rho/self.Ps*(self.betaHayashi+0.4))
#         print
         
        self.betaHayashi = (13.3/(np.sqrt(Amm*4./np.pi)**0.3))**2.*2.*self.rho/self.Ps*self.betaHayashi 
        
        
    def A(self, P):
        P = P-self.externalPressure
        return self.As*( 1.0 + np.log( P/self.Ps ) / self.betaHayashi )**2.0    
    def C(self, P):
        P = P-self.externalPressure
        return 2.0* self.As / self.betaHayashi * ( 1.0 + np.log( P/self.Ps ) / self.betaHayashi ) / P
    
    def A_Node(self, P, node):
        P = P[node]-self.externalPressure
        return self.As[node]*( 1.0 + np.log( P/self.Ps ) / self.betaHayashi[node] )**2.0
    def C_Node(self, P, node):
        P = P[node]-self.externalPressure
        return 2.0* self.As[node] / self.betaHayashi[node] * ( 1.0 + np.log( P/self.Ps ) / self.betaHayashi[node] ) / P
    
class Reymond(Compliance):
    """
    Compliance model found in Reymond et al. 2009
    
    special parameter:
    
    Cs   : compliance at reference pressure of compliance model
    """
    def __init__(self, rho, As):
        
        Compliance.__init__(self, rho, As)
        self.Cs             = None # compliance at reference pressure ps
        
        # defined constants from the paper
        self.a1 = 0.4
        self.b1 = 5.0
        self.PmaxC  = 20. *133.32
        self.Pwidth = 30. *133.32
    
    def initialize(self,complianceDataDict):
        """
        Initilalize compliance class with type specific variables
        and calculate set the marterial parameters
        """
        self.update(complianceDataDict)
        self.Cs              = self.Cs*np.ones_like(self.As)
            
    def A(self, P):
        a1 = self.a1
        b1 = self.b1
        PmaxC  = self.PmaxC
        Pwidth = self.Pwidth
        P = P-self.externalPressure
        intConstant = self.As - self.Cs*(a1*self.Ps + b1*Pwidth*np.arctan((self.Ps-PmaxC)/Pwidth))
        return self.Cs*(a1*(P) + b1*Pwidth*np.arctan(((P)-PmaxC)/Pwidth)) + intConstant
    
    def C(self, P):
        a1 = self.a1
        b1 = self.b1
        PmaxC  = self.PmaxC
        Pwidth = self.Pwidth
        P = P-self.externalPressure
        return self.Cs* (a1 + b1/(1.0+(((P)-PmaxC)/Pwidth)**2.0))
    
    def A_Node(self, P, node):
        a1 = self.a1
        b1 = self.b1
        PmaxC  = self.PmaxC
        Pwidth = self.Pwidth
        Pnode = P[node]-self.externalPressure
        intConstant = self.As[node] - self.Cs[node]*(a1*self.Ps + b1*Pwidth*np.arctan((self.Ps-PmaxC)/Pwidth))
        return self.Cs[node]*(a1*(Pnode) + b1*Pwidth*np.arctan(((Pnode)-PmaxC)/Pwidth)) + intConstant
    
    def C_Node(self, P, node):
        a1 = self.a1
        b1 = self.b1
        PmaxC  = self.PmaxC
        Pwidth = self.Pwidth
        Pnode = P[node]-self.externalPressure
        return self.Cs[node]* (a1 + b1/(1.0+((Pnode-PmaxC)/Pwidth)**2.0))
        
    def adaptMaterialCoefficient(self, alpha_c, materialCoefficient = None):
        """
        Function which updates the specific stiffness/compliance coefficient of
        the wall model by a coefficient alpha_c
        
        Args alpha_c (float): correction factor alpha_c to update the wall model coefficient
        """
        if type(alpha_c) in [float,int,np.float64]:
            
            
            if materialCoefficient == None:
                self.Cs = self.Cs/alpha_c
            else:
                self.Cs = materialCoefficient/alpha_c
            
        else:
            raise TypeError("wall model Reymond.adaptMaterialCoefficient() got alpha_c of wrong format, should be float/int")
      
    def getMaterialCoefficient(self):
        return self.Cs
    
    
    def estimateMaterialCoefficientEmpiricalAreaRel(self):
        '''
        modified with diameter wavespeed equation found in Reymond et al. 2009
        '''
        # calculate an averaged beta coefficient based on the diameter area relation ship
        Amm = self.As*1e6 
        c = (13.3/(np.sqrt(Amm*4./np.pi)**0.3))
        self.Cs = self.As/(c**2.*self.rho)
        