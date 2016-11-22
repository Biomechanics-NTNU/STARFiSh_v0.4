import sys, os
import numpy as np

cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')
import UtilityLib.classStarfishBaseObject as cSBO

class System(cSBO.StarfishBaseObject):
    
    def __init__(self,vessel,simplifyEigenvalues,riemannInvariantUnitBase,currentTimeStep,dt):
        """
        Constructor of System with the SystemEqations
        
        Input:
            vessel <classVessel>              : vessel for which the system Equations should be build
            simplifyEigenvalues <bool>        : TRUE  == simplified eigenvalues lambda1/2 = +- c
                                                FALSE == lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            riemannInvariantUnitBase <String> : 'Pressure' == calculate R,L based on riemann invariants with unit
                                                              of Pressure
                                                'Flow'     == calculate R,L based on riemann invariants with unit
                                                              of Flow
        Variables (most important):
           
            m12 <np.array> : 12 value of system Matrix M
            m21 <np.array> : 21 value of system Matrix M
            m22 <np.array> : 22 value of system Matrix M
            b2  <np.array> : 2 value of righthandside vector B
            LAMBDA   <nested list> = [<np.array>,<np.array>] : Eigenvalues at each boundary node 
            R        <nested list> = [<np.array>,<np.array>] : Right-Eigenvector matrix at each boundary node 
            L        <nested list> = [<np.array>,<np.array>] : Left-Eigenvector matrix at each boundary node 
            Z        <nested list> = [<np.array>,<np.array>] : Impedances (Z1, Z2) (forward/backward) at each boundary node 
        
        """        
        self.name = ' '.join(['system equations', str(vessel.Id)])
        
        # compliance properties
        self.C    = vessel.C
        self.c    = vessel.c
        
        self.C_nID = vessel.C_nID
        
        #Fluid properties
        self.my     = vessel.my
        self.rho    = vessel.rho
        self.gamma  = vessel.gamma
        
        # vessel net gravity
        self.netGravity = vessel.netGravity
        
        # properties for Riemann invariant calculation
        self.z = vessel.z
        self.dt = dt 
        self.currentTimeStep = currentTimeStep
                
        # system matrices
        self.m12    = None
        self.m21    = None
        self.m22    = None
        self.M      = None
        self.b2     = None
        self.B      = None
        self.LAMBDA = np.ones(2)     
        self.R      = np.ones((2,2)) 
        self.L      = np.ones((2,2))
        self.du     = np.empty(2)
        
        #calculate dlt
        self.dlt = (self.gamma+2.0)/(self.gamma+1.0)
        
        # eigen values
        if simplifyEigenvalues   == True:
            if riemannInvariantUnitBase == 'Flow':
                self.updateLARL = self.updateLARLSys0InvariantFlow
            elif riemannInvariantUnitBase == 'Pressure':
                self.updateLARL = self.updateLARLSys0InvariantPressure
            
        elif simplifyEigenvalues == False:
            if riemannInvariantUnitBase == 'Flow':
                self.updateLARL = self.updateLARLSys1InvariantFlow
            elif riemannInvariantUnitBase == 'Pressure':
                self.updateLARL = self.updateLARLSys1InvariantPressure
    
        self.Re = 3000.0
                    
    def updateSystem(self,P,Q,A,pi=np.pi):
        """
        Update all the system-equation and matrices
        
        Input:
            P <np.array> : current pressure values of the vessel
            Q <np.array> : current flow values of the vessel
            A <np.array> : current area values of the vessel
        
        """
        # current time step of solution
        n = self.currentTimeStep[0]
        # calculate needed values
        C = self.C(P) 
        c = self.c(A,C)
        v = Q/A
        
# #         Re = np.max(np.sqrt(A/np.pi)*2.0*v/self.my*self.rho)    
# #         if Re > self.Re: # > 3000.0
# #             print "WARNING approaching high reynoldsnumber for vessel", self.name, Re
# #             self.Re = Re
#             
        dlt = self.dlt
        
        m12 = (1. / C)
        m21 = (C * (c**2 - dlt * v**2)) # vinz expression #||# (A/self.rho)[1:-1] # paul expression  # (C*(c**2-self.dlt*v**2))[1:-1] # vinz expression
        m22 = (2 * dlt * v)  #vinz expression #||# self.dlt #paul expression #(2*self.dlt*v)[1:-1]
        # set up B matrix
        b2 = (-2.0 / self.rho * pi * v * (self.gamma+2) * self.my + A * self.netGravity[n]) #Correction with net gravity force
        
        ## set up LAMBDA, L and R matrices
        #self.updateLARL(P,Q,A,Ct=C,ct=c)
        
        return m12,m21,m22,b2
            
    def updateLARLSys0InvariantFlow(self,P,Q,A,position):
        """
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys0 :          uses simplified eigenvalues lambda = c (wave speed) 
                            as result Z1=Z2=Zc
            invariantFlow:  calculates R,L using riemannInvariants based on unit FLow
                            as in Paper of Paul Roger Leinan but with R2*(-1) correction to
                            have additive invariants for pressure
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        """
        n = self.currentTimeStep[0]
               
        C = self.C_nID(P,position)
        c = self.c(A[position],C)
                        
        # Z1 = Z2 = Zc
        Zc =  1.0/(C*c)
                        
        LAMBDA = self.LAMBDA
        
        LAMBDA[0] = c
        LAMBDA[1] = -c
        
        L = self.L
        R = self.R     
        
        ## left eigenvalue matrix
        # riemannInvariants with unit flow    
        L[0][0] =  0.5/(Zc)
        L[0][1] =  0.5
        L[1][0] =  0.5/(Zc)
        L[1][1] = -0.5

        ## right eigenvalue matrix
        # riemannInvariants with unit flow 
        R[0][0] =  Zc
        R[0][1] =  Zc
        R[1][1] = -1.0
                                
        ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
        errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
        if errorIdentity > 5.e-16:
            self.warning("SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity), noException= True)
        
        ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
        # calculate omegas
        du = self.du
        zi = self.z[position] - [-c,c][position] * self.dt
                                        
        du[0]  = np.interp(zi,self.z,P) - P[position]
        du[-1] = np.interp(zi,self.z,Q) - Q[position]
        
        domega = np.dot( L[1+position], du) + self.dt * L[1+position][1] * A[position] * self.netGravity[n]
                
        return L,R,LAMBDA,Zc,Zc,domega[0]
             
    def updateLARLSys0InvariantPressure(self,P,Q,A,position):
        """
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys0 :          uses simplified eigenvalues lambda = c (wave speed)
                            as result Z1=Z2=Zc
            invariantFlow:  calculates R,L using riemannInvariants based on unit Pressure
                            as derieved by Leif Rune Hellevik
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        """
        n = self.currentTimeStep[0]
        
        C = self.C_nID(P,position)
        c = self.c(A[position],C)
                        
        # Z1 = Z2 = Zc
        Zc =  1.0/(C*c)
                        
        LAMBDA = self.LAMBDA
        
        LAMBDA[0] = c
        LAMBDA[1] = -c
        
        L = self.L
        R = self.R
        
        ## left eigenvalue matrix
        # riemannInvariants with unit pressure
        L[0][1] =  Zc
        L[1][1] = -Zc

        ## right eigenvalue matrix
        # riemannInvariants with unit pressure
        R[0][0] =  0.5
        R[0][1] =  0.5
        R[1][0] =  0.5/(Zc)
        R[1][1] = -0.5/(Zc)
                                
        ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
        errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
        if errorIdentity > 5.e-16:
            self.warning("SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity), noException= True)
        
        ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
        # calculate omegas
        du = self.du
        
        zi = self.z[position] - [-c,c][position] * self.dt
                    
        du[0]  = np.interp(zi,self.z,P) - P[position]
        du[-1] = np.interp(zi,self.z,Q) - Q[position]
            
        domega = np.dot( L[1+position], du) + self.dt * L[1+position][1] * A[position] * self.netGravity[n]    
            
        return L,R,LAMBDA,Zc,Zc,domega[0]
            
            
    def updateLARLSys1InvariantFlow(self,P,Q,A,position,sqrt=np.sqrt):
        """
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys1 :          uses correct eigenvalues lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            invariantFlow:  calculates R,L using riemannInvariants based on unit FLow
                            as in Paper of Paul Roger Leinan but with R2*(-1) correction to
                            have additive invariants for pressure
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        """
        dlt = self.dlt  
        n = self.currentTimeStep[0]
        
        Aid = A[position]
        
        C = self.C_nID(P,position)
        c = self.c(Aid,C)
        
        v = Q[position]/Aid
        
        l1 = dlt*v+sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
        l2 = dlt*v-sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
        
        LAMBDA = self.LAMBDA
        LAMBDA[0] = l1
        LAMBDA[1] = l2
        
        Z1 =  1.0/(C*l1)
        Z2 = -1.0/(C*l2)
                                            
        L = self.L
        R = self.R
        
        ## left eigenvalue matrix
        # riemannInvariants with unit flow    
        L[0][0] =  1.0/(Z1+Z2)
        L[0][1] =  Z2/(Z1+Z2)
        L[1][0] =  1.0/(Z1+Z2)
        L[1][1] = -Z1/(Z1+Z2)

        ## right eigenvalue matrix
        # riemannInvariants with unit flow 
        R[0][0] =  Z1
        R[0][1] =  Z2
        R[1][1] = -1.0
  
        ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
        errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
        if errorIdentity > 5.e-16:
            self.warning("SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity), noException = True)
        
        ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
        # calculate omegas
        du = self.du
        zi = self.z[position] - [l2,l1][position] * self.dt
                    
        du[0]  = np.interp(zi,self.z,P) - P[position]
        du[-1] = np.interp(zi,self.z,Q) - Q[position]
        
        domega = np.dot( L[1+position], du) + self.dt * L[1+position][1] * Aid * self.netGravity[n]   
          
        return L,R,LAMBDA,Z1,Z2,domega[0]
              
    def updateLARLSys1InvariantPressure(self,P,Q,A,position,sqrt=np.sqrt):
        """
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys1 :          uses correct eigenvalues lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            invariantFlow:  calculates R,L using riemannInvariants based on unit Pressure
                            as derieved by Leif Rune Hellevik
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        """
        dlt = self.dlt  
        n = self.currentTimeStep[0]
            
        Aid = A[position]       
                     
        C = self.C_nID(P,position)
        c = self.c(Aid,C)
        
        v = Q[position]/Aid
        
        l1 = dlt*v+sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
        l2 = dlt*v-sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
                        
        LAMBDA = self.LAMBDA
        LAMBDA[0] = l1
        LAMBDA[1] = l2
        
        Z1 =  1.0/(C*l1)
        Z2 = -1.0/(C*l2)
                                
        L = self.L
        R = self.R
        
        ## left eigenvalue matrix
        # riemannInvariants with unit poressure
        L[0][1] =  Z2
        L[1][1] = -Z1

        ## right eigenvalue matrix
        # riemannInvariants with unit pressure
        R[0][0] =  Z1/(Z1+Z2)
        R[0][1] =  Z2/(Z1+Z2)
        R[1][0] =  1.0/(Z1+Z2)
        R[1][1] = -1.0/(Z1+Z2)
                
        ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
                        
        ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
        errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
        if errorIdentity > 5.e-16:
            self.warning("SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity), noException= True)
            
        ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
        # calculate omegas
        du = self.du
        zi = self.z[position] - [l2,l1][position] * self.dt
                    
        du[0]  = np.interp(zi,self.z,P) - P[position]
        du[-1] = np.interp(zi,self.z,Q) - Q[position]
        
        domega = np.dot( L[1+position], du) + self.dt * L[1+position][1] * Aid * self.netGravity[n]   
          
        return L,R,LAMBDA,Z1,Z2,domega[0]