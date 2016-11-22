import sys, os
import numpy as np
from numpy.linalg import solve
from scipy.optimize import fsolve

from pprint import pprint as pp

cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')
import UtilityLib.classStarfishBaseObject as cSBO

from copy import copy as copy       
 
class Link():
    """
    Link object represends the connection between two vessels
    """
    def __init__(self, mother, motherSys, 
                     daughter, daughterSys,
                     currentMemoryIndex, dt, rigidAreas, solvingScheme):
        self.type = 'Link'
        
        self.name = ' '.join(['Link',str(mother.Id),str(daughter.Id)])
        
        self.rho             = []
        self.my              = []
        self.systemEquations = []
        self.z               = []
        self.A_func          = []
        self.positions       = []
        self.names           = []
#       self.vz = []

        self.dt = dt
        
        self.currentMemoryIndex  = currentMemoryIndex
        
        # equations to solve in f solve
        self.fsolveFunction = None
        self.jacobiMatrix = None
        
        #initialize Vessels
        # mother vessel
        self.rho.append(mother.rho)
        self.my.append(mother.my)
        self.z.append(mother.z)
        self.systemEquations.append(motherSys)
        self.positions.append(-1)
        self.names.append(mother.Id)
        self.A_func.append(mother.A_nID)
#        self.vz.append(-1)
        # SolutionVariables
        self.P_mother = mother.Psol
        self.Q_mother = mother.Qsol
        self.A_mother = mother.Asol
               
        # daughter vessel
        self.rho.append(daughter.rho)
        self.my.append(daughter.my)
        self.z.append(daughter.z)
        self.systemEquations.append(daughterSys)
        self.positions.append(0)
        self.names.append(daughter.Id)
        self.A_func.append(daughter.A_nID)
#        self.vz.append(1)
        # SolutionVariables
        self.P_daughter = daughter.Psol
        self.Q_daughter = daughter.Qsol
        self.A_daughter = daughter.Asol
            
#        if   rigidAreas == '02': 
#            self.fsolveFunction = self.fsolveConnectionSys0
##            self.jacobiMatrix = self.jacobiMatrixSys0
#        elif rigidAreas == '2':  
#            self.fsolveFunction = self.fsolveConnectionSys1
#            self.jacobiMatrix = self.jacobiMatrixSys1
#       else: 
#            print"ERROR: classConnections: EquSys not properly defined! system exit"
#            exit()

        self.rigidAreas = rigidAreas
        #solvingScheme = "Stenosis"
        solvingScheme = "NonLinear"
        # Define the call function depending on the solving Scheme
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        elif solvingScheme == "NonLinear":
            self.__call__ = self.callNonLinear
        elif solvingScheme == "Stenosis":
            self.__call__ = self.callStenosisYoungAndTsai
        else:
            raise ImportError("Connections wrong solving scheme! {}".format(solvingScheme))
    
        ## benchamark Test variables
        self.nonLin = False
        self.updateL = False
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        
        P1 = self.P_mother[n]
        Q1 = self.Q_mother[n]
        A1 = self.A_mother[n]
        
        P2 = self.P_daughter[n]
        Q2 = self.Q_daughter[n]
        A2 = self.A_daughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        
        # update system equation and store L1
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_2 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
                
        ###### Linear approach
        denom = (R1_12*R2_21-R1_22*R2_11)
        # -1 reflectionCoeff mother->daugther
        alpha1 = -(R1_11*R2_21-R1_21*R2_11)/denom
        # +1 transmission daugther->mother
        alpha2 = -(R2_11*R2_22-R2_21*R2_12)/denom 
        # +1 transmission daugther->mother
        beta1 = -(R1_11*R1_22-R1_12*R1_21)/denom 
        # -1 reflectionCoeff daugther->mother
        beta2 = -(R1_12*R2_22-R2_12*R1_22)/denom
        
        #print 'cC153 alphas',alpha1,alpha2,beta1,beta2
                
        domega1_2 = alpha1 * domega1_1 + alpha2 * domega2_2
        domega2_1 = beta1  * domega1_1 + beta2  * domega2_2
                
        P1_new = P1o + R1_11*domega1_1 + R1_12*domega1_2
        Q1_new = Q1o + R1_21*domega1_1 + R1_22*domega1_2
    
        P2_new = P2o + R2_11*domega2_1 + R2_12*domega2_2
        Q2_new = Q2o + R2_21*domega2_1 + R2_22*domega2_2
        
        # apply calculated values to next time step
        self.P_mother[n+1][pos1]   = P1_new
        self.Q_mother[n+1][pos1]   = Q1_new
        self.P_daughter[n+1][pos2] = P2_new
        self.Q_daughter[n+1][pos2] = Q2_new
        
        if P1_new < 0 or P2_new < 0:
            raise ValueError("Connection: {} calculated negative pressure, P1_new = {}, P2_new = {}, at time {} (n {},dt {})".format(self.name,P1_new, P2_new,n*dt,n,dt))
            #print P1_new, P2_new
            #exit()
                
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)          
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
        # apply Areas
        self.A_mother[n+1][pos1]   = A1n
        self.A_daughter[n+1][pos2] = A2n
              
        # Error estimation
        try: sumQError = abs(abs(Q1_new)-abs(Q2_new))/abs(Q1_new)
        except: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name, ' Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        if sumQError > 1.e-5:
            raise ValueError("Connection: {} too high error, sumQError = {}, in conservation of mass at time {} (n {},dt {})".format(self.name,sumQError,n*dt,n,dt))
            #print sumQError
            #exit()
            
        # Linear presssure equation
        sumPError = abs(abs(P1_new)-abs(P2_new))/abs(P1_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        if sumPError > 1.e-5:
            raise ValueError("Connection: {} too high error, sumPError = {}, in conservation of pressure at time {} (n {},dt {})".format(self.name,sumPError,n*dt,n,dt))
            #print sumPError
            #exit()
                                               
        ## Non linear pressure equation
        sumPErrorNonLin = abs(abs(P1_new+1000*0.5*(Q1_new/A1n)**2)-abs(abs(P2_new+1000*0.5*(Q2_new/A2n)**2)))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin
        # TODO: This was set to test sumPError, not sumPErrorNonLin
        # TODO: When set correctly, program wouldn't run. Commented out.
#        if sumPErrorNonLin > 1.e-10:
#           raise ValueError("Connection: {} too high error, sumPErrorNonLin = {} in conservation of pressure at time {} (n {},dt {})".format(self.name,sumPErrorNonLin,n*dt,n,dt))
            #print sumPError
            #exit()
            
        

    def callNonLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        #if n == 1:
            #print "using nonlinear link model"
        P1 = self.P_mother[n]
        Q1 = self.Q_mother[n]
        A1 = self.A_mother[n]
        
        P2 = self.P_daughter[n]
        Q2 = self.Q_daughter[n]
        A2 = self.A_daughter[n]
        
        rho1 = self.rho[0]
        rho2 = self.rho[1]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        A1o = A1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        A2o = A2[pos2]
        # update system equation and store L1
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_2 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
                

        
        # apply calculated values to next time step
                
        domega1_2_init = -R1_11*domega1_1/R1_12 #result if pressure is constant- could be optimized
        domega2_1_init = -R2_12*domega2_2/R2_11 #result if pressure is constant
        
 
        epsilonvalues = np.array([1,1])
        epsilonlimit = 1e-10
        domega1_2_last = domega1_2_init
        domega2_1_last = domega2_1_init

        Niterations = 0
        Xold = np.array([domega1_2_last, domega2_1_last])

        while epsilonvalues[0]>1e-14 or epsilonvalues[1]>0.0001 :
            """iterative Newton Rahpson solver
                domega1_2, domega2_1, domega3_1 are the unknowns that are changing from
            each iteration. A1, A2, A3 are also changing, but the value from the previous iteration is used. (should actually be for the next timestep, but this should converge to the correct solution)
            R1_11, R1_12, R1_21, R1_22, ,R2_11, R2_12, R2_21, R2_22, ,R3_11, R3_12, R3_21, R3_22, domega1_1, domega2_2, domega3_2, Q1o, Q2o, Q2o 
            are constant for each timestep, and are taken from previous timestep. domega1_1, domega2_2, domega3_2 are the field domegas.
            """
            domega1_2_last = Xold[0]
            domega2_1_last = Xold[1]

            
            Q1discretize = Q1o + R1_21*domega1_1 + R1_22*domega1_2_last
            Q2discretize = Q2o + R2_21*domega2_1_last + R2_22*domega2_2

            
            P1discretize = P1o + R1_11*domega1_1 + R1_12*domega1_2_last
            P2discretize = P2o + R2_11*domega2_1_last + R2_12*domega2_2
            
            if self.rigidAreas == False:
                A1_last = self.A_func[0]([P1discretize],pos1)
                A2_last = self.A_func[1]([P2discretize],pos2)         
            else:
                A1_last = A1o
                A2_last = A2o
            
            f1 = Q1discretize - Q2discretize#R1_21*domega1_1 + R1_22*domega1_2_last  - R2_21*domega2_1_last - R2_22*domega2_2 - R3_21*domega3_1_last - R3_22*domega3_2
            
            f2 = P1discretize + 0.5*rho1*(Q1discretize/A1_last)**2 - P2discretize - 0.5*rho2*(Q2discretize/A2_last)**2
            

            
            F = np.array([f1, f2])
            """Inverse Jacobi elements: """
            a = R1_22
            b = -R2_21

            d = R1_12 + rho1*R1_22*(Q1discretize)/(A1_last**2)
            e = - R2_11 - rho2*R2_21*(Q2discretize)/(A2_last**2)
            
            Determinant = a*e -b*d
            
            J_inv = np.array([[ e, -b],
                          [ -d, a]]) / (Determinant)
            
            Xnew = Xold - np.dot(J_inv,F)

            epsilonvalues = np.abs(F)
            Niterations = Niterations + 1
            
            
            if Niterations > 30:
                
                print "Niterations excedded in link calculation in vessel, Niterations: ", self.names[0], Niterations
                print "f1,f2: ", f1, f2
                
                break
            Xold = Xnew
            #exit()
        
        
        Q1_new = Q1discretize
        Q2_new = Q2discretize
        
        
        P1_new = P1discretize
        P2_new = P2discretize
        
        self.P_mother[n+1][pos1]   = P1_new
        self.Q_mother[n+1][pos1]   = Q1_new
        self.P_daughter[n+1][pos2] = P2_new
        self.Q_daughter[n+1][pos2] = Q2_new
        
        if P1_new < 0 or P2_new < 0:
            raise ValueError("Connection: {} calculated negative pressure, P1_new = {}, P2_new = {}, at time {} (n {},dt {})".format(self.name,P1_new,P2_new,n*dt,n,dt))
            #print P1_new, P2_new
            #exit()
                 
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)          
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
        # apply Areas
        self.A_mother[n+1][pos1]   = A1n
        self.A_daughter[n+1][pos2] = A2n
        

    def callStenosisYoungAndTsai(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        #if n == 1:
            #print "using nonlinear link model"
        P1 = self.P_mother[n]
        Q1 = self.Q_mother[n]
        A1 = self.A_mother[n]
        
        P2 = self.P_daughter[n]
        Q2 = self.Q_daughter[n]
        A2 = self.A_daughter[n]
        
        rho1 = self.rho[0]
        rho2 = self.rho[1]
        
        my1 = self.my[0]
        my2 = self.my[1]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        A1o = A1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        A2o = A2[pos2]
        # update system equation and store L1
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_2 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        
        # apply calculated values to next time step
                
        domega1_2_init = -R1_11*domega1_1/R1_12 #result if pressure is constant- could be optimized
        domega2_1_init = -R2_12*domega2_2/R2_11 #result if pressure is constant
        
 
        epsilonvalues = np.array([1,1])
        epsilonlimit = 1e-10
        domega1_2_last = domega1_2_init
        domega2_1_last = domega2_1_init

        Niterations = 0
        Xold = np.array([domega1_2_last, domega2_1_last])
#         print "\n"
#         print "domega1_2_last, domega2_1_last, domega3_1_last, domega1_1, domega2_2, domega3_2 : ", domega1_2_last, domega2_1_last, domega3_1_last, domega1_1, domega2_2, domega3_2
#         print "\n"
        while epsilonvalues[0]>1e-14 or epsilonvalues[1]>0.0001 :
            """iterative Newton Rahpson solver
                domega1_2, domega2_1, domega3_1 are the unknowns that are changing from
            each iteration. A1, A2, A3 are also changing, but the value from the previous iteration is used. (should actually be for the next timestep, but this should converge to the correct solution)
            R1_11, R1_12, R1_21, R1_22, ,R2_11, R2_12, R2_21, R2_22, ,R3_11, R3_12, R3_21, R3_22, domega1_1, domega2_2, domega3_2, Q1o, Q2o, Q2o 
            are constant for each timestep, and are taken from previous timestep. domega1_1, domega2_2, domega3_2 are the field domegas.
            """
            domega1_2_last = Xold[0]
            domega2_1_last = Xold[1]

            
            Q1discretize = Q1o + R1_21*domega1_1 + R1_22*domega1_2_last
            Q2discretize = Q2o + R2_21*domega2_1_last + R2_22*domega2_2

            
            P1discretize = P1o + R1_11*domega1_1 + R1_12*domega1_2_last
            P2discretize = P2o + R2_11*domega2_1_last + R2_12*domega2_2
            
            if self.rigidAreas == False:
                A1_last = self.A_func[0]([P1discretize],pos1)
                A2_last = self.A_func[1]([P2discretize],pos2)         
            else:
                A1_last = A1o
                A2_last = A2o
            
            f1 = Q1discretize - Q2discretize#R1_21*domega1_1 + R1_22*domega1_2_last  - R2_21*domega2_1_last - R2_22*domega2_2 - R3_21*domega3_1_last - R3_22*domega3_2
            
            U1 = Q1discretize/A1_last
            U2 = Q2discretize/A2_last
            D1 = np.sqrt(4*A1_last/np.pi)
            Re1 = rho1*U1*D1/my1
            
            Astenosis = A1_last*0.56
            
            Kv = 500. # pouseille friction correction coeficient
            Kt = 0.9 # expansion coefficient
            if abs(Re1<0.01):
                deltaPstenosis = 0
            else:
                deltaPstenosis = rho1*U1**2*(Kv/Re1 + (Kt/2.)*(A1_last/Astenosis - 1)**2)
            
            
            f2 = P1discretize + 0.5*rho1*(Q1discretize/A1_last)**2 - P2discretize - 0.5*rho2*(Q2discretize/A2_last)**2 - deltaPstenosis
            

            
            F = np.array([f1, f2])
            """Inverse Jacobi elements: """
            a = R1_22
            b = -R2_21

            d = R1_12 + rho1*R1_22*(Q1discretize)/(A1_last**2)
            e = - R2_11 - rho2*R2_21*(Q2discretize)/(A2_last**2)
            
            Determinant = a*e -b*d
            
            J_inv = np.array([[ e, -b],
                          [ -d, a]]) / (Determinant)
            
            Xnew = Xold - np.dot(J_inv,F)

            epsilonvalues = np.abs(F)
            Niterations = Niterations + 1
            
            
            if Niterations > 30:
                
                print "Niterations excedded in link calculation in vessel, Niterations: ", self.names[0], Niterations
                print "f1,f2: ", f1, f2
                
                break
            Xold = Xnew
            #exit()

        
        
        Q1_new = Q1discretize
        Q2_new = Q2discretize
        
        
        P1_new = P1discretize
        P2_new = P2discretize
        
        self.P_mother[n+1][pos1]   = P1_new
        self.Q_mother[n+1][pos1]   = Q1_new
        self.P_daughter[n+1][pos2] = P2_new
        self.Q_daughter[n+1][pos2] = Q2_new
        
        if P1_new < 0 or P2_new < 0:
            raise ValueError("Connection: {} calculated negative pressure, P1_new = {}, P2_new = {}, at time {} (n {},dt {})".format(self.name,P1_new,P2_new,n*dt,n,dt))
            #print P1_new, P2_new
            #exit()
                 
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)          
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
        # apply Areas
        self.A_mother[n+1][pos1]   = A1n
        self.A_daughter[n+1][pos2] = A2n
              
#     def callMacCormackField1(self):
#         """
#         Call function for vessel-vessel connection
#         """        
#         dt = self.dt
#         n = self.currentMemoryIndex[0]
#         pos1 = self.positions[0]
#         pos2 = self.positions[1]
#         
#         #""" Predictor Step """positions
#         if self.step == "predictor":
#             P1 = self.P_leftMother[n]
#             Q1 = self.Q_leftMother[n]
#             A1 = self.A_leftMother[n]
#             
#             P2 = self.P_daughter[n]
#             Q2 = self.Q_daughter[n]
#             A2 = self.A_daughter[n]
#                  
#         #"""Corrector Step"""    A
#         elif self.step == "corrector":
#             P1 = self.P_mother_pre
#             Q1 = self.Q_mother_pre
#             A1 = self.A_mother_pre
#             
#             P2 = self.P_daughter_pre
#             Q2 = self.Q_daughter_pre
#             A2 = self.A_daughter_pre
#                  
#         P1o = P1[pos1]
#         Q1o = Q1[pos1]
#         P2o = P2[pos2]
#         Q2o = Q2[pos2]
#         
#         # update system equation and store L1o
#         self.systemEquations[0].updateLARL(P1,Q1,A1,idArray=[pos1],update='L') #
#         L1o = self.systemEquations[0].L[pos1][pos1+1]
#         # calculate domega1
#         z1 = self.z[0][pos1] - self.systemEquations[0].LAMBDA[pos1][0] * dt
#         du1 = np.array([np.interp(z1,self.z[0],P1)-P1o,np.interp(z1,self.z[0],Q1)-Q1o])
#         domega1 =  np.dot(L1o,du1)
#         
#         # update system equation and store L2o
#         self.systemEquations[1].updateLARL(P2,Q2,A2,idArray=[pos2],update='L') #
#         L2o = self.systemEquations[1].L[pos2][pos2+1]
#         # calculate domega2
#         z2 = self.z[1][pos2] - self.systemEquations[1].LAMBDA[pos2][1] * dt
#         du2 = np.array([np.interp(z2,self.z[1],P2)-P2o,np.interp(z2,self.z[1],Q2)-Q2o])
#         domega2 =  np.dot(L2o,du2)
#         
#         # setup solve function
#         args = [A1[pos1],A2[pos2],pos1,pos2,self.vz[0],self.vz[1],P1o,Q1o,P2o,Q2o,domega1,domega2,self.rho[0],self.rho[1],L1o,L2o,du1,du2] 
#         x = [P1o,Q1o,Q2o,P2o]
#         sol = fsolve(self.fsolveFunction ,x ,args = args, fprime = self.jacobiMatrix)
#         
#         #sol,infodict,a,b = fsolve(self.fsolveFunction ,x ,args = args, fprime = self.jacobiMatrix,full_output=True)
#         #print infodict['nfev'],infodict['njev']
#         
#         #print sol[0],sol[3]
#         
#         #""" Predictor Step """
#         if self.step == "predictor":
#             self.step = "corrector"   
#             # apply changed values
#             self.P_mother_pre[pos1]   = sol[0]
#             self.Q_mother_pre[pos1]   = sol[1]
#             self.P_daughter_pre[pos2] = sol[3]
#             self.Q_daughter_pre[pos2] = sol[2]
#             
#             if self.rigidAreas == False:
#                 # calculate new areas
#                 A1n = self.A_func[0]([sol[0]],pos1)
#                 A2n = self.A_func[1]([sol[3]],pos2)
#                 self.A_mother_pre[pos1]   = A1n
#                 self.A_daughter_pre[pos2] = A2n             
#             else:
#                 self.A_mother_pre[pos1]   = A1[pos1]
#                 self.A_daughter_pre[pos2] = A2[pos2]
#         
#         #"""Corrector Step"""    
#         elif self.step == "corrector":
#             self.step = "predictor"
#             
#             # apply changed values
#             self.P_leftMother[n+1][pos1]   = sol[0]
#             self.Q_leftMother[n+1][pos1]   = sol[1]
#             self.P_daughter[n+1][pos2] = sol[3]
#             self.Q_daughter[n+1][pos2] = sol[2]
#             
#             sumQError = abs(sol[1] - sol[2])
#             if sumQError > 0.0: 
#                 self.sumQErrorCount = self.sumQErrorCount+1
#             if sumQError > self.maxQError:
#                 self.maxQError  = sumQError
#             print 'Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
#             
#             if self.rigidAreas == False:
#                 # calculate new areas
#                 A1n = self.A_func[0]([sol[0]],pos1)
#                 A2n = self.A_func[1]([sol[3]],pos2)
#                 self.A_leftMother[n+1][pos1]   = A1n
#                 self.A_daughter[n+1][pos2] = A2n             
#             else:
#                 self.A_leftMother[n+1][pos1]   = A1[pos1]
#                 self.A_daughter[n+1][pos2] = A2[pos2]
#     
#     def fsolveConnectionSys0(self,x,args):
#         """
#         Residual Function with equations to solve for at the Link
#         Using constant areas, i.e. initial areas
#         
#         Input:     x = array [P1,Q1,Q2,P2]
#                    args = args with local variables
#         Returns array with residuals 
#         """
#         P1,Q1,Q2,P2 = x
#         A1,A2,pos1,pos2,vz1,vz2,P1o,Q1o,P2o,Q2o,domega1,domega2,rho1,rho2,L1,L2,du1,du2 = args
#         
#         du1[0] = P1 - P1o 
#         du1[1] = Q1 - Q1o
#         du2[0] = P2 - P2o
#         du2[1] = Q2 - Q2o
#         
#         
#         self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[pos1],update='L') 
#         self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[pos2],update='L') 
#         
#         #calculate residuals
#         res1 = vz1*Q1+vz2*Q2
#         res2 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz2*P2+vz2*rho2*0.5*(Q2/A2)**2.
#         res3 = np.dot(self.systemEquations[0].L[pos1][pos1+1],du1) - domega1
#         res4 = np.dot(self.systemEquations[1].L[pos2][pos2+1],du2) - domega2
#         
#         return [res3,res2,res1,res4]
#     
#     def jacobiMatrixSys0(self,x,args):
#         """
#         Returns the jabcobi matrix, bifurcation-functions and x; J = dF/dx
#         Using constant areas, i.e. initial areas
#         """
#         P1,Q1,Q2,P2 = x
#         A1,A2,pos1,pos2,vz1,vz2,P1o,Q1o,P2o,Q2o,domega1,domega2,rho1,rho2,L1,L2,du1,du2 = args
#         
#         return np.array([[L1[0], L1[1]            , 0                , 0    ],
#                          [vz1  , vz1*rho1*Q1/A1**2, vz2*rho2*Q2/A2**2, vz2  ],
#                          [0    , vz1              , vz2              , 0    ],
#                          [0    , 0                , L2[1]            , L2[0]]])
#          
#         
#     def fsolveConnectionSys1(self,x,args):
#         """
#         Residual Function with equations to solve for at the Link
#         Using recalculated areas depending on the new pressure values
#         
#         Input:     x = array [P1,Q1,Q2,P2]
#                    args = args with local variables
#         Returns array with residuals 
#         """
#         P1,Q1,Q2,P2 = x
#         A1,A2,pos1,pos2,vz1,vz2,P1o,Q1o,P2o,Q2o,domega1,domega2,rho1,rho2,L1,L2,du1,du2 = args
#                         
#         A1 = self.A_func[pos1]([P1],pos1)
#         A2 = self.A_func[pos2]([P2],pos2)
#         
#         du1[0] = P1 - P1o 
#         du1[1] = Q1 - Q1o
#         du2[0] = P2 - P2o
#         du2[1] = Q2 - Q2o 
#         
#         if self.updateL == True:      
#             self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[pos1],update='L')
#             self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[pos2],update='L')
#         
#         
#         #calculate residuals
#         res1 = vz1*Q1+vz2*Q2
#         
#         if self.nonLin == True:
#             # non linear pressure equation
#             res2 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz2*P2+vz2*rho2*0.5*(Q2/A2)**2.
#         else:
#             # linear pressure equation p1 - p2 = 0
#             res2 = vz1*P1+vz2*P2
#         
#         res3 = np.dot(self.systemEquations[0].L[pos1][pos1+1],du1) - domega1
#         res4 = np.dot(self.systemEquations[1].L[pos2][pos2+1],du2) - domega2
#         
#         return [res3,res2,res1,res4]
#  
#     def jacobiMatrixSys1(self,x, args):
#         """
#         Returns the jabcobi matrix, bifurcation-functions and x; J = dF/dx
#         Using recalculated areas depending on the new pressure values
#         """
#         P1,Q1,Q2,P2 = x
#         A1,A2,pos1,pos2,vz1,vz2,P1o,Q1o,P2o,Q2o,domega1,domega2,rho1,rho2,L1,L2,du1,du2 = args
#         
#         
#         A1 = self.A_func[pos1]([P1],pos1)
#         A2 = self.A_func[pos2]([P2],pos2)
#         
#         if self.nonLin == True:
#             # J non linear pressure equ.
#             J =    np.array([[L1[0], L1[1]            , 0                , 0    ],
#                              [vz1  , vz1*rho1*Q1/A1**2, vz2*rho2*Q2/A2**2, vz2  ],
#                              [0    , vz1              , vz2              , 0    ],
#                              [0    , 0                , L2[1]            , L2[0]]])
#         else:
#             # J linear pressure equation
#             J =    np.array([[L1[0], L1[1]            , 0                , 0    ],
#                              [vz1  , 0                , 0                , vz2  ],
#                              [0    , vz1              , vz2              , 0    ],
#                              [0    , 0                , L2[1]            , L2[0]]])
#             
#         return J
#     
#     
#    def fsolveMultibranchingSys0(self,x):
#        
#        #self.count = self.count+1
#        
#        out = [0 for i in x]
#        #outt = ['' for i in x]
#        for i in range(self.number):
#           
#            Ai = self.A[i]
#            out[0]=out[0] + self.vz[i]*x[i*2+1]
#            if i != 0: out[i]= (x[0]+(self.rho[0]/2.)*(x[1]/self.A[0])**2.) - (x[i*2]+(self.rho[i]/2.)*(x[i*2+1]/Ai)**2.)
#                                 
#            du_t = np.array([x[i*2],x[i*2+1]])-np.array([self.P[i],self.Q[i]])  
#            self.systemEquations[i].updateLARL([x[i*2]],[x[i*2+1]],[Ai],idArray=[self.positions[i]],update='L') 
#            
#            out[i+self.number]= np.dot(self.systemEquations[i].L[self.positions[i]][self.positions[i]+1],du_t)/self.domega[i] - 1
#            
#            #outt[0] = outt[0]+str('mass eq'+str(i)+'+')
#            #if i != 0: outt[i]= 'P0 - P'+str(i)
#            #outt[i+self.number] = 'domega '+str(i)
#            
#        #if self.count == 8:print outt
#        return out
    
    
class Bifurcation():
    
    def __init__(self, mother, motherSys,
                       leftDaughter, leftDaughterSys,
                       rightDaughter, rightDaughterSys, 
                       currentMemoryIndex, dt, rigidAreas, solvingScheme):
        # vessel variables initially set, constant through simulation
        self.type = 'Bifurcation'
        
        self.name = ' '.join(['Bifurcation',str(mother.Id),str(leftDaughter.Id),str(rightDaughter.Id)])
        
        #System Variables
        self.dt = dt
        self.currentMemoryIndex = currentMemoryIndex
        
        self.rho = []
        self.systemEquations = []
        self.z = []
        self.A_func = []
        self.positions =[]
#         self.vz = []
        self.names = []
        
#         # equations to solve in f solve
#         self.fsolveFunction = None
#         self.jacobiMatrix = None
        
        ###initialize
        ##mother branch
        self.rho.append(mother.rho)
        self.z.append(mother.z)
        self.systemEquations.append(motherSys)
        self.positions.append(-1)
        self.names.append(mother.Id)
        self.A_func.append(mother.A_nID)
#        self.vz.append(-1)
        #SolutionVariables
        self.P_leftMother = mother.Psol
        self.Q_leftMother = mother.Qsol
        self.A_leftMother = mother.Asol
              
        
        ##left daughter
        self.rho.append(leftDaughter.rho)
        self.z.append(leftDaughter.z)
        self.systemEquations.append(leftDaughterSys)
        self.positions.append(0)
        self.names.append(leftDaughter.Id)
        self.A_func.append(leftDaughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_leftDaughter = leftDaughter.Psol
        self.Q_leftDaughter = leftDaughter.Qsol
        self.A_leftDaughter = leftDaughter.Asol
        
        ##right daughter
        self.rho.append(rightDaughter.rho)
        self.z.append(rightDaughter.z)
        self.systemEquations.append(rightDaughterSys)
        self.positions.append(0)
        self.names.append(rightDaughter.Id)
        self.A_func.append(rightDaughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_rightDaughter = rightDaughter.Psol
        self.Q_rightDaughter = rightDaughter.Qsol
        self.A_rightDaughter = rightDaughter.Asol
        
#         if   rigidAreas == False:
#             self.fsolveFunction = self.fsolveBifurcationSys0
#             self.jacobiMatrix = self.jacobiMatrixBifSys0
#         elif rigidAreas == False:
#             self.fsolveFunction = self.fsolveBifurcationSys1
#             self.jacobiMatrix = self.jacobiMatrixBifSys1
#         else: print "ERROR classConnections: EquSys not properly defined!";exit()

        self.rigidAreas = rigidAreas
        solvingScheme = "NonLinear"
        # Define the call function depending on the solving Scheme
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        elif solvingScheme == "NonLinear":
            self.__call__ = self.callNonLinear
        else:
            raise ImportError("Connections wrong solving scheme! {}".format(solvingScheme))
        
        ## benchamark Test variables
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_leftDaughter[n]
        Q2 = self.Q_leftDaughter[n]
        A2 = self.A_leftDaughter[n]
        
        P3 = self.P_rightDaughter[n]
        Q3 = self.Q_rightDaughter[n]
        A3 = self.A_rightDaughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
        
        ## update LARL
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_2 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
        L,R3,LMBD,Z1,Z2,domega3_2 = self.systemEquations[2].updateLARL(P3,Q3,A3,pos3)
        
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        ###### Linear approach
        denom = R1_12 * R2_11 * R3_21 + R1_12 * R2_21 * R3_11 - R1_22 * R2_11 * R3_11 
        
        alpha1 = -(R1_11 * R2_11 * R3_21 + R1_11 * R2_21 * R3_11 - R1_21 * R2_11 * R3_11)/denom
        alpha2 = -(R2_11 * R2_22 * R3_11 - R2_12 * R2_21 * R3_11)/denom
        alpha3 = -(R2_11 * R3_22 * R3_11 - R2_11 * R3_12 * R3_21)/denom
        
        beta1 = -(R1_11 * R1_22 * R3_11 - R1_12 * R1_21 * R3_11)/denom
        beta2 = -(R1_12 * R2_12 * R3_21 + R1_12 * R2_22 * R3_11 - R1_22 * R2_12 * R3_11)/denom
        beta3 = -(R1_12 * R3_11 * R3_22 - R1_12 * R3_12 * R3_21)/denom
        
        gamma1 = -(R1_11 * R1_22 * R2_11 - R1_12 * R1_21 * R2_11)/denom
        gamma2 = -(R1_12 * R2_11 * R2_22 - R1_12 * R2_12 * R2_21)/denom
        gamma3 = -(R1_12 * R2_11 * R3_22 + R1_12 * R2_21 * R3_12 - R1_22 * R2_11 * R3_12)/denom
        
        domega1_2 = alpha1 * domega1_1  + alpha2 * domega2_2 + alpha3 * domega3_2 
        domega2_1 = beta1  * domega1_1  + beta2  * domega2_2 + beta3  * domega3_2 
        domega3_1 = gamma1 * domega1_1  + gamma2 * domega2_2 + gamma3 * domega3_2 
        
        P1_new = P1o + (R1_11*domega1_1 + R1_12*domega1_2)
        Q1_new = Q1o + (R1_21*domega1_1 + R1_22*domega1_2)
    
        P2_new = P2o + (R2_11*domega2_1 + R2_12*domega2_2)
        Q2_new = Q2o + (R2_21*domega2_1 + R2_22*domega2_2)
        
        P3_new = P3o + (R3_11*domega3_1 + R3_12*domega3_2)
        Q3_new = Q3o + (R3_21*domega3_1 + R3_22*domega3_2)
               
        
        # apply calculated values to next time step
        self.P_leftMother[n+1][pos1]    = P1_new
        self.Q_leftMother[n+1][pos1]    = Q1_new
        self.P_leftDaughter[n+1][pos2]  = P2_new
        self.Q_leftDaughter[n+1][pos2]  = Q2_new
        self.P_rightDaughter[n+1][pos3] = P3_new
        self.Q_rightDaughter[n+1][pos3] = Q3_new
        print "using linear bifurcation model"
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            raise ValueError("Connection: {} calculated negative pressure, P1_new = {}, P2_new = {}, P3_new = {}, at time {} (n {},dt {})".format(self.name,P1_new, P2_new, P3_new, n*dt,n,dt))
            #print P1_new, P2_new, P3_new
            #exit()
        
        
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]       = A1n
        self.A_leftDaughter[n+1][pos2]     = A2n       
        self.A_rightDaughter[n+1][pos3]    = A3n  
            
        ## non linear error        
        try: sumQError = abs(Q1_new-Q2_new-Q3_new)/abs(Q1_new)
        except Exception: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name,' \n Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        if sumQError > 1.e-5:
            raise ValueError("Connection: {} too high error, sumQError = {} in conservation of mass at time {} (n {},dt {})".format(self.name,sumQError,n*dt,n,dt))
            #print sumQError
            #exit()
        
        sumPError = abs(P1_new-P2_new)/abs(P1_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        #print self.name,' Error P lin    ',  sumPError, self.maxPError ,' - ', n, self.sumPErrorCount
        if sumPError > 1.e-10:
            raise ValueError("Connection: {} too high error, sumPError = {}, in conservation of pressure at time {} (n {},dt {}), exit system".format(self.name,sumPError,n*dt,n,dt))
            #print sumPError
            #exit()
        
        sumPErrorNonLin = abs(P1_new+500*(Q1_new/A1n)**2-(P2_new+500*(Q2_new/A2n)**2))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin


    def callNonLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        #if n == 1:
            #print "using nonlinear bifurcation model"
        #print "using nonlinear bifurcation model"
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_leftDaughter[n]
        Q2 = self.Q_leftDaughter[n]
        A2 = self.A_leftDaughter[n]
        
        P3 = self.P_rightDaughter[n]
        Q3 = self.Q_rightDaughter[n]
        A3 = self.A_rightDaughter[n]
        
        rho1 = self.rho[0]
        rho2 = self.rho[1]
        rho3 = self.rho[2]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        A1o = A1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        A2o = A2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
        A3o = A3[pos3]

        ## update LARL
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_2 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
        L,R3,LMBD,Z1,Z2,domega3_2 = self.systemEquations[2].updateLARL(P3,Q3,A3,pos3)
        
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        if n>0:
            P1o2 = self.P_leftMother[n-1][pos1]
            P2o2 = self.P_leftDaughter[n-1][pos2]
            P3o2 = self.P_rightDaughter[n-1][pos3]
            
            deltaP1_last = P1o - P1o2
            deltaP2_last = P2o - P2o2
            deltaP3_last = P3o - P3o2
            
            domega1_2_init = (deltaP1_last - R1_11*domega1_1)/R1_12 #same omega as previous timestep
            domega2_1_init = (deltaP2_last - R2_12*domega2_2)/R2_11 #same omega as previous timestep
            domega3_1_init = (deltaP3_last - R3_12*domega3_2)/R3_11 #same omega as previous timestep
        else:
            domega1_2_init = 0#-R1_11*domega1_1/R1_12 #result if pressure is constant- could be optimized
            domega2_1_init = 0#-R2_12*domega2_2/R2_11 #result if pressure is constant
            domega3_1_init = 0#-R3_12*domega3_2/R3_11 #result if pressure is constant

        epsilonvalues = np.array([1,1,1])
        epsilonlimit = 1e-10
        domega1_2_last = domega1_2_init
        domega2_1_last = domega2_1_init
        domega3_1_last = domega3_1_init
        A1_last = A1o
        A2_last = A2o
        A3_last = A3o
        Niterations = 0
        Xold = np.array([domega1_2_last, domega2_1_last, domega3_1_last])

        while epsilonvalues[0]>1e-14 or epsilonvalues[1]>0.001 or epsilonvalues[2]>0.001:
            """iterative Newton Rahpson solver
                domega1_2, domega2_1, domega3_1 are the unknowns that are changing from
            each iteration. A1, A2, A3 are also changing, but the value from the previous iteration is used. (should actually be for the next timestep, but this should converge to the correct solution)
            R1_11, R1_12, R1_21, R1_22, ,R2_11, R2_12, R2_21, R2_22, ,R3_11, R3_12, R3_21, R3_22, domega1_1, domega2_2, domega3_2, Q1o, Q2o, Q2o 
            are constant for each timestep, and are taken from previous timestep. domega1_1, domega2_2, domega3_2 are the field domegas.
            """
            domega1_2_last = Xold[0]
            domega2_1_last = Xold[1]
            domega3_1_last = Xold[2]
            
            Q1discretize = Q1o + R1_21*domega1_1 + R1_22*domega1_2_last
            Q2discretize = Q2o + R2_21*domega2_1_last + R2_22*domega2_2
            Q3discretize = Q3o + R3_21*domega3_1_last + R3_22*domega3_2
            
            P1discretize = P1o + R1_11*domega1_1 + R1_12*domega1_2_last
            P2discretize = P2o + R2_11*domega2_1_last + R2_12*domega2_2
            P3discretize = P3o + R3_11*domega3_1_last + R3_12*domega3_2
            

            if self.rigidAreas == False:
                try:
                    A1_last = self.A_func[0]([P1discretize],pos1)
                    A2_last = self.A_func[1]([P2discretize],pos2)
                    A3_last = self.A_func[2]([P3discretize],pos3)
                except FloatingPointError as E:
                    print "Floating Point error in Connection {}".format(self.name)
                    raise E
            else:
                A1_last = A1[pos1]
                A2_last = A2[pos2]
                A3_last = A3[pos3] 

            
            f1 = Q1discretize - Q2discretize - Q3discretize#R1_21*domega1_1 + R1_22*domega1_2_last  - R2_21*domega2_1_last - R2_22*domega2_2 - R3_21*domega3_1_last - R3_22*domega3_2
            
            f2 = P1discretize + 0.5*rho1*((Q1discretize/A1_last)**2) - P2discretize - 0.5*rho2*((Q2discretize/A2_last)**2)
            
            f3 = P1discretize + 0.5*rho1*((Q1discretize/A1_last)**2) - P3discretize - 0.5*rho3*((Q3discretize/A3_last)**2)
            
            F = np.array([f1, f2, f3])
            """Inverse Jacobi elements: """
            a = R1_22
            b = -R2_21
            c = -R3_21
            d = R1_12 + rho1*R1_22*(Q1discretize)/(A1_last**2)
            e = - R2_11 - rho2*R2_21*(Q2discretize)/(A2_last**2)
            f = d
            g = - R3_11 - rho3*R3_21*(Q3discretize)/(A3_last**2)
            
            Determinant = a*e*g -b*d*g -c*e*f
            
            J_inv = np.array([[ e*g, -b*g, -c*e  ],
                          [ -d*g, a*g-c*f, c*d ],
                          [  -e*f, b*f, a*e-b*d ]]) / (Determinant)
            
            Xnew = Xold - np.dot(J_inv,F)

            epsilonvalues = np.abs(F)
            Niterations = Niterations + 1
            
            if Niterations > 30:
                print "\n"
                print "Niterations excedded in Bifurcation calculation in vessel, Niterations: ", self.names[0], Niterations

                print "Xnew: ", Xnew
                print "Xold: ", Xold
                
                print "f1: ", f1
                print "f2: ", f2
                print "f3: ", f3
                print "epsilonvalues: ", epsilonvalues
                print "Q1discretize, Q1o: ", Q1discretize, Q1o
                print "Q2discretize, Q2o: ", Q2discretize, Q2o
                print "Q3discretize, Q3o: ", Q3discretize, Q3o
                print "P1discretize, P1o: ", P1discretize, P1o
                print "P2discretize, P2o: ", P2discretize, P2o
                print "P3discretize, P3o: ",P3discretize, P3o
                
                break
            Xold = Xnew
            #exit()
        
        
        Q1_new = Q1discretize
        Q2_new = Q2discretize
        Q3_new = Q3discretize
        
        P1_new = P1discretize
        P2_new = P2discretize
        P3_new = P3discretize
        
        # apply calculated values to next time step
        self.P_leftMother[n+1][pos1]    = P1_new
        self.Q_leftMother[n+1][pos1]    = Q1_new
        self.P_leftDaughter[n+1][pos2]  = P2_new
        self.Q_leftDaughter[n+1][pos2]  = Q2_new
        self.P_rightDaughter[n+1][pos3] = P3_new
        self.Q_rightDaughter[n+1][pos3] = Q3_new
        
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            print "ERROR: Connection: {} calculated negative pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P1_new, P2_new, P3_new
            print "Niterations: ", Niterations
            
            print "solving nonlinear/total pressure Bifurcation"
            print "domega1_2_init: ", domega1_2_init
            print "domega1_1: ", domega1_1
                    
            print "f1: ", f1
            print "f2: ", f2
            print "f3: ", f3
            print "epsilonvalues: ", epsilonvalues
            print "Q1discretize, Q1o: ", Q1discretize, Q1o
            print "Q2discretize, Q2o: ", Q2discretize, Q2o
            print "Q3discretize, Q3o: ", Q3discretize, Q3o
            print "P1discretize, P1o: ", P1discretize, P1o
            print "P2discretize, P2o: ", P2discretize, P2o
            print "P3discretize, P3o: ",P3discretize, P3o
            exit()
        
        
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]       = A1n
        self.A_leftDaughter[n+1][pos2]     = A2n       
        self.A_rightDaughter[n+1][pos3]    = A3n 
        #print self.name,' Error P non lin',  sumPErrorNonLin, self.maxPErrorNonLin ,' - ', n, self.sumPErrorNonLinCount
        
        
        
                
#     def callMacCormackField2(self):
#         """
#         Call function for a bifurcation
#         """  
#         #print self.counter
#         #self.counter = 1+ self.counter
#         
#         
#         dt = self.dt
#         n = self.currentMemoryIndex[0]  
#         pos1 = self.positions[0]
#         pos2 = self.positions[1]
#         pos3 = self.positions[2]
#         
#         #""" Predictor Step """positions
#         if self.step == "predictor":
#             P1 = self.P_leftMother[n]
#             Q1 = self.Q_leftMother[n]
#             A1 = self.A_leftMother[n]
#             
#             P2 = self.P_rightMother[n]
#             Q2 = self.Q_rightMother[n]
#             A2 = self.A_rightMother[n]
#             
#             P3 = self.P_daughter[n]
#             Q3 = self.Q_daughter[n]
#             A3 = self.A_daughter[n]
#                  
#         #"""Corrector Step"""    
#         elif self.step == "corrector":
#             P1 = self.P_mother_pre
#             Q1 = self.Q_mother_pre
#             A1 = self.A_mother_pre
#             
#             P2 = self.P_leftDaughter_pre
#             Q2 = self.Q_leftDaughter_pre
#             A2 = self.A_leftDaughter_pre
#             
#             P3 = self.P_rightDaughter_pre
#             Q3 = self.Q_rightDaughter_pre
#             A3 = self.A_rightDaughter_pre
#                  
#         P1o = P1[pos1]
#         Q1o = Q1[pos1]
#         P2o = P2[pos2]
#         Q2o = Q2[pos2]
#         P3o = P3[pos3]
#         Q3o = Q3[pos3]
#         
#         # update system equation and store L1o
#         self.systemEquations[0].updateLARL(P1,Q1,A1,idArray=[pos1],update='L') #
#         L1o = self.systemEquations[0].L[pos1][pos1+1]
#         # calculate domega1
#         z1 = self.z[0][pos1] - self.systemEquations[0].LAMBDA[pos1][0] * dt
#         du1 = np.array([np.interp(z1,self.z[0],P1)-P1o,np.interp(z1,self.z[0],Q1)-Q1o])
#         domega1 =  np.dot(L1o,du1)
#         
#         # update system equation and store L2o
#         self.systemEquations[1].updateLARL(P2,Q2,A2,idArray=[pos2],update='L') #
#         L2o = self.systemEquations[1].L[pos2][pos2+1]
#         # calculate domega2
#         z2 = self.z[1][pos2] - self.systemEquations[1].LAMBDA[pos2][1] * dt
#         du2 = np.array([np.interp(z2,self.z[1],P2)-P2o,np.interp(z2,self.z[1],Q2)-Q2o])
#         domega2 =  np.dot(L2o,du2)
#         
#         # update system equation and store L2o
#         self.systemEquations[2].updateLARL(P3,Q3,A3,idArray=[pos3],update='L') #
#         L3o = self.systemEquations[2].L[pos3][pos3+1]
#         # calculate domega2
#         z3 = self.z[2][pos3] - self.systemEquations[2].LAMBDA[pos3][1] * dt
#         du3 = np.array([np.interp(z3,self.z[2],P3)-P3o,np.interp(z3,self.z[2],Q3)-Q3o])
#         domega3 =  np.dot(L3o,du3)
#         
#         # setup solve function
#         args = [A1[pos1],A2[pos2],A3[pos3],pos1,pos2,pos3,self.vz[0],self.vz[1],self.vz[2],P1o,Q1o,P2o,Q2o,P3o,Q3o,domega1,domega2,domega3,self.rho[0],self.rho[1],self.rho[2],L1o,L2o,L3o,du1,du2,du3]      
#         x = [P2o,Q2o,P1o,Q1o,Q3o,P3o]
#         #sol = fsolve(self.fsolveFunction ,x ,args = args, fprime = self.jacobiMatrix)
#         
#         #error = sum([abs(i) for i in self.fsolveFunction(x,args)])
#         
#         #if error < 1.E-4:
#             #print error
#             #return [P1o,Q1o,A[0][pos1],P2o,Q2o,A[1][pos2],P3o,Q3o,A[2][pos3]]
#         # solve system
#         sol,infodict,a,b = fsolve(self.fsolveFunction ,x ,args = args, fprime = self.jacobiMatrix,full_output=True)
#         #print "cC",infodict['nfev'],infodict['njev'], b
#         #if infodict['nfev'] > 2: raw_input("")
#         #""" Predictor Step """
#         if self.step == "predictor":
#             self.step = "corrector"
#
#             # apply changed values
#             self.P_mother_pre[pos1]        = sol[2]
#             self.Q_mother_pre[pos1]        = sol[3]
#             self.P_leftDaughter_pre[pos2]  = sol[0]
#             self.Q_leftDaughter_pre[pos2]  = sol[1]
#             self.P_rightDaughter_pre[pos3] = sol[5]
#             self.Q_rightDaughter_pre[pos3] = sol[4]
#
#             if self.rigidAreas == False:
#                 # calculate new areas
#                 A1n = self.A_func[0]([sol[2]],pos1)
#                 A2n = self.A_func[1]([sol[0]],pos2)
#                 A3n = self.A_func[2]([sol[5]],pos3)
#                 
#                 self.A_mother_pre[pos1]        = A1n
#                 self.A_leftDaughter_pre[pos2]  = A2n
#                 self.A_rightDaughter_pre[pos3] = A3n  
#             else:
#                 self.A_mother_pre[pos1]        = A1[pos1]
#                 self.A_leftDaughter_pre[pos2]  = A2[pos2]
#                 self.A_rightDaughter_pre[pos3] = A3[pos3]
#         
#         #"""Corrector Step"""    
#         elif self.step == "corrector":
#             self.step = "predictor"
#     
#             # apply changed values
#             self.P_leftMother[n+1][pos1]        = sol[2]
#             self.Q_leftMother[n+1][pos1]        = sol[3]
#             self.P_rightMother[n+1][pos2]  = sol[0]
#             self.Q_rightMother[n+1][pos2]  = sol[1]
#             self.P_daughter[n+1][pos3] = sol[5]
#             self.Q_daughter[n+1][pos3] = sol[4]
#             
#             if self.rigidAreas == False:
#                 # calculate new areas
#                 A1n = self.A_func[0]([sol[2]],pos1)
#                 A2n = self.A_func[1]([sol[0]],pos2)
#                 A3n = self.A_func[2]([sol[5]],pos3)
#                 
#                 self.A_leftMother[n+1][pos1]        = A1n
#                 self.A_rightMother[n+1][pos2]  = A2n
#                 self.A_daughter[n+1][pos3] = A3n  
#             else:
#                 self.A_leftMother[n+1][pos1]        = A1[pos1]
#                 self.A_rightMother[n+1][pos2]  = A2[pos2]
#                 self.A_daughter[n+1][pos3] = A3[pos3]
#         
#             print   
#             sumQError = abs(abs(sol[3])-abs(sol[1])-abs(sol[4]))/abs(sol[3])
#             if sumQError > 0.0: 
#                 self.sumQErrorCount = self.sumQErrorCount+1
#             if sumQError > self.maxQError:
#                 self.maxQError  = sumQError
#             print 'Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
#             
#             sumPError = abs(abs(sol[2])-abs(sol[0]))/abs(sol[2])
#             if sumPError > 0.0: 
#                 self.sumPErrorCount = self.sumPErrorCount+1
#             if sumPError > self.maxPError:
#                 self.maxPError  = sumPError
#             print 'Error P lin    ',  sumPError, self.maxPError ,' - ', n, self.sumPErrorCount
#             
#             ## non linear error
#             sumPErrorNonLin = abs((sol[2]+0.5*(sol[3]/A1n)**2)-(sol[0]+0.5*(sol[1]/AL12n)**2))/abs(sol[2]+0.5*(sol[3]/A1n)**2)
#             if sumPErrorNonLin > 0.0: 
#                 self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
#             if sumPErrorNonLin > self.maxPErrorNonLin:
#                 self.maxPErrorNonLin  = sumPErrorNonLin
#             print 'Error P non lin',  sumPErrorNonLin, self.maxPErrorNonLin ,' - ', n, self.sumPErrorNonLinCount
#             
#     
#     def fsolveBifurcationSys0(self,x,args):
#         """
#         Residual Function with equations to solve for at the bifuraction
#         Using constant areas, i.e. initial areas
#         
#         Input:     x = array [P2,Q2,P1,Q1,Q3,P3]
#                    args = args with local variables
#         Returns array with residuals 
#         """
#         P2,Q2,P1,Q1,Q3,P3 = x
#         A1,A2,A3,pos1,pos2,pos3,vz1,vz2,vz3,P1o,Q1o,P2o,Q2o,P3o,Q3o,domega1,domega2,domega3,rho1,rho2,rho3,L1,L2,L3,du1,du2,du3 = args
#         
#         du1[0] = P1 - P1o    
#         du1[1] = Q1 - Q1o
#         du2[0] = P2 - P2o
#         du2[1] = Q2 - Q2o
#         du3[0] = P3 - P3o
#         du3[1] = Q3 - Q3o
#         
#         self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[pos1],update='L')
#         self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[pos2],update='L')
#         self.systemEquations[2].updateLARL([P3],[Q3],[A3],idArray=[pos3],update='L')
#                 
#         #calculate residuals
#         res1 = vz1*Q1+vz2*Q2+vz3*Q3
#         res2 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz2*P2+vz2*rho2*0.5*(Q2/A2)**2.
#         res3 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz3*P3+vz3*rho3*0.5*(Q3/A3)**2.    
#         res4 = np.dot(self.systemEquations[0].L[pos1][pos1+1],du1) - domega1
#         res5 = np.dot(self.systemEquations[1].L[pos2][pos2+1],du2) - domega2
#         res6 = np.dot(self.systemEquations[2].L[pos3][pos3+1],du3) - domega3
#         
#         return [res5,res2,res4,res1,res3,res6]
#     
#     def jacobiMatrixBifSys0(self,x, args):
#         """
#         Returns the jabcobi matrix, bifurcation-functions and x; J = dF/dx
#         Using constant areas, i.e. initial areas
#         """
#         P2,Q2,P1,Q1,Q3,P3 = x
#         A1,A2,A3,pos1,pos2,pos3,vz1,vz2,vz3,P1o,Q1o,P2o,Q2o,P3o,Q3o,domega1,domega2,domega3,rho1,rho2,rho3,L1,L2,L3,du1,du2,du3 = args
#                 
#         return np.array([[L2[0], L2[1]            , 0    , 0                , 0                , 0    ],
#                          [vz2  , vz2*rho2*Q2/A2**2, vz1  , vz1*rho1*Q1/A1**2, 0                , 0    ],
#                          [0    , 0                , L1[0], L1[1]            , 0                , 0    ],
#                          [0    , vz2              , 0    , vz1              , vz3              , 0    ],
#                          [0    , 0                , vz1  , vz1*rho1*Q1/A1**2, vz3*rho3*Q3/A3**2, vz3  ],
#                          [0    , 0                , 0    , 0                , L3[1]            , L3[0]]])
#     
#     def fsolveBifurcationSys1(self,x,args):
#         """
#         Residual Function with equations to solve for at the bifuraction
#         Using recalculated areas depending on the new pressure values
#         
#         Input:     x = array [P2,Q2,P1,Q1,Q3,P3]
#                    args = args with local variables
#         Returns array with residuals 
#         """
#         P2,Q2,P1,Q1,Q3,P3 = x
#         A1,A2,A3,pos1,pos2,pos3,vz1,vz2,vz3,P1o,Q1o,P2o,Q2o,P3o,Q3o,domega1,domega2,domega3,rho1,rho2,rho3,L1,L2,L3,du1,du2,du3 = args
#                         
#         A1 = self.A_func[0]([P1],0)
#         A2 = self.A_func[1]([P2],0)
#         A3 = self.A_func[2]([P3],0)
#         
#         du1[0] = P1 - P1o    
#         du1[1] = Q1 - Q1o
#         du2[0] = P2 - P2o
#         du2[1] = Q2 - Q2o
#         du3[0] = P3 - P3o
#         du3[1] = Q3 - Q3o
#         
#         if self.updateL == True:   
#             self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[pos1],update='L')
#             self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[pos2],update='L')
#             self.systemEquations[2].updateLARL([P3],[Q3],[A3],idArray=[pos3],update='L')
#         
#         
#         res1 = vz1*Q1+vz2*Q2+vz3*Q3
#         res4 = np.dot(self.systemEquations[0].L[pos1][pos1+1],du1) - domega1
#         res5 = np.dot(self.systemEquations[1].L[pos2][pos2+1],du2) - domega2
#         res6 = np.dot(self.systemEquations[2].L[pos3][pos3+1],du3) - domega3
#         
#         if self.nonLin == True:
#         #calculate residuals
#             res2 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz2*P2+vz2*rho2*0.5*(Q2/A2)**2.
#             res3 = vz1*P1+vz1*rho1*0.5*(Q1/A1)**2.+vz3*P3+vz3*rho3*0.5*(Q3/A3)**2.    
#         else:
#             res2 = vz1*P1+vz2*P2
#             res3 = vz1*P1+vz3*P3  
#         
#         return [res5,res2,res4,res1,res3,res6]
#     
#     def jacobiMatrixBifSys1(self,x, args):
#         """
#         Returns the jabcobi matrix, bifurcation-functions and x; J = dF/dx
#         Using recalculated areas depending on the new pressure values
#         """
#         P2,Q2,P1,Q1,Q3,P3 = x
#         A1,A2,A3,pos1,pos2,pos3,vz1,vz2,vz3,P1o,Q1o,P2o,Q2o,P3o,Q3o,domega1,domega2,domega3,rho1,rho2,rho3,L1,L2,L3,du1,du2,du3 = args
#         
#         A1 = self.A_func[0]([P1],0)
#         A2 = self.A_func[1]([P2],0)
#         A3 = self.A_func[2]([P3],0)
#         
#         
#         if self.nonLin == True:
#             jacobi = np.array([[L2[0], L2[1]            , 0    , 0                , 0                , 0    ],
#                              [vz2  , vz2*rho2*Q2/A2**2, vz1  , vz1*rho1*Q1/A1**2, 0                , 0    ],
#                              [0    , 0                , L1[0], L1[1]            , 0                , 0    ],
#                              [0    , vz2              , 0    , vz1              , vz3              , 0    ],
#                              [0    , 0                , vz1  , vz1*rho1*Q1/A1**2, vz3*rho3*Q3/A3**2, vz3  ],
#                              [0    , 0                , 0    , 0                , L3[1]            , L3[0]]])
#         
#         else:
#             jacobi = np.array([[L2[0], L2[1]          , 0    , 0                , 0                , 0    ],
#                              [vz2  , 0                , vz1  , 0                , 0                , 0    ],
#                              [0    , 0                , L1[0], L1[1]            , 0                , 0    ],
#                              [0    , vz2              , 0    , vz1              , vz3              , 0    ],
#                              [0    , 0                , vz1  , 0                , 0                , vz3  ],
#                              [0    , 0                , 0    , 0                , L3[1]            , L3[0]]])
#         #print max(np.linalg.eigvals(jacobi))
#         
#         return jacobi
#         
        
class Anastomosis():
    
    def __init__(self, leftMother, leftMotherSys,
                       rightMother, rightMotherSys,
                       daughter, daughterSys, 
                       currentMemoryIndex, dt, rigidAreas, solvingScheme):
        
        
        # vessel variables initially set, constant through simulation
        self.type = 'Anastomosis'
        
        self.name = ' '.join(['Anastomosis',str(leftMother.Id),str(rightMother.Id),str(daughter.Id)])
        
        #System Variables
        self.dt = dt
        self.currentMemoryIndex = currentMemoryIndex
        
        self.rho = []
        self.systemEquations = []
        self.z = []
        self.A_func = []
        self.positions =[]
#         self.vz = []
        self.names = []
        
#         # equations to solve in f solve
#         self.fsolveFunction = None
#         self.jacobiMatrix = None
        
        ###initialize
        ##left mother branch
        self.rho.append(leftMother.rho)
        self.z.append(leftMother.z)
        self.systemEquations.append(leftMotherSys)
        self.positions.append(-1)
        self.names.append(leftMother.Id)
        self.A_func.append(leftMother.A_nID)
        #SolutionVariables
        self.P_leftMother = leftMother.Psol
        self.Q_leftMother = leftMother.Qsol
        self.A_leftMother = leftMother.Asol
        
        ##left mother branch
        self.rho.append(rightMother.rho)
        self.z.append(rightMother.z)
        self.systemEquations.append(rightMotherSys)
        self.positions.append(-1)
        self.names.append(rightMother.Id)
        self.A_func.append(rightMother.A_nID)
        #SolutionVariables
        self.P_rightMother = rightMother.Psol
        self.Q_rightMother = rightMother.Qsol
        self.A_rightMother = rightMother.Asol
        
        ##right daughter branch
        self.rho.append(daughter.rho)
        self.z.append(daughter.z)
        self.systemEquations.append(daughterSys)
        self.positions.append(0)
        self.names.append(daughter.Id)
        self.A_func.append(daughter.A_nID)
#        self.vz.append(1)
        #SolutionVariables
        self.P_daughter = daughter.Psol
        self.Q_daughter = daughter.Qsol
        self.A_daughter = daughter.Asol

        self.rigidAreas = rigidAreas
    
        # Define the call function depending on the solving Scheme
        solvingScheme = "NonLinear"
        if solvingScheme == "Linear": 
            self.__call__ = self.callLinear
        elif solvingScheme == "NonLinear":
            print "classconnection 1652: using nonlinear anastomosis model: {0}".format(self.name) 
            self.__call__ = self.callNonLinear
        else:
            raise ValueError("Connections; wrong solving scheme! {}".format(solvingScheme))
        
        ## benchamark Test variables
        self.sumQErrorCount = 0
        self.maxQError = 0
        self.maxPErrorNonLin = 0 
        self.maxPError = 0
        self.sumPErrorCount = 0
        self.sumPErrorNonLinCount = 0
    
    def callLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_rightMother[n]
        Q2 = self.Q_rightMother[n]
        A2 = self.A_rightMother[n]
        
        P3 = self.P_daughter[n]
        Q3 = self.Q_daughter[n]
        A3 = self.A_daughter[n]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
                           
                
        # update LARL
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_1 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
        L,R3,LMBD,Z1,Z2,domega3_2 = self.systemEquations[2].updateLARL(P3,Q3,A3,pos3)
            
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        ###### Linear approach
        
        ####### change?!
        denom = R1_12 * R2_12 * R3_21 - R1_12 * R2_22 * R3_11 - R1_22 * R2_12 * R3_11 
        
        alpha1 = -( R1_11 * R2_12 * R3_21 - R1_11 * R2_22 * R3_11 - R1_21 * R2_12 * R3_11)/denom
        alpha2 = -( R2_11 * R2_22 * R3_11 - R2_12 * R2_21 * R3_11)/denom
        alpha3 = -( R2_12 * R3_22 * R3_11 - R2_12 * R3_12 * R3_21)/denom
        
        beta1 = -(R1_11 * R1_22 * R3_11 - R1_12 * R1_21 * R3_11)/denom
        beta2 = -(R1_12 * R2_11 * R3_21 - R1_12 * R2_21 * R3_11 - R1_22 * R2_11 * R3_11)/denom
        beta3 = -(R1_12 * R3_11 * R3_22 - R1_12 * R3_12 * R3_21)/denom
        
        gamma1 = -( R1_11 * R1_22 * R2_12 - R1_12 * R1_21 * R2_12)/denom
        gamma2 = -( R1_12 * R2_11 * R2_22 + R1_12 * R2_12 * R2_21)/denom
        gamma3 = -( R1_12 * R2_12 * R3_22 - R1_12 * R2_22 * R3_12 - R1_22 * R2_12 * R3_12)/denom
        
        ############        
        domega1_2 = (alpha1 * domega1_1  + alpha2 * domega2_1 + alpha3 * domega3_2 )
        domega2_2 = (beta1  * domega1_1  + beta2  * domega2_1 + beta3  * domega3_2 )
        domega3_1 = (gamma1 * domega1_1  + gamma2 * domega2_1 + gamma3 * domega3_2 )
        
        P1_new = P1o + (R1_11*domega1_1 + R1_12*domega1_2)
        Q1_new = Q1o + (R1_21*domega1_1 + R1_22*domega1_2)
    
        P2_new = P2o + (R2_11*domega2_1 + R2_12*domega2_2)
        Q2_new = Q2o + (R2_21*domega2_1 + R2_22*domega2_2)
        
        P3_new = P3o + (R3_11*domega3_1 + R3_12*domega3_2)
        Q3_new = Q3o + (R3_21*domega3_1 + R3_22*domega3_2)
               
        
        # apply calculated values to next time step solution
        self.P_leftMother[n+1][pos1]  = P1_new
        self.Q_leftMother[n+1][pos1]  = Q1_new
        self.P_rightMother[n+1][pos2] = P2_new
        self.Q_rightMother[n+1][pos2] = Q2_new
        self.P_daughter[n+1][pos3]    = P3_new
        self.Q_daughter[n+1][pos3]    = Q3_new
        
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            tmpstring = "Connection: {} calculated negative pressure at time {} (n {},dt {})".format(self.name,n*dt,n,dt)
            tmpstring = tmpstring + "\n the values were: P1_new = {}, P2_new = {}, P3_new = {}".format(P1_new, P2_new, P3_new)
            raise ValueError(tmpstring)
            #exit()
                
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]   = A1n
        self.A_rightMother[n+1][pos2]  = A2n       
        self.A_daughter[n+1][pos3]     = A3n  
            
        ## non linear error        
        try: sumQError = abs(Q1_new+Q2_new-Q3_new)#/abs(Q1_new)
        except: sumQError = 0.0
        if sumQError > 0.0: 
            self.sumQErrorCount = self.sumQErrorCount+1
        if sumQError > self.maxQError:
            self.maxQError  = sumQError
        #print self.name,' \n Error cons mass',  sumQError, self.maxQError ,' - ', n, self.sumQErrorCount
        if sumQError > 1.e-5:
            print "Warning: Connection: {} to high error in conservation of mass at time {} (n {},dt {})".format(self.name,n*dt,n,dt)
            print sumQError, ' <- Q1,Q2,Q3 ',Q1_new, Q2_new, Q3_new
            #exit()
        
        sumPError = abs(P1_new-P3_new)/abs(P3_new)
        if sumPError > 0.0: 
            self.sumPErrorCount = self.sumPErrorCount+1
        if sumPError > self.maxPError:
            self.maxPError  = sumPError
        #print self.name,' Error P lin    ',  sumPError, self.maxPError ,' - ', n, self.sumPErrorCount
        if sumPError > 1.e-5:
            print "ERROR: Connection: {} to high error in conservation of pressure at time {} (n {},dt {})".format(self.name,n*dt,n,dt)
            print sumPError, ' <- P1,P2,P3 ',P1_new, P2_new, P3_new
            #exit()
        
        sumPErrorNonLin = 1050./2.*(Q1_new/A1n)**2#abs(P1_new+500*(Q1_new/A1n)**2-(P2_new+500*(Q2_new/A2n)**2))/abs(P1_new+0.5*(Q1_new/A1n)**2)
        if sumPErrorNonLin > 0.0: 
            self.sumPErrorNonLinCount = self.sumPErrorNonLinCount+1
        if sumPErrorNonLin > self.maxPErrorNonLin:
            self.maxPErrorNonLin  = sumPErrorNonLin
        print self.name,' Error P non lin',  sumPErrorNonLin, self.maxPErrorNonLin ,' - ', n, self.sumPErrorNonLinCount
        
        print self.name,'dynamic Pressures',1050./2.*(Q1_new/A1n)**2,1050./2.*(Q2_new/A2n)**2,1050./2.*(Q3_new/A3n)**2, '--',1050./2.*(Q1_new/A1n)**2+-1050./2.*(Q3_new/A3n)**2
        
        
    def callNonLinear(self):
        """
        Call function for vessel-vessel connection
        """        
        dt = self.dt
        n = self.currentMemoryIndex[0]
        #if n == 1:
            #print "using nonlinear bifurcation model"
        #print "using nonlinear bifurcation model"
        pos1 = self.positions[0]
        pos2 = self.positions[1]
        pos3 = self.positions[2]
        
        P1 = self.P_leftMother[n]
        Q1 = self.Q_leftMother[n]
        A1 = self.A_leftMother[n]
        
        P2 = self.P_rightMother[n]
        Q2 = self.Q_rightMother[n]
        A2 = self.A_rightMother[n]
        
        P3 = self.P_daughter[n]
        Q3 = self.Q_daughter[n]
        A3 = self.A_daughter[n]
        
        rho1 = self.rho[0]
        rho2 = self.rho[1]
        rho3 = self.rho[2]
        
        P1o = P1[pos1]
        Q1o = Q1[pos1]
        A1o = A1[pos1]
        P2o = P2[pos2]
        Q2o = Q2[pos2]
        A2o = A2[pos2]
        P3o = P3[pos3]
        Q3o = Q3[pos3]
        A3o = A3[pos3]

        ## update LARL
        L,R1,LMBD,Z1,Z2,domega1_1 = self.systemEquations[0].updateLARL(P1,Q1,A1,pos1)
        L,R2,LMBD,Z1,Z2,domega2_1 = self.systemEquations[1].updateLARL(P2,Q2,A2,pos2)
        L,R3,LMBD,Z1,Z2,domega3_2 = self.systemEquations[2].updateLARL(P3,Q3,A3,pos3)
        
        # local R matrices
        R1_11 = R1[0][0]
        R1_12 = R1[0][1]
        R1_21 = R1[1][0]
        R1_22 = R1[1][1]
        
        R2_11 = R2[0][0]
        R2_12 = R2[0][1]
        R2_21 = R2[1][0]
        R2_22 = R2[1][1]
        
        R3_11 = R3[0][0]
        R3_12 = R3[0][1]
        R3_21 = R3[1][0]
        R3_22 = R3[1][1]
        
        if n>0:
            P1o2 = self.P_leftMother[n-1][pos1]
            P2o2 = self.P_rightMother[n-1][pos2]
            P3o2 = self.P_daughter[n-1][pos3]
            
            deltaP1_last = P1o - P1o2
            deltaP2_last = P2o - P2o2
            deltaP3_last = P3o - P3o2
            
            domega1_2_init = (deltaP1_last - R1_11*domega1_1)/R1_12 #same omega as previous timestep
            domega2_2_init = (deltaP2_last - R2_11*domega2_1)/R2_12 #same omega as previous timestep
            domega3_1_init = (deltaP3_last - R3_12*domega3_2)/R3_11 #same omega as previous timestep
        else:
            domega1_2_init = 0#-R1_11*domega1_1/R1_12 #result if pressure is constant- could be optimized
            domega2_2_init = 0#-R2_12*domega2_1/R2_11 #result if pressure is constant
            domega3_1_init = 0#-R3_12*domega3_2/R3_11 #result if pressure is constant

        epsilonvalues = np.array([1,1,1])
        epsilonlimit = 1e-10
        domega1_2_last = domega1_2_init
        domega2_2_last = domega2_2_init
        domega3_1_last = domega3_1_init
        A1_last = A1o
        A2_last = A2o
        A3_last = A3o
        Niterations = 0
        Xold = np.array([domega1_2_last, domega2_2_last, domega3_1_last])
#         print "\n"
#         print "domega1_2_last, domega2_1_last, domega3_1_last, domega1_1, domega2_1, domega3_2 : ", domega1_2_last, domega2_1_last, domega3_1_last, domega1_1, domega2_1, domega3_2
#         print "\n"
        while epsilonvalues[0]>1e-14 or epsilonvalues[1]>0.001 or epsilonvalues[2]>0.001:
            """iterative Newton Rahpson solver
                domega1_2, domega2_1, domega3_1 are the unknowns that are changing from
            each iteration. A1, A2, A3 are also changing, but the value from the previous iteration is used. (should actually be for the next timestep, but this should converge to the correct solution)
            R1_11, R1_12, R1_21, R1_22, ,R2_11, R2_12, R2_21, R2_22, ,R3_11, R3_12, R3_21, R3_22, domega1_1, domega2_1, domega3_2, Q1o, Q2o, Q2o 
            are constant for each timestep, and are taken from previous timestep. domega1_1, domega2_1, domega3_2 are the field domegas.
            """
            domega1_2_last = Xold[0]
            domega2_2_last = Xold[1]
            domega3_1_last = Xold[2]
            
            Q1discretize = Q1o + R1_21*domega1_1 + R1_22*domega1_2_last
            Q2discretize = Q2o + R2_21*domega2_1 + R2_22*domega2_2_last
            Q3discretize = Q3o + R3_21*domega3_1_last + R3_22*domega3_2
            
            P1discretize = P1o + R1_11*domega1_1 + R1_12*domega1_2_last
            P2discretize = P2o + R2_11*domega2_1 + R2_12*domega2_2_last
            P3discretize = P3o + R3_11*domega3_1_last + R3_12*domega3_2
            

            if self.rigidAreas == False:
                try:
                    A1_last = self.A_func[0]([P1discretize],pos1)
                    A2_last = self.A_func[1]([P2discretize],pos2)
                    A3_last = self.A_func[2]([P3discretize],pos3)
                except FloatingPointError as E:
                    print "Floating Point error in Connection {}".format(self.name)
                    raise E
            else:
                A1_last = A1[pos1]
                A2_last = A2[pos2]
                A3_last = A3[pos3] 
            
            f1 = Q1discretize + Q2discretize - Q3discretize#R1_21*domega1_1 + R1_22*domega1_2_last  - R2_21*domega2_1_last - R2_22*domega2_1 - R3_21*domega3_1_last - R3_22*domega3_2
            
            f2 = P1discretize + 0.5*rho1*((Q1discretize/A1_last)**2) - P3discretize - 0.5*rho3*((Q3discretize/A3_last)**2)
            
            f3 = P2discretize + 0.5*rho2*((Q2discretize/A2_last)**2) - P3discretize - 0.5*rho3*((Q3discretize/A3_last)**2)
            
            F = np.array([f1, f2, f3])
            """Inverse Jacobi elements: """
            a = R1_22
            b = R2_22
            c = -R3_21
            d = R1_12 + rho1*R1_22*(Q1discretize)/(A1_last**2)
            e = - R3_11 - rho3*R3_21*(Q3discretize)/(A3_last**2)
            f = R2_12 + rho2*R2_22*(Q2discretize)/(A2_last**2)
            
            
            Determinant = a*e*f + b*d*e -c*d*f
            
            J_inv = np.array([[ e*f, b*e - c*f, -b*e  ],
                          [ -d*e, -a*e, a*e - c*d ],
                          [  -d*f, a*f, b*d ]]) / (Determinant)
            
            Xnew = Xold - np.dot(J_inv,F)

            epsilonvalues = np.abs(F)
            Niterations = Niterations + 1
            
            if Niterations > 50:
                print "\n"
                print "Niterations excedded in anastomosis calculation in vessel, Niterations: ", self.names[0], Niterations

                print "Xnew: ", Xnew
                print "Xold: ", Xold
                
                print "f1: ", f1
                print "f2: ", f2
                print "f3: ", f3
                print "epsilonvalues: ", epsilonvalues
                print "Q1discretize, Q1o: ", Q1discretize, Q1o
                print "Q2discretize, Q2o: ", Q2discretize, Q2o
                print "Q3discretize, Q3o: ", Q3discretize, Q3o
                print "P1discretize, P1o: ", P1discretize, P1o
                print "P2discretize, P2o: ", P2discretize, P2o
                print "P3discretize, P3o: ",P3discretize, P3o
                #exit()
                break
            Xold = Xnew
            #exit()
        
        
        Q1_new = Q1discretize
        Q2_new = Q2discretize
        Q3_new = Q3discretize
        
        P1_new = P1discretize
        P2_new = P2discretize
        P3_new = P3discretize
        
        # apply calculated values to next time step
        self.P_leftMother[n+1][pos1]    = P1_new
        self.Q_leftMother[n+1][pos1]    = Q1_new
        self.P_rightMother[n+1][pos2]  = P2_new
        self.Q_rightMother[n+1][pos2]  = Q2_new
        self.P_daughter[n+1][pos3] = P3_new
        self.Q_daughter[n+1][pos3] = Q3_new
        
        if P1_new < 0 or P2_new < 0 or P3_new < 0:
            print "ERROR: Connection: {} calculated negative pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P1_new, P2_new, P3_new
            print "Niterations: ", Niterations
            
            print "solving nonlinear/total pressure anastomosis"
            print "domega1_2_init: ", domega1_2_init
            print "domega1_1: ", domega1_1
                    
            print "f1: ", f1
            print "f2: ", f2
            print "f3: ", f3
            print "epsilonvalues: ", epsilonvalues
            print "Q1discretize, Q1o: ", Q1discretize, Q1o
            print "Q2discretize, Q2o: ", Q2discretize, Q2o
            print "Q3discretize, Q3o: ", Q3discretize, Q3o
            print "P1discretize, P1o: ", P1discretize, P1o
            print "P2discretize, P2o: ", P2discretize, P2o
            print "P3discretize, P3o: ",P3discretize, P3o
            exit()
        
        
        # calculate new areas
        if self.rigidAreas == False:
            A1n = self.A_func[0]([P1_new],pos1)
            A2n = self.A_func[1]([P2_new],pos2)
            A3n = self.A_func[2]([P3_new],pos3)
        else:
            A1n = A1[pos1]
            A2n = A2[pos2]
            A3n = A3[pos3] 
           
        self.A_leftMother[n+1][pos1]       = A1n
        self.A_rightMother[n+1][pos2]     = A2n       
        self.A_daughter[n+1][pos3]    = A3n