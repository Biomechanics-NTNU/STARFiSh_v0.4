import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
#sys.path.append(cur+'/NetworkLib')

import numpy as np


class Field():
    
    def __init__(self, vessel, currentMemoryIndex, dt, systemEquation, rigidArea, solvingSchemeField = 'MacCormack_Flux'):
        """
        Constructor of Field object
        
        calculates the interior field of a vessel
        with a MackKormack Predictor-Corrector Schmea
        """
        
        self.name = ' '.join(['Field',str(vessel.Id)])
        
        self.vessel = vessel
        #System and Vessel Variables
        self.dz = vessel.dz
        self.dt = dt
        # current time index in memory
        self.currentMemoryIndex = currentMemoryIndex
        
        self.systemEquation = systemEquation
        self.AFunction = vessel.A
        self.rigidArea = rigidArea
        
        #SolutionVariables
        self.P = vessel.Psol
        self.P_pre = np.ones_like(vessel.Psol[0])
        self.Q = vessel.Qsol
        self.Q_pre = np.ones_like(vessel.Qsol[0])
        self.A = vessel.Asol
        self.A_pre = np.ones_like(vessel.Asol[0])
        
        self.step = "predictor"
                       
        if solvingSchemeField == 'MacCormack_Matrix':
            self.__call__ = self.MacCormackMatrix
        elif solvingSchemeField == "MacCormack_Flux":
            self.__call__ = self.MacCormackFlux
        else:
            raise ValueError('Classfields51: error, scheme for solving field not correct')
        
    def F(self,u,A,C,Aconst,Cconst):
    
        """Flux based on conservative form of governing mass and momentum equations"""
        gamma = self.vessel.gamma
        alpha = 1 #(gamma+2.0)/(gamma+1.0) #1
        flux=np.zeros_like(u)
        p=u[0,:]
        q=u[1,:]
        rho = self.vessel.rho
        flux[0,:] = q/Cconst
        flux[1,:] = Aconst*p/rho + alpha*q**2/A
        return flux 
    
    
    def MacCormackFlux(self):
        
        dt = self.dt
        #dz = self.dz[0] #currently only equidistant spacing is implemented and thus dz can be taken as float instead of vector
        
        # the current position in solution memory
        n = self.currentMemoryIndex[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
        C = self.vessel.C(P)
        netGravity = self.vessel.netGravity[n]
                
        #dt = self.dt
        dx = self.dz[0] #currently only equidistant spacing is implemented and thus dz can be taken as float instead of vector
        
        # the current position in solution memory
        #n = self.currentMemoryIndex[0]
        #print "using Flux based formulation in field calculation"
        #if n== 1:
            #print "using Flux based formulation in field calculation"
#             P = self.P[n]
#             Q = self.Q[n]
#             A = self.A[n]
        #C = self.vessel.C(P)
        
        rho = self.vessel.rho
        gamma = self.vessel.gamma
        alpha = (gamma+2.0)/(gamma+1.0) #flux corrector term associated with nonlinear integration  of velocityprofile in convective term (flat profile yields alpha=1, pousouille yields alpha =4/3)
        my = self.vessel.my
        
        u = np.zeros((2,len(P))) #matrix for storing pressure and flow to caalculate fluxes
        up = np.zeros((2,len(P))) # predictor values for pressure and flow
        u[0,:] = P
        u[1,:] = Q
        up=u.copy()
        
        Utemp1=u[1,:-1]/A[:-1]
        Aconst1 = A[:-1]
        Cconst1 = C[:-1]
        A1 = A[1:]
        A2 = A[:-1]
        C1 = C[1:]
        C2 = C[:-1]
        up[:,:-1] = u[:,:-1] - dt*(self.F(u[:,1:],A1,C1,Aconst1,Cconst1)-self.F(u[:,:-1],A2,C2,Aconst1,Cconst1))/dx 
        # TODO: make sure the areas used are correct
        A_grav = A[:-1]
        up[1,:-1] = up[1,:-1]-dt*2*(gamma+2)*my*Utemp1*np.pi/rho +  dt* A_grav * netGravity
        
        if self.rigidArea == True:
            A_p = A
            C_p = C
        else:
            try:
                A_p = self.vessel.A(up[0,:])
                C_p = self.vessel.C(up[0,:])
            except FloatingPointError as E:
                print "Floating Point error in Field {}".format(self.name)
                raise E
        
        A_p1 = A_p[1:]
        A_p2 = A_p[:-1]
        C_p1 = C_p[1:]
        C_p2 = C_p[:-1]
        Utemp2=up[1,1:]/A_p1
        Aconst2 = A_p[1:]
        Cconst2 = C_p[1:]
        u[:,1:] = .5*(u[:,1:]+up[:,1:] -  dt/dx*(self.F(up[:,1:],A_p1,C_p1,Aconst2,Cconst2)-self.F(up[:,:-1],A_p2,C_p2,Aconst2,Cconst2)))
        A_grav = A_p1
        u[1,1:] = u[1,1:]-0.5*dt*2*(gamma+2)*my*Utemp2*np.pi/rho + 0.5 * dt * A_grav * netGravity
        
        Pnewinterior = u[0,:]
        Qnewinterior = u[1,:]
        
        if self.rigidArea == True:
            Anewinterior = A
        else:
            Anewinterior = self.vessel.A(Pnewinterior)
        
        Pnewinterior = Pnewinterior[1:-1]
        Qnewinterior = Qnewinterior[1:-1]
        Anewinterior=Anewinterior[1:-1]
        
        self.P[n+1][1:-1] = Pnewinterior
        self.Q[n+1][1:-1] = Qnewinterior
        self.A[n+1][1:-1] = Anewinterior

        #TODO: Please explain this if statement in a comment.
        if (self.P[n+1] < 0).any():
            raise ValueError("ERROR: {} calculated negative pressure in corrector step at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt))
        
    def MacCormackMatrix(self):
        """
        Mac Cormack Predictor-Corrector
        """
        # solve vessel objects
        dt = self.dt
        
        # Any need to print this out?
        # print "using Matrix based formulation in field calculation"
        
        # the current position in solution memory
        n = self.currentMemoryIndex[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
        
        # set bc values of predictor step arrays
        self.P_pre[0]   = P[0]
        self.P_pre[-1]  = P[-1] 
        self.Q_pre[0]   = Q[0]
        self.Q_pre[-1]  = Q[-1] 
        self.A_pre[0]   = A[0]
        
        self.A_pre[-1]  = A[-1]
        
        #""" Predictor Step """            
        # update matrices               
        m12,m21,m22,b2 = self.systemEquation.updateSystem(P,Q,A)
         
        # create lokal variables for predictor varibales with
        P_pre = self.P_pre
        Q_pre = self.Q_pre 
        A_pre = self.A_pre
          
        dzPre = self.dz
        dzCor = self.dz[0:-1]
        
        # calculate derivatives forward # predict bc-predictor value at P[0] etc.
        dPdz  = (P[1::] - P[0:-1])/dzPre
        dQdz  = (Q[1::] - Q[0:-1])/dzPre
        dQ2A  = pow(Q,2.)/A
        dQ2dz = (dQ2A[1::] - dQ2A[0:-1])/dzPre
        
        # solve vessel
        P_pre[0:-1] = (P[0:-1] - (m12[0:-1]*dQdz)*dt)
        Q_pre[0:-1] = (Q[0:-1] - (m21[0:-1]*dPdz + m22[0:-1]*dQ2dz - b2[0:-1] )*dt)
         
        # check pressure solution
        if (P_pre < 0).any():
            raise ValueError("{} calculated negative pressure P_pre = {} in predictor step at time {} (n {},dt {}), exit system".format(self.name, P_pre,n*dt,n,dt))

        # solve area
        if self.rigidArea == True:
            A_pre = A
        else:
            A_pre[0:-1] = self.AFunction(P_pre)[0:-1]        
                                          
        #"""Corrector Step"""    
        # update matrices  
        m12,m21,m22,b2 = self.systemEquation.updateSystem(P_pre,Q_pre,A_pre)
         
        # calculate derivatives backward 
        dPdz  = (P_pre[1:-1] - P_pre[0:-2])/dzCor
        dQdz  = (Q_pre[1:-1] - Q_pre[0:-2])/dzCor
        dQ2A  = pow(Q_pre,2.)/A_pre
        dQ2dz = (dQ2A[1:-1] - dQ2A[0:-2])/dzCor
                 
        #solve vessel
        self.P[n+1][1:-1] = (P[1:-1] + P_pre[1:-1] - (m12[1:-1]*dQdz)*dt)/2.0
        self.Q[n+1][1:-1] = (Q[1:-1] + Q_pre[1:-1] - (m21[1:-1]*dPdz + m22[1:-1]*dQ2dz - b2[1:-1])*dt )/2.0
         
        # check pressure solution
        if (self.P[n+1] < 0).any():
            raise ValueError("{} calculated negative pressure self.P[n+1] = {} in corrector step at time {} (n {},dt {}), exit system".format(self.name,self.P[n+1],n*dt,n,dt))
         
        # solve area
        if self.rigidArea == True:
            self.A[n+1] = A_pre
        else:
            self.A[n+1][1:-1] = self.AFunction(self.P[n+1])[1:-1]
