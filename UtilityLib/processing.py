import numpy as np



    

def linearWaveSplitting(pressureArray,flowArray,areaArray,waveSpeedArray,rho,maxLength=500):
    """
    calculates the linear wave splitting for a given P,Q,A,c (np.arrays) and rho (float) 
    return values are: P_forward, P_backward, Q_forward, Q_backward 
    """

    ## calculate Zo and recast// delete last element
    Zo =  rho*waveSpeedArray/areaArray
    Zo = np.ones_like(Zo)*np.mean(Zo)
     
    ##calculateing dP and dQ
    dP = pressureArray[1::] - pressureArray[0:-1]      
    dQ = flowArray[1::] - flowArray[0:-1] 
    
    dP_div_Z = pressureArray[1::]/Zo[1::] - pressureArray[0:-1]/Zo[0:-1]
    dQ_multi_Z = flowArray[1::]*Zo[1::] - flowArray[0:-1]*Zo[0:-1]
    
    ## calculate dp_f, dp_b and dQ_f, dq_b     
    dp_f = (dP + dQ_multi_Z)/2.0 
    dp_b = (dP - dQ_multi_Z)/2.0
    dQ_f = (dQ + dP_div_Z)/2.0 
    dQ_b = (dQ - dP_div_Z)/2.0 
    
    p_f = np.cumsum(dp_f)
    p_b = np.cumsum(dp_b)
    q_f = np.cumsum(dQ_f)
    q_b = np.cumsum(dQ_b)   
    
    # alternative
    #pf = (pressureArray-pm+(flowArray-qm)*Zo)/2.0
        
    """
    ## old version 
#     
#     #cut (reshape) vectors if they are very long to avoid memory problems
#     lendp_f = len(dp_f)
#     if lendp_f > maxLength:
#         # state if len(dp_f) is a prime number
#         isPrime = False
#         # find divisor
#         divisor = 0
#         # exclude split into two parts if the result is still to memory consuming
#         startValue = 3
#         #find divisors
#         divs = np.linspace(startValue,int(lendp_f/2.0),int(lendp_f/2.0)+1-startValue)
#         check = 0
#         for div in divs:
#             if 	lendp_f%div == 0:
#                 divisor = div
#                 dlength = lendp_f/div
#                 if dlength < maxLength:
#                     break
#             elif lendp_f%div == 1 and check == 0:
#                 divisorPrimeCase = div
#                 check = 1
#        
#         
#         if divisor == 0:
#             isPrime = True
#             divisor = divisorPrimeCase
#             dlength = (lendp_f-1)/divisorPrimeCase
#         
#         print 'original signal length',lendp_f,' splitted in ',divisor,' parts a,', dlength, '+',int(isPrime),'(prime number)'
#         
#         
#         if isPrime:
#             #get last entry
#             dp_fPrime = dp_f[-1]
#             dp_bPrime = dp_b[-1]
#             dQ_fPrime = dQ_f[-1]
#             dQ_bPrime = dQ_b[-1]
#             #removing last entry from signals
#             dp_f = dp_f[0:len(dp_f)-1]
#             dp_b = dp_b[0:len(dp_b)-1]
#             dQ_f = dQ_f[0:len(dQ_f)-1]
#             dQ_b = dQ_b[0:len(dQ_b)-1]
# 
#         #reshape dp_f, dp_b and dQ_f, dq_b
#         dp_f.shape = (divisor,dlength)
#         dp_b.shape = (divisor,dlength)
#         dQ_f.shape = (divisor,dlength)
#         dQ_b.shape = (divisor,dlength)
# 
#         matrix = np.tril(np.ones((dlength,dlength)))
#         
#         
#         # set initial Values
#         dp_f[0][0] = dp_f[0][0] #+Pm
#         dQ_f[0][0] = dQ_f[0][0] #+Qm
#         
#         #perform integration/summation for all splitted compartments and add them again
#         # the first one
#         lastElement = dp_f[0][-1]
#         p_f = np.dot(matrix,dp_f[0])
#         p_b = np.dot(matrix,dp_b[0])
#         q_f = np.dot(matrix,dQ_f[0])
#         q_b = np.dot(matrix,dQ_b[0])
#         # the rest
#         for i in np.arange(1,divisor):  
#             
#             p_f = np.append(p_f, np.dot(matrix,dp_f[i])+p_f[-1])
#             p_b = np.append(p_b, np.dot(matrix,dp_b[i])+p_b[-1]) 
#             q_f = np.append(q_f, np.dot(matrix,dQ_f[i])+q_f[-1])
#             q_b = np.append(q_b, np.dot(matrix,dQ_b[i])+q_b[-1])
# 
#         if isPrime:
#             #append the last element if it was removed before
#             p_f = np.append(p_f, dp_fPrime+p_f[-1])
#             p_b = np.append(p_b, dp_bPrime+p_b[-1]) 
#             q_f = np.append(q_f, dQ_fPrime+q_f[-1])
#             q_b = np.append(q_b, dQ_bPrime+q_b[-1])
# 
#     else:
#         ## caculate and set d_f, p_b and Q_f, q_b by summation over time to the iven timpoint tp
#         
#         dp_f[0] = dp_f[0] #+Pm
#         dQ_f[0] = dQ_f[0] #+Qm
#         
#         matrix = np.tril(np.ones((len(dp_f),len(dp_f))))
#         p_f = np.dot(matrix,dp_f) 
#         p_b = np.dot(matrix,dp_b) 
#         q_f = np.dot(matrix,dQ_f) 
#         q_b = np.dot(matrix,dQ_b)   
    """
        
    return p_f,p_b,q_f,q_b

def nonLinearWaveSplitting(pressureArray,flowArray,areaArray,waveSpeedArray,comlianceArray,rho,gamma=2.0,maxLength=500):
    """
    calculates the non linear wave splitting for a given P,Q,A,c (np.arrays), rho (float, ) 
    return values are: P_forward, P_backward, Q_forward, Q_backward 
    """
    
    ## calculate velocity
    v = np.mean(flowArray/areaArray)
    
    dlt = (gamma+2)/(gamma+1)
    ## calculate eigenvalues lambda1/2
    lambda1 = dlt * v + np.sqrt( waveSpeedArray**2 + dlt *(dlt-1)*v**2)
    lambda2 = dlt * v - np.sqrt( waveSpeedArray**2 + dlt *(dlt-1)*v**2) 
    ### recast lambdas
    lambda1 = lambda1
    lambda2 = lambda2
    
    ##calculate Compliance
    #C1 = dA / dP
    #C = areaArray / (rho*waveSpeedArray**2)
    #C = C[0:len(C)-1]
    C = comlianceArray
    
    #print C1-C
    ## callculate Z+/Z-
    Z_plus  =  1.0 / (lambda1 * C)
    Z_minus = -1.0 / (lambda2 * C)
    
    
#     ##calculate differentials dP, dQ and dA
#     dP = pressureArray[1::] - pressureArray[0:-1]   
#     dQ = flowArray[1::]     - flowArray[0:-1]
#     dA = areaArray[1::]     - areaArray[0:-1]
#     
# # #     pm = np.mean(pressureArray)
# # #     p0 = pressureArray[0]
# # #     qm = np.mean(flowArray)
# # #     q0 = flowArray[0]
# # #     pressureArray = pressureArray - pm
# # #     flowArray = flowArray - qm
# # #     
# # #     pf =  (pressureArray+Z_minus*flowArray)*Z_plus/(Z_plus+Z_minus)
# # #     pb =  (pressureArray-Z_plus*flowArray)*Z_minus/(Z_plus+Z_minus)
# # #     qf =  (pressureArray+Z_minus*flowArray)/(Z_plus+Z_minus)
# # #     qb = -(pressureArray-Z_plus*flowArray)/(Z_plus+Z_minus)
#     
#     dQZminus = Z_minus[1::]*flowArray[1::] - Z_minus[0:-1]*flowArray[0:-1]
#     dQZplus  = Z_plus[1::] *flowArray[1::] - Z_plus[0:-1] *flowArray[0:-1]
#     # calculate average impedance 
#     Z_plus  = (Z_plus[0:-1]  + Z_plus[1::])/2.0
#     Z_minus = (Z_minus[0:-1] + Z_minus[1::])/2.0
#     
#     #Z_plus  = Z_plus[0:-1] 
#     #Z_minus = Z_minus[0:-1] 
#     
#     
#     ## calculate differnetial characteristics dw1 /dw2
#     dw1 = dP + dQZminus# Z_minus * dQ #
#     dw2 = dP - dQZplus# Z_plus * dQ  #
#     
#     ## calculate dp_f, dp_b and dQ_f, dq_b     
#     dp_f =   (dw1*Z_plus) / (Z_minus+Z_plus)
#     dp_b =  (dw2*Z_minus) / (Z_minus+Z_plus)
#     dQ_f =   dw1 / (Z_minus+Z_plus)
#     dQ_b =  -dw2 / (Z_minus+Z_plus)
    
    ### calculate impedance factors
    ## Z+ / (Z- + Z+)
    # backward
    ZpDivZsumB = Z_plus[0:-1]/(Z_minus[0:-1]+Z_plus[0:-1])
    # forward
    ZpDivZsumF = Z_plus[1::]/(Z_minus[1::]+Z_plus[1::])
    
    ## Z- Z+ / (Z- + Z+)
    # backward
    ZpZmDivZsumB = (Z_minus[0:-1]*Z_plus[0:-1])/(Z_minus[0:-1]+Z_plus[0:-1])
    # forward
    ZpZmDivZsumF = (Z_minus[1::]*Z_plus[1::])/(Z_minus[1::]+Z_plus[1::])
    
    ## Z- / (Z- + Z+)
    # backward
    ZmDivZsumB =  Z_minus[0:-1]/(Z_minus[0:-1]+Z_plus[0:-1])
    # forward
    ZmDivZsumF =  Z_minus[1::]/(Z_minus[1::]+Z_plus[1::])
    
    ## 1.0 / (Z- + Z+)
    # backward
    OneDivZsumB = 1.0/(Z_minus[0:-1]+Z_plus[0:-1])
    # forward 
    OneDivZsumF = 1.0/(Z_minus[1::]+Z_plus[1::])
    
    # p_f = [Z+ / (Z- + Z+) * dp ] + [Z- Z+ / (Z- + Z+) * dq]
    dp_f = (ZpDivZsumF * pressureArray[1::] - ZpDivZsumB * pressureArray[0:-1]) + (ZpZmDivZsumF * flowArray[1::] - ZpZmDivZsumB * flowArray[0:-1])
    # p_b = [Z- / (Z- + Z+) * dp ] - [Z- Z+ / (Z- + Z+) * dq]
    dp_b = ((ZmDivZsumF * pressureArray[1::] - ZmDivZsumB * pressureArray[0:-1]) - (ZpZmDivZsumF * flowArray[1::] - ZpZmDivZsumB * flowArray[0:-1]))
    # q_f = [Z- / (Z- + Z+) * dp ] - [Z- Z+ / (Z- + Z+) * dq]
    dq_f = ((OneDivZsumF * pressureArray[1::] - OneDivZsumB * pressureArray[0:-1]) + (ZmDivZsumF * flowArray[1::] - ZmDivZsumB * flowArray[0:-1]))
    # q_b = [Z- / (Z- + Z+) * dp ] - [Z- Z+ / (Z- + Z+) * dq]
    dq_b =-((OneDivZsumF * pressureArray[1::] - OneDivZsumB * pressureArray[0:-1]) - (ZpDivZsumF * flowArray[1::] - ZpDivZsumB * flowArray[0:-1]))
     
    ## sum delta arrays useing cumsum as it is faster
    p_f = np.cumsum(dp_f)
    p_b = np.cumsum(dp_b)
    q_f = np.cumsum(dq_f)
    q_b = np.cumsum(dq_b)
    
    return p_f,p_b,q_f,q_b
    #return pf[1::]+pm,pb[1::]+pm,qf[1::],qb[1::]
    #return pf[0:-1],pb[1::],qf[1::],qb[1::]
    """
    ## old version - using matrix method
    #cut (reshape) vectors if they are very long to avoid memory problems
    lendp_f = len(dp_f)
    if lendp_f > maxLength:
        # state if len(dp_f) is a prime number
        isPrime = False
        # find divisor
        divisor = 0
        # exclude split into two parts if the result is still to memory consuming
        startValue = 3
        #find divisors
        divs = np.linspace(startValue,int(lendp_f/2.0),int(lendp_f/2.0)+1-startValue)
        check = 0
        for div in divs:
            if lendp_f%div == 0:
                divisor = div
                dlength = lendp_f/div
                if dlength < maxLength:
                    break
            elif lendp_f%div == 1 and check == 0:
                divisorPrimeCase = div
                check = 1   
                
        if divisor == 0:
            isPrime = True
            divisor = divisorPrimeCase
            dlength = (lendp_f-1)/divisorPrimeCase

        print 'original signal length',lendp_f,' splitted in ',divisor,' parts a,', dlength, '+',int(isPrime),'(prime number)'

        if isPrime:
            #get last entry
            dp_fPrime = dp_f[-1]
            dp_bPrime = dp_b[-1]
            dQ_fPrime = dQ_f[-1]
            dQ_bPrime = dQ_b[-1]
            #removing last entry from signals
            dp_f = dp_f[0:len(dp_f)-1]
            dp_b = dp_b[0:len(dp_b)-1]
            dQ_f = dQ_f[0:len(dQ_f)-1]
            dQ_b = dQ_b[0:len(dQ_b)-1]

        #reshape dp_f, dp_b and dQ_f, dq_b
        dp_f.shape = (divisor,dlength)
        dp_b.shape = (divisor,dlength)
        dQ_f.shape = (divisor,dlength)
        dQ_b.shape = (divisor,dlength)

        matrix = np.tril(np.ones((dlength,dlength)))
        
        #perform integration/summation for all splitted compartments and add them again
        # the first one
        lastElement = dp_f[0][-1]
        p_f = np.dot(matrix,dp_f[0])
        p_b = np.dot(matrix,dp_b[0])
        q_f = np.dot(matrix,dQ_f[0])
        q_b = np.dot(matrix,dQ_b[0])
        # the rest
        for i in np.arange(1,divisor):  
            
            p_f = np.append(p_f, np.dot(matrix,dp_f[i])+p_f[-1])
            p_b = np.append(p_b, np.dot(matrix,dp_b[i])+p_b[-1]) 
            q_f = np.append(q_f, np.dot(matrix,dQ_f[i])+q_f[-1])
            q_b = np.append(q_b, np.dot(matrix,dQ_b[i])+q_b[-1])

        if isPrime:
            #append the last element if it was removed before
            p_f = np.append(p_f, dp_fPrime+p_f[-1])
            p_b = np.append(p_b, dp_bPrime+p_b[-1]) 
            q_f = np.append(q_f, dQ_fPrime+q_f[-1])
            q_b = np.append(q_b, dQ_bPrime+q_b[-1])

    else:
        ## caculate and set d_f, p_b and Q_f, q_b by summation over time to the iven timpoint tp
        matrix = np.tril(np.ones((len(dp_f),len(dp_f))))
        p_f = np.dot(matrix,dp_f)
        p_b = np.dot(matrix,dp_b)
        q_f = np.dot(matrix,dQ_f)
        q_b = np.dot(matrix,dQ_b)  
    """
    
    


def minMaxFunction(arrayToEvaluate,timeValues=np.array([]),delta=0.025, seperateMinMax = False):
    """finds mins and max of all arrays
    if timeValues == 0 then indices of min and max will be returned instead
    returns if  seperateMinMax = False:
        2 arrays with [y-values (min+max)] [ x-values (min+max)] as they occure in the x-direction
    returns if  seperateMinMax = True: 
        2 arrays with maxs[[max1 y,max1 x],[max2 y,max2 x]...] mins[...]] """

    maxtab = []
    mintab  = []
    minMaxAmp  = []
    minMaxTime = []

    lenAtE = len(arrayToEvaluate)
    x = np.arange(lenAtE)
    if len(timeValues) == 0:
        timeValues = x
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    
    # adjust delta as percentage of min-max range of signal to adapt to different signal magnitudes
    delta = (np.max(arrayToEvaluate)-np.min(arrayToEvaluate))*delta
    
    # determine to look for min or max
    zeroNaNinputdetection = True
    startValue = arrayToEvaluate[0]
    
    for num in arrayToEvaluate:
        if startValue > num and startValue > (num + delta/2.0):
            # check for max first
            lookForMax = False
            zeroNaNinputdetection=False
            break                
        elif startValue < num and startValue < (num - delta/2.0):
            # check for max first
            lookForMax = True
            zeroNaNinputdetection = False
            break
    if zeroNaNinputdetection: 
        print 'WARNING: minMax Function: zero or NaN input detected'    
        if   seperateMinMax == False:
            return minMaxAmp, minMaxTime
        else:
            return np.array(maxtab),np.array(mintab)

    # scroll through the signal
    for i in np.arange(lenAtE):
        curr = arrayToEvaluate[i]

        if curr > mx:
            mx = curr
            mxpos = x[i]
        if curr < mn:
            mn = curr
            mnpos = x[i]
            
        if lookForMax:
            if curr < mx-delta:
                # save maximum
                timeT = timeValues[mxpos]
                maxtab.append( [mx, timeT] )
                
                minMaxAmp.append(mx)
                minMaxTime.append(timeT)
                  
                mn = curr
                mnpos = x[i]
                lookForMax = False
        else:
            if curr > mn+delta:
                #save minimum
                timeT = timeValues[mnpos]
                mintab.append([mn, timeT])
                
                minMaxAmp.append(mn)
                minMaxTime.append(timeT)                       
                
                mx = curr
                mxpos = x[i]
                lookForMax = True
                
    if   seperateMinMax == False:
        return minMaxAmp, minMaxTime
    else:
        return np.array(maxtab),np.array(mintab)
    
    
def calculateWaveShoulderPoint(waveData, timeData, startPointT, endPointT = None):
    """
    This method calculates the point inside [startPointT,endPointT] of the wave data 
    where the connection line between it and starPointX has the steepest slope 
    """
    # find start index through time
    startIndex = None
    endIndex   = None
    startIndex = np.where(timeData == startPointT)
    if startIndex == None:
        print "ERROR calculateWaveShoulderPoint: could not find startIndex! return"; return
    if endPointT: endIndex   = np.where(timeData == endPointT)
    else: endPointT = len(timeData)
    print startIndex
    print endIndex
    slope = 0
    slopeIndex = None
    for i in xrange(startIndex[0],endIndex[0]):
        m = (waveData[i]- waveData[startPointT]) / (timeData[i] - timeData[startPointT])
        if m > slope:
            slope = m
            slopeIndex = i
    
    if slopeIndex:
        return [waveData[slopeIndex],timeData[slopeIndex]]
    
def calculateWaveShoulderPointTangent(waveData, timeData, startPointT, endPointT = None):
    """
    This method calculates the point inside [startPointX,endPointX] of the wave data 
    where the tangent has the steepest slope  
    """
    raise NotImplementedError("ERROR: Function renamed to calculatePointOfInflection: this function will be removed soon!")
    
def calculatePointOfInflection(waveData, timeData, startPointT, endPointT = None):
    """
    This method calculates the point inside [startPointX,endPointX] of the wave data 
    where the tangent has the steepest slope  
    """
    # find start index through time
    startIndex = None
    endIndex   = None
    startIndex = np.where(timeData == startPointT)
    if startIndex == None:
        print "ERROR calculateWaveShoulderPoint: could not find startIndex! return"; return
    if endPointT != None: 
        endIndex = np.where(timeData == endPointT)
    else:
        endIndex = [len(timeData)-1]
        
    slope = 0
    slopeIndex = None
    for i in xrange(startIndex[0],endIndex[0]):
        m = (waveData[i+1]- waveData[i-1]) / (timeData[i+1] - timeData[i-1])
        if m > slope:
            slope = m
            slopeIndex = i
    
    if slopeIndex:
        return [waveData[slopeIndex],timeData[slopeIndex]]
    
def calculateReflectionCoefficientBifurcations(vascularNetwork, solutionDataSet = None):
    """
    This method calculates the reflection coefficients of a vascularNetwork by
    traversing it as binary tree.
    Inputparameter:
        vascularNetwork (instance::classVascularNetwork())
        solutionDataSet = None (optional)
        If a solutionDataSet is given, the coefficients are calculated transient taking the 
        soultion data into account.
        Else, it the coefficients are calculated with the initial condictions.
    Retrun:  Rf_transient= {}  == { str(bifurcation) : [Rf_mother, Rf_leftDaughter, Rf_rightDaughter] }
    """
    
    
    ## firstTraverse the tree and create a list of all Bifurcations
    bifurcationList = []
    viz = []
    #find root
    root = vascularNetwork.root[0]        
    # add root to the viz vessels if root has daughters:
    if vascularNetwork.vessels[root].leftDaugther != None:
        viz.append(root)
    # loop through tree until all daughters are conected
    while len(viz) != 0:
        # get the mother vessel (already added) and add its daughters
        motherVessel = viz.pop(0)
        # find left daughter
        leftDaughter  = vascularNetwork.vessels[motherVessel].leftDaugther
        rightDaughter = vascularNetwork.vessels[motherVessel].rightDaugther
        #append the mother to the calc list
        curCalcList = [motherVessel]
        
        if leftDaughter != None:
            #append the left daughter to the calc list
            curCalcList.append(leftDaughter)
           
            if rightDaughter != None:
                #append the left daughter to the calc list
                curCalcList.append(rightDaughter)
            else:
                curCalcList.append(None)
            # check if leftDaughter has also daughters 
            if vascularNetwork.vessels[leftDaughter].leftDaugther != None:
                viz.append(leftDaughter)
                
            if rightDaughter != None:
                # check if rightDaughter has also daughters 
                if vascularNetwork.vessels[rightDaughter].leftDaugther != None:
                    viz.append(rightDaughter)
        bifurcationList.append(curCalcList)
          
    # static // initial condition
    if solutionDataSet == None:
        Rf_static = {}
        rho = vascularNetwork.globalFluid['rho'] # local rho!!
        print 'Rf at bifurcations for initial condition'
        for bifurcation in bifurcationList:
            mother = bifurcation[0]
            leftDaughter = bifurcation[1]
            rightDaughter = bifurcation[2]
            # bifurcation    
            Y_m = C*c# = (vascularNetwork.vessels[mother].r[-1]**2*np.pi) / (rho*vascularNetwork.vessels[mother].c[-1])
            Y_lD # =  (vascularNetwork.vessels[leftDaughter].r[0]**2*np.pi) / (vascularNetwork.globalFluid[rho]*vascularNetwork.vessels[leftDaughter].c[0])
            if rightDaughter != None:
                Y_rD # = (vascularNetwork.vessels[rightDaughter].r[0]**2*np.pi) / (vascularNetwork.globalFluid[rho]*vascularNetwork.vessels[rightDaughter].c[0])          
            else :
                Y_rD = 0.0
            Rf_static[bifurcation] = (Y_m-Y_lD-Y_rD) / (Y_m+Y_lD+Y_rD)
             
        return Rf_static
    else:
    # over time // after simulation
        print 'Rf at bifurcations transient'
        Rf_transient = {}
        rho = vascularNetwork.globalFluid['rho'] # change to local
        for bifurcation in bifurcationList:
            mother = bifurcation[0]
            leftDaughter = bifurcation[1]
            rightDaughter = bifurcation[2]
            
            motherAsol = solutionDataSet['Area'][mother][:,[-1]]
            motherCsol = solutionDataSet['WaveSpeed'][mother][:,[-1]]
            leftDaughterAsol = solutionDataSet['Area'][leftDaughter][:,[0]]
            leftDaughterCsol = solutionDataSet['WaveSpeed'][leftDaughter][:,[0]]
            
            Y_m  = (motherAsol) / (rho*motherCsol)
            Y_lD  =  (leftDaughterAsol) / (rho*leftDaughterCsol)
            
            if rightDaughter != None:
                rightDaughterAsol = solutionDataSet['Area'][rightDaughter][:,[0]]
                rightDaughterCsol = solutionDataSet['WaveSpeed'][rightDaughter][:,[0]]
        
                Y_rD  =  (rightDaughterAsol) / (rho*rightDaughterCsol)
            else :
                Y_rD = 0.0
                
            Rf_transient[str(bifurcation)] = [(Y_m-Y_lD-Y_rD) / (Y_m+Y_lD+Y_rD),(Y_lD-Y_m-Y_rD) / (Y_m+Y_lD+Y_rD),(Y_rD-Y_lD-Y_m) / (Y_m+Y_lD+Y_rD)]
            
        return Rf_transient
