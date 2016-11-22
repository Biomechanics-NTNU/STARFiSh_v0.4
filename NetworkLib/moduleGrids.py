import numpy as np 
    
def uniform(length, radiusA, radiusB, N):
    """
    Function to calculate a Uniform grid 
    with a single radius over the hole vessel length
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """         
    A0 = np.ones(N)*np.pi*radiusA**2.
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    
    return z,dz,A0

def cone(length, radiusA, radiusB, N):
    """
    Function to calculate a Cone grid 
    with a radiusA at the proximal end and radiusB at the distal end
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """         
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    A0 = np.pi*(radiusA + ((radiusB-radiusA)*np.linspace(0.0,1.0,N)))**2
    
    return z,dz,A0

def constriction(length, radiusA, radiusB, N):
    """
    Function to calculate a constricted grid 
    with a radiusA as proximal and distal radius and radiusB
    as constriction radius; the constriction is calculated via a cos-function
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """             
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    R0 = radiusA + (radiusB-radiusA)*(1.0 + np.cos(np.linspace(-1.,1.,N)*np.pi))/2.
    A0 = np.pi*R0**2
    print "WARNING: constricted vessel geometry was never verified and vaildated!"
    return z,dz,A0


def mtStig2016(length, radiusA, radiusB, N):
    """
    Function to calculate a gird based on measured values.
    More details of the values see masterthesis of Stig 2016.
    
    Args:  length (float): vessel-length
           N (int): number of grid points 
    Returns:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """  
    
    print "we use this model yes! mtStig2016 gird"
    exit()
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    
    zg = np.array([0.,1.,2.,3.])
    rg = np.array([10.,12.,13.,14.])
    
    #r = np.interp(z, xp, fp, left, right)
    A0 = rg**2.0 * np.pi
    
    return z,dz,A0