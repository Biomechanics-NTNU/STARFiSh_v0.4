import sys,os,io
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')


#sys.path.append(cur+'/../NetworkLib')
from NetworkLib.classVascularNetwork import VascularNetwork 

import numpy as np
from math  import *

import copy 
import time
import gc

from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import * 

import pyglet
from pyglet.gl import *
from pyglet import *
from pyglet.window import mouse, key
from pyglet.window.key import MOTION_UP,MOTION_DOWN

import matplotlib.animation as animation
from matplotlib._png import read_png
from pylab import *    


class WindowLUT(pyglet.window.Window):
    '''
    This class defines the a GUI window for the
    look up table
    '''
    def __init__(self):
                
        height = 300 
        width  = 150
                
        super(WindowLUT, self).__init__(resizable = False, 
                                        height= height, 
                                        width= width)
        self.set_size(width, height)
        
        
        self.batch = pyglet.graphics.Batch() 
        self.vertexList = None
        
        pyglet.gl.glClearColor(0.75,0.75,0.75,1)
        self.set_caption("LUT")
        
        self.scaleText = []
        for i in range(7):
            self.scaleText.append(pyglet.text.Label('125.00', anchor_x='right', anchor_y='center',color=(0, 0, 0, 255)))
            
        self.quantity = pyglet.text.Label('Pressure', anchor_x='left', anchor_y='bottom',color=(0, 0, 0, 255), bold=1)
        self.unit     = pyglet.text.Label('mmHg', anchor_x='right', anchor_y='bottom',color=(0, 0, 0, 255), font_size = 10 )
        self.unitFactor = 1./133.32
        
        self.createVerticesAndTextPosition()
                      
        self.applyScaleChange([0,133.32])
        
    def on_draw(self):
        '''
        draws the look uptable
        '''
        self.switch_to()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.batch.draw()
        self.drawText()
    
    def drawText(self):
        '''
        
        '''
        for text in self.scaleText:
            text.draw()
        self.quantity.draw()
        self.unit.draw()
        
    def on_close(self):
        return pyglet.event.EVENT_HANDLED
        
    def createVerticesAndTextPosition(self):
        '''
        calculate vertices of the square
        '''
        
        height = 300 
        width  = 150
        
        anchorX   = 10.
        anchorY   = 10.
        spaceTop  = height/8.5
        spaceLeft = width/2.5
        
        self.numberOfVertices = 255
                
        yCoordinates = np.repeat([np.linspace(anchorY,height-spaceTop-anchorY,self.numberOfVertices)],2,axis=0)
        xCoordinates = np.repeat(np.linspace(anchorX,width-spaceLeft-anchorX,2),self.numberOfVertices)
        
        vertices = np.vstack((xCoordinates,yCoordinates.ravel())).T
        
        indices = []
        for i in range(self.numberOfVertices-1):
            indices.extend([i,i+self.numberOfVertices,i+self.numberOfVertices+1,i+1])
        
        vertices = vertices.ravel()
        
        colors = np.repeat([np.array([0, 255, 0])], self.numberOfVertices*2, axis = 0).ravel()        
        
        
        self.vertexList = self.batch.add_indexed(len(vertices)/2,  pyglet.gl.GL_QUADS, None, 
                                                 indices,
                                                 ('v2f/static',vertices),
                                                 ('c3B/dynamic',colors))
        
        ### text position and vertices
        numberOfScale = 7
        textPositionX = width - anchorX
        textWidth = 5
        self.scaleX = np.repeat(textPositionX,numberOfScale)
        self.scaleY = np.linspace(anchorY+2*textWidth,height-spaceTop-anchorY-textWidth,numberOfScale)
        for text,x,y in zip(self.scaleText,self.scaleX,self.scaleY):
            text.x = x
            text.y = y
                
        self.quantity.y = height-height/16.-anchorY
        self.quantity.x = anchorX
        
        self.unit.y = height-height/16.-anchorY
        self.unit.x = width-anchorY
        
    def applyColorChange(self,colors):
        '''
        applies the new colors to the colortable square
        '''
        nColors = len(colors)-1
        colors = colors[np.linspace(0, nColors, self.numberOfVertices).astype(int)]
        colors = np.repeat([colors],2,axis=0).ravel()
        self.vertexList.colors = colors
            
    def applyScaleChange(self,minMax, backwardForward = False):
        '''
        applies the new numbers to look up table
        '''
        if backwardForward:
            for i in range(7):
                self.scaleText[i].text =  ''
            self.scaleText[0].text =  'back.'
            self.scaleText[6].text =  'forw.'
        else:
            textSet = [ '{:.2f}'.format(i) for i in np.linspace(minMax[0]*self.unitFactor,minMax[1]*self.unitFactor,7) ]
            for i in range(7):
                self.scaleText[i].text =  textSet[i]
        
        
    def applyQuantityChange(self,quantity):
        '''
        applies the new quantity and unit to the LUT
        '''
        self.quantity.text  = quantity
        self.unit.text      = {'Pressure':'mmHg','Flow':'ml/s'}[quantity]
        self.unitFactor     = {'Pressure':1./133.32,'Flow':1.e6}[quantity]
        
        
class LUT(object):
    
    def __init__(self, range):
        
        self.range          = range # dict with  range of the quantities [min, max]
        self.rangeWaveSplit = None
        
        self.colorTable      = None
        self.colorTableRange = [0,0] 
        self.colorTableWaves = None
        
        self.quantityLUT = 'Pressure'
        
        # checks which color table is used
        self.changeColorTableBool = 0
        # list of available color tables to scroll through
        self.colorTables = [self.colorTableBlueGreenRed,
                            self.colorTableBlueGreen,
                            self.colorTableBlueRed,
                            self.colorTableGreenRed]
                
        self.waveSplit            = False
        self.exponentialWaveSplit = False
        self.waveSplitExponent    = 1
        self.leftRightColoring    = False
        self.colorTableWaveSplit()
        
                        
        self.windowLUT = WindowLUT()
        # initialize color tables
        self.changeColorTable()
        self.windowLUT.applyScaleChange(self.range[self.quantityLUT])
        self.windowLUT.applyQuantityChange(self.quantityLUT)
        
    def changeColorTable(self):
        '''
        to be changed as it is bad programmed _:)
        '''
        
        if self.waveSplit == False:
            
            tableIterator = range(len(self.colorTables))
            tableIterator.pop(0)
            tableIterator.append(0)
            # set color table and range
            self.colorTable = self.colorTables[self.changeColorTableBool]()
            self.colorTableRange = [0.,len(self.colorTable)-1]
            # update color table walker
            self.changeColorTableBool = tableIterator[self.changeColorTableBool]
            # apply it to LUT window
            self.updateLUTWindow()
            
        else:
            #activate exponential function for the wavesplit look up table
            self.waveSplitExponent = [1,2,3,4,5,6,1][self.waveSplitExponent]
            print " wave split exponent is set to: ",self.waveSplitExponent
    def updateLUTWindow(self):
        '''
        updates the LUT window with all current quantity, range and colors
        '''        
        if self.waveSplit == False:
            self.windowLUT.applyQuantityChange(self.quantityLUT)
            self.windowLUT.applyColorChange(self.colorTable)
            self.windowLUT.applyScaleChange(self.range[self.quantityLUT])
        else:
            self.windowLUT.applyQuantityChange(self.quantityLUT)
            self.windowLUT.applyScaleChange([],True)
            self.windowLUT.applyColorChange(self.colorTableWaves)
    
    def changeLUTquantity(self,newQuantitiy):
        '''
        change the quantity to new quantity ['Pressure','Flow']
        '''
        self.quantityLUT = newQuantitiy
        self.updateLUTWindow()
                                      
    def enableWaveSplitLookUp(self):
        '''
        activate / deactivate wave split look up
        '''
        self.waveSplit = [True,False][self.waveSplit]
        self.updateLUTWindow()
        self.waveSplitExponent = 1
           
    def enableLeftRightColoring(self):
        self.leftRightColoring = [True,False][self.leftRightColoring]
                
    def getColors(self,inputArray):
        '''
        creates index array of inputArray
        and returns colors depending on these indices
        '''
        min = self.range[self.quantityLUT][0]
        max = self.range[self.quantityLUT][1]
        maxRange = self.colorTableRange[1]
        indexArray = (inputArray-min)/(max-min)*maxRange
        
        return self.colorTable[indexArray.astype(int)]
    
    def getColorsWaveSplit(self,forwardInput, backwardInput, rangeWaveSplit):
        '''
        creates index array of forward and backward Input arrays
        and returns colors depending on these indices
        '''
        
        leftRightColoringFactor = 1
        
        ## forward
        min1 = rangeWaveSplit[self.quantityLUT][0][0]
        max1 = rangeWaveSplit[self.quantityLUT][0][1] 
        ## backward
        min2 = rangeWaveSplit[self.quantityLUT][1][0]
        max2 = rangeWaveSplit[self.quantityLUT][1][1] 
        
        # exponent reaching from 1,2,3,4 increasing the function
        exponent = self.waveSplitExponent
        ## forward
        min1 = min1*abs(min1)**(exponent-1)
        max1 = max1*abs(max1)**(exponent-1)
        ## backrard
        min2 = min2*abs(min2)**(exponent-1)
        max2 = max2*abs(max2)**(exponent-1)
                   
        forwardInput  = forwardInput*abs(forwardInput)**(exponent-1)
        backwardInput = backwardInput*abs(backwardInput)**(exponent-1)
                                        
        ## coloring backward forward on left and right
        if self.leftRightColoring == True:
            maxRange = len(self.colorTableBlueRed())-1
            indexArray1 = (forwardInput-min1)/(max1-min1)*maxRange
            indexArray2 = (backwardInput-min2)/(max2-min2)*maxRange
            
            colors1 = self.colorTableBlueRed()[indexArray1.astype(int)]
            colors2 = self.colorTableBlueGreen()[indexArray2.astype(int)]
            
            colors = np.hstack([colors1,colors2]).reshape(len(colors2)*2,3)
            colors[colors<0] = 0
            
            leftRightColoringFactor = 0.5
            
        else: ## coloring the hole vessel
            maxRange = len(self.colorTableSplit1)-1
            indexArray1 = (forwardInput-min1)/(max1-min1)*maxRange
            indexArray2 = (backwardInput-min2)/(max2-min2)*maxRange
            
            colors1 = self.colorTableSplit1[indexArray1.astype(int)]
            colors2 = self.colorTableSplit2[indexArray2.astype(int)]
            colorsM = self.colorTableMeanColor[np.zeros_like(indexArray1).astype(int)]
            
            colors = colors1+colors2+colorsM
            colors[colors<0] = 0
                        
        return colors,leftRightColoringFactor  
    
    def colorTableBlueRed(self):
        '''
        Creates color table from blue to magenta to red        
        '''
        nColorsChanges = 2
        colors = np.zeros((255*nColorsChanges+1,3))
        #start with blue and go to red
        colors[0][2] = 255
        # 1. add red
        for i in range(256):
            colors[i][0] = i # red
            colors[i][2] = 255 # blue
        # 2. remove blue
        for i in range(256):
            # red
            colors[i+255][0] = 255
            colors[i+255][2] = 255-i # blue
                        
        return colors.astype(int)
        
    def colorTableGreenRed(self):
        '''
        creates color table from blue to cyan to green to yellow to red
        '''
        nColorsChanges = 2 # blue to green (2x) green to red (2x)
        colors = np.zeros((255*nColorsChanges+1,3))
        #start with green and go to red
        colors[0][1] = 255
        # 1. add red
        for i in range(256):
            colors[i][0] = i    # red
            colors[i][1] = 255   # green
        # 2. remove green
        for i in range(256):
            colors[i+255][0] = 255    # red
            colors[i+255][1] = 255-i  # green   
                
        return colors.astype(int)
        
    def colorTableBlueGreen(self):
        '''
        creates color table from blue to cyan to green
        '''        
        nColorsChanges = 2 # blue to green (2x) green to red (2x)
        colors = np.zeros((255*nColorsChanges+1,3))
        #start with blue and go to red
        colors[0][2] = 255
        # 1. add green
        for i in range(256):
            colors[i][1] = i # green
            colors[i][2] = 255 # blue
        # 2. remove blue and get only green
        for i in range(256):
            colors[i+255][1] = 255   # green
            colors[i+255][2] = 255-i # blue
            
        return colors.astype(int)
        
    def colorTableWaveSplit(self):
        ## create addatove colorTable for forward
        blueRed   = np.zeros((255+1,3))
        for i in range(256):
            blueRed[i][0] =  i # red
            blueRed[i][2] = -i # blue
        ## create addative colorTable for backward
        blueGreen = np.zeros((255+1,3))
        for i in range(256):
            blueGreen[i][1] =  i # green
            blueGreen[i][2] = -i # blue
        ## create standart blue
        standardBlue = np.array([[0,0,255]])
        
        self.colorTableMeanColor = standardBlue.astype(int)
        self.colorTableSplit1 = blueRed.astype(int)
        self.colorTableSplit2 = blueGreen.astype(int)
                
        self.colorTableWaves = np.append(np.flipud(self.colorTableBlueGreen()),self.colorTableBlueRed(), axis=0).astype(int)
        
    def colorTableBlueGreenRed(self):
        '''
        creates color table from blue to cyan to green
        '''
        colors = np.append(self.colorTableBlueGreen(), self.colorTableGreenRed(), axis=0).astype(int)
        
        return colors.astype(int)
                
    def update(self,updateDict):
        '''
        updates the LUT data using a dictionary in form of 
        updateDict = {'variableName': value}
        '''
        for key,value in updateDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'Warning: LUT.update(): wrong key: %s, could not set up LUT' %key
