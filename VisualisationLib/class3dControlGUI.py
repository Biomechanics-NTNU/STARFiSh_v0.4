import pyglet
from pyglet import *

from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

import sys,os,io
cur = os.path.dirname( os.path.realpath( __file__ ) )

class Control(pyglet.event.EventDispatcher):
    x = y = 0
    width = height = 10
    
    dx = 40 #  horizontal distance between buttons
    dy = 40 # vertical distance between the buttons 
    def __init__(self, parent):
        super(Control, self).__init__()
        self.parent = parent
        
    
    def hit_test(self, x, y):
        return (self.x < x < self.x + self.width and  
                self.y < y < self.y + self.height)
    
    def capture_events(self):
        self.parent.push_handlers(self)
    
    def release_events(self):
        self.parent.remove_handlers(self)

class Button(Control):
    
    def __init__(self,*args, **kwargs):
        super(Button, self).__init__(*args, **kwargs)
        
        self.mouseXY = []
        
        self.buttonState   = False
        self.troggleButton = False
        
        self.charged = False
        self.hovered = False
        
        self._description = pyglet.text.Label('', anchor_x='left', anchor_y='top',color=(0, 0, 0, 255), font_size = 10 )

        self.enableDescription  = False
        self.hasDescription     = False
        
    def draw(self):
        
        drawed = False
        
        if self.hovered and self.charged == False:
            self.draw_rect(self.x, self.y, self.width, self.height, [0.0,0.8,1.0])
            drawed = True
            
        if self.troggleButton and self.buttonState == True:
            self.draw_rect(self.x, self.y, self.width, self.height, [84,84,84])
            drawed = True
        
        if self.charged:
            self.draw_rect(self.x, self.y, self.width, self.height, [0.0,0.4,1.0])
            drawed = True
            
        if drawed == False:
            self.draw_rect(self.x, self.y, self.width, self.height, [0.75,0.75,0.75])
        self.draw_label()

    def drawHelp(self):
        if self.hovered and self.charged == False:
            if self.enableDescription:
                self.draw_description(*self.mouseXY)

    def draw_description(self, x, y):
                    
        x = x-2
        y = y-15
        
        descriptionWidth = self._description.content_width+20
        descriptionHeigth = self._description.content_height+4
                    
        if y- descriptionHeigth < 0 :
            y = y + descriptionHeigth + 18
                
        if x + descriptionWidth > self.parent.width:
            x = x - 4 - descriptionWidth
        
        descriptionColor = [0.95,0.95,0.85]
        glBegin(GL_QUADS)
        glColor3f(*descriptionColor)
        glVertex2f(x, y ) #up left
        glVertex2f(x + descriptionWidth, y ) #up right
        glVertex2f(x + descriptionWidth, y- descriptionHeigth) #down right
        glVertex2f(x, y -descriptionHeigth) #down left        
        glEnd()
        
        self._description.x = x+10
        self._description.y = y-2
        self._description.draw()    

    def on_mouse_press(self, x, y, button, modifiers):
        self.capture_events()
        self.charged = True
            
    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        self.charged = self.hit_test(x, y)
    
    def on_mouse_release(self, x, y, button, modifiers):
        self.state()
        self.release_events()
        if self.hit_test(x, y):
            self.dispatch_event('on_press')
        self.charged = False
        
    def state(self):
        if self.buttonState == True: self.buttonState = False
        elif self.buttonState == False : self.buttonState = True 
        
    def draw_rect(self, x, y, width, height, color):
        
        c = 2
        glBegin(GL_QUADS)
        glColor3f(0.5,0.5,0.5)
        glVertex2f(x+c, y-c)
        glVertex2f(x+c + width, y-c)
        glVertex2f(x+c + width, y + height-c)
        glVertex2f(x+c, y + height-c)
        glEnd()
        
        glBegin(GL_QUADS)
        glColor3f(*color)
        glVertex2f(x, y)
        glVertex2f(x + width, y)
        glVertex2f(x + width, y + height)
        glVertex2f(x, y + height)
        glEnd()
        
    def set_description(self,description):
        self._description.text = description
        self.hasDescription = True
    def displayDescription(self):
        self.enableDescription = True
    def displayNotDescription(self):
        self.enableDescription = False

    description = property(lambda self: self._description.text, set_description) 
        
Button.register_event_type('on_press')

class TextButton(Button):
    def __init__(self, *args, **kwargs):
        super(TextButton, self).__init__(*args, **kwargs)
        self._text = pyglet.text.Label('', anchor_x='center', anchor_y='center',color=(0, 0, 0, 255), font_size = 10 )

    def draw_label(self):
        
        self._text.x = self.x + self.width / 2
        self._text.y = self.y + self.height / 2
        self._text.draw() 
        
    def set_text(self, text):
        self._text.text = text
        
    text = property(lambda self: self._text.text, set_text)
    
class ImageTextButton(TextButton):
    def __init__(self, *args, **kwargs):
        
        super(ImageTextButton, self).__init__(*args, **kwargs)
        self._text = pyglet.text.Label('', anchor_x='left', anchor_y='center',color=(0, 0, 0, 255),font_size = 10 )

        self._image = None
        self._sprite = None
    
    def draw_label(self):
        if self._sprite: self._sprite.draw()
        self._text.x = self.x + 22
        self._text.y = self.y + self.height / 2
        self._text.draw()

    def set_text(self, text):
        self._text.text = text

    def loadImage(self,image):
        self._image = pyglet.image.load(image)
        self._sprite = pyglet.sprite.Sprite(self._image)
        self._sprite.set_position(self.x+2 , self.y + 2)
        self._sprite.scale =  (self.height-4.)/self._image.height
        
    text  = property(lambda self: self._text.text, set_text)
    image = property(lambda self: self._image, loadImage)

class Slider(Button):
    THUMB_WIDTH = 6
    THUMB_HEIGHT = 10
    GROOVE_HEIGHT = 2
    
    value = 0
    min   = 0
    max   = 0
        
    def draw(self):
        center_y = self.y + self.height / 2
        self.draw_rect(self.x, center_y - self.GROOVE_HEIGHT / 2, self.width, self.GROOVE_HEIGHT, [0.75,0.75,0.75])
        pos = self.x + self.value * self.width / (self.max - self.min)
        self.draw_rect(pos - self.THUMB_WIDTH / 2, center_y - self.THUMB_HEIGHT / 2, self.THUMB_WIDTH, self.THUMB_HEIGHT, [0.75,0.75,0.75])

    def draw_rect(self, x, y, width, height, color):
        
        c = 2
        glBegin(GL_QUADS)
        glColor3f(0.5,0.5,0.5)
        glVertex2f(x+c, y-c)
        glVertex2f(x+c + width, y-c)
        glVertex2f(x+c + width, y + height-c)
        glVertex2f(x+c, y + height-c)
        glEnd()
            
        glBegin(GL_QUADS)
        glColor3f(*color)
        glVertex2f(x, y)
        glVertex2f(x + width, y)
        glVertex2f(x + width, y + height)
        glVertex2f(x, y + height)
        glEnd()
    
        glBegin(GL_LINE_LOOP)
        glColor3f(0.0,0.0,0.0)
        glVertex2f(x, y-c) #down left
        glVertex2f(x+c + width, y-c) #down right
        glVertex2f(x+c + width, y + height) #up right
        glVertex2f(x, y + height) #up left
        glEnd()
    
    def coordinate_to_value(self, x):
        return float(x - self. x) / self.width * (self.max - self.min) + self.min
        
    def on_mouse_press(self, x, y, button, modifiers):
        self.value = self.coordinate_to_value(x)
        self.capture_events()
        self.dispatch_event('on_begin_scroll')
        self.dispatch_event('on_change', self.value)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        self.value = min(max(self.coordinate_to_value(x), self.min), self.max)
        self.dispatch_event('on_change', self.value)
        
    def on_mouse_release(self, x, y, button, modifiers):
        self.release_events()
        self.dispatch_event('on_end_scroll')

        
Slider.register_event_type('on_begin_scroll')
Slider.register_event_type('on_end_scroll')
Slider.register_event_type('on_change')


class Dash(object):
    
    def __init__(self, x, y, hight, width):
        # user
        self.x = x
        self.y = y
        self.heightD = hight 
        # default
        self.widthD = width

    def draw(self):
        glBegin(GL_QUADS)
        glColor3f(0.4,0.4,0.4)
        glVertex2f(self.x, self.y) #down left
        glVertex2f(self.x+ self.widthD, self.y) #down right
        glVertex2f(self.x+ self.widthD, self.y + self.heightD) #up right
        glVertex2f(self.x, self.y + self.heightD) #up left
        glEnd()


#########################################################################################################################################################
class ControlWindow(pyglet.window.Window):
    def __init__(self, parent):
        self.parent = parent
        
        width, height = 400, 150
        super(ControlWindow, self).__init__(caption="Control Window",resizable=False,height= height, width= width)
        self.set_size(width, height)
        
        pyglet.gl.glClearColor(0.75,0.75,0.75,1)
        buttonDistanceToWindow = 10
        buttonHeightWidth = 30
        L, l = 400-2*buttonDistanceToWindow,2
                
        distanceBetweenButtonsX = buttonHeightWidth+20
        distanceBetweenButtonsY = buttonHeightWidth+20
                
        d = 100
        self.pathToIcons = ''.join([cur,'/icons3dViz/'])
                        
        # buttons for coordinate axis
               
        #### Button group1 : Buttons for coordinale axis: four buttons <-------------------------------------start
        buttonGroupViewX = buttonDistanceToWindow
        buttonGroupViewY = buttonDistanceToWindow+distanceBetweenButtonsY
        #self.viewXY = ImageTextButton(self)
        self.viewXY = TextButton(self)
        self.viewXY.x = buttonDistanceToWindow
        self.viewXY.y = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.viewXY.height = buttonHeightWidth
        self.viewXY.width = buttonHeightWidth
        self.viewXY.on_press = self.onPressViewXY
        self.viewXY.text = 'XY'
        self.viewXY.description = 'View XY plane (1)'
        #self.Buttonxplus.image = ''.join([self.pathToIcons,IconName.png'])
        #self.viewXY.image = ''.join([self.pathToIcons,'xyi.png'])
        
        
        #self.viewYZ = ImageTextButton(self)
        self.viewYZ = TextButton(self)
        self.viewYZ.x = buttonDistanceToWindow+distanceBetweenButtonsX
        self.viewYZ.y =buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.viewYZ.height = buttonHeightWidth
        self.viewYZ.width = buttonHeightWidth
        self.viewYZ.on_press = self.onPressViewYZ
        self.viewYZ.text = 'YZ'
        self.viewYZ.description = 'View YZ plane (2)'
        #self.Buttonxminus.image = 'b.jpeg'height
        #self.viewYZ.image = ''.join([self.pathToIcons,'xyi.png'])
        
        #self.viewXZ = ImageTextButton(self)
        self.viewXZ = TextButton(self)
        self.viewXZ.x =  buttonDistanceToWindow+distanceBetweenButtonsX*2
        self.viewXZ.y = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.viewXZ.height = buttonHeightWidth
        self.viewXZ.width = buttonHeightWidth
        self.viewXZ.on_press = self.onPressViewXZ
        self.viewXZ.text = 'XZ'
        self.viewXZ.description = 'View XZ plane (3)'
        #self.Buttonyplus.image = 'b.jpeg'
        #self.viewXZ.image = ''.join([self.pathToIcons,'xyi.png'])
        
        #self.viewXYZ = ImageTextButton(self)
        self.viewXYZ = TextButton(self)
        self.viewXYZ.x =   buttonDistanceToWindow+distanceBetweenButtonsX*3
        self.viewXYZ.y = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.viewXYZ.height = buttonHeightWidth
        self.viewXYZ.width = buttonHeightWidth
        self.viewXYZ.on_press = self.onPressViewXYZ
        self.viewXYZ.text = 'XYZ'
        self.viewXYZ.description = 'Initial view (4)'
        #self.Buttonyminus.image = 'b.jpeg'
        #self.viewXYZ.image = ''.join([self.pathToIcons,'xyz.png'])
        
        
        self.simMovie = ImageTextButton(self)
        #self.simMovie = TextButton(self)
        self.simMovie.x =  buttonDistanceToWindow+distanceBetweenButtonsX*4
        self.simMovie.y = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.simMovie.height = buttonHeightWidth
        self.simMovie.width = buttonHeightWidth 
        self.simMovie.on_press = self.onPressMovie
        #self.simMovie.text = 'MV'
        self.simMovie.description = 'Create movie (M)'
        self.simMovie.image = ''.join([self.pathToIcons,'Movie1.png'])
        
        self.photo                  = ImageTextButton(self)
        self.photo.x                = buttonDistanceToWindow+distanceBetweenButtonsX*5
        self.photo.y                = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.photo.height           = buttonHeightWidth
        self.photo.width            = buttonHeightWidth
        self.photo.on_press         = self.onPressViewPhoto
        #self.photo.text = 'XT1'
        self.photo.description       = 'Save screen shot (S)'
        self.photo.image = '' .join([self.pathToIcons,'Camera.png'])
        
        self.backgroundColor                = TextButton(self)
        self.backgroundColor.x              = buttonDistanceToWindow+distanceBetweenButtonsX*6
        self.backgroundColor.y              = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.backgroundColor.height         = buttonHeightWidth
        self.backgroundColor.width          = buttonHeightWidth
        self.backgroundColor.on_press       = self.onPressBackgroundColor
        self.backgroundColor.text           = "B"
        self.backgroundColor.description    = 'change background color (B)'
        
        

        
        self.helpButton = ImageTextButton(self)
        #self.helpButton = TextButton(self)
        self.helpButton.x =buttonDistanceToWindow+distanceBetweenButtonsX*7
        self.helpButton.y = buttonDistanceToWindow+distanceBetweenButtonsY*4-d
        self.helpButton.height =buttonHeightWidth
        self.helpButton.width  = buttonHeightWidth
        self.helpButton.on_press = self.onPressedHelpButton
        #self.helpButton.text = '?'
        self.helpButton.description = 'Help'
        self.helpButton.image = "".join([self.pathToIcons,'q.png'])
        #### End of buttons group 1 <-------------------------------------------------------End
        
                
        ## Buttons groupe 2: 3 buttons. for play/pause, reset and moovie <-----------------------------------------start
        #self.simPlayPause = TextButton(self)
        self.simPlayPause = ImageTextButton(self)
        self.simPlayPause.x = buttonDistanceToWindow 
        self.simPlayPause.y = buttonDistanceToWindow+distanceBetweenButtonsY*3-d
        self.simPlayPause.height = buttonHeightWidth
        self.simPlayPause.width = buttonHeightWidth
        self.simPlayPause.on_press = self.onPressPlayPause
        #self.simPlayPause.text = 'PP'
        self.simPlayPause.description = 'Play/Pause (P)'
        self.simPlayPause.image = ''.join([self.pathToIcons,'Play1.png'])
               
               
        
               
        #self.simReset = TextButton(self)
        self.simReset = ImageTextButton(self)
        self.simReset.x = buttonDistanceToWindow+distanceBetweenButtonsX
        self.simReset.y = buttonDistanceToWindow+distanceBetweenButtonsY*3-d
        self.simReset.height =buttonHeightWidth
        self.simReset.width =buttonHeightWidth
        self.simReset.on_press = self.onPressReset
        #self.simReset.text = 'RT'
        self.simReset.description = 'Stop and reset (O)'
        self.simReset.image = ''.join([self.pathToIcons,'reset1.png'])
        
        #self.realTimeSim            = ImageTextButton(self)
        self.realTimeSim = TextButton(self)
        self.realTimeSim.x          =buttonDistanceToWindow+distanceBetweenButtonsX*3
        self.realTimeSim.y          = buttonDistanceToWindow+distanceBetweenButtonsY*3-d
        self.realTimeSim.height     = buttonHeightWidth # buttonHeightWidth is the horizontal distance betwen the buttons
        self.realTimeSim.width      = buttonHeightWidth
        self.realTimeSim.on_press   = self.onPressRealTimeSim
        self.realTimeSim.text = 'RTS'
        self.realTimeSim.description = 'Real time speed (R)'
        #self.realTimeSim.image = "" .join([self.pathToIcons,'Movie1.png'])
        
        #self.halfTimeSim = ImageTextButton(self)
        self.halfTimeSim = TextButton(self)
        self.halfTimeSim.x =buttonDistanceToWindow+distanceBetweenButtonsX*4
        self.halfTimeSim.y = buttonDistanceToWindow+distanceBetweenButtonsY*3-d
        self.halfTimeSim.height =buttonHeightWidth
        self.halfTimeSim.width = buttonHeightWidth
        self.halfTimeSim.on_press = self.onPressHalfTimeSim
        self.halfTimeSim.text = 'HTS'
        self.halfTimeSim.description = 'Half time speed (H)'
        #self.halfTimeSim.image = "".join([self.pathToIcons,'half1.png'])

        ### End buttons  group <-----------------------------------------------------------------------------End
        
        
        ## create sliderSpeed button
        self.sliderSpeed = Slider(self)
        self.sliderSpeed.x = buttonDistanceToWindow+distanceBetweenButtonsX*5
        self.sliderSpeed.y = buttonDistanceToWindow+distanceBetweenButtonsY*3+5-d
        self.sliderSpeed.width = 123
        #self.sliderSpeed.height = 20
        self.sliderSpeed.on_begin_scroll = self.sliderBegingScroll
        self.sliderSpeed.on_end_scroll     = self.sliderEndScroll
        self.sliderSpeed.on_change         = self.onChangeSpeedSlider
        self.sliderSpeed.value = 1
        self.sliderSpeed.min   = 0
        self.sliderSpeed.max   = 1
        self.sliderSpeed.description = "Change visualisation speed"
        ###
        
        ## button group split <------------------------------------------start
        self.split = TextButton(self)
        #self.split = ImageTextButton(self)
        self.split.x = buttonDistanceToWindow
        self.split.y = buttonDistanceToWindow+distanceBetweenButtonsY*2-d
        self.split.height = buttonHeightWidth
        self.split.width = buttonHeightWidth
        self.split.on_press = self.onPressSplit

        self.split.text = 'WS'
        self.split.description = 'Enable/Disable wave split (W)'
        self.split.troggleButton = True
        #self.split.image = "".join([self.pathToIcons,'fbi.png'])
        
        
        self.splitExp = TextButton(self)
        #self.splitExp = ImageTextButton(self)
        self.splitExp.x = buttonDistanceToWindow+distanceBetweenButtonsX
        self.splitExp.y =  buttonDistanceToWindow+distanceBetweenButtonsY*2-d
        self.splitExp.height = buttonHeightWidth
        self.splitExp.width = buttonHeightWidth 
        self.splitExp.on_press = self.onPressSplitExp
        self.splitExp.text = 'L/R'
        self.splitExp.description = 'WS: left right coloring (E)'
        #self.splitExp.image = "".join([self.pathToIcons,'efbi.png'])
        
        
        #self.changeQuantity = ImageTextButton(self)
        self.changeQuantity = TextButton(self)
        self.changeQuantity.x = buttonDistanceToWindow+distanceBetweenButtonsX*2
        self.changeQuantity.y =  buttonDistanceToWindow+distanceBetweenButtonsY*2-d
        self.changeQuantity.height = buttonHeightWidth
        self.changeQuantity.width = buttonHeightWidth 
        self.changeQuantity.on_press = self.onPressChangeQuantity
        self.changeQuantity.text = 'P'
        self.changeQuantity.description = 'Change LUT quantity (Q)'
        #self.changeQuantity.image = "".join([self.pathToIcons,'pi.png'])
        
        self.lookUpTable = TextButton(self)
        #self.lookUpTable = ImageTextButton(self)
        self.lookUpTable.x = buttonDistanceToWindow+distanceBetweenButtonsX*3
        self.lookUpTable.y =  buttonDistanceToWindow+distanceBetweenButtonsY*2-d
        self.lookUpTable.height = buttonHeightWidth
        self.lookUpTable.width = buttonHeightWidth 
        self.lookUpTable.on_press = self.onPressLookUpTable
        self.lookUpTable.text = 'LUT'
        self.lookUpTable.description = 'Change look up table (L)'
        #self.lookUpTable.image = "".join([self.pathToIcons,'luti.png'])
        
        
        self.wallMovement = TextButton(self)
        #self.wallMovement = ImageTextButton(self)
        self.wallMovement.x = buttonDistanceToWindow+distanceBetweenButtonsX*4
        self.wallMovement.y =  buttonDistanceToWindow+distanceBetweenButtonsY*2-d
        self.wallMovement.height = buttonHeightWidth
        self.wallMovement.width = buttonHeightWidth 
        self.wallMovement.on_press = self.onPressWallMovement
        self.wallMovement.text = 'A'
        self.wallMovement.description = 'Enable/disable wall movement (A)'
        self.wallMovement.troggleButton = True
        self.wallMovement.buttonState = True
        #self.wallMovement.image = "".join([self.pathToIcons,'ai.png'])
        
        ## create sliderSpeed button
        self.wallSlider = Slider(self)
        self.wallSlider.x = buttonDistanceToWindow+distanceBetweenButtonsX*5
        self.wallSlider.y = buttonDistanceToWindow+distanceBetweenButtonsY*2+5-d
        self.wallSlider.width = 123
        #self.sliderSpeed.height = 20
        self.wallSlider.on_begin_scroll   = self.sliderBegingScroll2
        self.wallSlider.on_end_scroll     = self.sliderEndScroll2
        self.wallSlider.on_change         = self.onChangeWallSlider
        self.wallSlider.value = 0
        self.wallSlider.min   = 0
        self.wallSlider.max   = 1
        self.wallSlider.description = "Change area amplification"
        
        #### <-------------------------------------------------------------------end
        
        s=5

        ##########<----------------------------------------------------------------------------------------------------------------------
        #self.dash1 = Dash(buttonDistanceToWindow, buttonDistanceToWindow+distanceBetweenButtonsY*4+buttonHeightWidth+s,    l, L)
        self.dash2 = Dash(buttonDistanceToWindow, buttonDistanceToWindow+distanceBetweenButtonsY*3+buttonHeightWidth+s-d,    l, L) #<-------------------------------------------------------------Dashes
        self.dash3 = Dash(buttonDistanceToWindow, buttonDistanceToWindow+distanceBetweenButtonsY*2+buttonHeightWidth+s-d,    l, L)
        self.dash4 = Dash(buttonDistanceToWindow,  buttonDistanceToWindow+distanceBetweenButtonsY+buttonHeightWidth+s-d,    l, L)
        #self.dash5 = Dash(0,distanceBetweenButtonsY+s,l,L)
        #self.dash6 = Dash(L,distanceBetweenButtonsY+s,buttonDistanceToWindow+distanceBetweenButtonsY*3+buttonHeightWidth+l,l)
        #self.dash7 = Dash(0,distanceBetweenButtonsY+s, buttonDistanceToWindow+distanceBetweenButtonsY*3+buttonHeightWidth+l,l)
        #self.dash2 = Dash(uttonDistanceToWindow+distanceBetweenButtonsX*2, buttonDistanceToWindow+distanceBetweenButtonsY*3+buttonHeightWidth,L, l) #<-------------------------------------------------------------Dashes
                
#--------------------------------------------------End
#         
        ## buttons group end <--------------------------------------------end
        self.controls = [   self.realTimeSim,
                            self.halfTimeSim,
                            self.sliderSpeed,
                            self.wallSlider,
                            self.viewXY,
                            self.viewYZ,
                            self.viewXZ,
                            self.viewXYZ,
                            self.backgroundColor,
                            self.simPlayPause,
                            self.simReset,
                            self.photo,
                            self.simMovie,
                            self.split,
                            self.splitExp,
                            self.helpButton,
                            self.changeQuantity,
                            self.lookUpTable,
                            self.wallMovement ]# , self.SCS]
        
        self.dashes = [self.dash2,self.dash3,self.dash4 ] #self.dash1,self.dash5,self.dash6,self.dash7]
        pyglet.gl.glClearColor(0.75,0.75,0.75,1)

        
    def on_draw(self):
        '''
        draws the GUI 
        '''
        self.switch_to()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        for dash in self.dashes:
            dash.draw()    
        for control in self.controls:
            control.draw()
        for control in self.controls:
            #TODO: Try Except Pass should be fixed
            try: control.drawHelp()
            except: pass
    
    def on_close(self):
        return pyglet.event.EVENT_HANDLED
       
    def on_mouse_drag(self,x,y,dx,dy,buttons,modifiers):
        self.changePicture = True

    def on_mouse_press(self, x, y, button, modifiers):
        self.capture_events()
        for control in self.controls:
            if control.hit_test(x, y):
                control.on_mouse_press(x, y, button, modifiers)
        
    def on_mouse_release(self, x, y, button, modifiers):
        self.release_events()    
    
    def on_mouse_motion(self, x, y, dx, dy):
        for control in self.controls:
            if control.hit_test(x, y):
                control.hovered = True                
                control.mouseXY = [x, y]
                
            else: control.hovered = False
                
    ## call functions
    
    ### simulation speed
    def onPressRealTimeSim(self):
        self.dispatch_event('onPressRealTimeSim')
        self.sliderSpeed.value = 1
                    
    def onPressHalfTimeSim(self):
        self.dispatch_event('onPressHalfTimeSim')
        self.sliderSpeed.value = 0.5
    
    def sliderBegingScroll(self):
        pass

    def sliderEndScroll(self):
        pass
        
    def onChangeSpeedSlider(self, value):
        self.dispatch_event('onChangeSpeedSlider', value)
        ## change text
        #self.speedText.text = ':2f'.format(value)
    
    def sliderBegingScroll2(self):
        pass

    def sliderEndScroll2(self):
        pass
        
    def onChangeWallSlider(self, value):
        self.dispatch_event('onChangeWallSlider', value)
        ## change text
        #self.speedText.text = ':2f'.format(value)
                 
    
    ### view
    def onPressViewXY(self):
        self.dispatch_event('onPressViewXY')
        
    def onPressViewYZ(self):
        self.dispatch_event('onPressViewYZ')
        
    def onPressViewXZ(self):
        self.dispatch_event('onPressViewXZ')
        
    def onPressViewXYZ(self):
        self.dispatch_event('onPressViewXYZ')
        
    ### visualisation control
    def onPressPlayPause(self):                                           #<---------------------------------------------------------------------------------
        if self.simPlayPause.buttonState == True:
            self.simPlayPause.image = ''.join([self.pathToIcons,'pause3.png'])
        else:
            self.simPlayPause.image = ''.join([self.pathToIcons,'Play1.png'])
        self.dispatch_event('onPressPlayPause')
        

    def onPressReset(self):
        self.dispatch_event('onPressReset')
        self.simPlayPause.image = ''.join([self.pathToIcons,'Play1.png'])
        self.simPlayPause.buttonState = False
        
    ### save movie screenshot
    def onPressMovie(self):
        if self.simMovie.buttonState == True:
            self.simMovie.image = ''.join([self.pathToIcons,'record.png'])
             
        else:
            self.simMovie.image = ''.join([self.pathToIcons,'Movie1.png'])
            
        self.dispatch_event('onPressMovie')
            
    def onPressViewPhoto(self):
        self.dispatch_event('onPressViewPhoto')

    ### open help window
    def onPressViewHELP(self):
        self.dispatch_event("onPressViewHELP")
                
    ### wall movment
    def onPressWallMovement(self):
        self.wallMovement.buttonState == True  
        self.dispatch_event("onPressWallMovement")
           
    ### open help button
    def onPressedHelpButton(self):
        for control in self.controls:
            if control.hasDescription:
                if self.helpButton.buttonState:
                    control.displayDescription()
                else: control.displayNotDescription()
                            
    ### wave split
    def onPressSplit(self):
        self.split.buttonState == True  
        self.dispatch_event('onPressSplit')
        
    def onPressSplitExp(self):
        self.dispatch_event('onPressSplitExp')
    
    ## LUT 
    def onPressChangeQuantity(self):
        quantity =""
        if self.changeQuantity.buttonState == True:
            self.changeQuantity.text = "Q"
            quantity = "Flow"
        
        else:
            self.changeQuantity.text = "P"
            quantity = "Pressure"
        self.dispatch_event("onPressChangeQuantity",quantity)
        
        
    def onPressLookUpTable(self):
        self.dispatch_event("onPressLookUpTable")
   
    ## windowing
    def onPressBackgroundColor(self):
        self.dispatch_event("onPressBackgroundColor")
  
    ### connecting to mother window with events
    def capture_events(self):
        self.parent.push_handlers(self)
    
    def release_events(self):
        self.parent.remove_handlers(self)
                
ControlWindow.register_event_type('onPressRealTimeSim')
ControlWindow.register_event_type('onPressHalfTimeSim')
ControlWindow.register_event_type('onChangeSpeedSlider')
ControlWindow.register_event_type('onChangeWallSlider')
ControlWindow.register_event_type('onPressViewXY')        
ControlWindow.register_event_type('onPressViewYZ')        
ControlWindow.register_event_type('onPressViewXZ')        
ControlWindow.register_event_type('onPressViewXYZ')        
ControlWindow.register_event_type('onPressPlayPause')
ControlWindow.register_event_type('onPressReset')        
ControlWindow.register_event_type('onPressMovie')            
ControlWindow.register_event_type('onPressViewPhoto')
ControlWindow.register_event_type('onPressSplit')
ControlWindow.register_event_type('onPressSplitExp')
ControlWindow.register_event_type('onPressChangeQuantity')
ControlWindow.register_event_type('onPressLookUpTable')
ControlWindow.register_event_type('onPressWallMovement')
ControlWindow.register_event_type('onPressBackgroundColor')


    