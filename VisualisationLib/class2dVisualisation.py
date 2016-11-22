#!/usr/bin/env python
import gtk
import gobject
import matplotlib
matplotlib.use(u'GTKAgg')
import matplotlib.pyplot as plt   
from matplotlib.figure import Figure
# uncomment to select /GTK/GTKAgg/GTKCairo
# from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
# from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle

import sys, os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

#sys.path.append(cur + '/../UtilityLib')
import UtilityLib.processing as mProc
#from processing import nonLinearWaveSplitting

import UtilityLib.moduleStartUp as moduleStartUp
import UtilityLib.moduleFilePathHandler as mFPH
import UtilityLib.moduleXML as mXML


import numpy as np
import pprint
from copy import deepcopy as copy

class Visualisation2DPlotWindowAdjustValues(gtk.Window):
    
    def __init__(self, parent, dictVariableName, actualValuesDict, inititalValuesDict, variableTypes = None, variableFixedChoice = None):
        '''
        creates a window to update data written in a dictionary:
        
        dictVariableName    := <str> name of the dictionary in the parent class (used to update dict)
        actualValuesDict    := <dict> containing the data to be updated form {'variableName': [ item, item ...]}
        inititalValuesDict  := <dict> as actualValuesDict with initial values to reset initial values
        variableTypes       := <dict> as actualValuesDict but indicating the type of the variable which is used with e.g. 'str'
        variableFixedChoice := <dict> as actualValuesDict prescribe a fixed option for some variables, e.g. {'variableName': [ ['-',':'], ..]}
        
        for each item of each variable in the dictionary an entry element is created,
        if no type is defined (i.e. variableTypes == None) the variable is converted to float
        if a variableFixedChoice is defined for the variable a comboBox is created instead of a entry
        
        '''       
        super(Visualisation2DPlotWindowAdjustValues, self).__init__()
        
        self.set_resizable(False)
        #self.set_position(gtk.WIN_POS_CENTER)
                                
        self.parentWindow        = parent
        self.dictVariableName    = dictVariableName
        self.inititalValuesDict  = copy(inititalValuesDict)
        self.actualValuesDict    = actualValuesDict
        self.variableTypes       = variableTypes
        self.variableFixedChoice = variableFixedChoice
        
        vbox = gtk.VBox(False, 1)
               
        buttonUpdate = gtk.Button("update")
        buttonUpdate.connect("clicked", self.on_clickedUpdate)
        buttonUpdate.set_size_request(120,30)
                
        buttonReset = gtk.Button("reset")
        buttonReset.connect("clicked", self.on_clickedReset)
        buttonReset.set_size_request(120,30)
        
        hbox = gtk.HBox(False, 10) 
        hbox.pack_start(buttonUpdate, fill=False, expand=True)  
        hbox.pack_start(buttonReset, fill=False, expand=True)   
        vbox.pack_start(hbox, expand=True, fill=False)
                
        self.entriesCombo = {}
        
        countHeight = 1
        countWidth = 1
        
        keys = actualValuesDict.keys()
        keys.sort()
        for key in keys:
            hbox = gtk.HBox(False, 10) 
            label = gtk.Label(key)
            hbox.pack_start(label, fill=False, expand=True)     
            self.entriesCombo[key] = []
            
            countHeight = countHeight +1
            
            values = actualValuesDict[key]
            countWidth  = max(countWidth,len(values))
            for index,value in enumerate(values):
                #try to set up a combo box if defined 
                # if not create a text entry field
                comboChoices = None
                try: comboChoices = self.variableFixedChoice[key][index]
                except: pass
                
                if comboChoices != None:
                    entry = gtk.combo_box_new_text() 
                    for choice in comboChoices:
                        entry.append_text(choice)
                    entry.set_active(comboChoices.index(value)) 
                    entry.set_size_request(120,30)
                    #cbType.connect("changed", self.on_changePlotType)
                else:
                    entry = gtk.Entry()
                    entry.set_size_request(120,20)
                    #entry.add_events(gtk.gdk.KEY_RELEASE_MASK)
                    entry.set_text(str(value))
                
                hbox.pack_start(entry, fill=False, expand=False)  
                self.entriesCombo[key].append(entry)   
                
            vbox.pack_start(hbox, expand=True, fill=False)
            
        width = 130 + 120 * (countWidth)
        height  = 30 + 30*countHeight
        self.set_size_request(width,height)
            
        self.add(vbox)
        self.show_all()
        
    def on_clickedUpdate(self, widget):
        '''
        update dictionary from data 
        written in entriesCombo and choosen from comboboxes
        '''
        for key, entries in self.entriesCombo.iteritems():
            for index,entry in enumerate(entries):
                
                comboChoices = None
                try: comboChoices = self.variableFixedChoice[key][index]
                except: pass
                
                if comboChoices != None: 
                    value = entry.get_active_text()                    
                    try: type = self.variableTypes[key][index]
                    except: type = 'str'
                    #TODO: Try Except Pass should be fixed
                    try: self.actualValuesDict[key][index] = eval(type)(value)
                    except: print """WARNING: could not convert "{}" to "{}" """.format(text,type)
                    
                else:
                    text = entry.get_text()
                    if text is not '':
                        try: type = self.variableTypes[key][index]
                        except: type = 'float'
                        #TODO: Try Except Pass should be fixed
                        try: self.actualValuesDict[key][index] = eval(type)(text)
                        except: print """WARNING: could not convert "{}" to "{}" """.format(text,type)
                            
        self.parentWindow.update({self.dictVariableName:self.actualValuesDict})
        self.parentWindow.updatePlotWindow()
        
    def on_clickedReset(self, widget):
        '''
        Reset entriesCombo and and comboBoxes to entriesCombo of the beginning
        '''
        for key, entries in self.entriesCombo.iteritems():
            for index,entry in enumerate(entries):
                
                value = self.inititalValuesDict[key][index]
                
                comboChoices = None
                #TODO: Try Except Pass should be fixed
                try: comboChoices = self.variableFixedChoice[key][index]
                except: pass
                
                if comboChoices != None: 
                    valueIndex = comboChoices.index(value)
                    entry.set_active(valueIndex) 
                else:
                    entry.set_text(str(value))
                
        self.parentWindow.update({self.dictVariableName: copy(self.inititalValuesDict)})
        self.parentWindow.updatePlotWindow()

class Visualisation2DPlotWindowGui(gtk.Window):
    def __init__(self):
        # self.set_title(' '.join(['2d Visualisation - ',self.vascularNetwork.name,' : ',self.vessel.name,'- dataNumber ',dataNumber]))
        super(Visualisation2DPlotWindowGui, self).__init__()
        
        self.axisX = 'Time'
        self.scale = gtk.HScale()
        self.scale.set_range(0, 10)
        self.scale.set_increments(1, 1)
        self.scale.set_digits(0)
        self.scale.set_size_request(400, 30)
        self.scale.connect('value-changed', self.on_changedSlider)
        self.scale.set_draw_value(False)        
        
        self.nodeLabel = gtk.Label('Node 0')
        
        # # image legend
        self.cBlegend = gtk.CheckButton("show legend")
        self.cBlegend.connect("toggled", self.on_changedLegend)
        self.cBlegend.set_size_request(120, 30)
                
        # medical units
#         self.medicalUnit = True
#         cBmedical = gtk.CheckButton("use medical units")
#         #cBmedical.connect("toggled", self.on_changedMedicalBool)
#         cBmedical.set_size_request(120,30)
#         cBmedical.set_active(1) 

        # show decriptions
        self.buttonDescription = gtk.CheckButton('show descriptions')
        self.buttonDescription.connect("toggled", self.on_changedDescriptions)
        self.buttonDescription.set_size_request(120, 30)
        # render movie button
        self.buttonRenderMovie = gtk.Button("render movie")
        # self.buttonRenderMovie.connect("clicked", self.on_clickedRenderMovie)
        self.buttonRenderMovie.set_size_request(120, 30)
        self.buttonRenderMovie.set_sensitive(False)

        #### second row of buttons
        # min max plots
        self.buttonMinMaxPoints = gtk.ToggleButton('MinMax points')
        self.buttonMinMaxPoints.connect("toggled", self.on_changedPlotMinMax)
        self.buttonMinMaxPoints.set_size_request(120, 30)
        
        buttonDeltas = gtk.ToggleButton("update Deltas")
        buttonDeltas.connect("toggled", self.on_changeDeltas)
        buttonDeltas.set_size_request(120, 30)        
        
        #checkInflation = gtk.CheckButton('Point of Inflection')
        # checkInflation.connect("toggled", self.on_changedDescription)
        #checkInflation.set_size_request(120, 30)
        
#         cBchangeMinMaxDelta = gtk.Button('change MinMax-deltas')
#         cBchangeMinMaxDelta.connect("clicked", self.on_clickedMinMaxDelta)
#         cBchangeMinMaxDelta.set_size_request(120,30)
#         cBchangeMinMaxDelta.set_sensitive(False)

        buttonLimits = gtk.ToggleButton("update limits")
        buttonLimits.connect("toggled", self.on_changeLimits)
        buttonLimits.set_size_request(120, 30)

        buttonLabels = gtk.ToggleButton("update labels")
        buttonLabels.connect("toggled", self.on_changeLabels)
        buttonLabels.set_size_request(120, 30)
        
        self.buttonLines = gtk.ToggleButton("update lines")
        self.buttonLines.connect("toggled", self.on_changeLines)
        self.buttonLines.set_size_request(120, 30)
        
        self.fig = Figure(figsize=(5, 4), dpi=100)
        # self.plot(0)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        width = 640
        height = 500 #690
        self.canvas.set_size_request(width, height)
        
        toolbar = NavigationToolbar(self.canvas, self)
        
        cbType = gtk.combo_box_new_text()
        cbType.append_text('Plot P,Q')
        cbType.append_text('Plot P,Q with wave split')
        cbType.append_text('Plot P,Q with wave split - mean values')
        cbType.append_text('Plot CFL, wave speed')
        cbType.append_text('Plot Area, Compliance')
        cbType.append_text('Plot netGravity')
        cbType.append_text('Plot P,Q ws - m v - centeroids')
        cbType.set_active(0) 
        cbType.connect("changed", self.on_changePlotType)
        
        self.cbXaxis = gtk.combo_box_new_text()
        self.cbXaxis.connect("changed", self.on_changePlotXaxis)
        self.cbXaxis.append_text('Plot over Time')
        self.cbXaxis.append_text('Plot over Space')
        self.cbXaxis.set_active(0) 
        self.cbXaxis.set_sensitive(False)

        vbox = gtk.VBox(False, 1)
        hboxLegend = gtk.HBox(False, 10) 
        hboxCheckboxes2 = gtk.HBox(False, 10) 
        hbox = gtk.HBox(False, 1)
        hbox2 = gtk.HBox(False, 1)
                
        # Checkbutton series
        hboxLegend.pack_start(self.cBlegend)
        hboxLegend.pack_start(buttonLabels)
        hboxLegend.pack_start(self.buttonDescription)
        
        # hboxLegend.pack_start(cBmedical)
        hboxLegend.pack_start(self.buttonRenderMovie)
        hboxCheckboxes2.pack_start(self.buttonMinMaxPoints)
        hboxCheckboxes2.pack_start(buttonDeltas)
        hboxCheckboxes2.pack_start(buttonLimits)
        hboxCheckboxes2.pack_start(self.buttonLines)  

        # align pictures canvas
        alignIm = gtk.Alignment(0, 1 , 1, 0)
        alignIm.add(self.canvas)
        # align node switcher scale
        hbox.pack_start(self.nodeLabel)
        hbox.pack_start(self.scale)
        alignHbox = gtk.Alignment(0, 1, 1, 0)
        alignHbox.add(hbox)
        # align combobox
        hbox2.pack_start(cbType)
        hbox2.pack_start(self.cbXaxis)
        alignCB = gtk.Alignment(0, 1, 1, 0)  
        alignCB.add(hbox2)
        # align navigation toolbox
        alignNT = gtk.Alignment(0, 1, 1, 0)  
        alignNT.add(toolbar)
        
        # put all together
        vbox.pack_start(hboxLegend)
        vbox.pack_start(hboxCheckboxes2)
        vbox.pack_start(alignCB)
        vbox.pack_start(alignNT)
        vbox.pack_start(alignIm)
        vbox.pack_start(alignHbox)
        
        self.add(vbox)         
        self.show_all()
            
    def on_changePlotXaxis(self, widget):
        pass
        
    def on_changedSlider(self, widget):
        pass
            
    def on_changeLimits(self, widget):
        pass
    
    def on_changedPlotMinMax(self,widget):
        pass
    
    def on_changeDeltas(self,widget):
        pass
    
    def on_changedLegend(self,widget):
        pass
    
    def on_changeLabels(self,widget):
        pass
    
    def on_changeLines(self,widget):
        pass
    
    def on_changedDescriptions(self,widget):
        pass
    
class Visualisation2DPlotWindow(Visualisation2DPlotWindowGui):
    def __init__(self, selectedNetworks, selectedVesselIds, selectedExternalData, selectedCaseNames):
                        
        self.plot = lambda :''
                
        # # variables for the slider
        self.sliderValue    = 0  # grid node / time point
        self.limitsWindow   = None
        self.limits         = {'Time':      [0, 10],
                               'Space':     [0, 10],
                               'gridNodes': [0, 5]}
        self.limitsInit     = {}
        
        self.deltasWindow   = None
        self.deltas         = {}
        self.deltasInit     = {}
        
        self.axis = {}
        self.lines = {}
        self.points = {}
        
        self.legend       = None
        self.labelsWindow = None
        self.labels       = {}
        self.labelsInit   = {}
        self.labelsTypes  = {}
        
        self.descriptions = None
        
        self.lineWindow            = None
        self.lineProperties        = {}
        self.linePropertiesInit    = {}
        self.linePropertiesChoices = {}
        self.linePropertiesTypes   = {}
        
        # # initialize super class
        super(Visualisation2DPlotWindow, self).__init__()
        
        self.selectedNetworks     = selectedNetworks
        self.selectedVesselIds    = selectedVesselIds
        self.selectedExternalData = selectedExternalData
        self.selectedCaseNames    = selectedCaseNames
        
        # # activate space/time changer if only one network is loaded
        if len(self.selectedNetworks) == 1:
            self.cbXaxis.set_sensitive(True)
                
        self.linewidth = 1.5  
        self.fontSizeLabel = 14
        
        self.unitPtext = '$mmHg$'
        self.unitFtext = '$ml/s$'
        
        self.createGraph()
        self.estimatePlotLimits()
        self.plot = self.updateLinesPQ
                
        # # update slider and so --> update the plot
        self.scale.set_range(*self.limits['gridNodes'])
        self.updatePlotWindow()
    
    def update(self, visualisationData):
        '''
        updates the vascularNetwork data using a dictionary in form of 
        visualisationData = {'dictVariableName': value}
        '''
        for key,value in visualisationData.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'WARNING vascularNetwork.update(): wrong key: %s, could not update vascularNetwork' %key 
    
       
    def createGraph(self):
        '''
        create graph with 2 subplots and all necessary lines 
        '''
        self.fig = plt.figure(figsize=(5, 4), dpi=100, edgecolor='k')
        self.fig.subplots_adjust(right=0.86)
        self.fig.subplots_adjust(left=0.17)
        self.fig.subplots_adjust(top=0.95)
        self.fig.subplots_adjust(bottom=0.15)
        self.fig.subplots_adjust(hspace=0.18)
        
        fontLegend = FontProperties()
        fontLegend.set_size(self.fontSizeLabel)
        
        from matplotlib import rc
        from matplotlib import rcParams
        
        rcParams['text.usetex'] = False
        rcParams['text.latex.unicode'] = True
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = self.fontSizeLabel
        rcParams['savefig.dpi'] = 300.
        
        ax1 = plt.subplot(2, 1, 1, frameon=True)
        
        ax12 = ax1.twinx()
        ax12.set_visible(False)
        
        plt.xticks(np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1], 2), ['', '']) 
        ax1.tick_params(axis='x', top='off', bottom='off')
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='y', right='off')
                
        #plt.yticks(np.linspace(ax12.get_xlim()[0], ax12.get_xlim()[1], 2), ['', '']) 
                
        ax2 = plt.subplot(2, 1, 2, frameon=True)
        ax2.spines['top'].set_visible(False)
        ax2.tick_params(axis='x', top='off')
        ax2.spines['right'].set_visible(False)
        ax2.tick_params(axis='y', right='off')
           
        ax2.set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
           
        ax22 = ax2.twinx()
        ax22.set_visible(False)
        
        #plt.yticks(np.linspace(ax22.get_xlim()[0], ax22.get_xlim()[1], 2), ['', '']) 
                           
        self.axis = {'axis1':ax1,'axis1Twin':ax12, 'axis2':ax2,'axis2Twin':ax22}
        
        # # add 3 lines for each network case to each subplot
        colors = ['b', 'r', 'm', 'g', 'c', 'k']
        linestyles = {'axis1': ['-'],'axis1Twin': ['--',':'], 'axis2': ['-'], 'axis2Twin':['--',':']}
                
        # lines and points are refered with the id(int) in the selectedNetworks list == 'caseId'
        # lines  = dict[caseId][axis][linestyle]
        # points = dict[caseId][axis][linestyle] ## here the linestyle define the according line only
        for caseId in xrange(len(self.selectedVesselIds)):
            # create dictionaies accroding to current defined axis
            self.lines[caseId]  = {key: {} for key in self.axis.keys()}
            self.points[caseId] = {key: {} for key in self.axis.keys()}
            # fill in lines according to defined axis, linesytles and color 
            for axis in self.axis.keys():
                for linestyle in linestyles[axis]:
                    
                    currentName = ' '.join([str(caseId),axis,linestyle])
                    if axis == 'axis1':
                        self.labels[currentName]     = [self.selectedCaseNames[caseId]]
                        self.labelsInit[currentName] = [self.selectedCaseNames[caseId]]
                    else:
                        self.labels[currentName]     = ['']
                        self.labelsInit[currentName] = ['']
                    self.labelsTypes[currentName]    = ['str'] 
                    self.deltas[currentName]         = [0.025]
                    self.deltasInit[currentName]     = [0.025]
                    
                    self.lineProperties[currentName]        = [colors[caseId], linestyle, self.linewidth]
                    self.linePropertiesInit[currentName]    = [colors[caseId], linestyle, self.linewidth]
                    self.linePropertiesTypes[currentName]   = ['str', 'str', 'float']
                    self.linePropertiesChoices[currentName] = [colors, ['-','--',':'], None]
                    
                    self.lines[caseId][axis][linestyle] = self.axis[axis].plot( -1, 0,
                                                                                color=colors[caseId],
                                                                                linestyle = linestyle,
                                                                                linewidth = self.linewidth,
                                                                                label = currentName)[0]
                    self.points[caseId][axis][linestyle] = self.axis[axis].plot( -1, 0,
                                                                                 color=colors[caseId],
                                                                                 linestyle = '',
                                                                                 marker='o',
                                                                                 linewidth=self.linewidth)[0]
        
        if self.selectedExternalData != None:
            self.lines['external'] = {}
            self.lines['external']['axis1'] = {'-':self.axis['axis1'].plot( -1, 0,
                                                                            color= 'k',
                                                                            linestyle = '-',
                                                                            linewidth = self.linewidth)[0]}
                                      
            self.lines['external']['axis2'] = {'-':self.axis['axis2'].plot( -1, 0,
                                                                            color= 'k',
                                                                            linestyle = '-',
                                                                            linewidth = self.linewidth)[0]}
                
            self.points['external']= {'axis1':{'-':self.axis['axis1'].plot( -1, 0,
                                                                            color= 'k',
                                                                            linestyle = '',
                                                                            marker='o',
                                                                            linewidth=self.linewidth)[0]},
                                      'axis2':{'-':self.axis['axis2'].plot( -1, 0,
                                                                            color= 'k',
                                                                            linestyle = '',
                                                                            marker='o',
                                                                            linewidth=self.linewidth)[0]}}
            self.deltas['external axis1 -'] = [0.025]
            self.deltas['external axis2 -'] = [0.025]
                                    
    def clearPlotWindow(self, lines = True, points = True, labels = False):
        '''
        set all line data to single point at (-1, 0) out of the view field
        '''
        for caseId,axisDict in self.lines.iteritems():
            for axis,lineDict in axisDict.iteritems():
                if lines:
                    for line in lineDict.itervalues():
                        line.set_data([-1], [0])
                    self.axis['axis1Twin'].set_visible(False)
                    self.axis['axis2Twin'].set_visible(False)
                    self.axis['axis1'].spines['right'].set_visible(False)
                    self.axis['axis2'].spines['right'].set_visible(False)
                if points:
                    try:
                        for point in self.points[caseId][axis].itervalues():
                            point.set_data([-1], [0])
                    except: pass
                if labels:
                    for linestyle in lineDict.iterkeys():
                        currentName = ' '.join([str(caseId),axis,linestyle])
                        self.labels[currentName] = '' 
                                   
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()             
        
    def updatePlotWindow(self):
        '''
        update the plot canvas in the GUI
        '''    
        self.plot()
        
        if self.buttonMinMaxPoints.get_active():
            self.updatePoints()
            
        if self.cBlegend.get_active():
            self.updateLegend()
            try:  self.legend.set_visible(True)
            except: pass
        else:
            try:  self.legend.set_visible(False)
            except: pass
        
        if self.buttonLines.get_active():
            self.updateLineProperties()    
                        
        if self.buttonDescription.get_active():
            
            self.updateDescriptions()
            self.descriptions.set_visible(True)
        else:
            #TODO: Try Except Pass should be fixed
            try: self.descriptions.set_visible(False)
            except: pass
                   
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()        
                
    def updateLegend(self):
        '''
        Stores the labels of the lines in the legend
        '''
        handles = []
        labels  = []
        for axis in self.axis.itervalues():
            ha,la = axis.get_legend_handles_labels()
            for currentAxis,currentLine in zip(ha,la):
                label = self.labels[currentLine][0]
                if label != '':
                    labels.append(label)
                    handles.append(currentAxis)
        #TODO: Try Except Pass should be fixed
        try:  
            self.legend.set_visible(False)
            self.legend.remove()   
        except: pass
        self.legend = self.fig.legend(handles,labels, loc=3, borderaxespad=0., frameon = False, fontsize = self.fontSizeLabel/4.*3.)
        
    def updateDescriptions(self):
        '''
        Stores the descriptions of the cases in the legend
        '''
        handles = []
        labels  = []
        for caseId,vascularNetwork in enumerate(self.selectedNetworks):
            currentColor = self.lines[caseId]['axis1']['-'].get_color()
            feakLine = Rectangle((-1, -1.005), -1.005, -1.005, fc=currentColor, fill=True, edgecolor='none', linewidth=0)
            handles.append(feakLine)
            labels.append(vascularNetwork.description)
        #TODO: Try Except Pass should be fixed
        try:  
            self.descriptions.set_visible(False)
            self.descriptions.remove()   
        except: pass
                
        self.descriptions = self.fig.legend(handles,labels, loc=4, borderaxespad=0., frameon = False, fontsize = self.fontSizeLabel/4.*3.)
        
    def updateLineProperties(self):
        '''
        updates line informations:
        
           color
           linestyle
           linewidth
        '''
        for caseId,axisDict in self.lines.iteritems():
            for axis,lineDict in axisDict.iteritems():
                for linestyle,line in lineDict.iteritems():
                    currentName = ' '.join([str(caseId),axis,linestyle])
                    color,style,width = self.lineProperties[currentName]
                    line.set_color(color)
                    line.set_linestyle(style)
                    line.set_linewidth(width)
                    
    def updateLinesPQ(self):
        '''
        create line plot with P in the first axis and Q the second axis
        '''
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis['axis2'].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
                
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    yData00 = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]] / 133.32
                    ydata10 = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]] * 1e6    
                    xData = vascularNetwork.tsol
                                                           
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)
        
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                    
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                    
                if self.selectedExternalData != None:
                    #TODO: Try Except Pass should be fixed
                    try:
                        xData   = self.selectedExternalData['SolutionData'][vesselId]['Time']+self.limits['external time shift']
                        try:
                            yData00 = self.selectedExternalData['SolutionData'][vesselId]['Pressure']
                            self.lines['external']['axis1']['-'].set_data(xData, yData00)
                        except: pass
                        try:
                            ydata10 = self.selectedExternalData['SolutionData'][vesselId]['Flow']
                            self.lines['external']['axis2']['-'].set_data(xData, ydata10)
                        except: pass
                    except: pass
            
            elif self.axisX == "Space":
                try: 
                    yData00 = vascularNetwork.vessels[vesselId].Psol[gridNode] / 133.32
                    ydata10 = vascularNetwork.vessels[vesselId].Qsol[gridNode] * 1e6    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)      
                                        
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis2'].set_xlim(self.limits['Space'])                                  
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                    
                if self.selectedExternalData != None:
                    self.lines['external']['axis1']['-'].set_data([-1], [0])
                    self.lines['external']['axis2']['-'].set_data([-1], [0])
        
        self.axis['axis1'].set_ylim(self.limits['P'])
        self.axis['axis2'].set_ylim(self.limits['Q'])
                  
    def updateLinesPQsplitLin(self):
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis['axis1Twin'].set_ylabel('Contribution ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis['axis2'].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        self.axis['axis2Twin'].set_ylabel('Contribution ' + self.unitFtext, fontsize=self.fontSizeLabel)
        
        
        
        # 2. update lines for P and Q over time for grid node 0
        
        self.axis['axis1Twin'].set_visible(True)
        self.axis['axis2Twin'].set_visible(True) 
        self.axis['axis1'].spines['right'].set_visible(True)
        self.axis['axis2'].spines['right'].set_visible(True)
        
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                 try:
                    Psol = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[:, [gridNode]]   
                    csol = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = vascularNetwork.tsol
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                 except:
                    self.lines[i]['axis1']['-'].set_data(-1,0)
                    self.lines[i]['axis2']['-'].set_data(-1,0)
            
            elif self.axisX == "Space":
                try:                     
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[gridNode]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[gridNode]   
                    csol = vascularNetwork.vessels[vesselId].csol[gridNode]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlim(self.limits['Space']) 
                                 
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                    
        
        self.axis['axis1'].set_ylim(self.limits['P'])
        self.axis['axis2'].set_ylim(self.limits['Q'])
        
        self.axis['axis1Twin'].set_ylim(self.limits['Pfb'])
        self.axis['axis2Twin'].set_ylim(self.limits['Qfb'])
        
    def updateLinesPQsplitSubMean(self):
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis['axis1Twin'].set_ylabel('')
        self.axis['axis2'].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        self.axis['axis2Twin'].set_ylabel('')
        
        # 2. update lines for P and Q over time for grid node 0
        
        self.axis['axis1Twin'].set_visible(True)
        self.axis['axis2Twin'].set_visible(True) 
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                 try:
                    Psol = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[:, [gridNode]]   
                    csol = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32 - Psol[0]/ 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   - Qsol[0]* 1.e6
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    print sum(Psol_b)/len(Psol_b)
                    
                    xData = vascularNetwork.tsol
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                 except:
                    self.lines[i]['axis1']['-'].set_data(-1,0)
                    self.lines[i]['axis2']['-'].set_data(-1,0)
            
            elif self.axisX == "Space":
                try:                     
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[gridNode]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[gridNode]   
                    csol = vascularNetwork.vessels[vesselId].csol[gridNode]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlim(self.limits['Space']) 
                                 
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                    
        
        self.axis['axis1'].set_ylim(self.limits['PfbLev'])
        self.axis['axis2'].set_ylim(self.limits['QfbLev'])
        
        self.axis['axis1Twin'].set_ylim(self.limits['PfbLev'])
        self.axis['axis2Twin'].set_ylim(self.limits['QfbLev'])
        self.axis['axis1Twin'].get_yaxis().set_ticks([])
        self.axis['axis2Twin'].get_yaxis().set_ticks([])
        
    def updateLinesPQsplitSubMeanNonLinear(self):
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Pressure ' + self.unitPtext, fontsize=self.fontSizeLabel)
        self.axis['axis1Twin'].set_ylabel('')
        self.axis['axis2'].set_ylabel('Flow ' + self.unitFtext, fontsize=self.fontSizeLabel)
        self.axis['axis2Twin'].set_ylabel('')
        
        # 2. update lines for P and Q over time for grid node 0
        
        self.axis['axis1Twin'].set_visible(True)
        self.axis['axis2Twin'].set_visible(True) 
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                 try:
                    Psol = vascularNetwork.vessels[vesselId].Psol[:, [gridNode]]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[:, [gridNode]]   
                    csol = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]  
                    
                    Csol = vascularNetwork.vessels[vesselId].C(Psol)[:, [gridNode]]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.nonLinearWaveSplitting(Psol, Qsol, Asol, csol, Csol, vascularNetwork.vessels[vesselId].rho)
                                  
                    yData00 = Psol / 133.32 - Psol[0]/ 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32 
                    yData10 = Qsol * 1.e6   - Qsol[0]* 1.e6
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    print sum(Psol_b)/len(Psol_b)
                    xData = vascularNetwork.tsol
#                     
#                     # calculate centeroid of Psol_f
#                     A = np.trapz(yData01)
#                     t = xData[1:]
#                     x = np.trapz(yData01*t)/A            
#                     y = np.trapz(yData01**2./2)/A
#                     print "DB cernteroids start"
#                     print "A", A
#                     print x,y
#                     print "DB cernteroids end"
#                     
#                     self.points[i]['axis1Twin']['--'].set_data(x, y)
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                 except:
                    self.lines[i]['axis1']['-'].set_data(-1,0)
                    self.lines[i]['axis2']['-'].set_data(-1,0)
            
            elif self.axisX == "Space":
                try:                     
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]
                    Qsol = vascularNetwork.vessels[vesselId].Qsol[gridNode]   
                    Asol = vascularNetwork.vessels[vesselId].Asol[gridNode]   
                    csol = vascularNetwork.vessels[vesselId].csol[gridNode]  
                    
                    Psol_f, Psol_b, Qsol_f, Qsol_b = mProc.linearWaveSplitting(Psol, Qsol, Asol, csol, vascularNetwork.vessels[vesselId].rho)
                    
                    yData00 = Psol / 133.32
                    yData01 = Psol_f / 133.32
                    yData02 = Psol_b / 133.32
                    yData10 = Qsol * 1.e6   
                    yData11 = Qsol_f * 1.e6  
                    yData12 = Qsol_b * 1.e6  
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[i]['axis1']['-'].set_data(xData,     yData00)
                    self.lines[i]['axis1Twin']['--'].set_data(xData[1:], yData01)
                    self.lines[i]['axis1Twin'][':'].set_data(xData[1:], yData02)
                    
                    self.lines[i]['axis2']['-'].set_data(xData,     yData10)
                    self.lines[i]['axis2Twin']['--'].set_data(xData[1:], yData11)
                    self.lines[i]['axis2Twin'][':'].set_data(xData[1:], yData12)
                                        
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlim(self.limits['Space']) 
                                 
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                    
        
        self.axis['axis1'].set_ylim(self.limits['PfbLev'])
        self.axis['axis2'].set_ylim(self.limits['QfbLev'])
        
        self.axis['axis1Twin'].set_ylim(self.limits['PfbLev'])
        self.axis['axis2Twin'].set_ylim(self.limits['QfbLev'])
        self.axis['axis1Twin'].get_yaxis().set_ticks([])
        self.axis['axis2Twin'].get_yaxis().set_ticks([])
    
    def updateLinesWaveCFL(self):
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Wave Speed $m/s$', fontsize=self.fontSizeLabel)
        self.axis['axis2'].set_ylabel('CFL $-$' + self.unitFtext, fontsize=self.fontSizeLabel)
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    # wave speed
                    yData00 = vascularNetwork.vessels[vesselId].csol[:, [gridNode]]
                    # CFL 
                    dz = vascularNetwork.vessels[vesselId].dz
                    ydata10 = yData00 * vascularNetwork.dt / sum(dz) * len(dz)
                    
                    xData = vascularNetwork.tsol
                    
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)
                    
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
            
            elif self.axisX == "Space":
                try: 
                    # wave speed
                    yData00 = vascularNetwork.vessels[vesselId].csol[gridNode]
                    # CFL 
                    dz = vascularNetwork.vessels[vesselId].dz
                    ydata10 = yData00 * vascularNetwork.dt / sum(dz) * len(dz)
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    print np.shape(xData)
                    
                    print yData00, xData
                    
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)      
                    
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlim(self.limits['Space'])     
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
        
        self.axis['axis1'].set_ylim(self.limits['c'])
        self.axis['axis2'].set_ylim(self.limits['CFL'])
        
    def updateLinesAreaComp(self):
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Area $mm^2$', fontsize=self.fontSizeLabel)
        self.axis['axis2'].set_ylabel('Compliance $mm^2/mmHg$', fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    yData00 = vascularNetwork.vessels[vesselId].Asol[:, gridNode] * 1000 * 1000
                    Psol = vascularNetwork.vessels[vesselId].Psol   
                    ydata10 = vascularNetwork.vessels[vesselId].C(Psol)[:, gridNode] * 1000 * 1000 / 133.32 
                    
                    xData = vascularNetwork.tsol
                    
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)
                    
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
            
            elif self.axisX == "Space":
                try: 
                    yData00 = vascularNetwork.vessels[vesselId].Asol[gridNode] * 1000 * 1000
                    Psol = vascularNetwork.vessels[vesselId].Psol[gridNode]   
                    ydata10 = vascularNetwork.vessels[vesselId].C(Psol) * 1000 * 1000 / 133.32   
                    
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                    
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    self.lines[i]['axis2']['-'].set_data(xData, ydata10)      
                    
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlim(self.limits['Space'])              
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                           
        self.axis['axis1'].set_ylim(self.limits['A'])
        self.axis['axis2'].set_ylim(self.limits['C'])
        
    def updateLinesGravity(self):
        '''
        creates a plot of gravity in the first sub plot
        the second is free 
        '''
        
        gridNode = self.sliderValue
        # 1. set axis label
        self.axis['axis1'].set_ylabel('Gravity $m s^{-2}$', fontsize=self.fontSizeLabel)
        self.axis['axis2'].set_ylabel('', fontsize=self.fontSizeLabel)
        # 2. update lines for P and Q over time for grid node 0
                
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)), self.selectedNetworks, self.selectedVesselIds):        
            
            if self.axisX == 'Time':               
                try:
                    yData00 = vascularNetwork.vessels[vesselId].netGravity[:]
                    if len(yData00)!= 1:
                        print "WARNING 2dVisualisation.updateLinesGravity(): either no motion or no correct solution data"
                    #ydata10 = vascularNetwork.vessels[vesselId].Qsol[:, [gridNode]] * 1e6    
                    xData = vascularNetwork.tsol
                                                             
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    #self.lines[i]['axis2']['-'].set_data(xData, ydata10)
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
        
                    self.axis['axis1'].set_xlim(self.limits['Time'])
                    self.axis['axis2'].set_xlabel('Time $s$', fontsize=self.fontSizeLabel)
                    self.axis['axis2'].set_xlim(self.limits['Time'])
                    
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                                
            elif self.axisX == "Space":
                try:
                    yData00 = vascularNetwork.vessels[vesselId].netGravity[:]
                    #ydata10 = vascularNetwork.vessels[vesselId].Qsol[gridNode] * 1e6 
                    xData = np.linspace(0, vascularNetwork.vessels[vesselId].length, len(yData00)) * 100.
                         
                    self.lines[i]['axis1']['-'].set_data(xData, yData00)
                    #self.lines[i]['axis2']['-'].set_data(xData, ydata10)      
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
                                             
                    self.axis['axis1'].set_xlim(self.limits['Space'])
                    self.axis['axis2'].set_xlabel('Space $cm$', fontsize=self.fontSizeLabel)
                    self.axis['axis2'].set_xlim(self.limits['Space'])             
                    
                except:
                    self.lines[i]['axis1']['-'].set_data([-1], [0])
                    self.lines[i]['axis2']['-'].set_data([-1], [0])
  
#                 if self.selectedExternalData != None:
#                     self.lines['external']['axis1']['-'].set_data([-1], [0])
#                     self.lines['external']['axis2']['-'].set_data([-1], [0])
        
        self.axis['axis1'].set_ylim(self.limits['G'])
        #self.axis['axis2'].set_ylim(self.limits['Q'])
        
        
    def updatePoints(self):
        '''
        plot min and max points of all current lines
        '''
        # 1 find out current visible lines
        for caseId,axisDict in self.lines.iteritems():
            for axis,lineDict in axisDict.iteritems():
                for linestyle,line in lineDict.iteritems():
                    
                    xData = line.get_xdata()
                    if len(xData) != 1:
                        # if toggled add points if "untoggled" remove points
                        # 2 find min and max points for all current lines
                        currentName = ' '.join([str(caseId),axis,linestyle])
                        minMaxPoints = mProc.minMaxFunction(line.get_ydata(), xData, delta = self.deltas[currentName][0])
                        # 3 plot all min and max points for current lines
                        self.points[caseId][axis][linestyle].set_data(minMaxPoints[1], minMaxPoints[0])
        
    def estimatePlotLimits(self):
        '''
        This function evaluates all limits for the plots
        It may take a while to calculate all
        '''     
        self.limits = {'P':         [1e50, -1e50],
                       'Pf':        [1e50, -1e50],
                       'Pb':        [1e50, -1e50],
                       'Pfb':       [1e50, -1e50],
                       'PfbLev':    [1e50, -1e50],
                       'Q':         [1e50, -1e50],
                       'Qf':        [1e50, -1e50],
                       'Qb':        [1e50, -1e50],
                       'Qfb':       [1e50, -1e50],
                       'QfbLev':    [1e50, -1e50],
                       'Time':      [0, -1e50],
                       'Space':     [0, -1e50],
                       'c':         [1e50, -1e50],
                       'CFL':       [ 0, 1.1],
                       'A':         [1e50, -1e50],
                       'C':         [1e50, -1e50],
                       'gridNodes': [0, -1e50],
                       'G':         [-0.25, 0.25]}
        
        for i, vascularNetwork, vesselId in zip(xrange(len(self.selectedVesselIds)),
                                                self.selectedNetworks, self.selectedVesselIds):        
            # pressure
            limit = 'P'
            Psol = vascularNetwork.vessels[vesselId].Psol 
            # sol = Psol / 133.32
            self.limits[limit] = [min([self.limits[limit][0], np.min(Psol)/133.32]),
                                  max([self.limits[limit][1], np.max(Psol)/133.32])]
            # flow
            limit = 'Q'
            Qsol = vascularNetwork.vessels[vesselId].Qsol 
            #sol = Qsol * 1.e6
            self.limits[limit] = [min([self.limits[limit][0], np.min(Qsol)*1.e6]),
                                  max([self.limits[limit][1], np.max(Qsol)*1.e6])]  
            # time
            limit = 'Time'
            sol = vascularNetwork.tsol
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                  max(self.limits[limit][1], np.max(sol))]
            # space
            limit = 'Space'
            sol = vascularNetwork.vessels[vesselId].length * 100.
            self.limits[limit] = [0,
                                  max(self.limits[limit][1], np.max(sol))]
            # wave speed
            limit = 'c'
            csol = vascularNetwork.vessels[vesselId].csol
            sol = csol
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                  max([self.limits[limit][1], np.max(sol)])]
            # CFL
            limit = 'CFL'
            sol = csol * vascularNetwork.dt / vascularNetwork.vessels[vesselId].dz[0] 
            self.limits[limit] = [0, max([self.limits[limit][1], np.max(sol)])]
            # area
            limit = 'A'
            Asol = vascularNetwork.vessels[vesselId].Asol 
            # sol = Asol * 1000 * 1000
            self.limits[limit] = [min([self.limits[limit][0], np.min(Asol)*1.e6]),
                                  max([self.limits[limit][1], np.max(Asol)*1.e6])]
            # compliance
            limit = 'C'
            sol = vascularNetwork.vessels[vesselId].C(Psol) * 1000 * 1000 / 133.32
            self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                  max([self.limits[limit][1], np.max(sol)])]
            # gridNodes
            limit = 'gridNodes'
            sol = vascularNetwork.vessels[vesselId].N - 1
            self.limits[limit] = [0, max([self.limits[limit][1], np.max(sol)])]
            
            # pressure
            limit = 'G'
            #TODO: Try Except Pass should be fixed
            try: 
                netG = vascularNetwork.vessels[vesselId].netGravity 
                self.limits[limit] = [min([self.limits[limit][0], np.min(netG)]),
                                      max([self.limits[limit][1], np.max(netG)])]
            except: pass
                        
            # pressure / flow,  forward backward
            for n in [0,-1]:#xrange(int(vascularNetwork.vessels[vesselId].N)):
                pf,pb,qf,qb =  mProc.linearWaveSplitting(Psol[:,n],Qsol[:,n],Asol[:,n],csol[:,n],vascularNetwork.vessels[vesselId].rho)
                 
                limit = 'Pf'
                sol = pf/133.32
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                      max([self.limits[limit][1], np.max(sol)])]
                 
                limit = 'Pb'
                sol = pb/133.32
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                      max([self.limits[limit][1], np.max(sol)])]
                 
                limit = 'Qf'
                sol = qf*1.e6
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                      max([self.limits[limit][1], np.max(sol)])]
                 
                limit = 'Qb'
                sol = qb*1.e6
                self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]),
                                      max([self.limits[limit][1], np.max(sol)])]
                 
                limit = 'Pfb'
                self.limits[limit] = [min([self.limits[limit][0],
                                           np.min(pf/133.32),
                                           np.min(pb/133.32)]),
                                      max([self.limits[limit][1],
                                           np.max(pf/133.32),
                                           np.max(pb/133.32)])]
                 
                limit = 'Qfb'
                self.limits[limit] = [min([self.limits[limit][0],
                                           np.min(qf*1.e6),
                                           np.min(qb*1.e6)]),
                                      max([self.limits[limit][1],
                                           np.max(qf*1.e6),
                                           np.max(qb*1.e6)])]
             
                limit = 'PfbLev'
                self.limits[limit] = [min([self.limits[limit][0],
                                           np.min((Psol-Psol[0]) / 133.32),
                                           np.min(pf/133.32),
                                           np.min(pb/133.32)]),
                                      max([self.limits[limit][1],
                                           np.max((Psol-Psol[0]) / 133.32),
                                           np.max(pf/133.32),
                                           np.max(pb/133.32)])]
                 
                limit = 'QfbLev'
                self.limits[limit] = [min([self.limits[limit][0],
                                           np.min((Qsol-Qsol[0]) * 1.e6),
                                           np.min(qf*1.e6),
                                           np.min(qb*1.e6)]),
                                      max([self.limits[limit][1],
                                           np.max((Qsol-Qsol[0]) * 1.e6),
                                           np.max(qf*1.e6),
                                           np.max(qb*1.e6)])]
                         
            if self.selectedExternalData != None:
                #TODO: Try Except Pass should be fixed
                try:
                    # pressure
                    limit = 'P'
                    sol = self.selectedExternalData['SolutionData'][vesselId]['Pressure']
                    self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                except: pass
                try:
                    # flow
                    limit = 'Q'
                    sol = self.selectedExternalData['SolutionData'][vesselId]['Flow']
                    self.limits[limit] = [min([self.limits[limit][0], np.min(sol)]), max([self.limits[limit][1], np.max(sol)])]
                except: pass
                self.limits['external time shift'] = [0]
                
        self.limitsInit = copy(self.limits)
        
        
    # callback function
    def on_changePlotXaxis(self, widget):   
        '''
        callback function for changed plot axis menu time/space
        
        where space plot is only possible if one network is loaded
        '''    
        
        cbIndex = widget.get_active()
        if cbIndex == 0:
            self.axisX = 'Time'
            self.nodeLabel.set_text('Node ' + str(self.sliderValue))
            self.buttonRenderMovie.set_sensitive(False)
            self.scale.set_range(*self.limits['gridNodes']) 
        
        elif cbIndex == 1:
            self.axisX = 'Space'
            NodeValue = str(self.selectedNetworks[0].tsol[self.sliderValue])
            self.nodeLabel.set_text(' '.join(['Time ', str(NodeValue)[:6]]))
            self.buttonRenderMovie.set_sensitive(True)
            self.scale.set_range(0, len(self.selectedNetworks[0].tsol) - 1)   
            
        self.updatePlotWindow()
        
    def on_changedSlider(self, widget):
        '''
        callback function for changed node/time slider
        '''
        self.sliderValue = int(widget.get_value())
        if self.axisX == 'Time':
            self.nodeLabel.set_text('Node ' + str(self.sliderValue))
        elif self.axisX == 'Space':
            NodeValue = str(self.selectedNetworks[0].tsol[self.sliderValue])
            self.nodeLabel.set_text(' '.join(['Time ', str(NodeValue)[:6]]))
        
        self.updatePlotWindow()
        
    def on_changePlotType(self, widget):
        '''
        call function of the combobox to choose type of plot
        '''
        self.clearPlotWindow()
        cbIndex = widget.get_active()
        if   cbIndex == 0:
            self.plot = self.updateLinesPQ
        elif cbIndex == 1:
            self.plot = self.updateLinesPQsplitLin
        elif cbIndex == 2:
            self.plot = self.updateLinesPQsplitSubMean
        elif cbIndex == 3:
            self.plot = self.updateLinesWaveCFL
        elif cbIndex == 4:
            self.plot = self.updateLinesAreaComp
        elif cbIndex == 5:
            self.plot = self.updateLinesGravity
        elif cbIndex == 6:
            self.plot = self.updateLinesPQsplitSubMeanNonLinear
        self.plot()
        
    def on_changedPlotMinMax(self, widget):
        if widget.get_active():
            self.updatePlotWindow()
        else:
            self.clearPlotWindow(lines = False)
            #TODO: Try Except Pass should be fixed
            try: self.deltasWindow.destroy()
            except: pass
            self.deltasWindow = None
                    
    def on_changeDeltas(self, widget):
        pass
        '''
        create widnow to change the deltas for the min max function
        '''        
        if widget.get_active() and self.buttonMinMaxPoints.get_active():
            self.deltasWindow = Visualisation2DPlotWindowAdjustValues(self,'deltas',self.deltas,self.deltasInit)    
        else:
            #TODO: Try Except Pass should be fixed
            try: self.deltasWindow.destroy()
            except: pass
            self.deltasWindow = None
            widget.set_active(False)
         
        self.updatePlotWindow()
        
    def on_changeLimits(self, widget):
        '''
        create window to change plot limits
        '''
        if widget.get_active():
            self.limitsWindow = Visualisation2DPlotWindowAdjustValues(self,'limits',self.limits,self.limitsInit)    
        else:
            self.limitsWindow.destroy()
            self.limitsWindow = None
        
        self.updatePlotWindow()
        
    def on_changedLegend(self,widget):
        '''
        activate legend and deactivate it
        '''
        self.updatePlotWindow()
            
    def on_changeLabels(self, widget):
        '''
        create window to change plot limits
        '''
        if widget.get_active():
            self.labelsWindow = Visualisation2DPlotWindowAdjustValues(self,'labels',self.labels,self.labelsInit,self.labelsTypes)    
        else:
            self.labelsWindow.destroy()
            self.labelsWindow = None
        
        self.updatePlotWindow()
        
    def on_changeLines(self, widget):
        '''
        create window to change plot limits
        '''
        if widget.get_active():
            self.lineWindow = Visualisation2DPlotWindowAdjustValues(self,
                                                                    'lineProperties',
                                                                    self.lineProperties,   
                                                                    self.linePropertiesInit,
                                                                    self.linePropertiesTypes,
                                                                    self.linePropertiesChoices)
        else:
            self.lineWindow.destroy()
            self.lineWindow = None
        
        self.updatePlotWindow()
        
    def on_changedDescriptions(self,widget):
        '''
        activate descriptions
        '''
        self.updatePlotWindow()
        
class Visualisation2DMainCase(object):
    def __init__(self, number):
        
        self.networkInfo = {"choose simulation case" : [[], '-'] }
        self.networkCases = ["choose simulation case"]
        self.currentNetwork = self.networkCases[0]
        self.currentVesselId = None
                
        # description of the dataSet
        self.networkDescription = gtk.Label(' here is the description that can be rather long')
        self.networkDescription.set_size_request(400, 35)
        # # Combo boxes
        # vessel chooser
        self.comboBoxVessels = gtk.combo_box_new_text()
        self.comboBoxVessels.append_text("choose vessel")
        self.comboBoxVessels.connect("changed", self.on_changedCBvessels) 
        self.comboBoxVessels.set_active(0)  
        # Network SolutionData chooser
        self.comboBoxNetworks = gtk.combo_box_new_text()
        self.comboBoxNetworks.append_text("choose simulation case")
        self.comboBoxNetworks.connect("changed", self.on_changedNetworkComboBox) 
        self.comboBoxNetworks.set_size_request(250, 35)
        self.comboBoxNetworks.set_active(0)    
        textComboBoxNetworks = gtk.Label(' '.join(['Case', str(number)]))
        textComboBoxNetworks.set_size_request(100, 35)             
        
        # Case
        self.vBoxCase = gtk.VBox(False, 10)
        hBoxButtons = gtk.HBox(False, 10)
        hBoxDecription = gtk.HBox(False, 10)

        spacingText25Box = gtk.Label('')
        spacingText25Box.set_size_request(5, 35)
        spacingText25BoxDesc = gtk.Label('')
        spacingText25BoxDesc.set_size_request(25, 35)

        hBoxButtons.pack_start(spacingText25Box, fill=False, expand=False)
        hBoxButtons.pack_start(textComboBoxNetworks, fill=False, expand=False)
        hBoxButtons.pack_start(self.comboBoxNetworks, fill=False, expand=False)
        hBoxButtons.pack_start(self.comboBoxVessels, fill=False , expand=False)
         
        hBoxDecription.pack_start(spacingText25BoxDesc, fill=False, expand=False)
        hBoxDecription.pack_start(self.networkDescription, fill=False, expand=False)

        separator = gtk.HSeparator()
        
        self.vBoxCase.pack_start(hBoxButtons, expand=True, fill=False)
        self.vBoxCase.pack_start(hBoxDecription, expand=True, fill=False)
        self.vBoxCase.pack_start(separator, expand=True, fill=False)
        
    def updateVesselComboBox(self):
        '''
        Update the comboBox entriesCombo with vessel ids
        '''
        self.comboBoxVessels.get_model().clear()
        self.comboBoxVessels.append_text("choose vessel")
        for vesselName in self.networkInfo[self.currentNetwork][0]:
            self.comboBoxVessels.append_text(vesselName)
        self.comboBoxVessels.set_active(0)
                        
    def updateNetworkComboBox(self, networkCases, networkInfo):
        '''
        update the comboBox entriesCombo if availabel network names 
        have changed
        '''
        self.networkCases = networkCases
        self.networkInfo = networkInfo
        
        self.comboBoxNetworks.get_model().clear()
        for casName in networkCases:
            self.comboBoxNetworks.append_text(casName)     
        
        if len(networkCases) > 1:
            self.comboBoxNetworks.set_active(1)
            self.currentNetwork = self.networkCases[1]
        else:
            self.comboBoxNetworks.set_active(0)
            self.currentNetwork = self.networkCases[0]       
        
    def on_changedCBvessels(self, widget):
        '''
        call function of the combobox to choose vessel of the network
        '''
        cbIndex = widget.get_active()
        if cbIndex <= 0 : self.currentVesselId = None
        else: self.currentVesselId = self.networkInfo[self.currentNetwork][0][cbIndex - 1                                                                     ]
        
    def on_changedNetworkComboBox(self, widget):
        '''
        call function of the combobox to choose solution data set number
        '''
        cbIndex = widget.get_active()
        #TODO: Try Except Pass should be fixed
        try: self.currentNetwork = self.networkCases[cbIndex]
        except: pass
        self.updateVesselComboBox()
        
        self.networkDescription.set_text(self.networkInfo[self.currentNetwork][1])

class Visualisation2DMainGUI(gtk.Window):
    '''
    Class defining the GUI of the visualisation main window
    '''
    def __init__(self):
        '''
        Initialize
        '''
        self.networkInfo = {"choose simulation case" : [[], '-'] }
        self.networkCases = ["choose simulation case"]
        
        super(Visualisation2DMainGUI, self).__init__()
        # variables for windows                
        self.set_size_request(650, 290)
        self.set_position(gtk.WIN_POS_CENTER)
        self.connect("destroy", gtk.main_quit)
        self.mainTitle = "2D Visualisation - "
        self.set_title(self.mainTitle)
        self.set_resizable(False)
                       
        # open Button
        self.buttonOpenSolutionData = gtk.Button("Open solutionData")
        self.buttonOpenSolutionData.connect("clicked", self.on_clickedLoad)
        self.buttonOpenSolutionData.set_size_request(200, 35)
        
        # open plot window
        self.buttonOpenPlotWindow = gtk.Button("Show plots")
        self.buttonOpenPlotWindow.connect("clicked", self.on_clickedPlots)
        self.buttonOpenPlotWindow.set_size_request(120, 35)
        
         # add new Case
        self.cases = []
        self.buttonAddCase = gtk.Button("Add case")
        self.buttonAddCase.connect("clicked", self.on_clickedAddCase)
        self.buttonAddCase.set_size_request(120, 35)
        
        # ExternalData
        self.extDataDescription = gtk.Label("-")
        self.extDataDescription.set_size_request(400, 35)
        # open Button
        self.buttonOpenExternalData = gtk.Button('Open external data for comparison')
        self.buttonOpenExternalData.connect("clicked", self.on_clickedLoadExternal)
        self.buttonOpenExternalData.set_size_request(250, 35)
        # enable check box
        self.buttonEnableExternalData = gtk.CheckButton("plot external data")
        self.buttonEnableExternalData.set_size_request(150, 35)
        self.buttonEnableExternalData.set_active(0)
        
        # alignment of the boxes
        self.vBox = gtk.VBox(False, 10)
        
        # Load And Plot buttons
        hBox1 = gtk.HBox(False, 10)
        
        spacingText25Box1 = gtk.Label('')
        spacingText25Box1.set_size_request(25, 35)
        
        hBox1.pack_start(spacingText25Box1, fill=False, expand=False)
        hBox1.pack_start(self.buttonOpenSolutionData, fill=False, expand=False)
        hBox1.pack_start(self.buttonOpenPlotWindow, fill=False, expand=False)
        separator1 = gtk.HSeparator()
        
        self.vBox.pack_start(hBox1, expand=False, fill=False)
        self.vBox.pack_start(separator1, expand=False, fill=False)
        
        # # add first case button
        newCase = Visualisation2DMainCase(len(self.cases) + 1)
        self.cases.append(newCase)
        self.vBox.pack_start(newCase.vBoxCase, expand=False, fill=False)
                
        # add more button
        hBoxAdd = gtk.HBox(False, 10)
        spacingText25Add = gtk.Label('')
        spacingText25Add.set_size_request(25, 35)
        hBoxAdd.pack_start(spacingText25Add, fill=False, expand=False)
        hBoxAdd.pack_start(self.buttonAddCase, fill=False, expand=False)
        separator2 = gtk.HSeparator()
        
        self.vBox.pack_start(hBoxAdd, expand=False, fill=False)
        self.vBox.pack_start(separator2, expand=False, fill=False)
                                
        # external Data Set
        
        spacingText25Box4text = gtk.Label('')
        spacingText25Box4text.set_size_request(15, 25)
        spacingText25Box4Desc = gtk.Label('')
        spacingText25Box4Desc.set_size_request(25, 35)
        
        hBoxExternalData = gtk.HBox(False, 10)
        hBoxExternalDecription = gtk.HBox(False, 10)
        hBoxExternalData.pack_start(spacingText25Box4text, fill=False, expand=False)
        hBoxExternalData.pack_start(self.buttonOpenExternalData, fill=False, expand=False)
        hBoxExternalData.pack_start(self.buttonEnableExternalData, fill=False, expand=False)        
        hBoxExternalDecription.pack_start(spacingText25Box4Desc, fill=False, expand=False)
        hBoxExternalDecription.pack_start(self.extDataDescription, fill=False, expand=False)
        
        # external Data Set
        self.vBox.pack_start(hBoxExternalData, expand=False, fill=False)
        self.vBox.pack_start(hBoxExternalDecription, expand=False, fill=False)       
                     
        self.add(self.vBox)
        self.show_all()
        
    def on_clickedAddCase(self, widget):
        '''
        add new simulation case to compare plots
        bounded to 5 cases now ..
        '''
        if len(self.cases) < 6:
            newCase = Visualisation2DMainCase(len(self.cases) + 1)
            newCase.updateNetworkComboBox(self.networkCases, self.networkInfo)
            self.cases.append(newCase)
            self.vBox.pack_start(newCase.vBoxCase, expand=False, fill=False)
            self.vBox.reorder_child(newCase.vBoxCase, len(self.cases) + 1)
            
            width, height = self.get_size()
            self.set_size_request(width, height + 102)
            
            self.show_all()
                
    def on_clickedLoad(self, widget):
        pass
    
    def on_clickedPlots(self, widget):
        pass
        
    def on_clickedLoadExternal(self):
        pass
    
class Visualisation2DMain(Visualisation2DMainGUI):
    '''
    Class for the Main GUI window
    '''
    def __init__(self, networkName=None, dataNumber=None, connect=None):
        super(Visualisation2DMain, self).__init__()

        self.externalData = None
        
        self.networks     = { "choose simulation case"   : None}
        self.networkCases = [ "choose simulation case"]
        self.networkInfo  = { "choose simulation case"   : [[], '-'] }
        
        self.loadVascularNetwork(networkName, dataNumber)

    def on_clickedPlots(self, widget):
        '''
        Open new plot window to plot information of all selected vessels of the selected cases  
        '''
        # # check out selected networks and ids
        selectedNetworks = []
        selectedVesselIds = []
        selectedCases = []
        selectedCaseNames = []
        loadVesselIDdict ={} 
        for case in self.cases:
            currentNetwork = case.currentNetwork
            currentVesselId = case.currentVesselId
            if currentNetwork != "choose simulation case":
                if currentVesselId:
                    currentVesselId = int(currentVesselId.split('-')[0])                    
                    # check if it is not already in cases
                    add = True
                    addKeyToLoadDict = True
                    if len(selectedCases) > 0:
                        for case in selectedCases:
                            if currentNetwork == case[0]:
                                addKeyToLoadDict = False
                                if currentVesselId == case[1]:
                                    add = False
                                    
                    if addKeyToLoadDict == True:
                        loadVesselIDdict[currentNetwork] = []
                    
                    if add == True:
                        loadVesselIDdict[currentNetwork].append(currentVesselId)
                        selectedNetworks.append(self.networks[currentNetwork])
                        selectedVesselIds.append(currentVesselId)
                        selectedCases.append([currentNetwork,currentVesselId])
                        selectedCaseNames.append(' '.join([' '.join(currentNetwork.split('_')),'vessel',str(currentVesselId)]))
                                        
        # # check selected external data
        selectedExternalData = None
        if self.buttonEnableExternalData.get_active() == True:
            selectedExternalData = self.externalData
            
        # # open plot window        
        if selectedNetworks != []:
            for networkID in loadVesselIDdict.iterkeys():
                vesselIds = loadVesselIDdict[networkID]
                self.networks[networkID].loadSolutionDataRange(vesselIds, values=["All"])
            Visualisation2DPlotWindow(selectedNetworks, selectedVesselIds, selectedExternalData, selectedCaseNames)

    def loadVascularNetwork(self, networkName, dataNumber, networkXmlFile = None, pathSolutionDataFilename = None):
        '''
        This function actually loads the vascular networks with the given names and ids
            networkName
        The network is added to the existing ones if is not existing
        '''
        # load vascular network
        vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber = dataNumber, networkXmlFile = networkXmlFile, pathSolutionDataFilename = pathSolutionDataFilename)
        vascularNetwork.linkSolutionData()
        
        # # save it and refresh GUi setup
        networkSolutionName = '_'.join([networkName, dataNumber])  
        # # add data name and corresponding network
        self.networks[networkSolutionName] = vascularNetwork
        # # names of the solution data sets availabel sorted keys of the stuff above
        self.networkCases.append(networkSolutionName)
        # # get vessel names
        vesselNames = []
        for vessel in vascularNetwork.vessels.itervalues():
            vesselNames.append('{:3}-{}'.format(vessel.Id, vessel.name))
        # # reference network showing network case name and [corresponding number of vessels, vesselNames, description]
        self.networkInfo[networkSolutionName] = [ vesselNames, vascularNetwork.description]
        
        # # update comboBoxes
        for case in self.cases:
            case.updateNetworkComboBox(self.networkCases, self.networkInfo)

    def on_clickedLoad(self, widget):
        '''
        Call function of the Open Solution Data Button
        '''
        fileFilter = gtk.FileFilter()
        fileFilter.set_name("SolutionData")
        fileFilter.add_pattern("*.hdf5")
        
        filenames = self.LoadDialog(fileFilter)
        
        for filename in filenames:
            if '_SolutionData_' in filename: # normal network cases
                networkName = filename.split('/')[-1].split('_SolutionData_')[0]
                dataNumber = filename.split('/')[-1].split('_SolutionData_')[1].split('.')[0]
                if '_'.join([networkName, dataNumber]) not in self.networkCases:
                    self.loadVascularNetwork(networkName, dataNumber)
            elif '_evaluation_' in filename:# polynomial chaos cases
                print filename
                networkFileName = ''.join([filename.split('/')[-1].split('.')[0],'.xml'])
                networkName = filename.split('/')[-1].split('_evaluation_')[0]
                evaluationNumber = filename.split('/')[-1].split('_evaluation_')[1].split('.')[0]
                networkXmlFile = ''.join(['/'.join(filename.split('/')[0:-2]),'/evaluationNetworkFiles/',networkFileName])
                
                print networkFileName
                print networkName
                print evaluationNumber
                print networkXmlFile
                if '_'.join([networkName, evaluationNumber]) not in self.networkCases:
                    self.loadVascularNetwork(networkName, evaluationNumber,networkXmlFile,filename)
                
                          
    def LoadDialog(self, fileFilter):
        '''
        Dialog window function: 
        displaying the file browser where solution data files can be selected
        to be loaded into the program.
        
        return: [filename,..] (abs-path to the files)
        '''
        dialog = gtk.FileChooserDialog("Open Solution Data File..",
                                     None,
                                     gtk.FILE_CHOOSER_ACTION_OPEN,
                                     (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                      gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        directory = mFPH.getDirectory('workingDirectory', '', '', mode= 'read')
        dialog.set_current_folder(directory)
        
        dialog.set_select_multiple(True)
        
        dialog.add_filter(fileFilter)
        
        filenames = []
        
        response = dialog.run()
        
        if response == gtk.RESPONSE_OK:
            filenames = dialog.get_filenames()
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()    
        return filenames
    

    def on_clickedLoadExternal(self, widget):
        fileFilter = gtk.FileFilter()
        fileFilter.add_pattern("*.v1dfExD")
        filenames = self.LoadDialog(fileFilter)

        #TODO: Try Except Pass should be fixed
        try:
            fileName = filenames[0]
            self.externalData = mFPH.loadExternalDataSet(fileName)
            self.extDataDescription.set_text(self.externalData['Description'])
        except: 
            pass


if __name__ == '__main__':
               
    optionsDict = moduleStartUp.parseOptions(['f', 'n', 'c'], visualisationOnly=True)
    
    networkName = optionsDict['networkName']
    dataNumber = optionsDict['dataNumber']
    connect = optionsDict['connect']

    Visualisation2DMain(networkName=networkName, dataNumber=dataNumber, connect=connect)
    gtk.main()
