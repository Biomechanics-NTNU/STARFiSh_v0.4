#!/usr/bin/env python
import gtk
import gobject

import matplotlib.pyplot as plt   

from matplotlib.figure import Figure
# uncomment to select /GTK/GTKAgg/GTKCairo
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')

#sys.path.append(cur+'/../UtilityLib')
import UtilityLib.moduleXML as mXML 
import UtilityLib.processing as mProc
import UtilityLib.moduleStartUp as moduleStartUp
#from UtilityLib.modulePickle import loadSolutionDataFile
#from UtilityLib.modulePickle import loadExternalDataSet
#deprecated

#sys.path.append(cur+'/../NetworkLib')
from NetworkLib.classVascularNetwork import VascularNetwork 

from matplotlib.font_manager import FontProperties

import itertools
import subprocess

from optparse import OptionParser
import cPickle
import numpy as np

from pprint import pprint as pp

from numpy import arange, sin, pi

class Visualisation2D(gtk.Window):
    '''
    2D visualisation for vascular1Dflow simulation results
    '''
    def __init__(self, networkName = None, dataSetNumber = None, connect = False):
        '''
        Initialize
        If a path to solution data file is given, the one is
        loaded automatically. 
        If non is given, the program starts without initial loaded data        
        '''
        super(Visualisation2D, self).__init__()
        
        self.set_size_request(500,350)
        self.set_position(gtk.WIN_POS_CENTER)
        self.connect("destroy", gtk.main_quit)
        self.mainTitle = "2D Visualisation - "
        self.set_title(self.mainTitle)
        self.set_resizable(False)
        
        self.solutionDataSets = None
        
        self.currentNetworkName = "None"
        self.currentVesselID = None
        self.currentVesselIDComparison = None
        self.oldNumberOfVessels = 0
        self.numberOfVessels = 0
        self.vascularNetwork = None
        self.dataNumbers = []
        self.currentDataNumber = None
        self.currentDataNumberComparison = None
        self.oldDataNumber = 1
        self.dataSetDescriptions = {}
        self.externalData = None 
        self.externalData =  {'PressureTime':None,'Pressure':None,'PressureUnit':'',
                             'FlowTime':None,'Flow':None,'FlowUnit':'',
                             'Description': ''}
        
        # open Button
        self.ButtonOpenSolutionData = gtk.Button("Open SolutionData")
        self.ButtonOpenSolutionData.connect("clicked", self.on_clickedLoad)
        self.ButtonOpenSolutionData.set_size_request(200,35)
        # open plot window
        self.ButtonOpenPlotWindow = gtk.Button("Show Plots")
        self.ButtonOpenPlotWindow.connect("clicked", self.on_clickedPlots)
        self.ButtonOpenPlotWindow.set_size_request(120,35)
        # vessel chooser
        self.CBvessels = gtk.combo_box_new_text()
        self.CBvessels.append_text("choose vessel")
        self.CBvessels.connect("changed", self.on_changedCBvessels) 
        self.CBvessels.set_active(0)  
        # SolutionData chooser
        self.CBsolData = gtk.combo_box_new_text()
        self.CBsolData.append_text(" - ")
        self.CBsolData.connect("changed", self.on_changedCBsolData) 
        self.CBsolData.set_size_request(90,35)
        self.CBsolData.set_active(0)    
        textCBsolData = gtk.Label('Case 1')
        textCBsolData.set_size_request(100,35)             
        # description of the dataSet
        self.dataSetDescription = gtk.Label('')
        self.dataSetDescription.set_size_request(400,35)
        
        # compare to text
        compareText = gtk.Label('Choose second simulation case for comparison:')
        compareTextAlignment = gtk.Alignment(0, 1, 0, 1)  
        compareTextAlignment.add(compareText)
        # buttons for comparison
        # vessel chooser
        self.CBvesselsComparison = gtk.combo_box_new_text()
        self.CBvesselsComparison.append_text("choose vessel")
        self.CBvesselsComparison.connect("changed", self.on_changedCBvesselsComparison) 
        self.CBvesselsComparison.set_active(0)
        # SolutionData chooser
        self.CBsolDataComparison = gtk.combo_box_new_text()
        self.CBsolDataComparison.append_text(" - ")
        self.CBsolDataComparison.set_size_request(90,35)
        self.CBsolDataComparison.connect("changed", self.on_changedCBsolDataComparison) 
        self.CBsolDataComparison.set_active(0)    
        textCBsolDataComparison = gtk.Label('Case 2')
        textCBsolDataComparison.set_size_request(100,35)
        # description of the dataSet
        self.dataSetDescriptionComparison = gtk.Label('')
        self.dataSetDescriptionComparison.set_size_request(400,35)
        
        # ExternalData
        self.extDataDescription = gtk.Label(self.externalData['Description'])
        self.extDataDescription.set_size_request(400,35)
        # open Button
        self.ButtonOpenExternalData = gtk.Button('Open external data for comparison')
        self.ButtonOpenExternalData.connect("clicked", self.on_clickedLoadExternal)
        self.ButtonOpenExternalData.set_size_request(250,35)
        
        spacingText200Box1 = gtk.Label('')
        spacingText200Box1.set_size_request(200,35)
        spacingText25Box1 = gtk.Label('')
        spacingText25Box1.set_size_request(25,35)
        spacingText25Box2 = gtk.Label('')
        spacingText25Box2.set_size_request(25,35)
        spacingText25Box2Desc = gtk.Label('')
        spacingText25Box2Desc.set_size_request(25,35)
        spacingText25Box3text = gtk.Label('')
        spacingText25Box3text.set_size_request(25,35)
        spacingText25Box3 = gtk.Label('')
        spacingText25Box3.set_size_request(25,35)
        spacingText25Box3Desc = gtk.Label('')
        spacingText25Box3Desc.set_size_request(25,35)
        spacingText25Box4text = gtk.Label('')
        spacingText25Box4text.set_size_request(25,35)
        spacingText25Box4Desc = gtk.Label('')
        spacingText25Box4Desc.set_size_request(25,35)
        
        separator1 = gtk.HSeparator()
        separator2 = gtk.HSeparator()
        separator3 = gtk.HSeparator()
        
        #alignment of the boxes
        vBox = gtk.VBox(False,10)
        hBox1 = gtk.HBox(False,10)
        hBox2 = gtk.HBox(False,10)
        hBox2Decription = gtk.HBox(False,10)
        hBox3Text = gtk.HBox(False,10)
        hBox3 = gtk.HBox(False,10)
        hBox3Decription = gtk.HBox(False,10)
        hBox4Button = gtk.HBox(False,10)
        hBox4Decription = gtk.HBox(False,10)
        
        hBox1.pack_start(spacingText25Box1,fill = False, expand = False)
        hBox1.pack_start(self.ButtonOpenSolutionData,fill = False, expand = False)
        #hBox1.pack_start(spacingText200Box1,fill = False, expand = False)
        hBox1.pack_start(self.ButtonOpenPlotWindow,fill = False, expand = False)
        
        #hBox2.pack_start(spacingText50,fill = False)
        hBox2.pack_start(spacingText25Box2,fill = False, expand = False)
        hBox2.pack_start(textCBsolData,fill = False, expand = False)
        hBox2.pack_start(self.CBsolData,fill = False, expand = False)
        hBox2.pack_start(self.CBvessels, fill = False , expand = False)
        
        hBox2Decription.pack_start(spacingText25Box2Desc,fill = False, expand = False)
        hBox2Decription.pack_start(self.dataSetDescription,fill = False, expand = False)
        
        hBox3Text.pack_start(spacingText25Box3text, fill = False, expand = False)
        hBox3Text.pack_start(compareTextAlignment, fill = False, expand = False)
        
        hBox3.pack_start(spacingText25Box3, fill = False, expand = False)
        hBox3.pack_start(textCBsolDataComparison,fill = False, expand = False)
        hBox3.pack_start(self.CBsolDataComparison, fill = False, expand = False)
        #hBox3.pack_start(self.CBvesselsComparison, fill = False, expand = False)
        
        hBox3Decription.pack_start(spacingText25Box3Desc,fill = False, expand = False)
        hBox3Decription.pack_start(self.dataSetDescriptionComparison,fill = False, expand = False)
        
        hBox4Button.pack_start(spacingText25Box4text, fill = False, expand = False)
        hBox4Button.pack_start(self.ButtonOpenExternalData, fill = False, expand = False)
        hBox4Decription.pack_start(spacingText25Box4Desc,fill = False, expand = False)
        hBox4Decription.pack_start(self.extDataDescription,fill = False, expand = False)
        
        vBox.pack_start(hBox1, expand = True, fill = False)
        vBox.pack_start(separator1, expand = True, fill = False)
        vBox.pack_start(hBox2, expand = True, fill = False)
        vBox.pack_start(hBox2Decription, expand = True, fill = False)
        vBox.pack_start(separator2, expand = True, fill = False)
        vBox.pack_start(hBox3Text, expand = True, fill = False)
        vBox.pack_start(hBox3, expand = True, fill = False)
        vBox.pack_start(hBox3Decription, expand = True, fill = False)
        vBox.pack_start(separator3, expand = True, fill = False)
        vBox.pack_start(hBox4Button, expand = True, fill = False)
        vBox.pack_start(hBox4Decription, expand = True, fill = False)
        
        self.connectTo3DViz = connect
        if self.connectTo3DViz == True:
            gobject.timeout_add(500,self.invokedFrom3dViz)
             
                
        self.add(vBox)
        self.show_all()
        
        if networkName != None:
            self.initializeSolutionData(networkName, dataSetNumber)
            title = ' '.join([self.mainTitle,self.currentNetworkName])
            self.set_title(title)
            
            vID = self.vascularNetwork.vessels.keys()
            for i in range(1,self.numberOfVessels+1):
                self.CBvessels.append_text(' '.join([str(vID[i-1]),self.vascularNetwork.vessels[vID[i-1]].name]))
                self.CBvesselsComparison.append_text(' '.join([str(vID[i-1]),self.vascularNetwork.vessels[vID[i-1]].name]))
                
            self.oldNumberOfVessels = self.numberOfVessels
            
            self.CBsolData.remove_text(0)
            
            for dN in self.dataNumbers:
                self.CBsolData.append_text(dN)
                self.CBsolDataComparison.append_text(dN)
            self.CBsolData.set_active(0) 
            self.oldDataNumber = len(self.dataNumbers)
    
    def checkList(self, iterator):
        '''
        function that checks if a list contains only identical items
        '''
        try:
            iterator = iter(iterator)
            first = next(iterator)
            return all(first == rest for rest in iterator)
        except StopIteration:
            return True
        
    def on_clickedLoad(self,widget):
        '''
        Call function of the Open Solution Data Button
        '''
        fileFilter = gtk.FileFilter()
        fileFilter.set_name("SolutionData")
        fileFilter.add_pattern("*.pickle")
        
        filenames = self.LoadDialog(fileFilter)
        # check networkNames consistency
        networkNames =[]
        dataSetNumber = []
        for filename in filenames:
            networkNames.append(filename.split('/')[-1].split('_SolutionData_')[0])
            dataSetNumber.append(filename.split('/')[-1].split('_SolutionData_')[1].split('.')[0])
        
        if networkNames != [] and self.checkList(networkNames) == True:
            
            self.initializeSolutionData(networkNames[0],dataSetNumber)
                        
            title = ' '.join([self.mainTitle,self.currentNetworkName])
            self.set_title(title)
            
            for i in range(self.oldNumberOfVessels,0,-1):
                self.CBvessels.remove_text(i)
                self.CBvesselsComparison.remove_text(i)
                
            vID = self.vascularNetwork.vessels.keys()
            for i in range(1,self.numberOfVessels+1):
                self.CBvessels.append_text(' '.join([str(vID[i-1]),self.vascularNetwork.vessels[vID[i-1]].name]))
                self.CBvesselsComparison.append_text(' '.join([str(vID[i-1]),self.vascularNetwork.vessels[vID[i-1]].name]))
            
                    
            self.oldNumberOfVessels = self.numberOfVessels
            self.CBvessels.set_active(0)
    
            for i in range(self.oldDataNumber-1,-1,-1):
                self.CBsolData.remove_text(i)
                self.CBsolDataComparison.remove_text(i+1)
                
            for dN in self.dataNumbers:
                self.CBsolData.append_text(dN)
                self.CBsolDataComparison.append_text(dN)
                
            self.CBsolDataComparison.set_active(0) 
            self.CBsolData.set_active(0)  
            self.oldDataNumber = len(self.dataNumbers)
            
    
    def on_changedCBvessels(self,widget):
        '''
        call function of the combobox to choose vessel of the network
        '''
        if self.vascularNetwork != None:
            cbIndex = widget.get_active()-1
            if cbIndex < 0: self.currentVesselID = None
            else: self.currentVesselID = self.vascularNetwork.vessels.keys()[cbIndex] 
        
    def on_changedCBsolData(self,widget):
        '''
        call function of the combobox to choose solution data set number
        '''
        
        if self.dataNumbers != []: 
            self.currentDataNumber =  self.dataNumbers[widget.get_active()]
            self.dataSetDescription.set_text(self.dataSetDescriptions[self.currentDataNumber])
    
    def on_changedCBvesselsComparison(self,widget):
        '''
        call function of the combobox to choose vessel of the network for comparison
        '''
        if self.vascularNetwork != None:
            cbIndex = widget.get_active()-1
            if cbIndex < 0: self.currentVesselIDComparison  = None
            else: self.currentVesselIDComparison = self.vascularNetwork.vessels.keys()[cbIndex] 
    
    def on_changedCBsolDataComparison(self,widget):
        '''
        call function of the combobox to choose solution data set number
        '''
        
        if self.dataNumbers != [] and widget.get_active() != 0: 
            self.currentDataNumberComparison =  self.dataNumbers[widget.get_active()-1]
            self.dataSetDescriptionComparison.set_text(self.dataSetDescriptions[self.currentDataNumberComparison])
        elif self.dataNumbers != []:
            self.currentDataNumberComparison = None
            self.dataSetDescriptionComparison.set_text('')
    
    def invokedFrom3dViz(self): 
        '''
        call function when 2dviz is envoked from 3d visualisation by clicking on a vessel
        '''
        if os.path.isfile(cur+'/tempV.viz'):
            
            FILE = open( cur+'/tempV.viz', 'r')
            inPut = FILE.read().split('-')
            vesselId = int(inPut[0])
            dataNumber = inPut[1]
            networkName = inPut[2]
            FILE.close()
            os.remove(cur+'/tempV.viz')
            
            if networkName == self.vascularNetwork.name:
                
                description = [self.dataSetDescription.get_text(), self.dataSetDescriptionComparison.get_text()] 
                
                Window2dPlots( vesselId,
                               self.vascularNetwork,
                               dataNumber, 
                               self.solutionDataSets[self.dataNumbers.index(dataNumber)],
                               None,
                               None,
                               description,
                               None)
                
                if self.currentDataNumber != dataNumber:
                    self.currentDataNumber = dataNumber
                    self.CBsolData.set_active(self.dataNumbers.index(dataNumber))
            else:
                print "Error: different Networks loaded!"
        return True    
        
    def on_clickedPlots(self, widget):
        '''
        call function of the Show Plot buttom, it starts a new GUI class in an extra window
        showing the plots of the vessel.
        '''
        
        if self.currentVesselID != None:
            
            dataComparison = None
            if self.currentDataNumberComparison != None:
                dataComparison = self.solutionDataSets[self.dataNumbers.index(self.currentDataNumberComparison)]
            
            description = [self.dataSetDescription.get_text(), self.dataSetDescriptionComparison.get_text()] 
            
            external = self.externalData
            try:
                if self.externalData['Pressure'] is None and self.externalData['Flow'] is None:
                    external = None
            except: pass
            
            Window2dPlots( self.currentVesselID,
                           self.vascularNetwork,
                           self.currentDataNumber, 
                           self.solutionDataSets[self.dataNumbers.index(self.currentDataNumber)],
                           self.currentVesselIDComparison,
                           dataComparison,
                           description,
                           external)
    
    def on_clickedLoadExternal(self, widget):
        fileFilter = gtk.FileFilter()
        fileFilter.add_pattern("*.v1dfExD")
        filenames = self.LoadDialog(fileFilter)
        
        self.externalData
        
        try:
            fileName = filenames[0]
            self.externalData = loadExternalDataSet(fileName)
            self.extDataDescription.set_text(self.externalData['Description'])
        except: 
            pass
        
        
        
    def LoadDialog(self,fileFilter):
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
        
        directory = ''.join([cur,'/../NetworkFiles/'])
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
            
    def initializeSolutionData(self,networkName, dataSetNumber):
        '''
        Function loads the solution data defined by the filenames
        and creates plots for perssure and flow which are
        stored on the heard disk in a temp-directory
        '''
        vascularNetwork, solutionDataSets, simulationCaseDescriptions = loadSolutionDataFile(networkName, dataSetNumber)
        
        self.solutionDataSets = solutionDataSets        
        self.vascularNetwork = vascularNetwork
        self.simulationCaseDescriptions = simulationCaseDescriptions
                
        self.currentNetworkName = networkName
        self.numberOfVessels = len(self.vascularNetwork.vessels.keys())
        
        self.dataNumbers = dataSetNumber
        
        self.dataSetDescriptions = {}
        for description,dataN in simulationCaseDescriptions:
            if len(description) > 73:
                description = description[:72]+'\n'+description[72:]
            self.dataSetDescriptions[dataN] = description

class Window2dPlots(gtk.Window):
    '''
    This class creates a window showing the plots of a vessel.
    At the bottom a scrollbar where the different nodes can be scrolled 
    through!    
    '''
    def __init__(self,vesselId,vascularNetwork,dataNumber,solutionData,vesselIDComparison,solutionDataComparison,description,externalData):
        super(Window2dPlots,self).__init__()
        
        self.set_size_request(680,840)
        self.set_position(gtk.WIN_POS_CENTER)
        
        #plotData
        self.vessel = vascularNetwork.vessels[vesselId]
        self.vascularNetwork = vascularNetwork
        self.xAxisMaxValue = {'Time':vascularNetwork.totalTime, 'Space':self.vessel.length*100.0}
        self.vesselId = vesselId
        self.solutionData = solutionData
        self.solutionDataComparison = solutionDataComparison
        self.compare = False
        if solutionDataComparison != None: self.compare = True
        self.N = self.vessel.N
        self.val = 0
        
        self.externalData = externalData

        self.description = description[0]
        self.descriptionCompare = description[1]
        
        self.plotNormal = False
        self.linearWaveSplit = False
        self.nonlinearWaveSplit = False
        self.plotCFL = False
        self.plotImpedance = False
        
        self.axisX = 'Time'
        
        self.set_title(' '.join(['2d Visualisation - ',self.vascularNetwork.name,' : ',self.vessel.name,'- dataNumber ',dataNumber]))
        
        self.scale = gtk.HScale()
        Psol = self.solutionData['Pressure'][self.vesselId][:,[0]]
        self.hScaleTimeArray = np.linspace(0,vascularNetwork.totalTime,len(Psol))
        
        self.hScaleRange = {'Space':len(Psol), 'Time':self.vessel.N}
        self.hScaleIncrements = {'Space':20, 'Time':10}
        self.scale.set_range(0,self.hScaleRange[self.axisX]-1)
        self.scale.set_increments(1,10)
        self.scale.set_digits(0)
        self.scale.set_size_request(400,30)
        self.scale.connect('value-changed',self.on_changedSlider)
        self.scale.set_draw_value(False)        
        
        self.nodeLabel = gtk.Label('Node 0')
        
        self.legend = False
        cBlegend = gtk.CheckButton("show legend")
        cBlegend.connect("toggled", self.on_changedLegendBool)
        cBlegend.set_size_request(120,30)
        
        self.plotTitle = False
        cBtitle = gtk.CheckButton("show title")
        cBtitle.connect("toggled", self.on_changedTitleBool)
        cBtitle.set_size_request(120,30)
        
        self.medicalUnit = True
        cBmedical = gtk.CheckButton("use medical units")
        cBmedical.connect("toggled", self.on_changedMedicalBool)
        cBmedical.set_size_request(120,30)
        
        self.showDescription = False
        cBdescription = gtk.CheckButton('show descriptions')
        cBdescription.connect("toggled", self.on_changedDescription)
        cBdescription.set_size_request(120,30)
        
        self.buttonRenderMovie = gtk.Button("render movie")
        self.buttonRenderMovie.connect("clicked", self.on_clickedRenderMovie)
        self.buttonRenderMovie.set_size_request(120,30)
        self.buttonRenderMovie.set_sensitive(False)

        #### second row of buttons
        self.plotMinMaxPoints = False
        cBplotMinMaxPoints = gtk.CheckButton('plot MinMax-points')
        cBplotMinMaxPoints.connect("toggled", self.on_changedPlotMinMax)
        cBplotMinMaxPoints.set_size_request(120,30)
        
        cBchangeMinMaxDelta = gtk.Button('change MinMax-deltas')
        cBchangeMinMaxDelta.connect("clicked", self.on_clickedMinMaxDelta)
        cBchangeMinMaxDelta.set_size_request(120,30)
        cBchangeMinMaxDelta.set_sensitive(False)
        
        self.fig = Figure(figsize=(5,4), dpi=100)
        self.plot(0)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.canvas.set_size_request(640,690)
        
        toolbar = NavigationToolbar(self.canvas, self)
        
        cBtitle.set_active(0) 
        cBmedical.set_active(1) 
        cBdescription.set_active(0) 
        
        cbType = gtk.combo_box_new_text()
        cbType.connect("changed",self.on_changePlotType)
        cbType.append_text('Plot P,Q')
        cbType.append_text('Plot P,Q with linear wavesplitting')
        cbType.append_text('Plot P,Q with non-linear wavesplitting')
        cbType.append_text('Plot CFL, wave speed')
        cbType.append_text('Plot Area, Compliance, Impedance')
        cbType.set_active(0) 
        
        cbXaxis = gtk.combo_box_new_text()
        cbXaxis.connect("changed",self.on_changePlotXaxis)
        cbXaxis.append_text('Plot over Time')
        cbXaxis.append_text('Plot over Space')
        cbXaxis.set_active(0) 
              

        vbox = gtk.VBox(False,1)
        hboxLegend = gtk.HBox(False,10) 
        hboxCheckboxes2 = gtk.HBox(False,10) 
        hbox = gtk.HBox(False,1)
        hbox2 = gtk.HBox(False,1)
        
        
        # Checkbutton series
        hboxLegend.pack_start(cBlegend)
        hboxLegend.pack_start(cBdescription)
        #hboxLegend.pack_start(cBtitle)
        hboxLegend.pack_start(cBmedical)
        hboxLegend.pack_start(self.buttonRenderMovie)
        hboxCheckboxes2.pack_start(cBplotMinMaxPoints)
        hboxCheckboxes2.pack_start(cBchangeMinMaxDelta)

        # align pictures canvas
        alignIm = gtk.Alignment(0,1 , 1, 0)
        alignIm.add(self.canvas)
        # align node switcher scale
        hbox.pack_start(self.nodeLabel)
        hbox.pack_start(self.scale)
        alignHbox = gtk.Alignment(0, 1, 1, 0)
        alignHbox.add(hbox)
        # align combobox
        hbox2.pack_start(cbType)
        hbox2.pack_start(cbXaxis)
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
    
    def on_changedLegendBool(self,widget):
        self.legend =  widget.get_active()
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
    
    def on_changedTitleBool(self,widget):
        self.plotTitle =  widget.get_active()
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
    
    def on_changedDescription(self,widget):
        self.showDescription =  widget.get_active()
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
        
    def on_changedMedicalBool(self,widget):
        self.medicalUnit =  widget.get_active()
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
    
    def on_changedPlotMinMax(self, widget):
        self.plotMinMaxPoints =  widget.get_active()
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
    
    def on_clickedMinMaxDelta(self, widget):
        pass
    
    def on_changedSlider(self,widget):

        self.val = int(widget.get_value())
        
        if self.axisX == 'Time':
            self.nodeLabel.set_text('Node '+str(self.val))
        elif self.axisX == 'Space':
            NodeValue = str(self.hScaleTimeArray[self.val])
            self.nodeLabel.set_text(' '.join(['Time',str(NodeValue)[:6],'s']))
                
        self.plot(self.val)
        
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
        
    def on_clickedRenderMovie(self,widget):
        print " Render movie data: {} plots to save".format(len(self.hScaleTimeArray))
        
        
        dirName = ''.join([cur,'/','movieTempData/'])
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        
        self.scale.set_value(2)
        for n in xrange(len(self.hScaleTimeArray)):
             
            fileName = ''.join([dirName,str(n).zfill(3),'.png'])
            self.scale.set_value(n)
            self.canvas.draw()
            self.canvas.queue_resize()
            self.canvas.print_figure(fileName, dpi=60, format='png')

        
    def on_changePlotType(self,widget):
        '''
        call function of the combobox to choose type of plot
        '''
        
        self.plotNormal = False
        self.linearWaveSplit = False
        self.nonlinearWaveSplit = False
        self.plotCFL = False
        self.plotImpedance = False
        
        cbIndex = widget.get_active()
        
        if   cbIndex == 0:
            self.plotNormal = True
        elif cbIndex == 1:
            self.plotNormal = True
            self.linearWaveSplit = True
        elif cbIndex == 2:
            self.plotNormal = True
            self.nonlinearWaveSplit = True
        elif cbIndex == 3:
            self.plotCFL = True
        elif cbIndex == 4:
            self.plotImpedance = True
            
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
        
    def on_changePlotXaxis(self,widget):
        
        cbIndex = widget.get_active()
        
        if cbIndex == 0:
            self.axisX = 'Time'
            self.buttonRenderMovie.set_sensitive(False)
        elif cbIndex == 1:
            self.axisX = 'Space'
            NodeValue = str(self.hScaleTimeArray[self.val])
            self.nodeLabel.set_text(' '.join(['Time',str(NodeValue)[:6],'s']))
            self.buttonRenderMovie.set_sensitive(True)
        
        self.scale.set_range(0,self.hScaleRange[self.axisX]-1)
        
        self.plot(self.val)
        self.canvas.figure = self.fig
        self.fig.set_canvas(self.canvas)
        self.canvas.queue_resize()
        
    def plot(self, currentNode):
        n = currentNode
        #font sizes
        fontSizeLegend = 14
        fontSizeGlobal = 14
        fontSizeLabel  = 14
        
        startPlotTime = 0
        
        self.fig = Figure(figsize=(5,4), dpi=100,edgecolor='k')
        
        from matplotlib import rc
        from matplotlib import rcParams
        
        rcParams['text.usetex'] = True
        rcParams['text.latex.unicode'] = True
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.size'] = fontSizeGlobal

        unitP = 133.32
        unitF = 1e6
        unitPtext = '[mmHg]'
        unitFtext = '[ml/s]'
        if self.medicalUnit == False:
            unitP = 1
            unitF = 1
            unitPtext = '[Pa]'
            unitFtext = '[m^3/s]'
                
        Pmax = np.max(self.solutionData['Pressure'][self.vesselId])/unitP
        Pmin = np.min(self.solutionData['Pressure'][self.vesselId])/unitP
        Qmax = np.max(self.solutionData['Flow'][self.vesselId])*unitF
        Qmin = np.min(self.solutionData['Flow'][self.vesselId])*unitF
        if self.compare == True:
            PmaxComp = np.max(self.solutionDataComparison['Pressure'][self.vesselId])/unitP
            PminComp = np.min(self.solutionDataComparison['Pressure'][self.vesselId])/unitP
            QmaxComp = np.max(self.solutionDataComparison['Flow'][self.vesselId])*unitF
            QminComp = np.min(self.solutionDataComparison['Flow'][self.vesselId])*unitF
        
        extDataP = None
        extDataQ = None
        
        if self.axisX == 'Time':
        
            Psol = self.solutionData['Pressure'][self.vesselId][:,[n]]
            Qsol = self.solutionData['Flow'][self.vesselId][:,[n]]
        
            xAxisValues = np.linspace(0,self.xAxisMaxValue[self.axisX],len(Psol))
        
            Asol = self.solutionData['Area'][self.vesselId][:,[n]]
            csol = self.solutionData['WaveSpeed'][self.vesselId][:,[n]]
        
        
            compliancearray = self.vessel.C(Psol)[:,[n]]
                
            xlabel = 'Time [t]'
            
            if self.externalData is not None:
                ## external Pressure Data
                extDataP = self.externalData['SolutionData'][self.vesselId]['Pressure']
                if self.medicalUnit == True and self.externalData['PressureUnit'] == 'Pa':
                    extDataP = extDataP/133.32
                elif self.medicalUnit == False and self.externalData['PressureUnit'] == 'mmHg':
                    extDataP = extDataP*133.32
                extDataPtime = self.externalData['PressureTime']
                
                ## external Flow data
                extDataQ = self.externalData['SolutionData'][self.vesselId]['Flow']
                if self.medicalUnit == True and self.externalData['FlowUnit'] == 'm^3 s-1':
                    extDataQ = extDataP*1e6
                elif self.medicalUnit == False and self.externalData['FlowUnit'] == 'ml s-1':
                    extDataQ = extDataP/1e6
                extDataQtime = self.externalData['FlowTime']
            
            if self.compare == True:
                PsolComp = self.solutionDataComparison['Pressure'][self.vesselId][:,[n]]
                QsolComp = self.solutionDataComparison['Flow'][self.vesselId][:,[n]]
                AsolComp = self.solutionDataComparison['Area'][self.vesselId][:,[n]]
                csolComp = self.solutionDataComparison['WaveSpeed'][self.vesselId][:,[n]]
                
                xAxisValuesComp = np.linspace(0,self.xAxisMaxValue[self.axisX],len(PsolComp))
        
        elif self.axisX == 'Space':
            Psol = self.solutionData['Pressure'][self.vesselId][n]
            Qsol = self.solutionData['Flow'][self.vesselId][n]
        
            xAxisValues = np.linspace(0,self.xAxisMaxValue[self.axisX],len(Psol))
            
            Asol = self.solutionData['Area'][self.vesselId][n]
            csol = self.solutionData['WaveSpeed'][self.vesselId][n]
            
            if self.compare == True:
                PsolComp = self.solutionDataComparison['Pressure'][self.vesselId][n]
                QsolComp = self.solutionDataComparison['Flow'][self.vesselId][n]
                AsolComp = self.solutionDataComparison['Area'][self.vesselId][n]
                csolComp = self.solutionDataComparison['WaveSpeed'][self.vesselId][n]
                
                xAxisValuesComp = np.linspace(0,self.xAxisMaxValue[self.axisX],len(PsolComp))
            
            xlabel = 'length [cm]'
        
        if self.plotNormal:
            
            #Psol,Psol,Qsol,Qsol#
            p_f,p_b,q_f,q_b = mProc.linearWaveSplitting(Psol,Qsol,Asol,csol,self.vessel.rho)
            if self.compare == True:
                #PsolComp,PsolComp,QsolComp,QsolComp #
                p_fComp,p_bComp,q_fComp,q_bComp = mProc.linearWaveSplitting(PsolComp,QsolComp,AsolComp,csolComp,self.vessel.rho)

            if self.nonlinearWaveSplit:
                
                
                
                p_f,p_b,q_f,q_b = mProc.nonLinearWaveSplitting(Psol,Qsol,Asol,csol,compliancearray,self.vessel.rho)
                if self.compare == True:
                    p_fComp,p_bComp,q_fComp,q_bComp = mProc.nonLinearWaveSplitting(PsolComp,QsolComp,AsolComp,csolComp,self.vessel.rho)
            
            minMaxPoints = {}
            if self.plotMinMaxPoints == True:
                print "Min Max Points"
                try:
                    minMaxPoints['Psol'],minMaxPoints['PsolT'] = mProc.minMaxFunction(Psol.ravel()/unitP,xAxisValues,delta=0.025)
                    print "Pressure: "
                    pp(minMaxPoints['Psol'])
                    pp(minMaxPoints['PsolT'])
                except: pass 
                try:
                    minMaxPoints['Pf'],minMaxPoints['PfT'] = mProc.minMaxFunction(p_f.ravel()/unitP,xAxisValues,delta=0.025)
                    print "Pressure forward: "
                    pp(minMaxPoints['Pf'])
                    pp(minMaxPoints['PfT'])
                except: pass 
                try:
                    minMaxPoints['Pb'],minMaxPoints['PbT'] = mProc.minMaxFunction(p_b.ravel()/unitP,xAxisValues,delta=0.025)
                    print "Pressure backward: "
                    pp(minMaxPoints['Pb'])
                    pp(minMaxPoints['PbT'])
                except: pass 
                try:
                    minMaxPoints['Qsol'],minMaxPoints['QsolT'] = mProc.minMaxFunction(Qsol.ravel()*unitF,xAxisValues,delta=0.025)
                    print "Flow: "
                    pp(minMaxPoints['Qsol'])
                    pp(minMaxPoints['QsolT'])
                except: pass 
                try: 
                    minMaxPoints['Qf'],minMaxPoints['QfT'] = mProc.minMaxFunction(q_f.ravel()*unitF,xAxisValues,delta=0.025)
                    print "Flow forward: "
                    pp(minMaxPoints['Qf'])
                    pp(minMaxPoints['QfT'])
                except: pass 
                try: 
                    minMaxPoints['Qb'],minMaxPoints['QbT'] = mProc.minMaxFunction(q_b.ravel()*unitF,xAxisValues,delta=0.025)
                    print "Flow backward: "
                    pp(minMaxPoints['Qb'])
                    pp(minMaxPoints['QbT'])
                except: pass 
                
            
            meanPressure = Psol[0]
            
            #Pmax = max(np.max(p_f/unitP),np.max((p_b+meanPressure)/unitP),Pmax)
            #Pmin = min(np.min(p_f/unitP),np.min((p_b+meanPressure)/unitP),Pmin)
            Qmax = max([np.max(q_f*unitF),np.max(q_b*unitF),Qmax])
            Qmin = min([np.min(q_f*unitF),np.min(q_b*unitF),Qmin])   
            
            #if self.compare == True:
            #    Pmax = max(np.max(p_fComp/unitP),np.max(p_bComp/unitP),PmaxComp,Pmax)
            #    Pmin = min(np.min(p_fComp/unitP),np.min(p_bComp/unitP),PminComp,Pmin)
            #    Qmax = max([Qmax,np.max(q_fComp*unitF),np.max(q_bComp*unitF),QmaxComp])
            #    Qmin = min([Qmin,np.min(q_fComp*unitF),np.min(q_bComp*unitF),QminComp]) 
            if extDataP is not None:
                Pmax = max(Pmax,np.max(extDataP))
                Pmin = min(Pmin,np.min(extDataP))
            if extDataQ is not None:
                Qmax = max(Qmax,np.max(extDataQ))
                Qmin = min(Qmin,np.min(extDataQ))
                
            PminmaxCorr = (Pmax-Pmin)*0.1
            QminmaxCorr = (Qmax-Qmin)*0.1
            Pmax = Pmax + PminmaxCorr
            Pmin = Pmin - PminmaxCorr
            Qmax = Qmax + QminmaxCorr
            Qmin = Qmin - QminmaxCorr
                        
            timeSplit  = np.linspace(0,self.xAxisMaxValue[self.axisX],len(p_f))
            if self.compare == True:
                timeSplitComp  = np.linspace(0,self.xAxisMaxValue[self.axisX],len(q_fComp))
                   
            self.fig.subplots_adjust(hspace  = 0.4)   
            self.fig.subplots_adjust(right   = 0.85)
            self.fig.subplots_adjust(top     = 0.98)
            self.fig.subplots_adjust(bottom  = 0.2)
            self.fig.subplots_adjust(hspace  = 0.5)
            self.fig.set_figwidth(8.27)
            self.fig.set_figheight((11.69/3)*2.5)
            fontLegend = FontProperties()
            fontLegend.set_size(fontSizeLegend)
            
            ax = self.fig.add_subplot(2,1,1)
            
            if extDataP is not None:
                ax.plot(extDataPtime,extDataP,color='k' ,linestyle = '-',label='ext D', linewidth = 1.) 
            
            if self.plotMinMaxPoints:
                try:ax.plot(minMaxPoints['PsolT'],minMaxPoints['Psol'],color='b' ,linestyle = '',marker='o',label='', linewidth = 1.) 
                except:pass
            if self.compare == False:
                ax.plot(xAxisValues,Psol/unitP,color='b' ,linestyle = '-',label='total', linewidth = 1.) 
            else:
                ax.plot(xAxisValues,Psol/unitP,color='b' ,linestyle = '-',label='total case 1', linewidth = 1.)
                ax.plot(xAxisValuesComp,PsolComp/unitP,color='#E59900' ,linestyle = '-',label='total case 2', linewidth = 1.) 
            #ax.plot(timeNormal,Psol,color='b' ,linestyle = '-',label='total', linewidth = 1.5) 
            if self.nonlinearWaveSplit or self.linearWaveSplit:
                ax12 = ax.twinx()
                if self.compare == False:
                    ax12.plot(timeSplit,p_f/unitP,color='g' ,linestyle = '--',label='forward', linewidth = 1.) 
                    ax12.plot(timeSplit,p_b/unitP,color='m' ,linestyle = ':',label='backward', linewidth = 1.)
                else:
                    ax12.plot(timeSplit,p_f/unitP,color='b' ,linestyle = '--',label='forward - case 1', linewidth = 1.) 
                    ax12.plot(timeSplit,p_b/unitP,color='b' ,linestyle = ':',label='backward - case 1', linewidth = 1.)
                    ax12.plot(timeSplitComp,p_fComp/unitP,color='#E59900' ,linestyle = '--',label='forward - case 2', linewidth = 1.) 
                    ax12.plot(timeSplitComp,p_bComp/unitP,color='#E59900' ,linestyle = ':',label='backward - case 2', linewidth = 1.)
                if self.plotMinMaxPoints:
                    try:ax12.plot(minMaxPoints['PfT'],minMaxPoints['Pf'],color='b' ,linestyle = ' ',marker='o',label='', linewidth = 1.) 
                    except: pass 
                    try:ax12.plot(minMaxPoints['PbT'],minMaxPoints['Pb'],color='b' ,linestyle = '',marker='o',label='', linewidth = 1.) 
                    except: pass
                #ax12.set_ylim(Pmin-meanPressure/unitP,Pmax-meanPressure/unitP)
                ax12.set_ylabel('Contributions '+unitPtext, fontsize = fontSizeLabel)
                
                if self.legend == True:
                    handles, labels = ax12.get_legend_handles_labels()
                    ax12.legend(flip(handles, 2), flip(labels, 2),loc='upper left', bbox_to_anchor=(0.27, -0.15),frameon=False, ncol=2,prop = fontLegend)
            
            
            #ax.ticklabel_format(style='sci', axis='y', scilimits=(1,3), useOffset=False)
            ax.set_ylabel('Pressure '+unitPtext, fontsize = fontSizeLabel)
            ax.set_xlabel(xlabel, fontsize = fontSizeLabel)
            ax.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            ax.set_ylim(Pmin,Pmax)
            
            ## change frame and ticks visablility
            ax.spines['top'].set_visible(False)
            ax.tick_params(axis='x',top='off')
            if not (self.nonlinearWaveSplit or self.linearWaveSplit): 
                ax.spines['right'].set_visible(False)
                ax.tick_params(axis='y',right='off')
            
            
            if self.legend == True:
                ax.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
            
            
            ax2 = self.fig.add_subplot(2,1,2)    
            
            if self.plotMinMaxPoints:
                try:ax2.plot(minMaxPoints['QsolT'],minMaxPoints['Qsol'],color='r' ,linestyle = '',marker='o',label='', linewidth = 1.) 
                except:pass
            
            if extDataQ is not None:
                ax2.plot(extDataQtime,extDataQ,color='k' ,linestyle = '-',label='ext D', linewidth = 1.) 
            if self.compare == False:
                ax2.plot(xAxisValues,Qsol*unitF,color='r' ,linestyle = '-',label='total', linewidth = 1.) 
                
            else:
                ax2.plot(xAxisValues,Qsol*unitF,color='r' ,linestyle = '-',label='total - case 1', linewidth = 1.) 
                ax2.plot(xAxisValuesComp,QsolComp*unitF,color='g' ,linestyle = '-',label='total - case 2', linewidth = 1.) 
            
            #ax2.plot(timeNormal,Qsol,color='r' ,linestyle = '-',label='total', linewidth = 1.5) 
            if self.nonlinearWaveSplit or self.linearWaveSplit:
                ax22 = ax2.twinx()
                if self.compare == False:
                    ax22.plot(timeSplit,q_f*unitF,color='g' ,linestyle = '--',label='forward', linewidth = 1.) 
                    ax22.plot(timeSplit,q_b*unitF,color='m' ,linestyle = ':',label='backward', linewidth = 1.)
                else:
                    ax22.plot(timeSplit,q_f*unitF,color='r' ,linestyle = '--',label='forward - case 1', linewidth = 1.) 
                    ax22.plot(timeSplit,q_b*unitF,color='r' ,linestyle = ':',label='backward - case 1', linewidth = 1.)
                    ax22.plot(timeSplitComp,q_fComp*unitF,color='g' ,linestyle = '--',label='forward - case 2', linewidth = 1.) 
                    ax22.plot(timeSplitComp,q_bComp*unitF,color='g' ,linestyle = ':',label='backward - case 2', linewidth = 1.)
                if self.plotMinMaxPoints:
                    try:ax22.plot(minMaxPoints['QfT'],minMaxPoints['Qf'],color='r' ,linestyle = ' ',marker='o',label='', linewidth = 1.) 
                    except: pass 
                    try:ax22.plot(minMaxPoints['QbT'],minMaxPoints['Qb'],color='r' ,linestyle = '',marker='o',label='', linewidth = 1.) 
                    except: pass
                ax22.set_ylim(Qmin,Qmax)
                ax22.set_ylabel('Contributions '+unitFtext, fontsize = fontSizeLabel)
               
                if self.legend == True:
                    handles, labels = ax22.get_legend_handles_labels()
                    ax22.legend(flip(handles, 2), flip(labels, 2),loc='upper left', bbox_to_anchor=(0.27, -0.15),frameon=False, ncol=2,prop = fontLegend)
                    
            #ax2.ticklabel_format(style='sci', axis='y', scilimits=(2,0), useOffset=False)
            ax2.set_ylabel('Flow '+unitFtext, fontsize = fontSizeLabel)
            ax2.set_xlabel(xlabel, fontsize = fontSizeLabel)
            ax2.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            ax2.set_ylim(Qmin,Qmax)
            
            ## change frame and ticks visablility
            ax2.spines['top'].set_visible(False)
            ax2.tick_params(axis='x',top='off')
            if not (self.nonlinearWaveSplit or self.linearWaveSplit): 
                ax2.spines['right'].set_visible(False)
                ax2.tick_params(axis='y',right='off')
            
            if self.legend == True:
                ax2.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
                        
            if self.showDescription == True:
                if self.compare == True:
                    self.fig.text(0.05,0.04, 'Simulation case 1: '+self.description, fontsize=10)
                    self.fig.text(0.05,0.01, 'Simulation case 2: '+self.descriptionCompare, fontsize=10)
                else:
                    self.fig.text(0.05,0.04, 'Simulation case: '+self.description, fontsize=10)
        
        
        elif self.plotCFL:
            self.fig.subplots_adjust(hspace=0.5)   
            self.fig.set_figwidth(8.27)
            self.fig.set_figheight((11.69/3)*2)
        
            CFL = csol*self.vascularNetwork.dt/sum(self.vessel.dz)*len(self.vessel.dz)
            maxYRange = np.max(CFL)
            minYRange = np.min(CFL)
            #if self.compare == True:
                #CFLComp = csolComp*self.vascularNetwork.simulationContext['dt']/sum(self.vessel.dz)*len(self.vessel.dz)
                #maxYRange = max(max(CFL),max(CFLComp))
                #minYRange = min(min(CFL),min(CFLComp))
            ax = self.fig.add_subplot(2,1,1) 
            ax.plot(xAxisValues,CFL, color='k',linestyle = '-',label='CFL', linewidth = 1.)
            #if self.compare == True:
                #ax.plot(xAxisValuesComp,CFLComp, color='b',linestyle = '-',label='CFL Comparison', linewidth = 1.)
            ax.ticklabel_format(style='plain', axis='y', scilimits=(0,0), useOffset=False)
            ax.set_ylabel('CFL [-]',fontsize = fontSizeLabel)
            ax.set_xlabel(xlabel,fontsize = fontSizeLabel)
            ax.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            ax.set_ylim(minYRange-0.1*minYRange,maxYRange+0.1*maxYRange)
            
            ##### wave Speed
            ax2 = self.fig.add_subplot(2,1,2) 
            ax2.plot(xAxisValues,csol, color='b',linestyle = '-',label='wave Speed', linewidth = 1.)
            
            maxYRangeWS = np.max(csol)
            minYRangeWS = np.min(csol)
            
            ax2.ticklabel_format(style='plain', axis='y', scilimits=(0,0), useOffset=False)
            ax2.set_ylabel('wave speed [m/s]',fontsize = fontSizeLabel)
            ax2.set_xlabel(xlabel,fontsize = fontSizeLabel)
            ax2.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            ax2.set_ylim(minYRangeWS-0.1*minYRangeWS,maxYRangeWS+0.1*maxYRangeWS)
            
            if self.legend == True:
                ax2.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
            
            if self.showDescription == True:
                    self.fig.text(0.05,0.04, 'Simulation case: '+self.description, fontsize=10)
            
        elif self.plotImpedance:
            
            self.fig.subplots_adjust(hspace=0.5)   
            self.fig.set_figwidth(8.27)
            self.fig.set_figheight((11.69/3)*2)
            
            ### Area
            ax = self.fig.add_subplot(3,1,1) 
            AsolPlot = Asol*1000*1000
            ax.plot(xAxisValues,AsolPlot, color='k',linestyle = '-',label='Area', linewidth = 1.)
            
            ax.set_ylabel('Area [mm^2]',fontsize = fontSizeLabel)
            ax.ticklabel_format(style='plain', axis='y', scilimits=(0,0), useOffset=False)
            ax.set_xlabel(xlabel,fontsize = fontSizeLabel)
            ax.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            maxYRangeA = np.max(AsolPlot)
            minYRangeA = np.min(AsolPlot)
            ax.set_ylim(minYRangeA-0.1*minYRangeA,maxYRangeA+0.1*maxYRangeA)
            
            ##### Compliance
            ax2 = self.fig.add_subplot(3,1,2)
            
            Compliance = self.vessel.C(self.solutionData['Pressure'][self.vesselId])
            
            if self.axisX == 'Time': compliance_i = Compliance[:,[n]]*1000*1000*133.322
            if self.axisX == 'Space': compliance_i = Compliance[n]*1000*1000*133.322
            
            ax2.plot(xAxisValues,compliance_i, color='k',linestyle = '-',label='Compliance', linewidth = 1.)
            ax2.set_ylabel('Compliance\n[mm mmHg-1]',fontsize = fontSizeLabel)
            ax2.set_xlabel(xlabel,fontsize = fontSizeLabel)
            ax2.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            maxYRangeC = np.max(compliance_i)
            minYRangeC = np.min(compliance_i)
            ax2.set_ylim(minYRangeC-0.1*minYRangeC,maxYRangeC+0.1*maxYRangeC)
            
            #### Impedance
            ax3 = self.fig.add_subplot(3,1,3)
            Z = 1.0/(compliance_i*csol*1000)
            ax3.plot(xAxisValues,Z, color='k',linestyle = '-',label='Impedance', linewidth = 1.)
            #ax3.set_ylabel('Impedance')
            ax3.set_ylabel('Impedance\n[mmHg s mm-1]',fontsize = fontSizeLabel)
            ax3.ticklabel_format(axis='y', scilimits=(3,3), useOffset=True)
            ax3.set_xlabel(xlabel,fontsize = fontSizeLabel)
            ax3.set_xlim(startPlotTime,self.xAxisMaxValue[self.axisX])
            maxYRangeZ = np.max(Z)
            minYRangeZ = np.min(Z)
            ax3.set_ylim(minYRangeZ-0.1*minYRangeZ,maxYRangeZ+0.1*maxYRangeZ)
            
            
            if self.legend == True:
                ax.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
                ax2.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
                ax3.legend(loc='upper left', bbox_to_anchor=(-0.03, -0.15),frameon=False, ncol=1,prop = fontLegend)
            
            if self.showDescription == True:
                self.fig.text(0.05,0.04, 'Simulation case: '+self.description, fontsize=10)
            
def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])        


if __name__ == '__main__':
               
    optionsDict = moduleStartUp.parseOptions(['f','n','c'], visualisationOnly = True)
    
    networkName  = optionsDict['networkName']
    dataNumber   = optionsDict['dataNumber']
    connect      = optionsDict['connect']
             
    if networkName == None:
        Visualisation2D()
        gtk.main()
        
    else:
        Visualisation2D(networkName = networkName, dataSetNumber = dataNumber, connect = connect)
        gtk.main()
