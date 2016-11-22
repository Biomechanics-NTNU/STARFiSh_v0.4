###### XML description version 4.3
##########################################################################################
## class VascularNetwork

# simulationContext
simulationContextElements = ['description',
                             'totalTime',
                             'CFL',
                             'timeSaveBegin',
                             'minSaveDt',
                             'maxMemory',
                             'gravitationalField',
                             'gravityConstant']

# initialisation controls
solverCalibrationElements =  [
                              'solvingSchemeField',
#                              'solvingSchemeConnections',
                              'rigidAreas',
                              'simplifyEigenvalues',
                              'riemannInvariantUnitBase',
                              'automaticGridAdaptation'  ]

# initialisation controls
initialisationControlsElements = ['initialsationMethod',
                                  'initMeanFlow',
                                  'initMeanPressure',
                                  'estimateWindkesselCompliance',
                                  'compPercentageWK3',
                                  'compPercentageTree',
                                  'compTotalSys']

#----------------------------------------------------------------------------------------#
# global fluid
globalFluidElements = ['my',
                       'rho',
                       'gamma']




## External Stimuli such as Head up tilt, Valsalva and Mueller maneuvers
headUpTiltElements = ['startAngle',
                       'stopAngle',
                       'startTime',
                       'duration']
valsalvaElements = []
hemmorrhageElements = []

externalStimulusElements = {'headUpTilt' : headUpTiltElements,
                    'valsalva'   : valsalvaElements,
                    'hemmorrhage' : hemmorrhageElements,
                    }

externalStimuliElementReference = {'externalStimulus':externalStimulusElements}


##
# TODO what's the deal with this?
vascularNetworkElements = ['simulationContext',
                          'solverCalibration',
                          'initialisationControls']

##########################################################################################
## Communicator objects

communicatorRealTimeViz  = ['comType', 'comId', 'vesselId', 'node', 'dn', 'quantitiesToPlot' ]
communicatorBaroreceptor = ['comType', 'comId', 'vesselId', 'node']


## need capital C in Communicator.. as it is also name of a class which is instatiated (as BCs)
communicatorReference = {'CommunicatorRealTimeViz' : communicatorRealTimeViz,
                         'CommunicatorBaroreceptor': communicatorBaroreceptor
                         }

##########################################################################################
## class BoundaryConditions

#----------------------------------------------------------------------------------------#
# Dictionary with boundary class references! Name_in_xml : class in classBoundaryConditions.py
bcTagsClassReferences = { # BoundaryConditions normal
                          'Flow-PhysiologicalFunction'       :'PhysiologicalFunction',
                          'Flow-Sinus'                       :'Sinus',
                          'Flow-Sinus2'                      :'Sinus2',
                          'Flow-ExpFunc'                     :'ExpFunc',
                          'Flow-CCAInflow'                   :'CCAInflow',
                          'Flow-AortaInflow'                 :'AortaInflow',
                          'Flow-AoBifInflow'                 :'AoBifInflow',
                          'Flow-ExperimentalInflow'          :'ExperimentalInflow',
                          'Flow-Adan55InflowFromfile'        :'Adan55InflowFromfile',
                          'Flow-RampMean'                    :'RampMean',
                          'Flow-Fourier'                     :'Fourier',
                          'Flow-PhysiologicalData'           :'PhysiologicalData',
                          'Flow-FromFile'                    :'FlowFromFile',
                          'Velocity-Gaussian'                :'expVelocity',
                          'Pressure-Sinus'                   :'Sinus',
                          'Pressure-Sinus2'                  :'Sinus2',
                          'Pressure-RampMean'                :'RampMean',
                          'ReflectionCoefficient'            :'ReflectionCoefficient',
                          'ReflectionCoefficientTimeVarying' :'ReflectionCoefficientTimeVarying',
                          'Resistance'                       :'Resistance',
                          'Windkessel-2ElementsDAE'          :'Windkessel2DAE',
                          'Windkessel-3ElementsDAE'          :'Windkessel3DAE',
                          'Windkessel-2Elements'             :'Windkessel2',
                          'Windkessel-3Elements'             :'Windkessel3',
                          'L-network'                        :'L_network',
                          'VaryingElastanceHeart'            :'VaryingElastance',
                          'VaryingElastanceSimple'           :'VaryingElastanceSimple',
                          'VaryingElastanceSimpleDAE'        :'VaryingElastanceSimpleDAE',
                          # BoundaryCondition names if 1 Vessel '_' == end-positionreference
                          '_Flow-PhysiologicalFunction'       :'PhysiologicalFunction',
                          '_Flow-Sinus'                       :'Sinus',
                          '_Flow-Sinus2'                      :'Sinus2',
                          '_Flow-ExpFunc'                     :'ExpFunc',
                          '_Flow-CCAInflow'                   :'CCAInflow',
                          '_Flow-AortaInflow'                 :'AortaInflow',
                          '_Flow-AoBifInflow'                 :'AoBifInflow',
                          '_Flow-ExperimentalInflow'          :'ExperimentalInflow',
                          '_Flow-Adan55InflowFromfile'        :'Adan55InflowFromfile',
                          '_Flow-RampMean'                    :'RampMean',
                          '_Flow-Fourier'                     :'Fourier',
                          '_Flow-PhysiologicalData'           :'PhysiologicalData',
                          '_Flow-FromFile'                    :'FlowFromFile',
                          '_Velocity-Gaussian'                :'expVelocity',
                          '_Pressure-Sinus'                   :'Sinus',
                          '_Pressure-Sinus2'                  :'Sinus2',
                          '_Pressure-RampMean'                :'RampMean',
                          '_ReflectionCoefficient'            :'ReflectionCoefficient',
                          '_ReflectionCoefficientTimeVarying' :'ReflectionCoefficientTimeVarying',
                          '_Resistance'                       :'Resistance',
                          '_Windkessel-2ElementsDAE'          :'Windkessel2DAE',
                          '_Windkessel-3ElementsDAE'          :'Windkessel3DAE',
                          '_Windkessel-2Elements'             :'Windkessel2',
                          '_Windkessel-3Elements'             :'Windkessel3',
                          '_L-network'                        :'L_network',
                          }

#----------------------------------------------------------------------------------------#
# Tag dictionary defines the variables a boundaryCondition needs as input via xml!
boundaryConditionElements = {
          # Type 1 attribute
          # BoundaryConditions normal
          'ReflectionCoefficient'       :['Rt'],
          'ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
          'Resistance'                  :['Rc'],
          'Windkessel-2ElementsDAE'     :['Rc','C'],
          'Windkessel-3ElementsDAE'     :['Rc','Rtotal','C','Z'],
          'Windkessel-2Elements'        :['Rc','C'],
          'Windkessel-3Elements'        :['Rc','Rtotal','C','Z'],
          'Flow-PhysiologicalFunction'  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Sinus'                  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Sinus2'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-ExpFunc'                :['runtimeEvaluation','prescribe'],
          'Flow-CCAInflow'              :['runtimeEvaluation','prescribe','freq','Npulse'],
          'Flow-AortaInflow'            :['runtimeEvaluation','prescribe','freq','Npulse'],
          'Flow-AoBifInflow'            :['runtimeEvaluation','prescribe','freq','Npulse'],
          'Flow-ExperimentalInflow'     :['runtimeEvaluation','prescribe','freq','Npulse'],
          'Flow-Adan55InflowFromfile'   :['runtimeEvaluation','prescribe','freq','Npulse'],
          'Flow-RampMean'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Velocity-Gaussian'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','gaussC','prescribe'],
          'Pressure-Sinus'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Pressure-Sinus2'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Pressure-RampMean'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Fourier'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
          'Flow-PhysiologicalData'      :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-FromFile'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
          'L-network'                   :['Z','C'],
          'VaryingElastanceHeart'       :['T', 'Emax', 'Emin', 'Tpeak', 'V0', 'K'],
          'VaryingElastanceSimple'      :['T', 'Emax', 'Emin', 'Tpeak', 'V0', 'K'],
          'VaryingElastanceSimpleDAE'   :['T', 'Emax', 'Emin', 'Tpeak', 'V0', 'residualName'],
          # BoundaryCondition names if 1 Vessel '_' == end-positionreference
          '_Flow-PhysiologicalFunction' :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Sinus'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Sinus2'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-ExpFunc'               :['runtimeEvaluation','prescribe'],
          '_Flow-CCAInflow'             :['runtimeEvaluation','prescribe','freq','Npulse'],
          '_Flow-AortaInflow'           :['runtimeEvaluation','prescribe','freq','Npulse'],
          '_Flow-AoBifInflow'           :['runtimeEvaluation','prescribe','freq','Npulse'],
          '_Flow-ExperimentalInflow'    :['runtimeEvaluation','prescribe','freq','Npulse'],
          '_Flow-Adan55InflowFromfile'  :['runtimeEvaluation','prescribe','freq','Npulse'],
          '_Flow-RampMean'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Fourier'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
          '_Flow-PhysiologicalData'     :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-FromFile'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
          '_Velocity-Gaussian'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','C','prescribe'],
          '_Pressure-Sinus'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Pressure-Sinus2'            :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Pressure-RampMean'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_ReflectionCoefficient'      :['Rt'],
          '_ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
          '_Resistance'                 :['Rc'],
          '_Windkessel-2ElementsDAE'    :['Rc','C'],
          '_Windkessel-2Elements'       :['Rc','C'],
          '_Windkessel-3ElementsDAE'    :['Rc','Rtotal','C','Z'],
          '_Windkessel-3Elements'       :['Rc','Rtotal','C','Z'],
          '_L-network'                  :['Z','C'],
          'None'                        :['']}
##########################################################################################
## class Vessel

vesselAttributes = ['Id',
                    'name']

vesselTopologyElements = ['leftDaughter',
                          'rightDaughter',
                          'angleYMother']

vesselGeometryElements = ['geometryType',
                          'length',
                          'radiusProximal',
                          'radiusDistal',
                          'N']

vesselComplianceElements = {'Laplace'     :['complianceType','constantCompliance','externalPressure','Ps','As','betaLaplace'],
                            'Laplace2'    :['complianceType','constantCompliance','externalPressure','Ps','As','wallThickness','youngModulus'],
                            'Exponential' :['complianceType','constantCompliance','externalPressure','Ps','As','betaExponential'],
                            'Hayashi'     :['complianceType','constantCompliance','externalPressure','Ps','As','betaHayashi'],
                            'HayashiEmpirical' :['complianceType','constantCompliance','externalPressure','Ps','As','betaHayashi'],
                            'Reymond'     :['complianceType','constantCompliance','externalPressure','Ps','As','Cs','PmaxC','Pwidth','a1','b1']
                            }

vesselFluidElements = ['applyGlobalFluid',
                       'my',
                       'rho',
                       'gamma']

vesselElementReference = {'topology'     :vesselTopologyElements,
                          'geometry'     :vesselGeometryElements,
                          'compliance'   :vesselComplianceElements,
                          'fluid'        :vesselFluidElements}

vesselElements = ['topology',
                  'geometry',
                  'compliance',
                  'fluid']

xmlElements  = [ 'simulationContext',
                 'solverCalibration',
                 'initialisationControls',
                 'globalFluid',
                 'externalStimuli',
                 'communicators',
                 'boundaryConditions',
                 'randomInputManager',
                 'measurementRoutine',
                 'vessels',
                 "venousPool"]

xmlElementsReference = {'simulationContext'     : simulationContextElements,
                        'solverCalibration'     : solverCalibrationElements,
                        'initialisationControls': initialisationControlsElements,
                        'boundaryConditions'    : boundaryConditionElements,
                        'globalFluid'           : globalFluidElements,
                        'externalStimuli'       : externalStimulusElements,
                        'communicators'         : communicatorReference, #TODO why is this Reference?
                        'randomInputManager'    : None,
                        'measurementRoutine'     : None,
                        'vessels'               : vesselElements,
                        'venousPool' : None }


