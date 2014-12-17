trigger_list = [
        'L0HadronDecision',
        'L0MuonDecision',
        'L0DiMuonDecision',
        
        'Hlt1SingleMuonNoIPDecision',
        'Hlt1SingleMuonHighPTDecision',
        'Hlt1DiMuonHighMassDecision',
        'Hlt1DiMuonLowMassDecision',
        'Hlt1TrackAllL0Decision',
        'Hlt1TrackMuonDecision',
        'Hlt1MuTrackDecision',
        'Hlt1GlobalDecision',

        'Hlt2TopoOSTF2BodyDecision',
        'Hlt2TopoOSTF3BodyDecision',
        'Hlt2TopoOSTF4BodyDecision',
        'Hlt2Topo2BodySimpleDecision',
        'Hlt2Topo3BodySimpleDecision',
        'Hlt2Topo4BodySimpleDecision',
        'Hlt2Topo2BodyBBDTDecision',
        'Hlt2Topo3BodyBBDTDecision',
        'Hlt2Topo4BodyBBDTDecision',
        'Hlt2TopoMu2BodyBBDTDecision',
        'Hlt2TopoMu3BodyBBDTDecision',
        'Hlt2MuonFromHLT1Decision',
        'Hlt2SingleMuonDecision',
        'Hlt2SingleMuonHighPTDecision',
        'Hlt2SingleMuonLowPTDecision',
        'Hlt2DisplVerticesLowMassSingleDecision',
        'Hlt2DisplVerticesHighMassSingleDecision',
        'Hlt2DisplVerticesDoubleDecision',
        'Hlt2DisplVerticesSinglePostScaledDecision',
        'Hlt2DisplVerticesHighFDSingleDecision',
        'Hlt2DisplVerticesSingleDownDecision',
        'Hlt2DiMuonDecision',
        'Hlt2DiMuonLowMassDecision',
        'Hlt2DiMuonBDecision',
         'Hlt2DiMuonDetachedDecision',
        'Hlt2MuTrackDecision'
        'Hlt2GlobalDecision'
 ]

##############################################################################
 # histos for Lc -> p mu mu ; 
#####################################################################################
from Gaudi.Configuration import *
from Configurables import GaudiSequencer, FilterDesktop, DeterministicPrescaler
from Configurables import CombineParticles
#from Configurables import LoKi__Hybrid__FilterCriterion as LoKiFilterCriterion
from Configurables import LoKi__Hybrid__PlotTool as PlotTool

from Configurables import DecayTreeTuple, TupleToolTrigger,  TupleToolDecay, TupleToolTISTOS
from DaVinci.Configuration import *
from Gaudi.Configuration import *
from Configurables import GaudiSequencer, DecayTreeTuple, TupleToolDecay, TupleToolTrigger, TupleToolTISTOS, FilterDesktop
from PhysSelPython.Wrappers import AutomaticData, DataOnDemand, Selection, SelectionSequence

from os import environ
from GaudiKernel.SystemOfUnits import *
from Gaudi.Configuration import *
from PhysSelPython.Wrappers import Selection, SelectionSequence, DataOnDemand, AutomaticData
from Configurables import EventTuple, TupleToolMCBackgroundInfo, BackgroundCategory, TupleToolDecay

from Configurables import CombineParticles
from Configurables import OfflineVertexFitter, DaVinci

from PhysSelPython.Wrappers import SelectionSequence

###########################################################################
###########################################################################
#DaVinci().DDDBtag = "Sim08-20130503-1"
#DaVinci().CondDBtag = "Sim08-20130503-1-vc-mu100" 
###########################################################################
###########################################################################

from StrippingConf.Configuration import StrippingConf, StrippingStream
from StrippingSettings.Utils import strippingConfiguration
from StrippingArchive.Utils import buildStreams, cloneLinesFromStream
from StrippingArchive import strippingArchive


eventNodeKiller = EventNodeKiller('Stripkiller')
eventNodeKiller.Nodes = [ '/Event/AllStreams', '/Event/Strip' ]

stripping='stripping20'
#get the configuration dictionary from the database
config  = strippingConfiguration(stripping)
#get the line builders from the archive
archive = strippingArchive(stripping)
streams = buildStreams(stripping = config, archive = archive)
MyStream = StrippingStream("MyStream")
MyLines = [ 'StrippingTau23MuTau2PMuMuLine' ]

for stream in streams:
    for line in stream.lines:
        if line.name() in MyLines:
                        print 'StrippingTau23MuTau2PMuMuLine'
                        MyStream.appendLines( [ line ] )


from Configurables import ProcStatusCheck
filterBadEvents = ProcStatusCheck()

sc = StrippingConf( Streams = [ MyStream ],
                                        MaxCandidates = 2000,
                                        AcceptBadEvents = False,
                                        BadEventSelection = filterBadEvents )


from Configurables import DecayTreeTuple

tuple= DecayTreeTuple('LC_OS')
tuple.Inputs = [ "Phys/Tau23MuTau2PMuMuLine/Particles" ]
tuple.UseLoKiDecayFinders = False
#tuple.IgnoreP2PVFromInputLocations = True # ignore all stored Particle -> PV relations
tuple.ReFitPVs = True # re-fit the PVs
tuple.Decay = '[Lambda_c+ -> ^p+ ^mu+ ^mu- ]cc'
tuple.ToolList +=  [       "TupleToolEventInfo"
                           , "TupleToolGeneration"
                           , "TupleToolMCTruth"
                           , "TupleToolMCBackgroundInfo"
                           , "TupleToolPrimaries"
                           , "TupleToolTrackInfo"
                           , "TupleToolPid"
                           , "TupleToolGeometry"
                           , "TupleToolKinematic"
#                           , "TupleToolMuonVariables"
                           , "TupleToolPropertime"
                           , "LoKi::Hybrid::TupleTool/LoKiTool"
                           , "TupleToolMuonIso_LC"
                           , "TupleToolEventInfo"
]
tuple.addTool( TupleToolTrigger, name='TupleToolTrigger' )      
tuple.TupleToolTrigger.Verbose = True                           
tuple.TupleToolTrigger.TriggerList = trigger_list               
                                                 
tuple.addTool( TupleToolTISTOS, name='TupleToolTISTOS' )     
tuple.TupleToolTISTOS.Verbose = True                         
tuple.TupleToolTISTOS.TriggerList = trigger_list             
tuple.ToolList += [ "TupleToolTISTOS" ]                      






#######################################33333
## Tuple for wrong sign


tuple2= DecayTreeTuple('LC_SS')
tuple2.Inputs = [ "Phys/Tau23MuTau2PMuMuLine/Particles" ]
#tuple2.IgnoreP2PVFromInputLocations = True # ignore all stored Particle -> PV relations 
tuple2.ReFitPVs = True # re-fit the PVs 
tuple2.UseLoKiDecayFinders = False                                                
tuple2.Decay = '[Lambda_c+ -> p~- mu+ mu+ ]cc'
tuple2.ToolList +=  [       "TupleToolEventInfo"
                           , "TupleToolGeneration"
                           , "TupleToolMCTruth"
                           , "TupleToolMCBackgroundInfo"
                           , "TupleToolPrimaries"
                           , "TupleToolTrackInfo"
                           , "TupleToolPid"
                           , "TupleToolGeometry"
                           , "TupleToolKinematic"
#                           , "TupleToolMuonVariables"
                           , "TupleToolPropertime"
                           , "LoKi::Hybrid::TupleTool/LoKiTool"
                           , "TupleToolMuonIso_LC"
                           , "TupleToolEventInfo"
]
tuple2.addTool( TupleToolTrigger, name='TupleToolTrigger' )
tuple2.TupleToolTrigger.Verbose = True
tuple2.TupleToolTrigger.TriggerList = trigger_list

tuple2.addTool( TupleToolTISTOS, name='TupleToolTISTOS' )
tuple2.TupleToolTISTOS.Verbose = True
tuple2.TupleToolTISTOS.TriggerList = trigger_list
tuple2.ToolList += [ "TupleToolTISTOS" ]

                                                                      
               


                          
from Configurables import DaVinci
##from DaVinci.Configuration import 
##############################################################################

DaVinci().EvtMax = -1
DaVinci().PrintFreq = 100 
DaVinci().DataType = "2011"
#DaVinci().Simulation   = True
DaVinci().Simulation   = False
DaVinci().InputType = "DST"
DaVinci().TupleFile = "Lc2pmumu.root"
DaVinci().appendToMainSequence( [ eventNodeKiller ] )
DaVinci().appendToMainSequence( [ sc.sequence() ] )
DaVinci().UserAlgorithms = [tuple, tuple2 ] 
