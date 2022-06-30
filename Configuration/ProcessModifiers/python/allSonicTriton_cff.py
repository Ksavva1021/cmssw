import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton
from Configuration.ProcessModifiers.particleNetSonicTriton_cff import particleNetSonicTriton
from Configuration.ProcessModifiers.particleNetPTSonicTriton_cff import particleNetPTSonicTriton

# collect all SonicTriton-related process modifiers here
allSonicTriton = cms.ModifierChain(enableSonicTriton,particleNetSonicTriton)
