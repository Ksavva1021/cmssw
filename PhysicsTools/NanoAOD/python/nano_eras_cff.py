from Configuration.Eras.Modifier_run2_egamma_2016_cff import run2_egamma_2016
from Configuration.Eras.Modifier_run2_egamma_2017_cff import run2_egamma_2017
from Configuration.Eras.Modifier_run2_egamma_2018_cff import run2_egamma_2018
from Configuration.Eras.Modifier_run2_jme_2016_cff import run2_jme_2016
from Configuration.Eras.Modifier_run2_jme_2017_cff import run2_jme_2017
from Configuration.Eras.Modifier_run2_muon_2016_cff import run2_muon_2016
from Configuration.Eras.Modifier_run2_muon_2018_cff import run2_muon_2018

from Configuration.Eras.Modifier_run2_HLTconditions_2016_cff import run2_HLTconditions_2016
from Configuration.Eras.Modifier_run2_HLTconditions_2017_cff import run2_HLTconditions_2017
from Configuration.Eras.Modifier_run2_HLTconditions_2018_cff import run2_HLTconditions_2018

from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_nanoAOD_92X_cff import run2_nanoAOD_92X
from Configuration.Eras.Modifier_run2_nanoAOD_94XMiniAODv1_cff import run2_nanoAOD_94XMiniAODv1
from Configuration.Eras.Modifier_run2_nanoAOD_94XMiniAODv2_cff import run2_nanoAOD_94XMiniAODv2
from Configuration.Eras.Modifier_run2_nanoAOD_94X2016_cff import run2_nanoAOD_94X2016
from Configuration.Eras.Modifier_run2_nanoAOD_102Xv1_cff import run2_nanoAOD_102Xv1
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv1_cff import run2_nanoAOD_106Xv1
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv2_cff import run2_nanoAOD_106Xv2
from Configuration.Eras.Modifier_tracker_apv_vfp30_2016_cff import tracker_apv_vfp30_2016

from Configuration.Eras.Modifier_run3_nanoAOD_devel_cff import run3_nanoAOD_devel

run2_nanoAOD_ANY = (
    run2_miniAOD_80XLegacy | run2_nanoAOD_92X | run2_nanoAOD_94X2016 | run2_nanoAOD_94XMiniAODv1 |  run2_nanoAOD_94XMiniAODv2 | run2_nanoAOD_102Xv1 | run2_nanoAOD_106Xv1 | run2_nanoAOD_106Xv2
)
run2_nanoAOD_preUL = (
    run2_miniAOD_80XLegacy | run2_nanoAOD_92X | run2_nanoAOD_94X2016 | run2_nanoAOD_94XMiniAODv1 |  run2_nanoAOD_94XMiniAODv2 | run2_nanoAOD_102Xv1 
) 
