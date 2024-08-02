# ZHTauTau

Examples for  *ee->ZH->MuMuTauTau* ,  *ee->ZH->QQTauTau*,  *ee->ZH->eeTauTau* from a tau prototype.

Trees produced with this fcc script: <a href="https://github.com/mcepeda/FCCAnalyses/tree/master/examples/FCCee/higgs/tautau/test/analysis_stage1_fromjets_win23.py">  FCCee/higgs/tautau/test/analysis_stage1_fromjets_win23.py </a>

The Tau Algo is based on finding tracks and photons within a Jet. See in  be seen in  <a href="https://github.com/mcepeda/FCCAnalyses/blob/master/analyzers/dataframe/src/myUtils.cc#L2559"> analyzers/dataframe/src/myUtils.cc </a>

Command to run fccanalysis:

```
fccanalysis run examples/FCCee/higgs/tautau/test/analysis_stage1_fromjets_win23.py --files-list /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_eeH_Htautau_ecm240/events_* --output "wzp6_ee_eeH_Htautau_ecm240.root" 
```

WARNING: this is work in progress :).  Both the taus and the selection.

<img src="MuMuTauTau_collmass_example.png" width=400> <img src="QQTauTau_collmass_example.png " width=400>

Trees at ciemat: "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/DelphesTrees/taus16Feb_eekt4Jets_idInJet/"

Example:

```
# To create the histograms:
python runMuMuTauTau.py 

# To make the plots:
python pintarMuMuTauTau.py
```



