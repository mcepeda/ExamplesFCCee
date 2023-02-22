#! /usr/bin/env python3
#
# 1) Setup ROOT with: 
#   "source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc8-opt/setup.sh
# 2) Run this script as:
#   "python -i [python_script]"
#

import os, sys, math
import ROOT
import multiprocessing

# Enable multi-threading (parallel processing, do not use more than 8 cores)
maxcpus = multiprocessing.cpu_count()
ROOT.EnableImplicitMT(min(maxcpus,8))
usecpus = max(ROOT.GetThreadPoolSize(),1)
print("ROOT: %s CPUs used out of %s available ..." % (usecpus, maxcpus))

# Archivos de entrada
# Los archivos estan en 'path'. Vamos a usar tres archivos, X_skimmed_reduced.root, donde X = eeHZ
# (sennal) +  eeZZ,  eeWW (fondos)
sampleName=["wzp6_ee_uuHorddH_Htautau_ecm240","wzp6_ee_bbH_Htautau_ecm240","wzp6_ee_ssH_Htautau_ecm240","wzp6_ee_ccH_Htautau_ecm240","wzp6_ee_qqH_ecm240","p8_ee_ZZ_ecm240", "p8_ee_WW_ecm240"]
nProcesses=len(sampleName)
path=' /nfs/cms/cepeda/FCC/taus16Feb_eekt4Jets_idInJet/'

# Cada rootfile tiene un 'tree' llamado events que vamos a
# convertir en un dataframe y guardar en una lista de dataframes. 

print ("Reading Trees")
df = {}
allCutsReport={}

for p in sampleName:
   if p=="p8_ee_ZZ_ecm240":
      names = ROOT.std.vector('string')()
      for i in range(0,14):
          names.push_back(path+'p8_ee_ZZ_ecm240_{}.root'.format(i))
      print (names)
      df[p] = ROOT.RDataFrame("events",(names))
   elif p=="p8_ee_WW_ecm240":
      names = ROOT.std.vector('string')()
      for i in range(0,19):
          names.push_back(path+'p8_ee_WW_ecm240_{}.root'.format(i))
      print (names)
      df[p] = ROOT.RDataFrame("events",(names))
   else:
          df[p] = ROOT.RDataFrame("events",(os.path.join(path,"{}.root".format(p))))

processes = list(df.keys())

ROOT.gInterpreter.Declare("""
   
   ROOT::VecOps::RVec<int> FindBest(ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>> Jets, double mass){ 

        double minDistance=10000;
        int lead=-1;
        int second=-1;

        for (int i=0; i<Jets.size(); i++){
                for (int j=i+1; j<Jets.size(); j++){
                        ROOT::Math::PxPyPzMVector DiJet=Jets[i]+Jets[j];
                        //std::cout<<DiJet.M()<<std::endl;
                        if ( fabs(DiJet.M()-mass)<minDistance) {
                                        minDistance=fabs(DiJet.M()-mass);
                                        lead=i;
                                        second=j; 
                        }
        }
        }

  ROOT::VecOps::RVec<int> idx;
  idx.push_back(lead);
  idx.push_back(second); 

  return idx;
  }     
"""
)


print(processes)

# Ahora hacemos la seleccion y definimos variables nuevas 

for i in range(nProcesses):
    p=processes[i]

    df[p] = df[p].Define("genpart_idx","Nonzero(GenPart_PDG)")
    df[p] = df[p].Define("selGenTauHad_idx","genpart_idx[abs(GenPart_PDG)==15&&GenPart_type>=0]") 
    df[p] = df[p].Define("GenTauHad_p4","Take(GenPart_VisP4,selGenTauHad_idx)")
    df[p] = df[p].Define("selGenTauMu_idx","genpart_idx[abs(GenPart_PDG)==15&&GenPart_type==-13]")
    df[p] = df[p].Define("GenTauMu_p4","Take(GenPart_VisP4,selGenTauMu_idx)")
    df[p] = df[p].Define("selGenTauEle_idx","genpart_idx[abs(GenPart_PDG)==15&&GenPart_type==-11]") 
    df[p] = df[p].Define("GenTauEle_p4","Take(GenPart_VisP4,selGenTauEle_idx)")
    df[p] = df[p].Define("NGenTauHad","(Int_t)GenTauHad_p4.size()")
    df[p] = df[p].Define("NGenTauMu","(Int_t)GenTauMu_p4.size()")
    df[p] = df[p].Define("NGenTauEle","(Int_t)GenTauEle_p4.size()")

#    df[p].Display({"NGenTauHad","NGenTauMu","NGenTauEle"},10).Print()
    df[p] = df[p].Define("TauTauEvent","NGenTauHad")
#    df[p] = df[p].Filter("NGenTauHad==2","Only Tau_h Tau_h")  # Si quieres filtrar solo  HTauTau

#    df[p] = df[p].Filter("NMuon==0","No Muons")
    df[p] = df[p].Define("muon_idx","Nonzero(Muon_pt)")
    df[p] = df[p].Define("selmuon_idx","muon_idx[(Muon_pt>10)]") 
    df[p] = df[p].Define("selMuons_pt","Take(Muon_pt,selmuon_idx)")
    df[p] = df[p].Define("NSelMuon","(Int_t)selMuons_pt.size()")
#    df[p] = df[p].Filter("NElectron==0","No Electrons")
    df[p] = df[p].Define("Electron_idx","Nonzero(Electron_pt)")
    df[p] = df[p].Define("selElectron_idx","Electron_idx[(Electron_pt>10)]") 
    df[p] = df[p].Define("selElectrons_pt","Take(Electron_pt,selElectron_idx)")
    df[p] = df[p].Define("NSelElectron","(Int_t)selElectrons_pt.size()")

    df[p] = df[p].Filter("NSelMuon==0","No Muons")
    df[p] = df[p].Filter("NSelElectron==0","No Electrons")

    # Ahora vamos a buscar los Taus. Para hacer esto mejor podriamos comprobar el orden
    # Que funciona mejor, ordenar por pt? Por aislamiento? Cuantos Taus hay? 
    df[p] = df[p].Filter("NTauFromJet == 2", "Events with exactly two taus")
    df[p] = df[p].Define("Taus_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>  > > ( TauFromJet_px,TauFromJet_py,TauFromJet_pz,TauFromJet_mass)")
#    df[p] = df[p].Define("bestHiggs","FindBest(Taus_p4,125)")
    df[p] = df[p].Define("Lead","0")#bestHiggs[0]")
    df[p] = df[p].Define("Second","1")#bestHiggs[1]")
    df[p] = df[p].Define("TauLead_p4","Taus_p4[Lead]")
    df[p] = df[p].Define("TauSecond_p4","Taus_p4[Second]")

    df[p] = df[p].Filter("TauFromJet_pt[Lead]>10 && TauFromJet_pt[Second]>10", "Tau Pt>10 GeV")
    df[p] = df[p].Filter("TauFromJet_type[Lead]>=0 && TauFromJet_type[Second]>=0",   "Identified Taus")

    # Demasiado WW. Limpiemos la muestra:
#    df[p] = df[p].Filter("TauFromJet_type[Lead]<=10 && TauFromJet_type[Second]<=10",   "No 3 + Photon")
    df[p] = df[p].Filter("TauFromJet_mass[Lead]<2 && TauFromJet_mass[Second]<2",   "Tau mass < 2 GeV")

    df[p] = df[p].Define("DThetaTaus","abs(TauLead_p4.Theta()-TauSecond_p4.Theta())")
    df[p] = df[p].Define("DPhiTaus","abs(TauLead_p4.Phi()-TauSecond_p4.Phi())")

    df[p] = df[p].Filter("TauFromJet_charge[Lead] != TauFromJet_charge[Second]", "Taus with opposite charge")

    # Reconstruimos el DiTau, Sumando los cuadrimomentos (p4)
    df[p] = df[p].Define("DiTau_p4","Taus_p4[Lead]+Taus_p4[Second]")

    df[p] = df[p].Define("TauLead_type","TauFromJet_type[Lead]")
    df[p] = df[p].Define("TauSecond_type","TauFromJet_type[Second]")
    df[p] = df[p].Define("TauLead_mass","TauFromJet_mass[Lead]")
    df[p] = df[p].Define("TauSecond_mass","TauFromJet_mass[Second]")

    df[p] = df[p].Define("DiTau_vis_mass", "DiTau_p4.M()")
    df[p] = df[p].Define("DiTau_Pt", "DiTau_p4.Pt()")

    # Minimo de masa para el DiTau
    df[p] = df[p].Filter("DiTau_vis_mass>40", "DiTauMass > 40 GeV")

    # Define missing momentum
    df[p] = df[p].Define('Missing_e', ' sqrt( Missing_px*Missing_px+Missing_py*Missing_py+Missing_pz*Missing_pz )')

    # El momento 'que falta', Missing Pt y su angulo cos(Theta_miss)
    df[p] = df[p].Define('Missing_p4', "ROOT::Math::PxPyPzEVector(Missing_px,Missing_py,Missing_pz,Missing_e)")
    df[p] = df[p].Define('Missing_costheta', 'abs(cos(Missing_p4.Theta()))')


    #Aplicamos collinear approximation
    df[p] = df[p].Define("p12","TauLead_p4.Py()*TauSecond_p4.Px()-TauLead_p4.Px()*TauSecond_p4.Py()")
    df[p] = df[p].Define("r0","abs((Missing_p4.Py()*TauLead_p4.Px()-Missing_p4.Px()*TauLead_p4.Py())/p12)")
    df[p] = df[p].Define("f0","1/(1+r0)")
    df[p] = df[p].Define("r1","abs((Missing_p4.Py()*TauSecond_p4.Px()-Missing_p4.Px()*TauSecond_p4.Py())/p12)")
    df[p] = df[p].Define("f1","1/(1+r1)")
    df[p] = df[p].Define("DiTau_coll_mass","DiTau_vis_mass/sqrt(f0*f1)")

    df[p] = df[p].Define("jet_idx","Nonzero(JetReRun_pt)")
    df[p] = df[p].Define("selJets_idx","jet_idx[(JetReRun_pt>15)&&JetReRun_tauID<0]") #&&JetReRun_nmu==0&&JetReRun_nel==0&&JetReRun_tauID<0)]")
    df[p] = df[p].Define("selJets_pt","Take(JetReRun_pt,selJets_idx)")
    df[p] = df[p].Define("selJets_px","Take(JetReRun_px,selJets_idx)")
    df[p] = df[p].Define("selJets_py","Take(JetReRun_py,selJets_idx)")
    df[p] = df[p].Define("selJets_pz","Take(JetReRun_pz,selJets_idx)")
    df[p] = df[p].Define("selJets_energy","Take(JetReRun_energy,selJets_idx)")
    df[p] = df[p].Define("NGoodJets","(Int_t)selJets_pt.size()")
    df[p] = df[p].Filter("NGoodJets >= 2","Two jets")
    df[p] =  df[p].Define("Jets_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>>(selJets_px,selJets_py,selJets_pz,selJets_energy)")

#    df[p].Display({"JetReRun_pt" ,"JetReRun_nmu","JetReRun_nel","JetReRun_tauID","jet_idx","selJets_idx","selJets_pt","NGoodJets"},10).Print()
   # df[p].Display({"NGoodJets","selJets_pt","selJets_energy"},10).Print()

    df[p] = df[p].Define("DiJet_p4","Jets_p4[0]+Jets_p4[1]")
    df[p] = df[p].Define("DiJet_mass", "DiJet_p4.M()")
    df[p] = df[p].Define("DiJet_Pt", "DiJet_p4.Pt()")

    df[p] = df[p].Filter("selJets_pt[0]>20&&selJets_pt[1]>20","Jet Pt>20")
#    df[p] = df[p].Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge")
    df[p] = df[p].Filter("DiJet_mass>80 && DiJet_mass<100","80 < Mreco(diJet) < 100 GeV")
    df[p] = df[p].Filter("DiJet_Pt>20 && DiJet_Pt<70" ,"20 < Pt(DiJet) < 70 GeV")

    df[p] = df[p].Define("DThetaJets","abs(Jets_p4[0].Theta()-Jets_p4[1].Theta())")
    df[p] = df[p].Define("DPhiJets","abs(Jets_p4[0].Phi()-Jets_p4[1].Phi())")

    # Reconstruimos el recoil, que corresponde a la masa del Higgs:
    df[p] = df[p].Define("p4total","ROOT::Math::PxPyPzEVector(0.,0.,0.,240.)")
    df[p] = df[p].Define("recoil","(p4total-DiJet_p4).M()")

    # Filtramos en la masa del Recoil 
#    df[p] = df[p].Filter("recoil>100", "Recoil>100")
    df[p] = df[p].Filter("recoil>120 && recoil<140", "120<Recoil<140 GeV")
#    df[p] = df[p].Filter("recoil<100", "Recoil<100")


    # Por que tengo mas sucesos de qqH que de uuH/ddH+bbH+ccH+ssH pasandome los cortes? 
    # Que tipo de desintegraciones del Higgs selecciono? Y del Z?
    df[p] = df[p].Define("higgsDaugh_idx","genpart_idx[ GenPart_momPDG==25]") 
    df[p] = df[p].Define("higgsDaugh_PDG","Take(GenPart_PDG,higgsDaugh_idx)")

    df[p] = df[p].Define("zDaugh_idx","genpart_idx[ GenPart_momPDG==23 && GenPart_PDG!=23]")
    df[p] = df[p].Define("zDaugh_PDG","Take(GenPart_PDG,zDaugh_idx)")

#    df[p].Display({"GenPart_PDG","GenPart_momPDG"},100).Print()
#    df[p].Display({"higgsDaugh_PDG"}).Print()

    print("CutFlow for process", sampleName[i])
    allCutsReport[p]=df[p].Report()
    allCutsReport[p].Print()	


hMass={}
hRecoil={}
hDiTauVisMass={}
hDiTauCollMass={}
hLeadJetPt={}
hSecondJetPt={}
hLeadJetTheta={}
hSecondJetTheta={}
hLeadJetPhi={}
hSecondJetPhi={}
hZPt={}
hDiTauPt={}
hRecoilPt={}
hZRapidity={}
hRecoilRapidity={}
hZTheta={}
hRecoilTheta={}
hNMuons={}
hNElectrons={}
hNTaus={}
hCosRecoilTheta={}
hCosMissEtTheta={}
hLeadTauPt={}
hSecondTauPt={}
hLeadTauTheta={}
hSecondTauTheta={}
hLeadTauPhi={}
hSecondTauPhi={}
hLeadTauMass={}
hSecondTauMass={}
hLeadTauType={}
hSecondTauType={}
hLeadTauIso={}
hSecondTauIso={}
hDiTauCollMass_Fine={}
hNTauFromJets={}
hNJets={}
hNGoodJets={}

hDPhiTaus={}
hDPhiJets={}
hDThetaTaus={}
hDThetaJets={}

hJetPt={}
hJetTauID={}
hJetNConst={}
hJetNChargedHad={}
hJetNNeutralHad={}
hJetNPhotons={}

fileCutFlow={}
fileEff={}

hTauTauEvent={}
hHiggsDaugh_PDG={}
hZDaugh_PDG={}


# Rellenamos histogramas: 

for i in range(nProcesses):
   p = processes[i]

   hTauTauEvent[p] = df[p].Histo1D(("IsTauTauGen_{}".format(p), "Check Gen TauTau; Tau_h Tau_h Check; N_{Events}",2,0,1),"TauTauEvent")
   hHiggsDaugh_PDG[p] = df[p].Histo1D(("HiggsGenDaughter_{}".format(p), "Check Higgs Decay; Higgs Daughter; N_{Events}",51,-25,25),"higgsDaugh_PDG")
   hZDaugh_PDG[p] = df[p].Histo1D(("ZGenDaughter_{}".format(p), "Check Z Decay; Z Daughter; N_{Events}",51,-25,25),"zDaugh_PDG")


   hMass[p] = df[p].Histo1D(("Z_mass_{}".format(p), "DiJet mass;m_{jj} (GeV);N_{Events}",40,80, 120), "DiJet_mass")
   hRecoil[p] = df[p].Histo1D(("Recoil_mass_{}".format(p), "Z leptonic recoil mass; m_{recoil} (GeV);N_{Events}",60,120, 180), "recoil")

   hDiTauVisMass[p] = df[p].Histo1D(("DiTau_vis_mass_{}".format(p), "DiTau Vis mass;m_{#tau#tau} (GeV);N_{Events}",150,20,170), "DiTau_vis_mass")
   hDiTauCollMass[p] = df[p].Histo1D(("DiTau_coll_mass_{}".format(p), "DiTau Coll mass;collinear m_{#tau#tau} (GeV);N_{Events}",140,60,200), "DiTau_coll_mass")

   hDiTauCollMass_Fine[p] = df[p].Histo1D(("DiTau_coll_mass_Fine_{}".format(p), "DiTau Coll mass;coll m_{#tau#tau} (GeV);N_{Events}",100,100,150), "DiTau_coll_mass")

   hZPt[p] = df[p].Histo1D(("Z_Pt_{}".format(p), "DiJet Pt; Z Pt (GeV);N_{Events}",60,0,120), "DiJet_Pt")
   hRecoilPt[p] = df[p].Define("recoilPt","(p4total-DiJet_p4).Pt()").Histo1D(("Recoil_Pt_{}".format(p), "Recoil Pt; Z recoil Pt (GeV);N_{Events}",60,0,120), "recoilPt")
   hDiTauPt[p] = df[p].Histo1D(("DiTau_Pt_{}".format(p), "DiTau Pt; DiTau Pt (GeV);N_{Events}",60,0,120), "DiTau_Pt")

   hZTheta[p] = df[p].Define("Z_Theta","DiJet_p4.Theta()").Histo1D(("Z_theta_{}".format(p), "DiJet Theta; Z #theta;N_{Events}",50,0,3.16), "Z_Theta")

   hRecoilTheta[p] = df[p].Define("recoil_Theta","(p4total-DiJet_p4).Theta()").Histo1D(("Recoil_theta_{}".format(p), "Recoil Theta; Recoil #theta;N_{Events}",50,0,3.16), "recoil_Theta")

   hDPhiTaus[p]  =  df[p].Histo1D( ("dPhiTaus_{}".format(p),"DPhiTaus;#Delta #phi Taus;N_{Events}",100,0,6.28),"DPhiTaus")
   hDPhiJets[p]  =  df[p].Histo1D( ("dPhiJets_{}".format(p),"DPhiJets;#Delta #phi Jets;N_{Events}",100,0,6.28),"DPhiJets")
   hDThetaTaus[p]  = df[p].Histo1D(("dThetaTaus_{}".format(p),"DThetaTaus;#Delta #Theta Taus;N_{Events}",100,0,6.28),"DThetaTaus")
   hDThetaJets[p]  = df[p].Histo1D(("dThetaJets_{}".format(p),"DThetaJets;#Delta #Theta Jets;N_{Events}",100,0,6.28),"DThetaJets")



   hLeadJetPt[p] = df[p].Define("LeadJet_Pt","Jets_p4[0].Pt()").Histo1D(("LeadJet_Pt_{}".format(p), "Lead Jet Pt; Lead jet Pt (GeV);N_{Events}",60,0,120), "LeadJet_Pt")
   hLeadJetTheta[p] = df[p].Define("LeadJet_Theta","Jets_p4[0].Theta()").Histo1D(("LeadJet_Theta_{}".format(p), "Lead Jet #theta; Lead jet #theta;N_{Events}",50,0,3.16), "LeadJet_Theta")
   hLeadJetPhi[p] = df[p].Define("LeadJet_Phi","Jets_p4[0].Phi()").Histo1D(("LeadJet_Phi_{}".format(p), "Lead Jet #phi; Lead jet #phi;N_{Events}",50,-3.16,3.16), "LeadJet_Phi")

   hSecondJetPt[p] = df[p].Define("SecondJet_Pt","Jets_p4[1].Pt()").Histo1D(("SecondJet_Pt_{}".format(p), "Second Jet Pt; Second jet Pt (GeV);N_{Events}",60,0,120), "SecondJet_Pt")
   hSecondJetTheta[p] = df[p].Define("SecondJet_Theta","Jets_p4[1].Theta()").Histo1D(("SecondJet_Theta_{}".format(p), "Second Jet #theta; Second jet #theta;N_{Events}",50,0,3.16), "SecondJet_Theta")
   hSecondJetPhi[p] = df[p].Define("SecondJet_Phi","Jets_p4[1].Phi()").Histo1D(("SecondJet_Phi_{}".format(p), "Second Jet #phi; Second jet #phi;N_{Events}",50,-3.16,3.16), "SecondJet_Phi")

   hLeadTauPt[p] = df[p].Define("LeadTau_Pt","TauLead_p4.Pt()").Histo1D(("LeadTau_Pt_{}".format(p), "Lead Tau Pt; Lead #tau Pt (GeV);N_{Events}",60,0,120), "LeadTau_Pt")
   hLeadTauTheta[p] = df[p].Define("LeadTau_Theta","TauLead_p4.Theta()").Histo1D(("LeadTau_Theta_{}".format(p), "Lead Tau #theta; Lead #tau #theta;N_{Events}",50,0,3.16), "LeadTau_Theta")
   hLeadTauPhi[p] = df[p].Define("LeadTau_Phi","TauLead_p4.Phi()").Histo1D(("LeadTau_Phi_{}".format(p), "Lead Tau #phi; Lead #tau #phi;N_{Events}",50,-3.16,3.16), "LeadTau_Phi")

   hSecondTauPt[p] = df[p].Define("SecondTau_Pt","TauSecond_p4.Pt()").Histo1D(("SecondTau_Pt_{}".format(p), "Second Tau Pt; Second #tau Pt (GeV);N_{Events}",60,0,120), "SecondTau_Pt")
   hSecondTauTheta[p] = df[p].Define("SecondTau_Theta","TauSecond_p4.Theta()").Histo1D(("SecondTau_Theta_{}".format(p), "Second Tau #theta; Second #tau #theta;N_{Events}",50,0,3.16), "SecondTau_Theta")
   hSecondTauPhi[p] = df[p].Define("SecondTau_Phi","TauSecond_p4.Phi()").Histo1D(("SecondTau_Phi_{}".format(p), "Second Tau #phi; Second #tau #phi;N_{Events}",50,-3.16,3.16), "SecondTau_Phi")

   hLeadTauType[p] = df[p].Histo1D(("LeadTau_Type_{}".format(p), "Lead Tau Type; Lead #tau Type;N_{Events}",15,0, 15), "TauLead_type")
   hSecondTauType[p] = df[p].Histo1D(("SecondTau_Type_{}".format(p), "Second Tau Type; Second #tau Type;N_{Events}",15,0, 15), "TauSecond_type")
   hLeadTauMass[p] = df[p].Histo1D(("LeadTau_Mass_{}".format(p), "Lead Tau Mass; Lead #tau Mass;N_{Events}",100,0, 2.5), "TauLead_mass")
   hSecondTauMass[p] = df[p].Histo1D(("SecondTau_Mass_{}".format(p), "Second Tau Mass; Second #tau Mass;N_{Events}",100,0, 2.5), "TauSecond_mass")
#   hLeadTauIso[p] = df[p].Histo1D(("LeadTau_Iso_{}".format(p), "Lead Tau Iso; Lead #tau Iso;N_{Events}",100,0, 2), "TauLead_iso")
#   hSecondTauIso[p] = df[p].Histo1D(("SecondTau_Iso_{}".format(p), "Second Tau Iso; Second #tau Iso;N_{Events}",100,0, 2), "TauSecond_iso")

   hZRapidity[p] = df[p].Define("Z_y","DiJet_p4.Rapidity()").Histo1D(("Z_y_{}".format(p), "DiJet Rapidity; Z y;N_{Events}",50,-2,2), "Z_y")
   hRecoilRapidity[p] = df[p].Define("recoil_y","(p4total-DiJet_p4).Rapidity()").Histo1D(("Recoil_y_{}".format(p),"Recoil Rapidity; Recoil y;N_{Events}",50,-2,2), "recoil_y")

   hNMuons[p] = df[p].Histo1D(("NMuons_{}".format(p), "NMuons;N_{#mu};N_{Events}",5,0, 5), "NMuon")
   hNElectrons[p] = df[p].Histo1D(("NElectrons_{}".format(p), "NElectrons;N_{e};N_{Events}",5,0, 5),"NElectron")
   hNTaus[p] = df[p].Histo1D(("NTaus_{}".format(p), "NTaus;N_{#tau};N_{Events}",5,0, 5), "NTau")

   hCosRecoilTheta[p] = df[p].Define("cos_recoil_Theta","abs(cos((p4total-DiJet_p4).Theta()))").Histo1D(("cos_Recoil_theta_{}".format(p), "cos (Recoil Theta); cos(#theta_{recoil});N_{Events}",100,0,1), "cos_recoil_Theta")

   hCosMissEtTheta[p] = df[p].Histo1D( ("cos_missing_theta_{}".format(p),"cos (Missing Momentum Theta); cos (#theta_{Miss}); N_{Events}",100,0,1), "Missing_costheta")

   hNTauFromJets[p] = df[p].Histo1D(("NTauFromJets_{}".format(p), "NTauFromJets;N_{#tau};N_{Events}",10,0, 10), "NTauFromJet")

   hNGoodJets[p] = df[p].Histo1D(("NGoodJets_{}".format(p), "NGoodJets;N_{GoodJets};N_{Events}",10,0, 10), "NGoodJets")

   hNJets[p] = df[p].Histo1D(("NJets_{}".format(p), "NJets;N_{jets};N_{Events}",10,0, 10), "NJetReRun")
   hJetPt[p] = df[p].Histo1D(("JetPt_{}".format(p), "JetPt; JetPt (GeV);N_{Events}",110,10, 120), "JetReRun_pt")
   hJetNConst[p] = df[p].Histo1D(("JetNConst_{}".format(p), "JetNConst; JetNConst;N_{Events}",10,0,10), "JetReRun_nconst")
   hJetTauID[p] = df[p].Histo1D(("JetTauID_{}".format(p), "JetTauID; JetTauID;N_{Events}",31,-15, 15), "JetReRun_tauID")
   hJetNChargedHad[p] = df[p].Histo1D(("JetNChargedHad_{}".format(p), "JetNChargedHad; JetNChargedHad;N_{Events}",10,0,10), "JetReRun_nchargedhad")
   hJetNNeutralHad[p] = df[p].Histo1D(("JetNNeutralHad_{}".format(p), "JetNNeutralHad; JetNNeutralHad;N_{Events}",10,0,10), "JetReRun_nneutralhad")
   hJetNPhotons[p] = df[p].Histo1D(("JetNPhotons_{}".format(p), "JetNPhotons; JetNPhotons;N_{Events}",10,0,10), "JetReRun_nphoton")


# Salva los histogramas en un archivo root para pintarlos despues 

out = ROOT.TFile("histos.root","RECREATE")
out.cd()

for p in processes:
   hMass[p].Write()
   hRecoil[p].Write()
   hDiTauVisMass[p].Write()
   hDiTauCollMass[p].Write()
   hLeadJetPt[p].Write()
   hSecondJetPt[p].Write()
   hLeadJetTheta[p].Write()
   hSecondJetTheta[p].Write()
   hLeadJetPhi[p].Write()
   hSecondJetPhi[p].Write()
   hLeadTauPt[p].Write()
   hSecondTauPt[p].Write()
   hLeadTauTheta[p].Write()
   hSecondTauTheta[p].Write()
   hLeadTauPhi[p].Write()
   hSecondTauPhi[p].Write()
   hZPt[p].Write()
   hDiTauPt[p].Write()
   hRecoilPt[p].Write()
   hZRapidity[p].Write()
   hRecoilRapidity[p].Write()
   hZTheta[p].Write()
   hRecoilTheta[p].Write()
   hNMuons[p].Write()
   hNTaus[p].Write()
   hNElectrons[p].Write()
   hCosRecoilTheta[p].Write() 
   hCosMissEtTheta[p].Write()
   hLeadTauType[p].Write()
   hSecondTauType[p].Write()
#   hLeadTauIso[p].Write()
#   hSecondTauIso[p].Write()
   hLeadTauMass[p].Write()
   hSecondTauMass[p].Write()    
   hDiTauCollMass_Fine[p].Write()
   hJetTauID[p].Write()
   hNJets[p].Write()
   hNGoodJets[p].Write()
   hNTauFromJets[p].Write()
   hJetPt[p].Write()
   hJetNConst[p].Write()
   hJetNChargedHad[p].Write()
   hJetNPhotons[p].Write()
   hJetNNeutralHad[p].Write()
   hDThetaTaus[p].Write()
   hDThetaJets[p].Write()
   hDPhiTaus[p].Write()
   hDPhiJets[p].Write()

   hTauTauEvent[p].Write()
   hHiggsDaugh_PDG[p].Write()
   hZDaugh_PDG[p].Write()

   fileCutFlow[p]=ROOT.TH1D("CutFlow_"+p,"",15,0,15)
   fileEff[p]=ROOT.TH1D("Eff_"+p,"",15,0,15)

   events={}
   binI=0
   for cutInfo in allCutsReport[p]:
      events[cutInfo.GetName()]=cutInfo.GetAll()
      print (cutInfo.GetName(),cutInfo.GetAll(),cutInfo.GetPass(),cutInfo.GetEff())
      fileCutFlow[p].SetBinContent(binI,cutInfo.GetAll())
      fileEff[p].SetBinContent(binI,cutInfo.GetEff())
      binI+=1

   fileCutFlow[p].Write()
   fileEff[p].Write()
 

