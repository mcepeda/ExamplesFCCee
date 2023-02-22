#! /usr/bin/env python3
#
# 1) Setup ROOT with: 
#  source /cvmfs/fcc.cern.ch/sw/latest/setup.sh 
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
#print("ROOT: %s CPUs used out of %s available ..." % (usecpus, maxcpus))

# Archivos de entrada
# Los archivos estan en 'path'. Vamos a usar tres archivos, X_skimmed_reduced.root, donde X = eeHZ
# (sennal) +  eeZZ,  eeWW (fondos)
sampleName=["wzp6_ee_mumuH_Htautau_ecm240","wzp6_ee_mumuH_ecm240","wzp6_ee_qqH_ecm240","wzp6_ee_tautauH_ecm240","p8_ee_ZZ_ecm240","p8_ee_WW_ecm240"]
nProcesses=len(sampleName)
#path = "/afs/ciemat.es/user/a/alcaraz/public/FCCee/"
#path="files/"

# Archivos de entrada con el ultimo algoritmo de Taus que tengo (a fecha de 22/02/23)
path='/nfs/cms/cepeda/FCC/taus16Feb_eekt4Jets_idInJet/'


# Cada rootfile tiene un 'tree' llamado events que vamos a
# convertir en un dataframe y guardar en una lista de dataframes. 

#df = {}
#for p in sampleName: 
#	df[p] = ROOT.RDataFrame("events",(os.path.join(path,"{}.root".format(p))))

print ("Reading Trees")
df = {}
allCutsReport={}

for p in sampleName:
   # Las muestras de ZZ y WW son muy grandes, estan repartidas en ficheros menores 
   # (No he procesado toda la estadistica disponible de momento)
   if p=="p8_ee_ZZ_ecm240":
      names = ROOT.std.vector('string')()
      for i in range(0,10):
          names.push_back(path+p+'_{}.root'.format(i))
      df[p] = ROOT.RDataFrame("events",(names))
   elif p=="p8_ee_WW_ecm240":
      names = ROOT.std.vector('string')()
      for i in range(0,5):
          names.push_back(path+p+'_{}.root'.format(i))
      df[p] = ROOT.RDataFrame("events",(names))
   else:
    # Las muestras de Higgs van por separado
          df[p] = ROOT.RDataFrame("events",(os.path.join(path,"{}.root".format(p))))

processes = list(df.keys())


# Ejemplo de funcion en cc, que busca los indices de los mejores candidatos a una
# resonancia (Z, Higgs...)
ROOT.gInterpreter.Declare("""
   
   ROOT::VecOps::RVec<int> FindBest(ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>> P4vector, double mass){ 

        double minDistance=10000;
        int lead=-1;
        int second=-1;

        for (int i=0; i<P4vector.size(); i++){
                for (int j=i+1; j<P4vector.size(); j++){
                        ROOT::Math::PxPyPzMVector DiObj=P4vector[i]+P4vector[j];
                        if ( fabs(DiObj.M()-mass)<minDistance) {
                                        minDistance=fabs(DiObj.M()-mass);
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

    # Para calcular aceptancias en la señal
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

#    df[p] = df[p].Filter("NGenTauHad==2","Only Tau_h Tau_h")  # Para calcular la
#    aceptancia/efficiencia *solo* para sucesos HTauTau hadronico

    # Para simplificar la seleccion, de momento vamos a mirar sucesos con exactamente dos
    # muones 
    df[p] = df[p].Filter("NMuon >= 2", "Events with at least two muons")

    df[p] = df[p].Filter("Muon_pt[0]>10&&Muon_pt[1]>10","Muon Pt>10")
    df[p] = df[p].Define("Muon_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>>(Muon_px,Muon_py,Muon_pz,Muon_mass)")
    df[p] = df[p].Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge")

    # Reconstruimos el Z! Sumando los cuadrimomentos (p4)
    df[p] = df[p].Define("DiMuon_p4","Muon_p4[0]+Muon_p4[1]")

    # A estos p4 le podemos pedir directamente las variables que nos interesan, como la
    # masa:
    df[p] = df[p].Define("Dimuon_mass", "DiMuon_p4.M()")
    df[p] = df[p].Define("Dimuon_Pt", "DiMuon_p4.Pt()")

    # Filtramos en la masa del Z:
    df[p] = df[p].Filter("Dimuon_mass>86 && Dimuon_mass<96","86 < Mreco(dimuon) < 96 GeV")

    # Filtramos en el pt del Z:
#    df[p] = df[p].Filter("Dimuon_Pt>20 && Dimuon_Pt<70","20<Dimuon PT <70 GeV")

    # Reconstruimos el recoil, que corresponde a la masa del Higgs:
    df[p] = df[p].Define("p4total","ROOT::Math::PxPyPzEVector(0.,0.,0.,240.)")
    df[p] = df[p].Define("recoil","(p4total-DiMuon_p4).M()")

    # Filtramos en la masa del Recoil 
#    df[p] = df[p].Filter("recoil>120", "Recoil>120")
    df[p] = df[p].Filter("recoil>120 && recoil<140", "120<M_{H}<140 GeV")
#    df[p] = df[p].Filter("recoil<100", "Recoil<100")

    # Define missing momentum
    df[p] = df[p].Define('Missing_e', ' sqrt( Missing_px*Missing_px+Missing_py*Missing_py+Missing_pz*Missing_pz )')

    # El momento 'que falta', Missing Pt y su angulo cos(Theta_miss)
    df[p] = df[p].Define('Missing_p4', "ROOT::Math::PxPyPzEVector(Missing_px,Missing_py,Missing_pz,Missing_e)")
    df[p] = df[p].Define('Missing_costheta', 'abs(cos(Missing_p4.Theta()))')

    # Filtramos en cos(theta_miss) - este de momento no lo vamos a activar
    df[p] = df[p].Filter("Missing_costheta<0.98", "cos(theta_miss)<0.98")

    # Ahora vamos a buscar los Taus. Para hacer esto mejor podriamos comprobar el orden
    # Que funciona mejor, ordenar por pt? Por aislamiento? Cuantos Taus hay? 
    df[p] = df[p].Filter("NTauFromJet == 2", "Events with exactly two taus")
    df[p] = df[p].Define("Taus_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>  > > ( TauFromJet_px,TauFromJet_py,TauFromJet_pz,TauFromJet_mass)")
#    df[p] = df[p].Define("bestHiggs","FindBest(Taus_p4,125)")  # Me estoy quedando con
#    solo dos taus, no hace falta esto
    df[p] = df[p].Define("Lead","0")#bestHiggs[0]")
    df[p] = df[p].Define("Second","1")#bestHiggs[1]")
    df[p] = df[p].Define("TauLead_p4","Taus_p4[Lead]")
    df[p] = df[p].Define("TauSecond_p4","Taus_p4[Second]")

    # Filtro en el pt de los taus
    df[p] = df[p].Filter("TauFromJet_pt[Lead]>10 && TauFromJet_pt[Second]>10", "Tau Pt>10 GeV")
    # Carga opuesta
    df[p] = df[p].Filter("TauFromJet_charge[Lead] != TauFromJet_charge[Second]", "Taus with opposite charge")
 
    # Identificacion de taus: el algoritmo ahora mismo parte de un jet y busca en sus
    # constituyentes patrones compatibles con h+ ó h+h-h+ (y sus conjuados en carga), más
    # fotones. Los posibles tipos son:
    # -13, -11 -> muones, electrones (no deberian aparecer en esta muestra)
    # -1 -> jets con mas trazas (por ejemplo
    # 0 -> h+ o h-
    # 1 -> h+ + 1 fotón cercano (dTheta=0.15 con respecto a la traza)
    # 2, 3,4,...  -> h+ + 2,3,4... fotones cercanos
    # 10 -> h+h-h+ o h-h+h-
    # 11 -> hhh + 1 fotón cercano (dTheta=0.15 con respecto a la traza)
    # 12, 13, 14,...  -> hhh + 2,3,4... fotones cercanos
    # Con este corte, nos quedamos solo con los buenos candidatos a tau:
    df[p] = df[p].Filter("TauFromJet_type[Lead]>=0 && TauFromJet_type[Second]>=0",   "Identified Taus")

    # Reconstruimos el DiTau, Sumando los cuadrimomentos (p4)
    df[p] = df[p].Define("DiTau_p4","Taus_p4[Lead]+Taus_p4[Second]")

    df[p] = df[p].Define("TauLead_type","TauFromJet_type[Lead]")
    df[p] = df[p].Define("TauSecond_type","TauFromJet_type[Second]")
    df[p] = df[p].Define("TauLead_mass","TauFromJet_mass[Lead]")
    df[p] = df[p].Define("TauSecond_mass","TauFromJet_mass[Second]")

    df[p] = df[p].Define("DiTau_vis_mass", "DiTau_p4.M()")
    df[p] = df[p].Define("DiTau_Pt", "DiTau_p4.Pt()")

    # Minimo de masa para el DiTau
    df[p] = df[p].Filter("DiTau_vis_mass>20", "DiTauMass > 20 GeV")

    #Aplicamos collinear approximation
    df[p] = df[p].Define("p12","TauLead_p4.Py()*TauSecond_p4.Px()-TauLead_p4.Px()*TauSecond_p4.Py()")
    df[p] = df[p].Define("r0","abs((Missing_p4.Py()*TauLead_p4.Px()-Missing_p4.Px()*TauLead_p4.Py())/p12)")
    df[p] = df[p].Define("f0","1/(1+r0)")
    df[p] = df[p].Define("r1","abs((Missing_p4.Py()*TauSecond_p4.Px()-Missing_p4.Px()*TauSecond_p4.Py())/p12)")
    df[p] = df[p].Define("f1","1/(1+r1)")
    df[p] = df[p].Define("DiTau_coll_mass","DiTau_vis_mass/sqrt(f0*f1)")

    print("CutFlow for process", sampleName[i])
    allCutsReport[p]=df[p].Report()
    allCutsReport[p].Print()	


hMass={}
hRecoil={}
hDiTauVisMass={}
hDiTauCollMass={}
hLeadMuonPt={}
hSecondMuonPt={}
hLeadMuonTheta={}
hSecondMuonTheta={}
hLeadMuonPhi={}
hSecondMuonPhi={}
hZPt={}
hDiTauPt={}
hRecoilPt={}
hZRapidity={}
hRecoilRapidity={}
hZTheta={}
hRecoilTheta={}
hNMuons={}
hNElectrons={}
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
hJetPt={}
hJetTauID={}
hJetNConst={}
hJetNChargedHad={}
hJetNNeutralHad={}
hJetNPhotons={}

fileCutFlow={}
fileEff={}

# Rellenamos histogramas: 

for i in range(nProcesses):
   p = processes[i]
   hMass[p] = df[p].Histo1D(("Z_mass_{}".format(p), "Dimuon mass;m_{#mu#mu} (GeV);N_{Events}",40,80, 100), "Dimuon_mass")
   hRecoil[p] = df[p].Histo1D(("Recoil_mass_{}".format(p), "Z leptonic recoil mass; m_{recoil} (GeV);N_{Events}",80,120, 140), "recoil")

   hDiTauVisMass[p] = df[p].Histo1D(("DiTau_vis_mass_{}".format(p), "DiTau Vis mass;m_{#tau#tau} (GeV);N_{Events}",150,20,170), "DiTau_vis_mass")
   hDiTauCollMass[p] = df[p].Histo1D(("DiTau_coll_mass_{}".format(p), "DiTau Coll mass;collinear m_{#tau#tau} (GeV);N_{Events}",140,60,200), "DiTau_coll_mass")

   hDiTauCollMass_Fine[p] = df[p].Histo1D(("DiTau_coll_mass_Fine_{}".format(p), "DiTau Coll mass;coll m_{#tau#tau} (GeV);N_{Events}",100,100,150), "DiTau_coll_mass")

   hZPt[p] = df[p].Histo1D(("Z_Pt_{}".format(p), "Dimuon Pt; Z Pt (GeV);N_{Events}",60,0,120), "Dimuon_Pt")
   hRecoilPt[p] = df[p].Define("recoilPt","(p4total-DiMuon_p4).Pt()").Histo1D(("Recoil_Pt_{}".format(p), "Recoil Pt; Z recoil Pt (GeV);N_{Events}",60,0,120), "recoilPt")
   hDiTauPt[p] = df[p].Histo1D(("DiTau_Pt_{}".format(p), "DiTau Pt; DiTau Pt (GeV);N_{Events}",60,0,120), "DiTau_Pt")

   hZTheta[p] = df[p].Define("Z_Theta","DiMuon_p4.Theta()").Histo1D(("Z_theta_{}".format(p), "Dimuon Theta; Z #theta;N_{Events}",50,0,3.16), "Z_Theta")

   hRecoilTheta[p] = df[p].Define("recoil_Theta","(p4total-DiMuon_p4).Theta()").Histo1D(("Recoil_theta_{}".format(p), "Recoil Theta; Recoil #theta;N_{Events}",50,0,3.16), "recoil_Theta")

   hLeadMuonPt[p] = df[p].Define("LeadMuon_Pt","Muon_p4[0].Pt()").Histo1D(("LeadMuon_Pt_{}".format(p), "Lead Muon Pt; Lead #mu Pt (GeV);N_{Events}",60,0,120), "LeadMuon_Pt")
   hLeadMuonTheta[p] = df[p].Define("LeadMuon_Theta","Muon_p4[0].Theta()").Histo1D(("LeadMuon_Theta_{}".format(p), "Lead Muon #theta; Lead #mu #theta;N_{Events}",50,0,3.16), "LeadMuon_Theta")
   hLeadMuonPhi[p] = df[p].Define("LeadMuon_Phi","Muon_p4[0].Phi()").Histo1D(("LeadMuon_Phi_{}".format(p), "Lead Muon #phi; Lead #mu #phi;N_{Events}",50,-3.16,3.16), "LeadMuon_Phi")

   hSecondMuonPt[p] = df[p].Define("SecondMuon_Pt","Muon_p4[1].Pt()").Histo1D(("SecondMuon_Pt_{}".format(p), "Second Muon Pt; Second #mu Pt (GeV);N_{Events}",60,0,120), "SecondMuon_Pt")
   hSecondMuonTheta[p] = df[p].Define("SecondMuon_Theta","Muon_p4[1].Theta()").Histo1D(("SecondMuon_Theta_{}".format(p), "Second Muon #theta; Second #mu #theta;N_{Events}",50,0,3.16), "SecondMuon_Theta")
   hSecondMuonPhi[p] = df[p].Define("SecondMuon_Phi","Muon_p4[1].Phi()").Histo1D(("SecondMuon_Phi_{}".format(p), "Second Muon #phi; Second #mu #phi;N_{Events}",50,-3.16,3.16), "SecondMuon_Phi")

   hLeadTauPt[p] = df[p].Define("LeadTau_Pt","TauLead_p4.Pt()").Histo1D(("LeadTau_Pt_{}".format(p), "Lead Tau Pt; Lead #tau Pt (GeV);N_{Events}",60,0,120), "LeadTau_Pt")
   hLeadTauTheta[p] = df[p].Define("LeadTau_Theta","TauLead_p4.Theta()").Histo1D(("LeadTau_Theta_{}".format(p), "Lead Tau #theta; Lead #tau #theta;N_{Events}",50,0,3.16), "LeadTau_Theta")
   hLeadTauPhi[p] = df[p].Define("LeadTau_Phi","TauLead_p4.Phi()").Histo1D(("LeadTau_Phi_{}".format(p), "Lead Tau #phi; Lead #tau #phi;N_{Events}",50,-3.16,3.16), "LeadTau_Phi")

   hSecondTauPt[p] = df[p].Define("SecondTau_Pt","TauSecond_p4.Pt()").Histo1D(("SecondTau_Pt_{}".format(p), "Second Tau Pt; Second #tau Pt (GeV);N_{Events}",60,0,120), "SecondTau_Pt")
   hSecondTauTheta[p] = df[p].Define("SecondTau_Theta","TauSecond_p4.Theta()").Histo1D(("SecondTau_Theta_{}".format(p), "Second Tau #theta; Second #tau #theta;N_{Events}",50,0,3.16), "SecondTau_Theta")
   hSecondTauPhi[p] = df[p].Define("SecondTau_Phi","TauSecond_p4.Phi()").Histo1D(("SecondTau_Phi_{}".format(p), "Second Tau #phi; Second #tau #phi;N_{Events}",50,-3.16,3.16), "SecondTau_Phi")

   hLeadTauType[p] = df[p].Histo1D(("LeadTau_Type_{}".format(p), "Lead Tau Type; Lead #tau Type;N_{Events}",15,0, 15), "TauLead_type")
   hSecondTauType[p] = df[p].Histo1D(("SecondTau_Type_{}".format(p), "Second Tau Type; Second #tau Type;N_{Events}",15,0, 15), "TauSecond_type")
   hLeadTauMass[p] = df[p].Histo1D(("LeadTau_Mass_{}".format(p), "Lead Tau Mass; Lead #tau Mass;N_{Events}",100,0, 2), "TauLead_mass")
   hSecondTauMass[p] = df[p].Histo1D(("SecondTau_Mass_{}".format(p), "Second Tau Mass; Second #tau Mass;N_{Events}",100,0, 2), "TauSecond_mass")
#   hLeadTauIso[p] = df[p].Histo1D(("LeadTau_Iso_{}".format(p), "Lead Tau Iso; Lead #tau Iso;N_{Events}",100,0, 2), "TauLead_iso")
#   hSecondTauIso[p] = df[p].Histo1D(("SecondTau_Iso_{}".format(p), "Second Tau Iso; Second #tau Iso;N_{Events}",100,0, 2), "TauSecond_iso")

   hZRapidity[p] = df[p].Define("Z_y","DiMuon_p4.Rapidity()").Histo1D(("Z_y_{}".format(p), "Dimuon Rapidity; Z y;N_{Events}",50,-2,2), "Z_y")
   hRecoilRapidity[p] = df[p].Define("recoil_y","(p4total-DiMuon_p4).Rapidity()").Histo1D(("Recoil_y_{}".format(p),"Recoil Rapidity; Recoil y;N_{Events}",50,-2,2), "recoil_y")

   hNMuons[p] = df[p].Histo1D(("NMuons_{}".format(p), "NMuons;N_{#mu};N_{Events}",5,0, 5), "NMuon")
   hNElectrons[p] = df[p].Histo1D(("NElectrons_{}".format(p), "NElectrons;N_{e};N_{Events}",5,0, 5),"NElectron")

   hCosRecoilTheta[p] = df[p].Define("cos_recoil_Theta","abs(cos((p4total-DiMuon_p4).Theta()))").Histo1D(("cos_Recoil_theta_{}".format(p), "cos (Recoil Theta); cos(#theta_{recoil});N_{Events}",100,0,1), "cos_recoil_Theta")

   hCosMissEtTheta[p] = df[p].Histo1D( ("cos_missing_theta_{}".format(p),"cos (Missing Momentum Theta); cos (#theta_{Miss}); N_{Events}",100,0,1), "Missing_costheta")

   hNTauFromJets[p] = df[p].Histo1D(("NTauFromJets_{}".format(p), "NTauFromJets;N_{#tau};N_{Events}",10,0, 10), "NTauFromJet")

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
   hLeadMuonPt[p].Write()
   hSecondMuonPt[p].Write()
   hLeadMuonTheta[p].Write()
   hSecondMuonTheta[p].Write()
   hLeadMuonPhi[p].Write()
   hSecondMuonPhi[p].Write()
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
   hNElectrons[p].Write()
   hCosRecoilTheta[p].Write() 
   hCosMissEtTheta[p].Write()
   hLeadTauType[p].Write()
   hSecondTauType[p].Write()
   hLeadTauMass[p].Write()
   hSecondTauMass[p].Write()    
   hDiTauCollMass_Fine[p].Write()
   hJetTauID[p].Write()
   hNJets[p].Write()
   hNTauFromJets[p].Write()
   hJetPt[p].Write()
   hJetNConst[p].Write()
   hJetNChargedHad[p].Write()
   hJetNPhotons[p].Write()
   hJetNNeutralHad[p].Write()


   fileCutFlow[p]=ROOT.TH1D("CutFlow_"+p,"",15,0,15)
   fileEff[p]=ROOT.TH1D("Eff_"+p,"",15,0,15)

   events={}
   binI=0
 
   # Guarda en el archivo final el numero total de sucesos para normalizar:
   for cutInfo in allCutsReport[p]:
      events[cutInfo.GetName()]=cutInfo.GetAll()
      fileCutFlow[p].SetBinContent(binI,cutInfo.GetAll())
      fileEff[p].SetBinContent(binI,cutInfo.GetEff())
      binI+=1

   fileCutFlow[p].Write()
   fileEff[p].Write()
 

