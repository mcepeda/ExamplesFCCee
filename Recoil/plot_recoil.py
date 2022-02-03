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
#maxcpus = multiprocessing.cpu_count()
#ROOT.EnableImplicitMT(min(maxcpus,8))
#usecpus = max(ROOT.GetThreadPoolSize(),1)
#print("ROOT: %s CPUs used out of %s available ..." % (usecpus, maxcpus))

# Archivos de entrada
# Los archivos estan en 'path'. Vamos a usar tres archivos, X_skimmed_reduced.root, donde X = eeHZ
# (sennal) +  eeZZ,  eeWW (fondos)
sampleName=["eeHZ","eeZZ","eeWW"]
nProcesses=len(sampleName)
path = "/afs/ciemat.es/user/a/alcaraz/public/FCCee/"
 
# Cada rootfile tiene un 'tree' llamado events que vamos a
# convertir en un dataframe y guardar en una lista de dataframes. 

df = {}
for p in sampleName: 
	df[p] = ROOT.RDataFrame("events",(os.path.join(path,"{}_skimmed_reduced.root".format(p))))

processes = list(df.keys())

print(processes)

# Ahora hacemos la seleccion y definimos variables nuevas 

for i in range(nProcesses):
    p=processes[i]

    # Para simplificar la seleccion, de momento vamos a mirar sucesos con exactamente dos
    # muones 
    df[p] = df[p].Filter("NMuon >= 2", "Events with at least two muons")
    df[p] = df[p].Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge")

    # 1) Metodo 1: variable a variable
    # Vamos a definir mas variables cinematicas. 
    # Dada la geometria del detector es comodo trabajar en cilindricas 
    # en vez de cartesianas: vamos a definir el momento del muon en el plano transverso
    df[p] = df[p].Define("Muon_pt","sqrt( Muon_px*Muon_px+Muon_py*Muon_py ) ")
    # Tambien queremos la energia del muon. E^2=M^2+P^2, aqui c=1. La masa es muy pequenna,
    # podriamos descartarla 
    df[p] = df[p].Define("Muon_E","sqrt( Muon_px*Muon_px+ Muon_py*Muon_py+ Muon_pz*Muon_pz+Muon_mass*Muon_mass)")

    df[p] = df[p].Filter("Muon_pt[0]>10&&Muon_pt[1]>10","Muon Pt>10")

    # A partir de los dos muones [0] y [1], vamos a reconstruir la masa invariante del par
    # de muones.  
#    df[p] = df[p].Define("Dimuon_p_long", "sqrt( (Muon_px[0]+Muon_px[1])*(Muon_px[0]+Muon_px[1]) + (Muon_py[0]+Muon_py[1])*(Muon_py[0]+Muon_py[1]) + (Muon_pz[0]+Muon_pz[1])*(Muon_pz[0]+Muon_pz[1]))")
#    df[p] = df[p].Define("Dimuon_E_long", "Muon_E[0]+Muon_E[1]")
#    df[p] = df[p].Define("Dimuon_mass_long", "sqrt( Dimuon_E_long*Dimuon_E_long - Dimuon_p_long*Dimuon_p_long ) ")

    # 2) Metodo 2: cuadrimomentos
    # La vida es un poco mas comoda si directamente definimos cuadrivectores de Lorentz  
    df[p] = df[p].Define("Muon0_p4","ROOT::Math::PxPyPzEVector(Muon_px[0],Muon_py[0],Muon_pz[0],Muon_E[0])")
    df[p] = df[p].Define("Muon1_p4","ROOT::Math::PxPyPzEVector(Muon_px[1],Muon_py[1],Muon_pz[1],Muon_E[1])")

    df[p] = df[p].Define("Muon0_pt","Muon0_p4.Pt()")
    df[p] = df[p].Define("Muon1_pt","Muon0_p4.Pt()")
    df[p] = df[p].Define("Muon0_eta","Muon0_p4.Eta()")
    df[p] = df[p].Define("Muon1_eta","Muon0_p4.Eta()")
    df[p] = df[p].Define("Muon0_phi","Muon0_p4.Phi()")
    df[p] = df[p].Define("Muon1_phi","Muon0_p4.Phi()")

    # Podrias definir df[p] = df[p].Define("Muon0_pt","Muon0_p4.Pt()")  

    # Reconstruimos el Z! Sumando los cuadrimomentos (p4)
    df[p] = df[p].Define("DiMuon_p4","Muon0_p4+Muon1_p4")

    # A estos p4 le podemos pedir directamente las variables que nos interesan, como la
    # masa:
    df[p] = df[p].Define("Dimuon_mass", "DiMuon_p4.M()")
    df[p] = df[p].Define("Dimuon_Pt", "DiMuon_p4.Pt()")

    # Filtramos en la masa del Z:
    df[p] = df[p].Filter("Dimuon_mass>80 && Dimuon_mass<100","80 < Mreco(dimuon) < 100 GeV")

    # Reconstruimos el recoil, que corresponde a la masa del Higgs:
    df[p] = df[p].Define("p4total","ROOT::Math::PxPyPzEVector(0.,0.,0.,240.)")
    df[p] = df[p].Define("recoil","(p4total-DiMuon_p4).M()")

#    df[p].Display({"Muon_pt","Muon_charge","Dimuon_mass","recoil"},30 ).Print()

    print("CutFlow for process", sampleName[i])
    df[p].Report().Print()	


hMass={}
hRecoil={}
hLeadMuonPt={}
hSecondMuonPt={}
hLeadMuonTheta={}
hSecondMuonTheta={}
hLeadMuonPhi={}
hSecondMuonPhi={}
hZPt={}
hRecoilPt={}
hZRapidity={}
hRecoilRapidity={}
hZTheta={}
hRecoilTheta={}
hNMuons={}
hNElectrons={}

# Rellenamos histogramas: 

for i in range(nProcesses):
   p = processes[i]
   hMass[p] = df[p].Histo1D(("Z_mass_{}".format(p), "Dimuon mass;m_{#mu#mu} (GeV);N_{Events}",80,80, 100), "Dimuon_mass")
   hRecoil[p] = df[p].Histo1D(("Recoil_mass_{}".format(p), "Z leptonic recoil mass; m_{recoil} (GeV);N_{Events}",200,120, 140), "recoil")

   hZPt[p] = df[p].Histo1D(("Z_Pt_{}".format(p), "Dimuon Pt; Z Pt (GeV);N_{Events}",120,0, 120), "Dimuon_Pt")
   hRecoilPt[p] = df[p].Define("recoilPt","(p4total-DiMuon_p4).Pt()").Histo1D(("Recoil_Pt_{}".format(p), "Recoil Pt; Z recoil Pt (GeV);N_{Events}",100,0,120), "recoilPt")

   hZTheta[p] = df[p].Define("Z_Theta","DiMuon_p4.Theta()").Histo1D(("Z_theta_{}".format(p), "Dimuon Theta; Z #theta;N_{Events}",100,0,3.16), "Z_Theta")
   hRecoilTheta[p] = df[p].Define("recoil_Theta","(p4total-DiMuon_p4).Theta()").Histo1D(("Recoil_theta_{}".format(p), "Recoil Theta; Recoil #theta;N_{Events}",100,0,3.16), "recoil_Theta")

   hLeadMuonPt[p] = df[p].Define("LeadMuon_Pt","Muon0_p4.Pt()").Histo1D(("LeadMuon_Pt_{}".format(p), "Lead Muon Pt; Lead #mu Pt (GeV);N_{Events}",120,0, 120), "LeadMuon_Pt")
   hLeadMuonTheta[p] = df[p].Define("LeadMuon_Theta","Muon0_p4.Theta()").Histo1D(("LeadMuon_Theta_{}".format(p), "Lead Muon #theta; Lead #mu #theta;N_{Events}",100,0,3.16), "LeadMuon_Theta")
   hLeadMuonPhi[p] = df[p].Define("LeadMuon_Phi","Muon0_p4.Phi()").Histo1D(("LeadMuon_Phi_{}".format(p), "Lead Muon #phi; Lead #mu #phi;N_{Events}",100,-3.16,3.16), "LeadMuon_Phi")

   hSecondMuonPt[p] = df[p].Define("SecondMuon_Pt","Muon0_p4.Pt()").Histo1D(("SecondMuon_Pt_{}".format(p), "Second Muon Pt; Second #mu Pt (GeV);N_{Events}",120,0, 120), "SecondMuon_Pt")
   hSecondMuonTheta[p] = df[p].Define("SecondMuon_Theta","Muon0_p4.Theta()").Histo1D(("SecondMuon_Theta_{}".format(p), "Second Muon #theta; Second #mu #theta;N_{Events}",100,0,3.16), "SecondMuon_Theta")
   hSecondMuonPhi[p] = df[p].Define("SecondMuon_Phi","Muon0_p4.Phi()").Histo1D(("SecondMuon_Phi_{}".format(p), "Second Muon #phi; Second #mu #phi;N_{Events}",100,-3.16,3.16), "SecondMuon_Phi")

   hZRapidity[p] = df[p].Define("Z_y","DiMuon_p4.Rapidity()").Histo1D(("Z_y_{}".format(p), "Dimuon Rapidity; Z y;N_{Events}",100,-2,2), "Z_y")
   hRecoilRapidity[p] = df[p].Define("recoil_y","(p4total-DiMuon_p4).Rapidity()").Histo1D(("Recoil_y_{}".format(p),"Recoil Rapidity; Recoil y;N_{Events}",100,-2,2), "recoil_y")

   hNMuons[p] = df[p].Histo1D(("NMuons_{}".format(p), "NMuons;N_{#mu};N_{Events}",5,0, 5), "NMuon")
   hNElectrons[p] = df[p].Histo1D(("NElectrons_{}".format(p), "NElectrons;N_{e};N_{Events}",5,0, 5),"NElectron")




# Salva los histogramas en un archivo root para pintarlos despues 

out = ROOT.TFile("histos.root","RECREATE")
out.cd()

for p in processes:
   hMass[p].Write()
   hRecoil[p].Write()
   hLeadMuonPt[p].Write()
   hSecondMuonPt[p].Write()
   hLeadMuonTheta[p].Write()
   hSecondMuonTheta[p].Write()
   hLeadMuonPhi[p].Write()
   hSecondMuonPhi[p].Write()
   hZPt[p].Write()
   hRecoilPt[p].Write()
   hZRapidity[p].Write()
   hRecoilRapidity[p].Write()
   hZTheta[p].Write()
   hRecoilTheta[p].Write()
   hNMuons[p].Write()
   hNElectrons[p].Write()
 


