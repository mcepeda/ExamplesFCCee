import ROOT
import tdrstyle

# CMS Style (estilo para los plots)
tdrstyle.setTDRStyle()

# Carga el archivo con los histogramas

file=ROOT.TFile("histos.root","READONLY")

#histos=["Z_mass","Recoil_mass","LeadMuon_Pt","LeadMuon_Theta","LeadMuon_Phi",\
#	"SecondMuon_Pt","SecondMuon_Theta","SecondMuon_Phi","Z_theta","Recoil_theta",\
#	"Z_Pt","Recoil_Pt","Z_y","Recoil_y","NMuons","NElectrons","cos_Recoil_theta",\
#	"cos_missing_theta"] # para pintar todos a la vez

histos=["Recoil_mass","DiTau_coll_mass","CutFlow"] # para ir rapido, solo un plot

#histos=["Z_mass","Recoil_mass","DiTau_vis_mass","DiTau_coll_mass","LeadJet_Pt","SecondJet_Pt","LeadJet_Theta","SecondJet_Theta",
#	"LeadJet_Phi","SecondJet_Phi","LeadTau_Pt","SecondTau_Pt","LeadTau_Theta","SecondTau_Theta","LeadTau_Phi","SecondTau_Phi",
#        "Z_Pt","DiTau_Pt","Recoil_Pt","Z_y","Recoil_y","Z_theta","Recoil_theta","NMuons","NTaus","NElectrons","cos_missing_theta",
#	"LeadTau_Type","LeadTau_Mass","SecondTau_Type","SecondTau_Mass",
#	"NJets","JetTauID","JetPt","JetNChargedHad","JetNConst","JetNPhotons","JetNNeutralHad","NTauFromJets","CutFlow","dPhiTaus","dPhiJets","dThetaTaus","dThetaJets","IsTauTauGen" ]

# las muestras : 

sampleName=["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240",
"wzp6_ee_uuHorddH_Htautau_ecm240",'wzp6_ee_bbH_Htautau_ecm240', 'wzp6_ee_ssH_Htautau_ecm240', 'wzp6_ee_ccH_Htautau_ecm240',
"wzp6_ee_qqH_ecm240"]
legendName=["WW","ZZ","uuH/ddH_Htautau","bbH_Htautau","ssH_Htautau","ccH_Htautau","qqH"]
nProcesses=len(sampleName)

# Los factores de normalizacion de cada proceso vienen dados por la seccion eficaz:
# eeHZ: 0.201868 pb
# eeZZ: 1.35899  pb, total win32 is 56,162,093
# eeWW: 16.4385  pb, total win32 is 58,228,784
#mumuH_Hbb - 0.00394, 300k 
#mumuH_Htautau - 0.0004243    , 400k
#tautauH_Htautau - 0.0004235, 400k 
#qqH_Htautau - 0.003346, 200k  -- careful this one is only u/d !
#bbH_Htautau - 0.00188, 400k
#ccH_Htautau - 0.001464,400k
#ssh_Htautau - 0.001879, 400k 
#qqH - 5,400,000, 0.13635
#mumuH - 0.0067643, 1,200,000
#tautauH - 0.0067518, 1,200,000
xsection=[16.4385,1.35899,0.003346,0.00188,0.001879, 0.001464,0.13635]
# Comentario a mi misma: para que esto fuese mas elegante, seria mejor haber hecho un mapa
# o un diccionario (muestra-> seccion eficaz) 
# Referencia:
# http://fcc-physics-events.web.cern.ch/fcc-physics-events/FCCee/winter2023/Delphesevents_IDEA.php

# Luminosidad acumulada: 5 ab-1 = 5000 fb-1 = 500000 pb-1
luminosity = 5e6

# Numero de sucesos MC producidos originamente.
#totalNumberOfEvents=[10000000,1000000,1200000,300000,400000]#1200000]
totalNumberOfEvents=[10000000,1000000,1200000,1200000,5400000,200000,2000000]
#qqH_Htautau 200000
for a in range(0,nProcesses):
     cutflowName="CutFlow_"+sampleName[a]
     cfHisto=file.Get(cutflowName)      
     totalNumberOfEvents[a]=cfHisto.GetBinContent(0)
     print (sampleName[a],xsection[a],totalNumberOfEvents[a])

# La normalizacion de los histogramas seguira :
#  N_Sucesos = xsection * luminosidad
#  Peso      = xsection * luminosidad / SucesosGenerados en total 

# Color de los histogramas
color=[ROOT.kBlue,ROOT.kGreen+2,ROOT.kRed,ROOT.kPink+2,ROOT.kViolet,ROOT.kMagenta,ROOT.kBlack]

# Definimos una funcion para dar estilo y normalizar los histogramas  
def StyleHisto(sample,variab,label,histColor,xsec,totalEvents,suffix=""):
     histo =file.Get(variab+"_"+sample)		
     histo.SetName("h"+variab+"_"+label+suffix)
     # Estilo:
     histo.SetLineColor(histColor)
     if "ccH" not in label  and "ssH" not in label and "uuH" not in label and "bbH" not in label and "qqH" not in label:
     	histo.SetFillColor(histColor)
     histo.SetLineWidth(3)

     # Aqui se aplica la normalizacion: 
     histo.Scale(xsec*luminosity/totalEvents)

     return histo

# Loop sobre los histogramas guardados en 'histos.root'

for histoName in histos:
   # Definimos un 'stack' de histogramas apilados:
   h={}
   hStack    = ROOT.THStack()
   
   print ("YIELDS (from %s)" %histoName)
   for i in range(nProcesses):
         h[ legendName[i] ] = StyleHisto(sampleName[i],histoName,legendName[i],color[i],xsection[i],totalNumberOfEvents[i])
         print ("... %s %2d" %(sampleName[i],h[ legendName[i] ].Integral() ) )
#         if legendName[i]!="qqH_Htautau" or legendName[i]!="qqH_Hbb": # para pintar la señal por separado
         if "H_Htautau"  not in legendName[i]: 
           hStack.Add( h[ legendName[i] ] )


   # Ahora pintamos los histogramas:
   c1 = ROOT.TCanvas("c1","HZ analysis canvas", 650,600)
   ROOT.gStyle.SetOptTitle(True)
  
   if histoName=="CutFlow":
      c1.SetLogy()
 
   c1.SetTicks(1,1)
   c1.SetLeftMargin(0.15)
   c1.SetRightMargin(0.08)
   c1.SetTopMargin(0.06)
  
   hStack.Draw("hist")
   for i in range(nProcesses):
     if "H_Htautau"  in legendName[i]: 
      h[legendName[i]].Draw("hist,sames")	# Opcion para pintar la segnal por separado
 
   ylabel="events"
   xlabel=h[legendName[0]].GetXaxis().GetTitle() #"DiMuon Mass [GeV]"
   
   hStack.GetXaxis().SetTitle(xlabel)
   hStack.GetYaxis().SetTitle(ylabel)
   
   hStack.GetYaxis().SetTitleOffset(1.4)
   hStack.GetXaxis().SetTitleOffset(1.2)

   maxY=hStack.GetMaximum()
   if h["qqH"].GetMaximum()> maxY:
      maxY=h["qqH"].GetMaximum() 
   hStack.SetMaximum(maxY*1.2)
 
   # Leyenda
   leg = ROOT.TLegend(0.60,0.60,0.9,0.93)
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   for i in range(nProcesses):
   	entry=leg.AddEntry(h[ legendName[i] ],legendName[i],"lf")
   leg.Draw()
   
   # Etiquetas explicando las condiciones del proceso
   text = "#sqrt{{s}} = 240 GeV,   L = {:.0f} ab^{{-1}}".format(luminosity/1e6)
   channel = 'e^{+}e^{-} #rightarrow ZH #rightarrow jj + #tau^{+}_{h}#tau^{-}_{h}'
   
   Text = ROOT.TLatex()
   Text.SetNDC()
   Text.SetTextAlign(31);
   Text.SetTextSize(0.04)
   Text.DrawLatex(0.47, 0.95, channel)
   
   Text.SetTextAlign(31);
   Text.SetTextSize(0.04)
   Text.DrawLatex(0.93, 0.95, text)
   
   # Guardamos la grafica 
   c1.SaveAs("Plot_"+histoName+".png")
