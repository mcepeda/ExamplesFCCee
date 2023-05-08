import ROOT
import tdrstyle

# CMS Style (estilo para los plots)
tdrstyle.setTDRStyle()

# Carga el archivo con los histogramas

file=ROOT.TFile("histos.root","READONLY")

#histos=["Z_mass","Z_Pt","Recoil_mass","DiTau_vis_mass","DiTau_coll_mass","CutFlow"] 

# muchos mas...
histos=["Z_mass","Recoil_mass","DiTau_vis_mass","DiTau_coll_mass","LeadMuon_Pt","SecondMuon_Pt","LeadMuon_Theta","SecondMuon_Theta",
	"LeadMuon_Phi","SecondMuon_Phi","LeadTau_Pt","ElectronLead_Pt","LeadTau_Theta","ElectronLead_Theta","LeadTau_Phi","ElectronLead_Phi",
        "Z_Pt","DiTau_Pt","Recoil_Pt","Z_y","Recoil_y","Z_theta","Recoil_theta","NMuons","NElectrons","cos_missing_theta",
	"LeadTau_Type","LeadTau_Mass",
	"NJets","JetTauID","JetPt","JetNChargedHad","JetNConst","JetNPhotons","JetNNeutralHad","NTauFromJets","CutFlow","dPhiTaus","dThetaTaus"]

# las muestras : 
sampleName=["p8_ee_ZZ_ecm240", "p8_ee_WW_ecm240","wzp6_ee_mumuH_Hbb_ecm240","wzp6_ee_mumuH_HWW_ecm240","wzp6_ee_mumuH_Htautau_ecm240"]
legendName=["ZZ","WW","mumuH_Hbb","mumuH_HWW","mumuH_Htautau"]
nProcesses=len(sampleName)

# Los factores de normalizacion de cada proceso vienen dados por la seccion eficaz:
# A) Muestras Spring21
# eeHZ: 0.201868 pb
# eeZZ: 1.35899  pb, total win32 is 56,162,093
# eeWW: 16.4385  pb, total win32 is 58,228,784

# B) Muestras Winter23
# p8_ee_ZZ_ecm240: 1.35899  pb
# p8_ee_WW_ecm240: 16.4385  pb
#mumuH_Hbb - 0.00394
#mumuH_Htautau - 0.0004243   
#tautauH_Htautau - 0.0004235
#uuH/ddH_Htautau - 0.003346  
# ojo esta muestra se llama originalmente qqH_Htautau: mal nombre, es solo q=u,d
#qqH -  0.13635  # esta si que es Zqq q=u,d,c,s,b
#mumuH - 0.0067643
#tautauH - 0.0067518
# Referencia:
# http://fcc-physics-events.web.cern.ch/fcc-physics-events/FCCee/winter2023/Delphesevents_IDEA.php

xsection=[1.35899,16.4385,0.00394,0.001456,0.0004243]
# Comentario a mi misma: para que esto sea mas elegante y cause menos errores tontos, seria
# muchisimo mejor hacer un mapa o un diccionario (muestra-> seccion eficaz) 

# Luminosidad acumulada: 5 ab-1 = 5000 fb-1 = 500000 pb-1
luminosity = 5e6

# Numero de sucesos MC producidos originamente.
totalNumberOfEvents=[0]*nProcesses
for a in range(0,nProcesses):
     cutflowName="CutFlow_"+sampleName[a]
     print(cutflowName)
     cfHisto=file.Get(cutflowName)      
     totalNumberOfEvents[a]=cfHisto.GetBinContent(0)
     print (sampleName[a],xsection[a],totalNumberOfEvents[a])

# La normalizacion de los histogramas seguira :
#  N_Sucesos = xsection * luminosidad
#  Peso      = xsection * luminosidad / SucesosGenerados en total 

# Color de los histogramas. Esto tambien deberia ir a un JSON/Diccionario
color=[ROOT.kGreen+2,ROOT.kBlue,ROOT.kViolet,ROOT.kBlue+2,ROOT.kPink]

# Definimos una funcion para dar estilo y normalizar los histogramas  
def StyleHisto(sample,variab,label,histColor,xsec,totalEvents,suffix=""):
     histo =file.Get(variab+"_"+sample)		
     histo.SetName("h"+variab+"_"+label+suffix)
     # Estilo:
     histo.SetLineColor(histColor)
     if "tautauH" not in label and "qqH" not in label and "mumuH" not in label:
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
#         if legendName[i]!="mumuH_Htautau" : # annade solo uno entre mumuH y mumuH_Htautau 
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
#   h["mumuH_Htautau"].Draw("hist,sames")	# Opcion para pintar la segnal por separado
 
   ylabel="events"
   xlabel=h[legendName[0]].GetXaxis().GetTitle() #"DiMuon Mass [GeV]"
   
   hStack.GetXaxis().SetTitle(xlabel)
   hStack.GetYaxis().SetTitle(ylabel)
   
   hStack.GetYaxis().SetTitleOffset(1.4)
   hStack.GetXaxis().SetTitleOffset(1.2)

   maxY=hStack.GetMaximum()
   if h["mumuH_Htautau"].GetMaximum()> maxY:
      maxY=h["mumuH_Htautau"].GetMaximum() 
   hStack.SetMaximum(maxY*1.2)
 
   # Leyenda
   leg = ROOT.TLegend(0.60,0.60,0.93,0.93)
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   for i in range(nProcesses):
   	entry=leg.AddEntry(h[ legendName[i] ],legendName[i],"lf")
   leg.Draw()
   
   # Etiquetas explicando las condiciones del proceso
   text = "#sqrt{{s}} = 240 GeV,   L = {:.0f} ab^{{-1}}".format(luminosity/1e6)
   channel = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + #tau^{+}_{e}#tau^{-}_{h}'
   
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
