import ROOT
import tdrstyle

# CMS Style (estilo para los plots)
tdrstyle.setTDRStyle()

# Carga el archivo con los histogramas

file=ROOT.TFile("histos.root","READONLY")

histos=["Z_mass","Recoil_mass","LeadMuon_Pt","LeadMuon_Theta","LeadMuon_Phi",\
	"SecondMuon_Pt","SecondMuon_Theta","SecondMuon_Phi","Z_theta","Recoil_theta",\
	"Z_Pt","Recoil_Pt","Z_y","Recoil_y"]


# Loop sobre los histogramas guardados en 'histos.root'

for histoName in histos:

   sampleName=["WW","ZZ","HZ"]
   nProcesses=len(sampleName)
   
   # Los factores de normalizacion de cada proceso vienen dados por la seccion eficaz:
   # eeHZ: 0.201868 pb
   # eeZZ: 1.35899  pb
   # eeWW: 16.4385  pb
   xsection=[16.4385,1.35899,0.201868]
   
   # Luminosidad acumulada: 5 ab-1 = 5000 fb-1 = 500000 pb-1
   luminosity = 5e6
   
   # Numero de sucesos MC producidos originamente: 
   totalNumberOfEvents=10000000

   # La normalizacion de los histogramas seguira :
   #  N_Sucesos = xsection * luminosidad
   #  Peso      = xsection * luminosidad / SucesosGenerados en total 
   
   # Color de los histogramas
   color=[ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kRed]

   # Definimos una funcion para dar estilo y normalizar los histogramas  
   def StyleHisto(sample,variab,histColor,xsec):
        histo =file.Get(variab+"_ee"+sample)		
	
	# Estilo:
        histo.SetLineColor(histColor)
        if sample!="HZ":
        	histo.SetFillColor(histColor)
        histo.SetLineWidth(3)

	# Aqui se aplica la normalizacion: 
        histo.Scale(xsec*luminosity/totalNumberOfEvents)

        return histo


   # Definimos un 'stack' de histogramas apilados:
   h={}
   hStack    = ROOT.THStack()
   
   for i in range(nProcesses):
         h[ sampleName[i] ] = StyleHisto(sampleName[i],histoName,color[i],xsection[i])
         hStack.Add( h[ sampleName[i] ] )


   # Ahora pintamos los histogramas:
   c1 = ROOT.TCanvas("c1","HZ analysis canvas", 650,600)
   
#   ROOT.gStyle.SetPadGridX(True)
#   ROOT.gStyle.SetPadGridY(True)
#   ROOT.gStyle.SetGridStyle(3)
#   ROOT.gStyle.SetHistLineWidth(3)
#   ROOT.gStyle.SetHistMinimumZero()
   ROOT.gStyle.SetOptTitle(True)
   
   c1.SetTicks(1,1)
   c1.SetLeftMargin(0.15)
   c1.SetRightMargin(0.08)
   c1.SetTopMargin(0.06)
  
   hStack.Draw("hist")
   
   ylabel="Events"
   xlabel=h[sampleName[0]].GetXaxis().GetTitle() #"DiMuon Mass [GeV]"
   
   hStack.GetXaxis().SetTitle(xlabel)
   hStack.GetYaxis().SetTitle(ylabel)
   
   hStack.GetYaxis().SetTitleOffset(1.8)
   hStack.GetXaxis().SetTitleOffset(1.2)

   maxY=hStack.GetMaximum() 
   hStack.SetMaximum(maxY*1.2)
 
   # Leyenda
   leg = ROOT.TLegend(0.80,0.70,0.93,0.93)
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   for i in range(nProcesses):
   	entry=leg.AddEntry(h[ sampleName[i] ],sampleName[i],"lf")
   leg.Draw()
   
   # Etiquetas explicando las condiciones del proceso
   text = "#sqrt{{s}} = 240 GeV,   L = {:.0f} ab^{{-1}}".format(luminosity/1e6)
   channel = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + X'
   
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
