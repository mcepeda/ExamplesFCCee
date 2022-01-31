#! /usr/bin/env python3
#
# 1) Setup ROOT with: 
#   "source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc8-opt/setup.sh"
# 2) Run this script as:
#   "python -i [python_script]"
#

import os, sys, math
import ROOT
import multiprocessing

import tdrstyle

# CMS Style (estilo para los plots)
tdrstyle.setTDRStyle()

# Proceso que vamos a estudiar
sampleName="eeHZ"

# Archivo de entrada
file = "/afs/ciemat.es/user/a/alcaraz/public/FCCee/"+sampleName+"_skimmed_reduced.root"
    
# Este rootfile tiene un 'tree' llamado events que vamos a
# convertir en un dataframe para analizarlo: 
print("\nProcessing file '%s'..." % (file))
df = ROOT.RDataFrame("events",file)

# Para ver que variables hay en el tree: imprimamos unos cuantos sucesos
df.Display("").Print()

# Vamos a ver cuantos muones hay en cada suceso:
# Dibujemos las variables con Histo1D( ("nombre","titulo;ejeX;ejeY", bines, minX,maxX) ,
# "columna/rama")
hNMuons = df.Histo1D(("hNMuons", "Muones Por Suceso ;N_{#mu};N_{Events}",10,0,10), "NMuon")

# Podemos filtrar la muestra para seleccionar parte de
# los sucesos (df.Filter( Seleccion, Explicacion) 
# Por ejemplo, vamos a fijarnos solo en los sucesos con dos muones: 
dfMuons = df.Filter("NMuon == 2", "Events with exactly two muons")

# Cuantos sucesos sobreviven a este corte?
report = dfMuons.Report()
report.Print()

# Tambien podemos añadir variables. 
# Dada la geometria del detector es comodo trabajar en cilindricas 
# en vez de cartesianas: vamos a definir el momento del muon en el plano transverso y
# añadirlo al dataframe: 
dfMuons = dfMuons.Define("Muon_pt","sqrt( pow(Muon_px,2)+pow(Muon_py,2) ) ")

# Vamos a volver a imprimir sucesos: ahora puedes ver que está ahi el Pt
dfMuons.Display({"Muon_px","Muon_py","Muon_pz","Muon_pt","Muon_charge"},10 ).Print()


# Dibujemos las variables con Histo1D( ("nombre","titulo;ejeX;ejeY", bines, minX,maxX) , "columna")
hAllMuonsPt = dfMuons.Histo1D(("hAllMuonPt", "MuonPt ; p_{T} (#mu) (GeV);N_{Events}", 100,0,100), "Muon_pt")
hAllMuonsCharge = dfMuons.Histo1D(("hAllMuonCharge", "MuonCharge ; Muon or AntiMuon?;N_{Events}",5,-2,2), "Muon_charge")


# Por ultimo pintamos las graficas:
c1 = ROOT.TCanvas("c1","HZ analysis canvas", 1024, 768)

hists = [ hNMuons, hAllMuonsPt, hAllMuonsCharge ]

for i,h in enumerate(hists):
    h.Draw("hist")
    output_pdf = h.GetName()+"_reco.pdf"
    c1.SaveAs(output_pdf)
    input("File '%s' already saved; press ENTER in this window when ready to close thisca canvas ... " % (output_pdf))
