#!/usr/bin/env python

# Importamos las clases de ROOT: 
from ROOT import *

# Cargamos el archivo root:
f = TFile('/afs/ciemat.es/user/a/alcaraz/public/FCCee/eeHZ_skimmed_reduced.root')

gStyle.SetOptStat(0)

# Cargamos el tree indicando el directorio y el nombre del tree:
tree = f.Get("events")

# Cuantos sucesos hay en el tree?
entries=tree.GetEntries()
print( "El tree tiene %d sucesos" %entries)

hMuon_px=TH1F("hMuon_px","Px del Muon", 100,-50, 50)
hNMuon  =TH1F("hNMuon","Numero de Muones", 10,0,10)

for event in tree: # loop sobre los sucesos del tree 

	hNMuon.Fill(event.NMuon) # rellenamos el histograma de numero de muones

	if event.NMuon!=2 : 
		continue 

	for i in range(0,event.NMuon): # loop sobre los muones en este suceso
        	hMuon_px.Fill(event.Muon_px.at(i)) # rellenamos el histograma 

# Definimos un "lienzo" para pintar los histogramas
canvasMuPx=TCanvas("canvas","Momento de los Muones",800,600)  # (nombre, titulo, dimensionx, dimensiony)
canvasMuPx.cd() # carga el lienzo

#Estilo del Plot (Color, Etiquetas de los ejes, etc)
hMuon_px.SetLineWidth(2)
hMuon_px.SetXTitle("Muon Px [GeV]")
hMuon_px.SetYTitle("Sucesos / GeV")

# Pinta los histogramas
hMuon_px.Draw("hist")

#Leyenda del histograma. Importante para saber que estás pintando.
legend=TLegend(0.5,0.7,0.90,0.90) # (posInicio_x, posInicio_y, posFin_x, posFin_y)
legend.AddEntry(hMuon_px,"Muones","l")
legend.Draw()

# Guarda el canvas:
canvasMuPx.SaveAs("canvasMuPx.png")

canvasNMu=TCanvas("canvas","Numero de Muones",800,600)  # (nombre, titulo, dimensionx, dimensiony)
canvasNMu.cd() # carga el lienzo

#Estilo del Plot (Color, Etiquetas de los ejes, etc)
hNMuon.SetLineWidth(2)
hNMuon.SetXTitle("Numero de Muons")
hNMuon.SetYTitle("Sucesos / GeV")

# Pinta los histogramas
hNMuon.Draw("hist")

#Leyenda del histograma. Necesario cuando tienes mas de una grafica.
legend=TLegend(0.5,0.7,0.90,0.90) # (posInicio_x, posInicio_y, posFin_x, posFin_y)
legend.AddEntry(hNMuon,"Muones","l")
legend.Draw()

#Archivo con histogramas y canvas para estudiarlo luego:
out = TFile('histos.root',"RECREATE")
out.cd()
hMuon_px.Write()
hNMuon.Write()

