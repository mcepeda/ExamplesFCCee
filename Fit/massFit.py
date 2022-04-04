#####################################################################################
#   Slightly modified from the original                                             #
#   Macro to fit the higgs recoil mass out of FlatNtuples produced by FCCAnalysis   #
#   Author: Clement Helsens (clement.helsens@cern.ch)                               #
#   Date: November 2020                						    #
#   Edit: Andrea Barron and Maria Cepeda (April 2022)                               #
#####################################################################################

# Note! This uses a RooFit shape that is very recent. You need root 6.24 
# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/setup.sh
# (instead of LCG 99) to be able to access this function


import ROOT as r
import sys

if len(sys.argv)!=4 and len(sys.argv)!=6:
    print ('usage:   python massFit.py BASEDIR HISTONAME SELECTION BINLOW=120 BINHIGH=140')
    print ('example: python massFit.py ./ Recoil_mass NORM')
    print ('example: python massFit.py ./ Recoil_mass NORM 122 128')
    sys.exit(3)

basedir=sys.argv[1]
hname=sys.argv[2]
selection=sys.argv[3]
binlow=120
binhigh=140

if len(sys.argv)==6:
    binlow=float(sys.argv[4])
    binhigh=float(sys.argv[5])

if basedir[-1]!='/':
    basedir+='/'

tf=r.TFile(basedir+"/"+hname+"_"+selection+".root","READONLY")
histosig=tf.Get("h"+hname+"_eeHZ")
histobg1=tf.Get("h"+hname+"_eeZZ")
histobg2=tf.Get("h"+hname+"_eeWW")

binA=histosig.FindBin(binlow)
binB=histosig.FindBin(binhigh)

# print original yields:
print ("Expected yields at 5 ab-1")
print ("eeHZ  events: %4.1f" %histosig.Integral(binA,binB))
print ("eeZZ  events: %4.1f" %histobg1.Integral(binA,binB))
print ("eeWW  events: %4.1f" %histobg2.Integral(binA,binB))
print ("ZZ+WW events: %4.1f" %(histobg1.Integral(binA,binB)+histobg2.Integral(binA,binB)))
 
# If you want to fit signal + background:
#histosig.Add(histobg1)
#histosig.Add(histobg2)

#bins for the fit
x = r.RooRealVar("recoil", "M_{recoil} [GeV]", binlow, binhigh)

# data is breit-wigner convoluted with a gaussian, taken from histogram
dhData = r.RooDataHist("dhData", "dhData", r.RooArgList(x), r.RooFit.Import(histosig))

# Normal CB
cbmean  = r.RooRealVar("cbmean", "cbmean" , 125.0,  binlow, binhigh)
cbsigma = r.RooRealVar("cbsigma", "cbsigma" , 0.4, 0.0, 0.6)
n     = r.RooRealVar("n","n", 0.9,0.,100.0)
alpha   = r.RooRealVar("alpha","alpha", -1.2,-5.,-0.0001)
cball   = r.RooCBShape("cball", "crystal ball", x, cbmean, cbsigma, alpha, n)

# Gaussian
gmean  = r.RooRealVar("gmean", "gmean" , 125.0,  binlow, binhigh)
gsigma = r.RooRealVar("gsigma", "gsigma" , 0.8, 0.0, 2)
gauss = r.RooGaussian("gauss","gauss",x,gmean,gsigma)

# Asymmetric, gaussian core CB
cbmean  = r.RooRealVar("cbmean", "cbmean" , 125.0,  binlow, binhigh)
cbsigmaL = r.RooRealVar("cbsigmaL", "cbsigmaL" , 0.2, 0.0, 0.6)
cbsigmaR = r.RooRealVar("cbsigmaR", "cbsigmaR" , 0.2, 0.0, 0.6)
nL     = r.RooRealVar("nL","nL", 0.9,0.,100.0)
nR     = r.RooRealVar("nR","nR", 0.9,0.,100.0)
alphaL   = r.RooRealVar("alphaL","alphaL", 1.2,0,10)
alphaR   = r.RooRealVar("alphaR","alphaR", 1.2,0,10)
cballAssym   = r.RooCrystalBall("cballAssym", "crystal ball", x, cbmean, cbsigmaL, cbsigmaR,alphaL, nL,alphaR,nR)

lam = r.RooRealVar("lam","lam",-1e-3,-1,-1e-10)
bkg = r.RooExponential("bkg","bkg",x,lam)

nsig  = r.RooRealVar("nsig", "number of signal events", 20000, 0., 10000000)

# Combine signal PDFs if you want
#sigfrac1  = r.RooRealVar("sigfrac1", "cb1" , 0.8,0,1) 
##sigfrac2  = r.RooRealVar("sigfrac2", "cb2" , 0.2,0,1)
#signal = r.RooAddPdf("sig","Signal",r.RooArgList(cballAssym,gauss),r.RooArgList(sigfrac1))

#nbkg  = r.RooRealVar("nbkg", "number of background events", 20000, 10000, 50000)

#model = r.RooAddPdf("model","CBAssym+bkg",r.RooArgList(bkg,cballAssym),r.RooArgList(nbkg,nsig))
#model = r.RooAddPdf("model","CBassym+gauss+bkg",r.RooArgList(signal,bkg),r.RooArgList(nsig,nbkg))

# If you only want to fit the signal -> it makes sense to start like this :)
model = r.RooAddPdf("model","CB",r.RooArgList(cballAssym),r.RooArgList(nsig))
#model = r.RooAddPdf("model","CB",r.RooArgList(signal),r.RooArgList(nsig))

# To fit signal + background
#nbkg  = r.RooRealVar("nbkg", "number of background events", 20000, 10000, 50000)

#model = r.RooAddPdf("model","CBAssym+bkg",r.RooArgList(bkg,cballAssym),r.RooArgList(nbkg,nsig))
#model = r.RooAddPdf("model","CBassym+gauss+bkg",r.RooArgList(signal,bkg),r.RooArgList(nsig,nbkg))


result = model.fitTo(dhData,r.RooFit.Save(),r.RooFit.NumCPU(8,0),r.RooFit.Extended(True),r.RooFit.Optimize(False),r.RooFit.Offset(True),r.RooFit.Minimizer("Minuit2","migrad"),r.RooFit.Strategy(2))

frame = x.frame(r.RooFit.Title("recoil "), )
dhData.plotOn(frame,r.RooFit.Name("data"))
model.plotOn(frame,r.RooFit.Name("model"))

#ras_bkg = r.RooArgSet(bkg)
#model.plotOn(frame, r.RooFit.Components(ras_bkg), r.RooFit.LineStyle(r.kDashed), r.RooFit.Name("ras_bkg"))

ras_sig = r.RooArgSet(cballAssym)
model.plotOn(frame, r.RooFit.LineColor(r.kRed), r.RooFit.Components(ras_sig), r.RooFit.LineStyle(r.kDashed), r.RooFit.Name("ras_sig"))

c = r.TCanvas()
frame.Draw()

leg = r.TLegend(0.6,0.7,0.89,0.89)
leg.AddEntry(frame.findObject("data"),"FCC IDEA Delphes","ep")
leg.AddEntry(frame.findObject("model"),"S+B fit","l")
#leg.AddEntry(frame.findObject("ras_sig"),"Signal","l")
#leg.AddEntry(frame.findObject("ras_bkg"),"background","l")

leg.Draw()

r.gPad.SaveAs("fitResult.pdf")
r.gPad.SaveAs("fitResult.png")


# These were for the old parameters (it does not fit well)
#fCB  = r.TF1("fCB","crystalball")
# These are the results with the new samples!
#fCB.SetParameter(0, 1.)
#fCB.SetParameter(1, 125.)
#fCB.SetParameter(2, 0.5)
#fCB.SetParameter(3, -1.)
#fCB.SetParameter(4, 1.)
#fCB.SetRange(120.,128.)


#c = r.TCanvas()
#c.SetLogy()
#histosig.Draw()
#histosig.Fit(fCB)
