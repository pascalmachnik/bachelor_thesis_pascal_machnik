import ROOT as r
import argparse
import numpy as np

# argparse setup
parser = argparse.ArgumentParser()
parser.add_argument("--input_MC", action="store", nargs="+", dest="input_MC", type=str)
parser.add_argument("--input_data", action="store", nargs="+", dest="input_data", type=str)
parser.add_argument("--track_type", action="store", dest="track_type", type=str, choices=["ll", "dd"], required=True)
args = parser.parse_args()

# Enable multi-threading in ROOT
r.EnableImplicitMT(8)

# set global style configurations
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gStyle.SetOptStat(0)
r.gStyle.SetLegendFont(132)
r.gStyle.SetLegendTextSize(0.06)

# Create RDataFrame for input files
df_MC = r.RDataFrame("DecayTree", args.input_MC, ["Lambdab_LOKI_MASS_LConstr", "BDTG"])  # MC signal
df_data = r.RDataFrame("DecayTree", args.input_data, ["Lambdab_LOKI_MASS_LConstr", "BDTG"])  # data

# Create snapshot of the dataframes
df_MC.Snapshot("DecayTree", f"/ceph/users/pmachnik/temp/{args.track_type}_MC_reduced.root", ["Lambdab_LOKI_MASS_LConstr", "BDTG"])
df_data.Snapshot("DecayTree", f"/ceph/users/pmachnik/temp/{args.track_type}_Data_reduced.root", ["Lambdab_LOKI_MASS_LConstr", "BDTG"])

# Load datasets 
treeDataAll = r.TChain("DecayTree") #data
treeDataAllMC = r.TChain("DecayTree") #MC

treeDataAll.Add(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_reduced.root") #data
treeDataAllMC.Add(f"/ceph/users/pmachnik/temp/{args.track_type}_MC_reduced.root") #MC

BDTval= 0.8 #starting scan point 

# Plotting setup
c1 = r.TCanvas("c1","c1",600,900)
c1.Divide(1,3)
c2 = r.TCanvas("c2","c2",600,900)
c2.Divide(1,3)

#it's 40 scan points from from -1. to 1.; 0.99999 is used instead of 1 to avoid ROOT from creating an extra bin after 1.
hist = r.TH1F("hist","Punzi cut optimisation",40,BDTval,0.99999)#-1.0#40
hist2 = r.TH1F("hist2","Cut optimisation 2",40,BDTval,0.99999)#-1.0#40
histb = r.TH1F("histb","Background yield",40,BDTval,0.99999)
hists = r.TH1F("hists","Signal yield",40,BDTval,0.99999)
hist.GetXaxis().SetTitle("BDT cut value")
hist.GetYaxis().SetTitle("P")
hist2.GetXaxis().SetTitle("BDT cut value")
hist2.GetYaxis().SetTitle("P")
histb.GetXaxis().SetTitle("BDT cut value")
histb.GetYaxis().SetTitle("# bkg events")
hists.GetXaxis().SetTitle("BDT cut value")
hists.GetYaxis().SetTitle("signal eff")

while (BDTval<0.9995):
    #a string which defines the rare fit range and the BDT cut
    nstr = "BDTG>"+str(BDTval)
    cutBDT = r.TCut(nstr)

    f_new = r.TFile(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_rare.root","recreate")
    treeData = treeDataAll.CopyTree(nstr)
    
    #this must be the MC of your rare mode, preselected (see above)
    dfsiga = r.RDataFrame("DecayTree",args.input_MC) ###
    dfsiga_sel = dfsiga.Filter(cutBDT.GetTitle())
    dfsiga_sel.Snapshot("treesig", f"/ceph/users/pmachnik/temp/{args.track_type}_preselMC_rare_L.root", {"Lambdab_LOKI_MASS_LConstr"})

    f_news = r.TFile(f"/ceph/users/pmachnik/temp/{args.track_type}_preselMC_rare_L.root")
    treesig = f_news.Get("treesig")
    nmc = treesig.GetEntries()
    print("Number MC events=", nmc)

    #fit
    Lb_M = r.RooRealVar("Lambdab_LOKI_MASS_LConstr", "Lambdab_LOKI_MASS_LConstr", 5620, 5300, 5940)
    r.RooAbsData.setDefaultStorageType(r.RooAbsData.Vector)
    datasetMC = r.RooDataSet("n2", "n2", treesig, r.RooArgSet(Lb_M))
    print("MC loaded")

    dataset = r.RooDataSet("n","n",treeData,r.RooArgSet(Lb_M))#data
    print("data loaded")

    #first, fit the signal MC
    #Linear sum of double sided CB Model
    mean_Lb = r.RooRealVar("mean_Lb", "mu", 5620, 5560, 6670)
    sigma_Lb = r.RooRealVar("sigma_Lb","sigma_Lb",6.,3.5,20.)
    AlphaL_Lb = r.RooRealVar("AlphaL_Lb", "AlphaL_Lb",3., 0.05, 9)
    nCBL_Lb = r.RooRealVar("nCBL_Lb","nCBL_Lb",1., 0.05, 250.)
    AlphaR_Lb = r.RooRealVar("AlphaR_Lb", "AlphaR_Lb",3.,0.05,20)
    nCBR_Lb = r.RooRealVar("nCBR_Lb","nCBR_Lb",1., 0.05, 250.)
    frac_sigma2 = r.RooRealVar("frac_sigma2","frac_sigma2",1.83,1.1,5.)
    sigma_Lb2 = r.RooFormulaVar("sigma_Lb2", "@0*@1", r.RooArgList(sigma_Lb,frac_sigma2))
    # sigma_Lb2 = r.RooRealVar("sigma_Lb2","sigma_Lb2",6.,3.5,12.)
    AlphaL_Lb2 = r.RooRealVar("AlphaL_Lb2", "AlphaL_Lb2",3.,0.05,30)
    nCBL_Lb2 = r.RooRealVar("nCBL_Lb2","nCBL_Lb2",1., 0.05, 250.)
    AlphaR_Lb2 = r.RooRealVar("AlphaR_Lb2", "AlphaR_Lb2",3.,0.05,20)
    nCBR_Lb2 = r.RooRealVar("nCBR_Lb2","nCBR_Lb2",1., 0.05, 250.)
    CB_pdf_Lb1 = r.RooCrystalBall("CB_pdf_Lb1", "bifurCB", Lb_M, mean_Lb, sigma_Lb, AlphaL_Lb, nCBL_Lb, AlphaR_Lb, nCBR_Lb)
    CB_pdf_Lb2 = r.RooCrystalBall("CB_pdf_Lb2", "bifurCB", Lb_M, mean_Lb, sigma_Lb2, AlphaL_Lb2, nCBL_Lb2, AlphaR_Lb2, nCBR_Lb2)
    frac_Lb = r.RooRealVar("frac_Lb","frac_Lb",0.5,0.2,0.8)

    # create signal PDF
    CB_pdf_Lb = r.RooAddPdf("CB_pdf_Lb","CB_pdf_Lb", r.RooArgList(CB_pdf_Lb1,CB_pdf_Lb2),r.RooArgList(frac_Lb))
    CB_pdf_Lb.fitTo(datasetMC)
    print("first fit done")

    sigma_new = (sigma_Lb.getValV()*1) #*1.X to account for resolution in the data
    sigma_Lb.setVal(sigma_new)
    sigma_Lb.setConstant(r.kTRUE)
    frac_sigma2.setConstant(r.kTRUE)
    frac_Lb.setConstant(r.kTRUE)
    AlphaL_Lb.setConstant(r.kTRUE)
    nCBL_Lb.setConstant(r.kTRUE)
    AlphaR_Lb.setConstant(r.kTRUE)
    nCBR_Lb.setConstant(r.kTRUE)
    AlphaL_Lb2.setConstant(r.kTRUE)
    nCBL_Lb2.setConstant(r.kTRUE)
    AlphaR_Lb2.setConstant(r.kTRUE)
    nCBR_Lb2.setConstant(r.kTRUE)

    # create signal PDF for data
    signal_pdf_Lb1 = r.RooCrystalBall("CB_pdf_Lb1", "bifurCB", Lb_M, mean_Lb, sigma_Lb, AlphaL_Lb, nCBL_Lb, AlphaR_Lb, nCBR_Lb)
    signal_pdf_Lb2 = r.RooCrystalBall("CB_pdf_Lb2", "bifurCB", Lb_M, mean_Lb, sigma_Lb2, AlphaL_Lb2, nCBL_Lb2, AlphaR_Lb2, nCBR_Lb2)
    signal = r.RooAddPdf("signal_pdf_Lb","signal_pdf_Lb", r.RooArgList(signal_pdf_Lb1,signal_pdf_Lb2),r.RooArgList(frac_Lb))

    #background pdf
    start_slope = -0.001
    slope = r.RooRealVar("slope","exp slope",start_slope,-1,0.03)
    exp_pdf = r.RooExponential("exp_pdf", "exp_pdf", Lb_M, slope)

    #adapt to dataset
    start_yield_comb = 4000
    start_yield_signal =500
    yield_comb = r.RooRealVar("yield_comb","yield",start_yield_comb,0,100000) #background
    yield_signal = r.RooRealVar("yield_signal","yield",start_yield_signal,0,6000) #signal
    shapes = r.RooArgList(exp_pdf, signal)
    yields = r.RooArgList(yield_comb, yield_signal)

    #total PDF
    total = r.RooAddPdf("total","total",shapes,yields)
    total.fitTo(dataset,r.RooFit.NumCPU(10), r.RooFit.Extended(1), r.RooFit.Offset(True))

    # print all fits for cross-check
    nbins = 46#25
    frame = Lb_M.frame(r.RooFit.Title("Constraint double sided CB fit to \Lambda_{b} mass"), r.RooFit.Bins(nbins))
    dataset.plotOn(frame, r.RooFit.MarkerSize(0.9), r.RooFit.Name("datahistogram"))
    total.plotOn(frame, r.RooFit.Name("Signal"), r.RooFit.Components('signal_pdf_Lb'), r.RooFit.LineColor(r.kRed))
    total.plotOn(frame, r.RooFit.Name("Comb_BKG"), r.RooFit.Components('comb'), r.RooFit.LineColor(r.kGreen-3), r.RooFit.LineStyle(r.kDashed), r.RooFit.DrawOption('F'), r.RooFit.FillColor(r.kGreen-3), r.RooFit.MoveToBack())
    total.plotOn(frame, r.RooFit.Name("PDF"), r.RooFit.LineColor(r.kBlue))
    total.paramOn(frame, r.RooFit.Format("NELU", r.RooFit.AutoPrecision(2)), r.RooFit.Layout(0.6, 0.9, 0.6)) #0.82
    frame.getAttText().SetTextSize(0.028)

    # plot pulls
    framePull = Lb_M.frame()
    pH = frame.pullHist()
    framePull.addPlotable(pH, "P")

    # adjust stuff for the pull plot
    framePull.GetXaxis().SetTitle("m(\Lambda \mu^{#minus}\mu^{#plus}) [MeV/c^{2}]")#\Psi(2S)(\mu^{#minus}\mu^{#plus})
    framePull.GetYaxis().SetTitle("pull")
    framePull.GetXaxis().SetLabelSize(.15)
    framePull.GetXaxis().SetTitleSize(.15)
    framePull.GetYaxis().SetTitleSize(.15)
    framePull.GetYaxis().SetTitleOffset(.4)
    framePull.GetYaxis().SetNdivisions(306, r.kTRUE)
    framePull.GetYaxis().SetLabelSize(.135)
    framePull.GetYaxis().SetRangeUser(-6.5,6.5)
    framePull.GetXaxis().SetTitleOffset(1.15)

    # plotting foo
    frame.SetTitle("Constraint double sided CB fit to \Lambda_{b}^{0} mass")
    frame.GetXaxis().SetTitle("m(\Lambda \mu^{#minus}\mu^{#plus}) [MeV/c^{2}]")#\Psi(2S)(\mu^{#minus}\mu^{#plus})
    frame.GetYaxis().SetTitle("events")

    # create legend and add all the individual components
    Legend = r.TLegend(0.65,0.65,0.92,0.85)#(0.75, 0.7, 0.9, 0.85)
    Legend.SetFillColor(0)

    #Legend.SetLineColor(r.kWhite)
    Legend.AddEntry(frame.getObject(1), 'Data','LEP')
    Legend.AddEntry(frame.getObject(3), 'PDF', 'L')
    Legend.AddEntry(frame.getObject(2), 'Signal','L')
    Legend.AddEntry(frame.getObject(0), 'Comb BKG','L')
    #LHCb preliminary
    lhcbName2mc = r.TPaveText(r.gStyle.GetPadLeftMargin() + 0.045,0.89 - r.gStyle.GetPadTopMargin(), r.gStyle.GetPadLeftMargin() + 0.31, 0.97 - r.gStyle.GetPadTopMargin(),"BRNDC") #0.31->0.3
    lhcbName2mc.AddText("LHCb preliminary")
    lhcbName2mc.SetFillColor(0)
    lhcbName2mc.SetBorderSize(0)

    # plotting foo part 2
    can = r.TCanvas("massFit", "massFit", 600, 600)
    pad1 = r.TPad("pad1","pad1",0,0.3,1,1)
    pad2 = r.TPad("pad2","pad2",0,0,1,0.3)
    pad1.SetBottomMargin(0.02)#20)
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.45)
    pad1.Draw()
    pad2.Draw()
    pad1.cd()
    frame.Draw()
    lhcbName2mc.Draw("same")
    pad2.cd()
    framePull.Draw()
    can.cd()
    BDTval = round(BDTval, 4)  # round BDTval to 4 decimal places
    can.SaveAs('../plots/FOM/R2/{}/Data_rare_dCB_Fit_{}BDT.pdf'.format(args.track_type, BDTval))
    can.Close()

    #Silence RooFit to prevent it printing the fit result:
    r.RooMsgService.instance().setSilentMode(r.kTRUE)

    #And suppress INFO messages for plot range yields:
    r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.WARNING)
    yield_comb.setConstant(r.kTRUE)
    yield_signal.setConstant(r.kTRUE)
    slope.setConstant(r.kTRUE)
    sigma_Lb.setConstant(r.kTRUE)

    #signal efficiency at a given BDT cut
    S = 1.*nmc/(treeDataAllMC.GetEntries())
    sl = slope.getValV()
    width = sigma_Lb.getValV()

    #get normalisation constant of the exponential, where 5300 and 5940 are fit range edges
    coef = (yield_comb.getValV())/((np.exp(sl*5940)-np.exp(sl*5300))/sl)

    #extrapolate to get bkg yield inside the 3sigma window around Lambdab mass (5620) LHCb
    B = coef * ((np.exp(sl*(5620+width*3))-np.exp(sl*(5620-width*3)))/sl)

    #Punzi definition
    Punzi = S/(5/2+np.sqrt(B))
    S1 = 1360 #full Run2 #adapt to dataset
    sign2 = S1*S/(np.sqrt(S1*S+B))
    print("cut=", BDTval)
    print("signal=", S)
    print("bkg=", B)
    print("Punzi significance=", Punzi)
    print("significance2=", sign2)
    print(BDTval)
    print(Punzi)
    print(sign2)
    print("///**************///***************///***************///*************///")

    hist.Fill(BDTval,Punzi)
    hist2.Fill(BDTval,sign2)
    hists.Fill(BDTval,S)
    histb.Fill(BDTval,B)
    #0.0125 is the step needed to get 160 scan points from -1.0 to 1., 0.025 for 80 etc
    BDTval=BDTval+0.01 #0.05#0.04 for 50 points 
    if (BDTval>0.9995):
        BDTval = 0.9995

c1.cd(1)
hist.SetMarkerStyle(8)
hist.Draw("Phist")

c1.cd(2)
histb.SetMarkerStyle(8)
histb.Draw("Phist")

c1.cd(3)
hists.SetMarkerStyle(8)
hists.Draw("Phist")

c2.cd(1)
hist2.SetMarkerStyle(8)
hist2.Draw("Phist")

c2.cd(2)
histb.SetMarkerStyle(8)
histb.Draw("Phist")

c2.cd(3)
hists.SetMarkerStyle(8)
hists.Draw("Phist")
c1.SaveAs(f"../plots/FOM/R2/TMVA_BDT_rare_Punzi_{args.track_type}.pdf")
c2.SaveAs(f"../plots/FOM/R2/TMVA_BDT_rare_{args.track_type}.pdf")
out = r.TFile(f"/ceph/users/pmachnik/Lambdab/FOM/TMVA_BDT_rare_FOM_{args.track_type}.root", "RECREATE")
hist.Write()
hist2.Write()
hists.Write()
histb.Write()
out.Close()

# remove temporary files
import os
if os.path.exists(f"/ceph/users/pmachnik/temp/{args.track_type}_preselMC_rare_L.root"):
    os.remove(f"/ceph/users/pmachnik/temp/{args.track_type}_preselMC_rare_L.root")
if os.path.exists(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_rare.root"):
    os.remove(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_rare.root")
if os.path.exists(f"/ceph/users/pmachnik/temp/{args.track_type}_MC_reduced.root"):
    os.remove(f"/ceph/users/pmachnik/temp/{args.track_type}_MC_reduced.root")
if os.path.exists(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_reduced.root"):
    os.remove(f"/ceph/users/pmachnik/temp/{args.track_type}_Data_reduced.root")