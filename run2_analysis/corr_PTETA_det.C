#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooAddModel.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TTree.h"
#include "TKDTreeBinning.h"
#include "TStyle.h"
#include "TROOT.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace RooFit;
using namespace RooStats ;


TH2Poly *getHistogram2DPoly(TTree *tree, TH2F *hHisto, string varX, Float_t vMinX, Float_t vMaxX, string varY, Float_t vMinY, Float_t vMaxY, Int_t nBins, string weight = "", TCut cut = "1", double frac = -1, string option = "") {
  cout << "start" <<endl;
  Int_t nDim = 2;
  string name = "hHisto_" + varX + "_" + varY;
  cout << "step 2" <<endl;

  if (varX.find("PT") != string::npos) {
    varX  += "*1e-4";
    vMinX *= 1e-4;
    vMaxX *= 1e-4;
  }
  if (varY.find("PT") != string::npos) {
    varY  += "*1e-4";
    vMinY *= 1e-4;
    vMaxY *= 1e-4;
  }
  string var = varY + ":" + varX;
  cout << "step 3" <<endl;
  Int_t nEvt = hHisto->Integral();

  cout << "step 4" <<endl;
  cout << nEvt  <<endl;
  cout << nDim <<endl;
  cout << nBins <<endl;

  while ((nEvt + nDim) % nBins != 0)
    nEvt--;
  cout<<"new nEvt = "<<nEvt<<endl;

  if ((nEvt + nDim) < nBins)
    return 0;
  TCanvas *cDatas1TMP = new TCanvas("cDatas1TMP","tmp",800,600);
  cout << tree->GetEntries() << endl;
  tree->SetEstimate(tree->GetEntries()*100);//0);
  cout<<"plotting 1 ... "<<endl;
  tree->Draw("(Lambdab_PT*1e-4):Lambdab_ETA","Lambdab_ETA>1.5 && Lambdab_ETA<7. && Lambdab_PT<75000. && Lambdab_PT>0.");
  cout<<"must have plotted 1 ... "<<endl;
  Double_t *vX = tree->GetV2();
  Double_t *vY = tree->GetV1();

  Int_t size = 2 * (nEvt + nDim);//const int
  cout << size <<endl;
  cout << endl;
  cout << "Filling TKDTreeBinning..." << endl;
  cout << endl;

  cout << "step 5" <<endl;
  Double_t *data = new Double_t[size];
  Int_t c = 0;
  for (Int_t i = 0; i < nEvt; ++i) {
    data[i] = vX[i];
    c++;
  }
  data[c] = vMinX;
  c++;
  data[c] = vMaxX;
  c++;
  for (Int_t i = 0; i < nEvt; ++i) {
    data[i + (nEvt + nDim)] = vY[i];
    c++;
  }
  data[c] = vMinY;
  c++;
  data[c] = vMaxY;
  c++;

  nEvt += nDim;
  cout << endl;
  cout << "Creating TKDTreeBinning..." << endl;
  cout << endl;
  cout << "NDim = " << nDim << endl;
  cout << "NEvt = " << nEvt << endl;
  cout << "NBin = " << nBins << endl;
  cout << endl;

  TKDTreeBinning *kdBins = new TKDTreeBinning(nEvt, nDim, data, nBins);

  const Double_t* binsMinEdges = kdBins->GetBinsMinEdges();
  const Double_t* binsMaxEdges = kdBins->GetBinsMaxEdges();

  cout << kdBins->GetDataMin(0) << endl;
  cout << kdBins->GetDataMax(0) << endl;
  cout << kdBins->GetDataMin(1) << endl;
  cout << kdBins->GetDataMax(1) << endl;

  TH2Poly *h2Poly = new TH2Poly("h2Poly", "", kdBins->GetDataMin(0), kdBins->GetDataMax(0), kdBins->GetDataMin(1), kdBins->GetDataMax(1));

  h2Poly->SetXTitle(varX.c_str());
  h2Poly->SetYTitle(varY.c_str());

  double content = 0;
  double volume  = 0;
  double density = 0;

  for (int i = 0; i < nBins+1; ++i) {///here=
    int edgeDim = i * nDim;
    Int_t number =  h2Poly->AddBin(binsMinEdges[edgeDim], binsMinEdges[edgeDim + 1], binsMaxEdges[edgeDim], binsMaxEdges[edgeDim + 1]);
    h2Poly->SetBinContent(number, kdBins->GetBinDensity(i));
    density += kdBins->GetBinDensity(i);
  }

  content /= (double) nBins;
  volume  /= (double) nBins;
  density /= (double) nBins;

  cout << "Range X = " << h2Poly->GetXaxis()->GetXmin() << " , " << h2Poly->GetXaxis()->GetXmax() << endl;
  cout << "Range Y = " << h2Poly->GetYaxis()->GetXmin() << " , " << h2Poly->GetYaxis()->GetXmax() << endl;
  cout << endl;
  cout << "Min density = " << kdBins->GetBinMinDensity() << endl;
  cout << "Max density = " << kdBins->GetBinMaxDensity() << endl;
  cout << endl;
  cout << "<content> = " << content << endl;
  cout << "<volume>  = " << volume << endl;
  cout << "<density> = " << density << endl;
  cout << endl;
  return h2Poly;
}

void corr_PTETA_det(){
  gROOT->ProcessLine(".x ./lhcbStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  RooRealVar signal_yield_sw("signal_yield_sw","yield",-10.,10.);
  RooRealVar Lambdab_PT("Lambdab_PT","Lambdab_PT",0.,75000.);
  RooRealVar Lambdab_ETA("Lambdab_ETA","Lambdab_ETA",1.5,7.);
  //=====================Data=============================//
  // TFile* file = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Selection/LHCb_Jpsi_data_withSWeight.root");
  // TFile* file = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Selection/R1_Jpsi_data_withSWeight.root");
  TFile* file = new TFile("/ceph/users/pmachnik/Lambdab/Jpsi/Jpsi_sWeights.root");
  TTree* tree1 = (TTree*)file->Get("DecayTree");

  TCanvas *cDatas1 = new TCanvas("cDatas1","sPlot (signal)",800,600);
  TH2F* hDatasEP = new TH2F("hDatasEP", "sPlot" ,25,1.5,7,25,0,75000) ;
  hDatasEP->Sumw2();
  tree1->Draw("Lambdab_PT:Lambdab_ETA>>hDatasEP","signal_yield_sw");
  Double_t normEP = hDatasEP->Integral();
  cout<<normEP<<endl;
  hDatasEP->Scale(1./normEP);
  hDatasEP->Draw("COLZ");

  TChain *tree11 = new TChain("tree11");
  // tree11->Add("/ceph/users/jnicolini/scripts/Lambdab/Samples/PhD_Sel_MC/Lb2Lmm_JpsiMC_R1_TM_sel_trig_Jpsi_test.root/DecayTree");
  tree11->Add("/ceph/users/pmachnik/Lambdab/MC/weights/PTETA_pre.root/DecayTree");

  TCanvas *cDatas11 = new TCanvas("cDatas11","MC",800,600);
  TH2F* hDatas11EP = new TH2F("hDatas11EP", " " ,25,1.5,7.,25,0.,75000.) ;
  hDatas11EP->Sumw2();
  // tree11->Draw("Lambdab_PT:Lambdab_ETA>>hDatas11EP","((year==2011)*0.5 + (year==2012)*1)");
  tree11->Draw("Lambdab_PT:Lambdab_ETA>>hDatas11EP","((year==2015)*0.51 + (year==2016)*0.84 + (year==2017)*0.73 + (year==2018)*1)");

  Double_t norm11EP = hDatas11EP->Integral();
  cout<<norm11EP<<endl;

  const Int_t nbins = 250;//200; for R1

  TH2Poly *h2Poly     = NULL;
  TH2Poly *h2PolyData = NULL;

  h2Poly = (TH2Poly*)  getHistogram2DPoly(tree11,hDatas11EP, "Lambdab_ETA", 1.5, 7., "Lambdab_PT", 0., 75000., nbins , "", "1",  -1,  "");

  TCanvas *cPoly = new TCanvas("cPoly","Template histo",800,600);
  h2Poly->Draw("COLZ");

  /// K-FOLDING
  Int_t nfolds = 1;//10;
  TH2Poly *hFold[nfolds];
  int i=1;
  //for (int i = 0; i < nfolds; i++){
    cout<<"Hello powder"<<endl;
  delete gROOT->FindObject("hHisto1");
  delete gROOT->FindObject("hHisto2");
  delete gROOT->FindObject("hDataspdEP");
  delete gROOT->FindObject("hHistoD");

  stringstream str;
    str << i;
    string ss = str.str();
  const char* cs11 = "";//"kfold!=";
  string cs1a = cs11 + ss;
  const char *cs1 = cs1a.c_str();
  cout<<cs1<<endl;

  TCanvas *cPoly1 = new TCanvas("cPoly1","Histo 1",800,600);
  TH2Poly *hHisto1 = (TH2Poly*) h2Poly->Clone("hHisto1");
  hHisto1->ClearBinContents();
  tree1->Draw("(Lambdab_PT*1e-4):Lambdab_ETA>>hHisto1","signal_yield_sw");

  Double_t norm1 = hHisto1->Integral();
  hHisto1->Scale(1./norm1);
  hHisto1->Draw("COLZ");

  TCanvas *cPoly2 = new TCanvas("cPoly2","Histo 2",800,600);
  TH2Poly *hHisto2 = (TH2Poly*) h2Poly->Clone("hHisto2");
  //hHisto2->Sumw2();
  hHisto2->ClearBinContents();
  // tree11->Draw("(Lambdab_PT*1e-4):Lambdab_ETA>>hHisto2","((year==2011)*0.5 + (year==2012)*1)");
  tree11->Draw("(Lambdab_PT*1e-4):Lambdab_ETA>>hHisto2","((year==2015)*0.51 + (year==2016)*0.84 + (year==2017)*0.73 + (year==2018)*1)");

  Double_t norm2 = hHisto2->Integral();
  hHisto2->Scale(1./norm2);
  hHisto2->Draw("COLZ");
  cout<<"Step 10"<<endl;

  TCanvas *cFinal = new TCanvas("cFinal","weights",800,600);
  gPad ->SetRightMargin(0.15);
  TH2Poly *hHistoD = (TH2Poly*) h2Poly->Clone("hHistoD");
  hHistoD->ClearBinContents();
  hHistoD->Divide(hHisto1,hHisto2,1,1,"B");
  hHistoD->GetYaxis()->SetTitle("p_{T}(#Lambda_{b}^{0}) [10 GeV/c]");
  hHistoD->GetXaxis()->SetTitle("#eta(#Lambda_{b}^{0})");
  hHistoD->Draw("COLZ");

  const char* ss11 = "fold";
  const char* ss21 =  ";Lambdab_PT;Weight";
  string scs11 = ss11 + ss + ss21;

  const char *Title = scs11.c_str();

  const char* nom = "hFold";
  string prenom = nom + ss;
  const char *Prenom = prenom.c_str();
  hFold[i] = (TH2Poly*) hHistoD->Clone(Prenom);

  // string filenamepdf = "PTETAweights_det_larger_final_R1.pdf";
  string filenamepdf = "../plots/weights/PTETAweights.pdf";
  const char *Fname = filenamepdf.c_str();

  cFinal->SaveAs(Fname);

  //random checks
  Int_t nBins = hHistoD->GetNumberOfBins();
  TList *bin_list=hHistoD->GetBins();
  Double_t spdWeight88 = hHistoD->GetBinContent(0);
  cout<<"weight= "<<spdWeight88<<endl;
  Double_t spdWeight881 = hHistoD->GetBinContent(1);
  cout<<"weight= "<<spdWeight881<<endl;
  Double_t spdWeight882 = hHistoD->GetBinContent(2);
  cout<<"weight= "<<spdWeight882<<endl;
  Double_t spdWeight883 = hHistoD->GetBinContent(nbins-1);
  cout<<"weight= "<<spdWeight883<<endl;
  Double_t spdWeight884 = hHistoD->GetBinContent(nbins-2);
  cout<<"weight= "<<spdWeight884<<endl;

  // TFile out("PTETAweights_det_larger_final_R1.root", "RECREATE");
  TFile out("/ceph/users/pmachnik/Lambdab/MC/weights/PTETAweights.root", "RECREATE");
  hFold[1]->Write();
  out.Close();

};
