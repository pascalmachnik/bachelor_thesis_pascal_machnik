#include <string>
#include <sstream>

void addL0DiMuonweights_Run2(TString tuple) {

  TFile* newFile = new TFile(tuple, "UPDATE");
  TTree *newTree = (TTree*)newFile->Get("DecayTree");
  //TTree *newTree = (TTree*)newFile->Get("tree");

  Double_t L1_PT,L2_PT;
  Int_t year;
  newTree->SetBranchAddress("L1_PT",&L1_PT);
  newTree->SetBranchAddress("L2_PT",&L2_PT);
  newTree->SetBranchAddress("year",&year);
  Int_t totalEntries = newTree->GetEntries();

  Double_t L0Weight = 0.;//interpolated
  Double_t L0Weight_b = 0.;//binned for cross-check
  Double_t L0Weight_err = 0.;
  //create new branches
  Double_t L0Weight_18, L0Weight_17, L0Weight_16, L0Weight_15;
  Double_t L0Weight_b_18, L0Weight_b_17, L0Weight_b_16, L0Weight_b_15;
  Double_t L0Weight_err_18, L0Weight_err_17, L0Weight_err_16, L0Weight_err_15;

  TBranch *newbranch = newTree->Branch("L0DiMuonWeight_Run2_v1", &L0Weight, "L0DiMuonWeight_Run2_v1/D");
  TBranch *newbranch_b = newTree->Branch("L0DiMuonWeight_Run2_v1_b", &L0Weight_b, "L0DiMuonWeight_Run2_v1_b/D");
  TBranch *newbranch_err = newTree->Branch("L0DiMuonWeight_Run2_v1_err", &L0Weight_err, "L0DiMuonWeight_Run2_v1_err/D");

  //read the files with weights
  TFile* weightsFile_18 = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Calibration/2023_fixes/2018/L0DiMuonweights_18_2023fixes.root","READ");
  TFile* weightsFile_17 = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Calibration/2023_fixes/2017/L0DiMuonweights_17_2023fixes.root","READ");
  TFile* weightsFile_16 = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Calibration/2023_fixes/2016/L0DiMuonweights_16_2023fixes.root","READ");
  TFile* weightsFile_15 = new TFile("/ceph/users/jnicolini/scripts/Lambdab/Calibration/2023_fixes/2015/L0DiMuonweights_15_2023fixes.root","READ");

  Double_t totWeight = 0.;

  //now, have fun!
  //read one histo with weights per year
  TH1F* hWeights_18 = (TH1F*)weightsFile_18->Get("hFold0");
  TH1F* hWeights_17 = (TH1F*)weightsFile_17->Get("hFold0");
  TH1F* hWeights_16 = (TH1F*)weightsFile_16->Get("hFold0");
  TH1F* hWeights_15 = (TH1F*)weightsFile_15->Get("hFold0");
  Int_t nBins = hWeights_18->GetNbinsX(); // all them have the same binning!
  cout<<"nBins = "<<nBins<<endl;

  //now, loop over events
  for(Int_t i=0; i<totalEntries; i++){
    newTree->GetEntry(i);
    Double_t L1L2PT =L1_PT *L2_PT;
    L0Weight = 0.;  // all set to 0. if out of range
    L0Weight_b = 0.;
    //interpolated correction
    L0Weight_18 = hWeights_18->Interpolate(L1L2PT);
L0Weight_17 = hWeights_17->Interpolate(L1L2PT);
L0Weight_16 = hWeights_16->Interpolate(L1L2PT);
L0Weight_15 = hWeights_15->Interpolate(L1L2PT);
L0Weight = (year == 2018) * L0Weight_18 + (year == 2017) * L0Weight_17 + (year == 2016) * L0Weight_16 + (year == 2015) * L0Weight_15;

//binned for a cross-check
    Int_t nbin_18 = hWeights_18->FindBin(L1L2PT);
    L0Weight_b_18 = hWeights_18->GetBinContent(nbin_18);
    L0Weight_err_18 = hWeights_18->GetBinError(nbin_18);
    Int_t nbin_17 = hWeights_17->FindBin(L1L2PT);
    L0Weight_b_17 = hWeights_17->GetBinContent(nbin_17);
    L0Weight_err_17 = hWeights_17->GetBinError(nbin_17);
    Int_t nbin_16 = hWeights_16->FindBin(L1L2PT);
    L0Weight_b_16 = hWeights_16->GetBinContent(nbin_16);
    L0Weight_err_16 = hWeights_16->GetBinError(nbin_16);
    Int_t nbin_15 = hWeights_15->FindBin(L1L2PT);
    L0Weight_b_15 = hWeights_15->GetBinContent(nbin_15);
    L0Weight_err_15 = hWeights_15->GetBinError(nbin_15);
L0Weight_b = (year == 2018) * L0Weight_b_18 + (year == 2017) * L0Weight_b_17 + (year == 2016) * L0Weight_b_16 + (year == 2015) * L0Weight_b_15;
L0Weight_err = (year == 2018) * L0Weight_err_18 + (year == 2017) * L0Weight_err_17 + (year == 2016) * L0Weight_err_16 + (year == 2015) * L0Weight_err_15;

    totWeight += L0Weight;//check
    //fill new branches
    newbranch->Fill();
    newbranch_b->Fill();
    newbranch_err->Fill();
  }
  cout << "totWeight/nEvts = " << totWeight << " / " << totalEntries << " = " << totWeight/((Double_t)totalEntries);


  newFile->Write();

}
