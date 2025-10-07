#include "TH2Poly.h"

void addPTETAweights_det(TString tuple, Int_t kfold=0) {

  TFile* newFile = new TFile(tuple, "UPDATE");
  TTree *newTree = (TTree*)newFile->Get("DecayTree");

  Double_t Lambdab_PT, Lambdab_ETA;
  Double_t Lambdab_PT_rec, Lambdab_ETA_rec;
  Double_t Lambdab_PT_tr;
  Double_t Lambdab_ETA_tr, Lambdab_PZ_tr;
  Int_t Lambdab_BKGCAT; //Int_t

  newTree->SetBranchAddress("Lambdab_PT",&Lambdab_PT_rec);
  newTree->SetBranchAddress("Lambdab_ETA",&Lambdab_ETA_rec);
  newTree->SetBranchAddress("Lambdab_TRUEPT",&Lambdab_PT_tr);
  newTree->SetBranchAddress("Lambdab_TRUEP_Z",&Lambdab_PZ_tr);
  newTree->SetBranchAddress("Lambdab_BKGCAT",&Lambdab_BKGCAT);

  Int_t totalEntries = newTree->GetEntries();

  Double_t kinWeight = 0.;//binned only. Interpolate() does not work for TH2Poly.
  Double_t kinWeightMax = 0;
  Double_t kinWeightMin = 0;
  Double_t totWeight = 0.;

  TBranch *newbranch = newTree->Branch("kinWeight_Ap24", &kinWeight, "kinWeight_Ap24/D");
  // TFile* weightsFile = new TFile("PTETAweights_det_larger_final_R1.root","READ");
  TFile* weightsFile = new TFile("/ceph/users/pmachnik/Lambdab/MC/weights/PTETAweights.root","READ");
  //now, have fun!
  //if no kfolding
  if (kfold==0){
    //read one histo with weights
    TH2Poly* hWeights = (TH2Poly*)weightsFile->Get("hFold1");
    Int_t nBins = hWeights->GetNumberOfBins();
    cout<<"nBins = "<<nBins<<endl;
    TList *bin_list=hWeights->GetBins();

    //now, loop over events
    for(Int_t i=0; i<totalEntries; i++){
      newTree->GetEntry(i);
      if ((Lambdab_PT_tr == 0) || (Lambdab_BKGCAT==60)){
	      Lambdab_PT = Lambdab_PT_rec;
	      Lambdab_ETA = Lambdab_ETA_rec;
      }
      else {
	      Lambdab_PT = Lambdab_PT_tr;
	      Lambdab_ETA=TMath::ASinH(Lambdab_PZ_tr/Lambdab_PT_tr);
      }
      Lambdab_PT=Lambdab_PT*0.0001;
      //Lambdab_PT = Lambdab_PT*(0.0001);
      kinWeight = 0.;
      for(Int_t j=0;j< nBins; j++) {
	      TH2PolyBin *bin= (TH2PolyBin*)bin_list->At(j);
	      if (bin->IsInside(Lambdab_ETA,Lambdab_PT)){
	        kinWeight = hWeights->GetBinContent(j+1);
          break; //we found the bin, no need to continue
	        //Explanation.
	        //There is some weird thing with the binning that boundaries of the bin are declared from 0, while bin values are declared from 1. That's why we do j+1 when reading weights. This have been tested on a simple example of 3 bins, and shown to work properly.
	      }
	      if (Lambdab_PT>50000*0.0001){
	        if (bin->IsInside(Lambdab_ETA,49999*0.0001)){
            kinWeight = hWeights->GetBinContent(j+1);
          }
	      }
      }
      totWeight += kinWeight;
      newbranch->Fill();
    }
  }
  //what if we are wise and do k-folding
  else if (kfold==1){
    TBranch *newbranch_max = newTree->Branch("kinWeight_Nov18_true3_max", &kinWeightMax, "Weight_max_of_folds/D");
    TBranch *newbranch_min = newTree->Branch("kinWeight_Nov18_true3_min", &kinWeightMin, "Weight_min_of_folds/D");
    Double_t nFold_d;
    newTree->SetBranchAddress("kfold",&nFold_d);
    TH2Poly* hWeights[10];
    Int_t nBins[10];
    TList *bin_list[10];

    //Load all the 10 histograms first
    for(Int_t i=0; i<10; i++){
      stringstream str;
      str << i;
      string ss = str.str();
      const char* nom = "hFold";
      string prenom = nom + ss;
      const char *Prenom = prenom.c_str();
      hWeights[i] = (TH2Poly*)weightsFile->Get(Prenom);
      nBins[i] = hWeights[i]->GetNumberOfBins();
      cout<<"fold number "<<i<< ", nBins = "<<nBins[i]<<endl;
      bin_list[i]=hWeights[i]->GetBins();
      //      cout<<"fold number "<<i<< ", bin list:  "<<bin_list[i]<<endl;
    }
    totWeight = 0.;
    //Now, loop over events
    for(Int_t i=0; i<totalEntries; i++){
      newTree->GetEntry(i);
      //Lambdab_PT = Lambdab_PT*(0.0001);
      if ((Lambdab_PT_tr == 0) || (Lambdab_BKGCAT==60)){
	      Lambdab_PT = Lambdab_PT_rec;
	      Lambdab_ETA = Lambdab_ETA_rec;
      }
	    else {
        Lambdab_PT = Lambdab_PT_tr;
        Lambdab_ETA=TMath::ASinH(Lambdab_PZ_tr/Lambdab_PT_tr);
	    }
      Lambdab_PT=Lambdab_PT*0.0001;
      kinWeight = 0.;
      Int_t nfold = (Int_t)nFold_d;

      for(Int_t j=0;j< nBins[nfold]; j++) {
	      TH2PolyBin *bin= (TH2PolyBin*)bin_list[nfold]->At(j);
	      if (bin->IsInside(Lambdab_ETA,Lambdab_PT)){
	        kinWeight = hWeights[nfold]->GetBinContent(j+1);
          kinWeightMax = std::max({hWeights[0]->GetBinContent(j+1),hWeights[1]->GetBinContent(j+1),hWeights[2]->GetBinContent(j+1),hWeights[3]->GetBinContent(j+1),hWeights[4]->GetBinContent(j+1),hWeights[5]->GetBinContent(j+1),hWeights[6]->GetBinContent(j+1),hWeights[7]->GetBinContent(j+1),hWeights[8]->GetBinContent(j+1),hWeights[9]->GetBinContent(j+1)});
          kinWeightMin = std::min({hWeights[0]->GetBinContent(j+1),hWeights[1]->GetBinContent(j+1),hWeights[2]->GetBinContent(j+1),hWeights[3]->GetBinContent(j+1),hWeights[4]->GetBinContent(j+1),hWeights[5]->GetBinContent(j+1),hWeights[6]->GetBinContent(j+1),hWeights[7]->GetBinContent(j+1),hWeights[8]->GetBinContent(j+1),hWeights[9]->GetBinContent(j+1)});
	        //Explanation.
	        //There is some weird thing with the binning that boundaries of the bin are declared from 0, while bin values are declared from 1. That's why we do j+1 when reading weights. This have been tested on a simple example of 3 bins, and shown to work properly.
	      }
	      if (Lambdab_PT>50000*0.0001){
	        if (bin->IsInside(Lambdab_ETA,49999*0.0001)){
	          kinWeight = hWeights[nfold]->GetBinContent(j+1);
            kinWeightMax = std::max({hWeights[0]->GetBinContent(j+1),hWeights[1]->GetBinContent(j+1),hWeights[2]->GetBinContent(j+1),hWeights[3]->GetBinContent(j+1),hWeights[4]->GetBinContent(j+1),hWeights[5]->GetBinContent(j+1),hWeights[6]->GetBinContent(j+1),hWeights[7]->GetBinContent(j+1),hWeights[8]->GetBinContent(j+1),hWeights[9]->GetBinContent(j+1)});
            kinWeightMin = std::min({hWeights[0]->GetBinContent(j+1),hWeights[1]->GetBinContent(j+1),hWeights[2]->GetBinContent(j+1),hWeights[3]->GetBinContent(j+1),hWeights[4]->GetBinContent(j+1),hWeights[5]->GetBinContent(j+1),hWeights[6]->GetBinContent(j+1),hWeights[7]->GetBinContent(j+1),hWeights[8]->GetBinContent(j+1),hWeights[9]->GetBinContent(j+1)}); 
	        }
	      }
      }
      totWeight += kinWeight;
      newbranch->Fill();
      newbranch_max->Fill();
      newbranch_min->Fill();
    }
  }
  else {
    cout << "Second parameter should be 0 (no kfold) or 1 (kfold)" << endl;
    return;
  }

  cout << "totWeight/nEvts = " << totWeight << " / " << totalEntries << " = " << totWeight/((Double_t)totalEntries) << endl;
  newFile->Write();

  // Plot histogram of the new weights
  newFile->cd();
  TH1D* hWeightDist = new TH1D("hWeightDist", "Distribution of kinWeight_Ap24;kinWeight_Ap24;Events", 100, 0, 5);
  newTree->Draw("kinWeight_Ap24 >> hWeightDist");
  TCanvas* c1 = new TCanvas("c1", "Weight Distribution", 800, 600);
  hWeightDist->Draw();
  c1->SaveAs("../plots/weights/kinWeight_Ap24_distribution.pdf");
}
