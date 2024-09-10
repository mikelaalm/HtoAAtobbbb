#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.h"


using namespace TMVA;

void TMVAClassification()
{
  TFile *outputFile = TFile::Open("TMVA_M30.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outputFile, "!V");

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("MVAnalysis");

  TFile *input_signal = TFile::Open("histos_M30.root");
  TFile *input_DY4 = TFile::Open("histos_DY_4.root");
  TFile *input_DY3 = TFile::Open("histos_DY_3.root");
  TFile *input_DY2 = TFile::Open("histos_DY_2.root");
  TFile *input_DY1 = TFile::Open("histos_DY_1.root");
  TFile *input_TTDileptonic = TFile::Open("histos_TT_Dileptonic.root");
  TFile *input_TTSemileptonic = TFile::Open("histos_TT_Semileptonic.root");

  TTree *sigTree = (TTree*)input_signal->Get("t1");
  TTree *bgTree_DY4 = (TTree*)input_DY4->Get("t1");
  TTree *bgTree_DY3 = (TTree*)input_DY3->Get("t1");
  TTree *bgTree_DY2 = (TTree*)input_DY2->Get("t1");
  TTree *bgTree_DY1 = (TTree*)input_DY1->Get("t1");
  TTree *bgTree_TTDileptonic = (TTree*)input_TTDileptonic->Get("t1");
  TTree *bgTree_TTSemileptonic = (TTree*)input_TTSemileptonic->Get("t1");
  
  float bg_weight_DY4(1.0);
  float bg_weight_DY3(1.0);
  float bg_weight_DY2(1.0);
  float bg_weight_DY1(1.0);
  float bg_weight_TTDileptonic(1.0);
  float bg_weight_TTSemileptonic(1.0);
  float sig_weight(1.0);


  TLeaf *xpos_bg_DY4 = bgTree_DY4->GetLeaf("weight"); xpos_bg_DY4->GetBranch()->GetEntry(1);
  bg_weight_DY4 = xpos_bg_DY4->GetValue();

  TLeaf *xpos_bg_DY3 = bgTree_DY3->GetLeaf("weight"); xpos_bg_DY3->GetBranch()->GetEntry(1);
  bg_weight_DY3 = xpos_bg_DY3->GetValue();

  TLeaf *xpos_bg_DY2 = bgTree_DY2->GetLeaf("weight"); xpos_bg_DY2->GetBranch()->GetEntry(1);
  bg_weight_DY2 = xpos_bg_DY2->GetValue();

  TLeaf *xpos_bg_DY1 = bgTree_DY1->GetLeaf("weight"); xpos_bg_DY1->GetBranch()->GetEntry(1);
  bg_weight_DY1 = xpos_bg_DY1->GetValue();
  
  TLeaf *xpos_bg_TTDileptonic = bgTree_TTDileptonic->GetLeaf("weight"); xpos_bg_TTDileptonic->GetBranch()->GetEntry(1);
  bg_weight_TTDileptonic = xpos_bg_TTDileptonic->GetValue();

  TLeaf *xpos_bg_TTSemileptonic = bgTree_TTSemileptonic->GetLeaf("weight"); xpos_bg_TTSemileptonic->GetBranch()->GetEntry(1);
  bg_weight_TTSemileptonic = xpos_bg_TTSemileptonic->GetValue();

  TLeaf *xpos_sig = sigTree->GetLeaf("weight"); xpos_sig->GetBranch()->GetEntry(1);
  sig_weight = xpos_sig->GetValue();

  dataloader->AddSignalTree(sigTree, sig_weight);
  dataloader->AddBackgroundTree(bgTree_DY4, bg_weight_DY4);
  dataloader->AddBackgroundTree(bgTree_DY3, bg_weight_DY3);
  dataloader->AddBackgroundTree(bgTree_DY2, bg_weight_DY2);
  dataloader->AddBackgroundTree(bgTree_DY1, bg_weight_DY1);
  dataloader->AddBackgroundTree(bgTree_TTDileptonic, bg_weight_TTDileptonic);
  dataloader->AddBackgroundTree(bgTree_TTSemileptonic, bg_weight_TTSemileptonic);
  
  dataloader->AddVariable("m_H", 'F');
  dataloader->AddVariable("pt_H", 'F');
  dataloader->AddVariable("pt_b1", 'F');
  dataloader->AddVariable("dR_min", 'F');
  dataloader->AddVariable("pt_2b", 'F');
  dataloader->AddVariable("dM_A1_A2", 'F');
  dataloader->AddVariable("HT", 'F');
  dataloader->AddVariable("pt_Z", 'F');
  dataloader->AddVariable("dR_ll", 'F');
  dataloader->AddVariable("dPhi_ZH", 'F');
  dataloader->AddVariable("n_jets_after_cuts", 'F');
  dataloader->AddVariable("btag_3", 'F');
  dataloader->AddVariable("met_pt", 'F');

 
			       
  dataloader->PrepareTrainingAndTestTree("" ,"SplitMode=Random:!V");
					
  
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT","!H:!V:NTrees=1000:MinNodeSize=5.0%:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
  

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  Double_t AUC = factory->GetROCIntegral(dataloader, "BDT");

  std::cout << "Area under ROC curve: " << AUC << std::endl;
  
  outputFile->Close();
  delete factory;
}
