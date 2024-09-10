#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"

#include "TMath.h"
#include "TLorentzVector.h"


using namespace std;

//defining the getDeltaR function:

double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
  {
    double delta_phi;
    double delta_eta;

    delta_phi = vec_1.Phi() - vec_2.Phi();
    delta_eta = vec_1.Eta() - vec_2.Eta();

    return std::sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
    }

struct JetAndBtag{
  TLorentzVector jet;
  float btag;
};

bool sortBtag(const JetAndBtag &jet_i, const JetAndBtag &jet_j){
  return jet_i.btag > jet_j.btag;
}


bool jet_matched(TLorentzVector jet, std::vector<TLorentzVector> vec_bb){
  bool match = false;
  float dRmin = 999;
  
  for (int i = 0; i < vec_bb.size(); i++){
    float dR = getDeltaR(jet, vec_bb[i]);
    if(dR < dRmin) dRmin = dR;
  }
  if(dRmin < 0.4) match = true;
  return match;
}

bool jet_matched_lowm(TLorentzVector jet, std::vector<TLorentzVector> vec_bb){
  bool match_lowm = false;

  float dR1 = getDeltaR(jet, vec_bb[0]);
  float dR2 = getDeltaR(jet, vec_bb[1]);

  if(dR1 < 0.4 && dR2 < 0.4) match_lowm = true;
  return match_lowm;
}

void MyClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MyClass.C
//      root> MyClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;
  bool verbose(false);
  bool signal(false);

  
  TString filename=std::string(f->GetName());
  cout << "input file name: " << filename << endl;

  //if(filename.Contains("_M")) signal = true;
  if(filename == "analysis_M12.root" || filename == "analysis_M15.root" || filename == "analysis_M20.root" || filename == "analysis_M25.root" || filename == "analysis_M30.root" || filename == "analysis_M60.root") signal = true;
  TString njet;
  TString mass_point;
  
  if (filename.Contains("DY_0")) njet = "_0";
  else if (filename.Contains("DY_1")) njet = "_1";
  else if (filename.Contains("DY_2")) njet = "_2";
  else if (filename.Contains("DY_3")) njet = "_3";
  else if (filename.Contains("DY_4")) njet = "_4";

  if (filename.Contains("M12")) mass_point = "_M12";
  else if (filename.Contains("M15")) mass_point = "_M15";
  else if (filename.Contains("M20")) mass_point = "_M20";
  else if (filename.Contains("M25")) mass_point = "_M25";
  else if (filename.Contains("M30")) mass_point = "_M30";
  else if (filename.Contains("M60")) mass_point = "_M60";
				       
			     
  TString fname="histos.root";
   if(signal) {
     fname=TString::Format("histos%s_lowm.root", mass_point.Data());
   } else if (filename.Contains("DY")) {
     fname=TString::Format("histos_DY%s_lowm.root", njet.Data());
   } else if (filename.Contains("Dileptonic")) {
     fname="histos_TT_Dileptonic_lowm.root";
   } else {
     fname="histos_TT_Semileptonic_lowm.root";
   }
   
     

   TFile fout(fname.Data(),"RECREATE");

   cout << "output file name: " << fname.Data() << endl;

  
  //create histograms for the Higgs boson

  TH1F *h_ptH  = new TH1F("h_ptH" ," ; Higgs p_{T} [GeV] ; Events",100,0.,600.);
  TH1F *h_etaH = new TH1F("h_etaH"," ; Higgs #eta  ; Events",100,-5.,5.);
  TH1F *h_phiH = new TH1F("h_phiH"," ; Higgs #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mH   = new TH1F("h_mH"  ," ; Higgs mass [GeV]  ; Events",100,120.,130.);
 

  //create histograms for the Z boson

  TH1F *h_ptZ  = new TH1F("h_ptZ", " ; Z boson p_{T} ; Events ",100,0.,500.);
  TH1F *h_etaZ = new TH1F("h_etaZ"," ; Z boson #eta  ; Events",100,-5.,5.);
  TH1F *h_phiZ = new TH1F("h_phiZ"," ; Z boson #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mZ   = new TH1F("h_mZ"  ," ; Z boson mass  ; Events",100,80.,100.);

  //create histograms for the first A particle (imc=6)

  TH1F *h_ptA1  = new TH1F("h_ptA1" ," ; A_{1} p_{T} ; Events",100,0.,500.);
  TH1F *h_etaA1 = new TH1F("h_etaA1"," ; A_{1} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phiA1 = new TH1F("h_phiA1"," ; A_{1} #phi  ; Events",100,-TMath::Pi(),TMath::Pi()); 
  TH1F *h_mA1   = new TH1F("h_mA1"  ," ; A_{1} mass  ; Events",100,10.,70.);

 
  //create histograms for the second A particle (imc=7)

  TH1F *h_ptA2  = new TH1F("h_ptA2" ," ; A_{2} p_{T} ; Events",100,0.,500.);
  TH1F *h_etaA2 = new TH1F("h_etaA2"," ; A_{2} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phiA2 = new TH1F("h_phiA2"," ; A_{2} #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mA2   = new TH1F("h_mA2"  ," ; A_{2} mass  ; Events",100,10.,70.);

   
  //create histograms for the first b quark (imc=8 & mc_momidx=6)

  TH1F *h_ptb11  = new TH1F("h_ptb11" , " ; b_{11} p_{T} ; Events",100,0.,500.);
  TH1F *h_etab11 = new TH1F("h_etab11", " ; b_{11} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phib11 = new TH1F("h_phib11", " ; b_{11} #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mb11   = new TH1F("h_mb11"  , " ; b_{11} mass  ; Events",100,0.,20.);
  

  //create histograms for the second b quark (imc=9 & mc_momidx=6)

  TH1F *h_ptb12  = new TH1F("h_ptb12" , " ; b_{12} p_{T} ; Events",100,0.,500.);
  TH1F *h_etab12 = new TH1F("h_etab12", " ; b_{12} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phib12 = new TH1F("h_phib12", " ; b_{12} #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mb12   = new TH1F("h_mb12"  , " ; b_{12} mass  ; Events",100,0.,20.);

  //create histograms for the third b quark (imc=10 & mc_momidx=7)
  
  TH1F *h_ptb21  = new TH1F("h_ptb21" , " ; b_{21} p_{T} ; Events",100,0.,500.);
  TH1F *h_etab21 = new TH1F("h_etab21", " ; b_{21} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phib21 = new TH1F("h_phib21", " ; b_{21} #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mb21   = new TH1F("h_mb21"  , " ; b_{21} mass  ; Events",100,0.,20.);

  //create histograms for the fourth b quark (imc=11 & mc_momidx=7)

  TH1F *h_ptb22  = new TH1F("h_ptb22" , " ; b_{22} p_{T} ; Events",100,0.,500.);
  TH1F *h_etab22 = new TH1F("h_etab22", " ; b_{22} #eta  ; Events",100,-5.,-5.);
  TH1F *h_phib22 = new TH1F("h_phib22", " ; b_{22} #phi  ; Events",100,-TMath::Pi(),TMath::Pi());
  TH1F *h_mb22   = new TH1F("h_mb22"  , " ; b_{22} mass  ; Events",100,0.,20.);

  //pT(2b)-m(2b) plane (truth)
  
  TH2F *h_pt2b_m2b_truth_vector = new TH2F("h_pt2b_m2b_truth_vector", " ; m(2b) from A_{1} ; p_{T}(2b) from A_{1}",100, 0., 400., 100, 0., 500.);
  TH1F *h_pt2b_truth_vector = new TH1F("h_pt2b_truth_vector", " ; p_{T}(2b) from A_{1}; Events", 100, 0., 500.);
  TH2F *h_pt2b_m2b_truth_vector_other = new TH2F("h_pt2b_m2b_truth_vector_other", " ; m(2b) from A_{2} ; p_{T}(2b) from A_{2}",100, 0., 400., 100, 0., 500.);
  TH1F *h_pt2b_truth_vector_other = new TH1F("h_pt2b_truth_vector_other", " ; p_{T}(2b) from A_{2} ; Events", 100, 0., 500.);

  //create histograms for the antilepton

  TH1F *h_pt_l1  = new TH1F("h_pt_l1" ," ; l_{1} p_{T} [GeV]; Events",100,0.,500.);
  TH1F *h_eta_l1 = new TH1F("h_eta_l1"," ; l_{1} #eta ; Events",100,-5.,5.);
  TH1F *h_phi_l1 = new TH1F("h_phi_l1"," ; l_{1} #phi ; Events",100,-TMath::Pi(),TMath::Pi());
  //TH1F *h_m_l1 = new TH1F("h_m_l1"  ," ; l_{1} mass [GeV] ; Events",100,0.,1.);
 

  //create histograms for the lepton

  TH1F *h_pt_l2  = new TH1F("h_pt_l2" ," ; l_{2} p_{T} [GeV]; Events",100,0.,500.);
  TH1F *h_eta_l2 = new TH1F("h_eta_l2"," ; l_{2} #eta ; Events",100,-5.,5.);
  TH1F *h_phi_l2 = new TH1F("h_phi_l2"," ; l_{2} #phi ; Events",100,-TMath::Pi(),TMath::Pi());
  //TH1F *h_m_l2 = new TH1F("h_m_l2"  ," ; l_{2} mass [GeV] ; Events",100,0.,1.);
 
  
  //create histograms for deltaR

  TH1F *h_dRA    = new TH1F("h_dRA", " ; #DeltaR(A,A)          ; Events", 100, 0.,7.); // A particle pair
  TH1F *h_dR1    = new TH1F("h_dR1", " ; #DeltaR(bb) from A_{1} ; Events", 50, 0.,7.); // first  b quark pair: b11 and b12
  TH1F *h_dR2    = new TH1F("h_dR2", " ; #DeltaR(bb) from A_{2} ; Events", 50, 0.,7.); // second b quark pair: b21 and b22
  
  
  //create a histogram for the scalar and vector  sum of the pT of the q/g's (hadronic system)

  TH1F *h_HT_G = new TH1F("h_HT_G", " ; H_{T} [GeV] ; Events", 100, 0. ,600.);
  TH1F *h_ptb_v_sum = new TH1F("h_ptb_v_sum"," ; p_{T} bbb(b) [GeV] ; Events", 100, 0., 600.);

  //leptonic
  
  TH1F *h_pt_Z_G = new TH1F("h_pt_Z_G", " ; p_{T} (2l) [GeV] ; Events", 100., 0., 600.);
  
  //create a histogram for the absolute delta Phi phi.Had-phi.Lep

  TH1F *h_dPhi_ZH_G   = new TH1F("h_dPhi_ZH_G", " ; |#Delta#phi|(H,Z) ; Events", 100, 0., TMath::Pi());
  TH1F *h_phi_H_G = new TH1F("h_phi_H_G", " ; phi had ; Events", 100, -TMath::Pi(), TMath::Pi());
  TH1F *h_phi_Z_G = new TH1F("h_phi_Z_G", " ; phi lep ; Events", 100, -TMath::Pi(), TMath::Pi());
 
  //create histograms for dM_min and dM_max
  
  TH1F *h_dM_min = new TH1F("h_dM_min", " ; #DeltaM,min (b pairs) ; Events", 100, 0., 0.01); 
  TH1F *h_dM_max = new TH1F("h_dM_max", " ; #DeltaM,max (b pairs) ; Events", 100, 0., 500.);
  
  //create histograms for m_Had and m_Lep
  
  TH1F *h_m_H_G = new TH1F("h_m_H_G"," ; m_{bbb(b)} [GeV] ; Events", 100, 0., 800.);
  TH1F *h_m_Z_G = new TH1F("h_m_Z_G"," ; m(2l) [GeV] ; Events", 100, 0., 300.);

  //create pT histograms for each b quark after sorting
  
  TH1F *h_pt_1_truth   = new TH1F("h_pt_1_truth"  , " ; p_{T} of the leading b quark (truth) ; Events"   , 100, 0., 300.); 
  TH1F *h_pt_2_truth   = new TH1F("h_pt_2_truth"  , " ; b_{2} p_{T}/truth ; Events"   , 100, 0., 300.);
  TH1F *h_pt_3_truth   = new TH1F("h_pt_3_truth"  , " ; b_{3} p_{T}/truth ; Events"   , 100, 0., 300.); 
  TH1F *h_pt_4_truth   = new TH1F("h_pt_4_truth"  , " ; b_{4} p_{T}/truth ; Events"   , 100, 0., 300.);

  TH1F *h_pt_1_global   = new TH1F("h_pt_1_global"  , " ; (q/g)_{1} p_{T} [GeV] ; Events"   , 100, 0., 400.); 
  TH1F *h_pt_2_global   = new TH1F("h_pt_2_global"  , " ; (q/g)_{2} p_{T} [GeV] ; Events"   , 100, 0., 400.);
  TH1F *h_pt_3_global   = new TH1F("h_pt_3_global"  , " ; (q/g)_{3} p_{T} [GeV] ; Events"   , 100, 0., 400.); 
  TH1F *h_pt_4_global   = new TH1F("h_pt_4_global"  , " ; (q/g)_{4} p_{T} [GeV] ; Events"   , 100, 0., 400.);


  //create pT  histograms for leading lepton and sub-leading lepton

  TH1F *h_pt_leadL    = new TH1F("h_pt_leadL",    " ; l_{1} p_{T} [GeV]  ; Events", 100, 0., 300.);
  TH1F *h_pt_subleadL = new TH1F("h_pt_subleadL", " ; l_{2} p_{T} [GeV]  ; Events", 100, 0., 300.);
  TH1F *h_dR_ll_G = new TH1F("h_dR_ll_G", " ; #DeltaR(ll) ; Events", 100, 0., 7.);
			       

  //create histograms for low dR analysis
 
  TH2F *h_pt2b_m2b_mindR_vector = new TH2F("h_pt2b_m2b_mindR_vector"," ; m(A_{1}) [GeV] ; p_{T}(A_{1}) [GeV]",100, 0., 200., 100, 0., 600.);
  TH2F *h_pt2b_m2b_other_vector  = new TH2F("h_pt2b_m2b_other_vector", " ; m(A_{2}) [GeV]; p_{T}(A_{2}) [GeV]",100, 0., 200., 100, 0., 600.);
  
  TH2F *h_m1_m2          = new TH2F("h_m1_m2",     " ; m(A_{1}); m(A_{2})", 100, 0., 400., 100, 0., 400.);
  TH1F *h_m2b_mindR      = new TH1F("h_m2b_mindR", " ; m(A_{1}) [GeV]; Events", 100, 0., 200.);
  TH1F *h_m2b_other      = new TH1F("h_m2b_other", " ; m(A_{2}) [GeV] ; Events", 100, 0., 200.);
				    
  TH2F *h_pT1_pT2_vector = new TH2F("h_pT1_pT2_vector", " ; p_{T}(A_{1});  p_{T}(A_{2})", 100, 0., 500., 100, 0., 500.);
  TH2F *h_dRmin_dRother   = new TH2F("h_dRmin_dRother", " ; #DeltaR_{min}; #DeltaR_{other}", 50, 0., 7., 50, 0., 7.);

  TH1F *h_dR_min_G = new TH1F("h_dR_min_G" , " ; #DeltaR_{min}  ; Events", 100, 0., 7.);
  TH1F *h_dR_other_G = new TH1F("h_dR_other_G", " ; #DeltaR_{other} ; Events", 100, 0., 7.);

  TH1F *h_pT_min_vector  = new TH1F("h_pT_min_vector"  ," ; p_{T}(A_{1}) [GeV]; Events", 100, 0., 600.);
  TH1F *h_pT_other_vector= new TH1F("h_pT_other_vector"," ; p_{T}(A_{2}) [GeV]; Events", 100, 0., 600.);

  //b quarks multiplicty
  TH1F *h_nb = new TH1F("h_nb", " ; number of quarks/gluons ; Events", 10, 0., 10.);

  //new mass variables 
  TH1F *h_dM_A1_A2_gen_level = new TH1F("h_dM_A1_A2_gen_level", " ; m(A_{1}) - m(A_{2}) ; Events", 100, -400., 200.);

  //DETECTOR LEVEL ANALYSIS

  //jet mulitplicity
  TH1F *h_n_jets = new TH1F("h_n_jets", " ; number of jets ; Events", 10, 0., 12.);
  TH1F *h_n_jets_after_cuts = new TH1F("h_n_jets_after_cuts", " ; number of jets ; Events", 10, 0., 12.);
  TH1F *h_n_bjets = new TH1F("h_n_bjets", " ; number of b jets ; Events", 10, 0., 12.);
  TH1F *h_n_bjets_after_cuts = new TH1F("h_n_bjets_after_cuts", " ; number of b jets ; Events", 10, 0., 12.);

  TH1F *h_btag_discriminator = new TH1F("h_btag_discriminator", " ; b - tag discriminator ; Events", 100, 0., 1.); 

  //hadronic and leptonic system
  TH1F *h_pt_H   = new TH1F("h_pt_H" , " ; p_{T} bbbj [GeV]  ; Events"  , 50, 0., 600.);
  TH1F *h_m_H    = new TH1F("h_m_H"  , " ; m_{bbbj} [GeV]     ; Events"  , 40, 0.,800.);
  TH1F *h_phi_H = new TH1F("h_phi_H",  " ; #phi (4b) ; Events"  , 25, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_H = new TH1F("h_eta_H",  " ; #eta (4b) ; Events"  , 25, -5., 5.);
  
  TH1F *h_pt_Z  = new TH1F("h_pt_Z"  , " ; p_{T}(2l) [GeV]; Events"  , 50, 0., 600.);
  TH1F *h_m_Z   = new TH1F("h_m_Z"   , " ; m(2l) [GeV]    ; Events"  , 50, 0., 300.);
  TH1F *h_phi_Z = new TH1F("h_phi_Z" , " ; #phi (2l); Events"  , 25, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_Z = new TH1F("h_eta_Z" , " ; #eta (2l); Events"  , 25, -5., 5.);
  TH1F *h_mll_before_cut = new TH1F("h_mll_before_cut", " ; m_{Z} ; Events", 50, 0., 300.);
  
  TH1F *h_dPhi_ZH    = new TH1F("h_dPhi_ZH"    , " ; |#Delta#phi|(H,Z); Events", 25, 0., TMath::Pi());
  
  //scalar sum of pT of jets
  TH1F *h_HT = new TH1F("h_HT" , " ; H_{T} [GeV]; Events", 33, 0., 600.);
  
  //jet kinematics - 4 highest b tag discriminator jets
  TH1F *h_pt_jet1  = new TH1F("h_pt_jet1",  " ;  b_{1}, p_{T} [GeV] ; Events", 100, 0., 400.);
  TH1F *h_eta_jet1 = new TH1F("h_eta_jet1", " ; b_{1}, #eta; Events", 100, -5., 5.);
  TH1F *h_phi_jet1 = new TH1F("h_phi_jet1", " ; b_{1}, #phi; Events", 100,-TMath::Pi(),TMath::Pi());
  TH1F *h_m_jet1   = new TH1F("h_m_jet1",   " ; b_{1}, mass; Events", 100, 0., 20.);
  
  TH1F *h_pt_jet2  = new TH1F("h_pt_jet2",  " ; b_{2}, p_{T} [GeV]; Events", 100, 0., 400.);
  TH1F *h_eta_jet2 = new TH1F("h_eta_jet2", " ; b_{2}, #eta; Events", 100, -5., 5.);
  TH1F *h_phi_jet2 = new TH1F("h_phi_jet2", " ; b_{2}, #phi; Events", 100,-TMath::Pi(),TMath::Pi());
  TH1F *h_m_jet2   = new TH1F("h_m_jet2",   " ; b_{2}, mass; Events", 100, 0., 20.);
  
  TH1F *h_pt_jet3  = new TH1F("h_pt_jet3",  " ; b_{3}, p_{T} [GeV] ; Events", 100, 0., 400.);
  TH1F *h_eta_jet3 = new TH1F("h_eta_jet3", " ; b_{3}, #eta; Events", 100, -5., 5.);
  TH1F *h_phi_jet3 = new TH1F("h_phi_jet3", " ; b_{3}, #phi; Events", 100,-TMath::Pi(),TMath::Pi());
  TH1F *h_m_jet3   = new TH1F("h_m_jet3",   " ; b_{3}, mass; Events", 100, 0., 20.);

  TH1F *h_pt_jet4  = new TH1F("h_pt_jet4",  " ; b_{4}, p_{T} [GeV] ; Events", 100, 0., 400.);
  TH1F *h_eta_jet4 = new TH1F("h_eta_jet4", " ; b_{4}, #eta; Events", 100, -5., 5.);
  TH1F *h_phi_jet4 = new TH1F("h_phi_jet4", " ; b_{4}, #phi; Events", 100,-TMath::Pi(),TMath::Pi());
  TH1F *h_m_jet4   = new TH1F("h_m_jet4",   " ; b_{4}, mass; Events", 100, 0., 20.);

  TH1F *h_pt_b1 = new TH1F("h_pt_b1", " ;  b_{1}, p_{T} [GeV] ; Events", 50, 0., 400.);
  TH1F *h_pt_2b = new TH1F("h_pt_2b", " ;  p_{T} (2b) [GeV]   ; Events", 50, 0., 400.);
  
  //2 most energetic leptons
  TH1F *h_pt_lepton1 = new TH1F("h_pt_lepton1", " ; lepton 1, p_{T} [GeV] ; Events", 50, 0., 300.);
  TH1F *h_pt_lepton2 = new TH1F("h_pt_lepton2", " ; lepton 2, p_{T} [GeV]; Events", 50, 0., 300.);

  //minimum dR analysis
  TH2F *h_pt2jets_m2jets_mindR = new TH2F("h_pt2jets_m2jets_mindR", " ; m(A_{1}) [GeV] ; p_{T}(A_{1}) [GeV]", 30, 0., 200., 50, 0., 600.);
  TH2F *h_pt2jets_m2jets_other = new TH2F("h_pt2jets_m2jets_other", " ; m(A_{2}) [GeV] ; p_{T}(A_{2}) [GeV]", 30, 0., 200., 50, 0., 600.);
  
  TH1F *h_dR_min = new TH1F("h_dR_min", " ; #DeltaR_{min} ; Events", 50, 0., 5.);
  TH1F *h_dR_other = new TH1F("h_dR_other", " ; #DeltaR_{other} ; Events", 50, 0., 5.);

  TH1F *h_pt_2jets_mindR = new TH1F("h_pt_2jets_mindR", " ; p_{T} (A_{1}) [GeV] ; Events", 50, 0., 600.);
  TH1F *h_pt_2jets_other = new TH1F("h_pt_2jets_other", " ; p_{T} (A_{2}) [GeV] ; Events", 50, 0., 600.);

  TH1F *h_m_2jets_mindR = new TH1F("h_m_2jets_mindR", " ; m(A_{1}) [GeV] ; Events", 30, 0., 200.);
  TH1F *h_m_2jets_other = new TH1F("h_m_2jets_other", " ; m(A_{2}) [GeV] ; Events", 30, 0., 200.);

  //dR(ll)
  TH1F *h_dR_ll = new TH1F("h_dR_ll", "; #DeltaR(ll) ; Events", 50, 0., 7.);
			   

  //b tag discriminator histograms
  TH1F *h_btag1 = new TH1F("h_btag1", " ; b - tag (1) ; Events", 50, 0., 1.);
  TH1F *h_btag2 = new TH1F("h_btag2", " ; b - tag (2) ; Events", 50, 0., 1.);
  TH1F *h_btag3 = new TH1F("h_btag3", " ; b - tag (3) ; Events", 50, 0., 1.);
  TH1F *h_btag4 = new TH1F("h_btag4", " ; b - tag (4) ; Events", 50, 0., 1.);

  //mass variables
  TH1F *h_dM_A1_A2 = new TH1F("h_dM_A1_A2", " ; m(A_{1}) - m(A_{2}) ; Events", 50, -400., 200.);

  //MET
  TH1F *h_met_pt = new TH1F("h_met_pt", " ; missing E_{T} ; Events", 50, 0., 500.);

  //deltaR best matched pair
  TH1F *h_best_matched_b_dR = new TH1F("h_best_matched_b_dR", " ; #DeltaR ; Events", 50, 0., 2.);
  TH1F *h_deltaR_correct_pair = new TH1F("h_deltaR_correct_pair", " ; #DeltaR ; Events", 50, 0., 5.);

  TH2F *h_ptjet_pt2b = new TH2F("h_ptjet_pt2b", " ; p_{T} jet [GeV] ; p_{T} (2b), truth [GeV]", 50, 0., 300., 50., 0., 300.);

  
  
  TTree t1("t1", "A Tree with variables for the MVA");
  

  float m_H;
  float pt_H;
  float pt_b1;
  // float dR_min;
  float pt_2b;
  float dM_A1_A2;
  float HT;
  float pt_Z;
  float dR_ll;
  float dPhi_ZH;
  int n_jets_after_cuts;
  float weight(1.);
  float btag_1;
  float btag_2;
  float met_pt_;
  
  t1.Branch("m_H", &m_H, "m_H/F");
  t1.Branch("pt_H", &pt_H, "pt_H/F");
  t1.Branch("pt_b1", &pt_b1, "pt_b1/F");
  // t1.Branch("dR_min", &dR_min, "dR_min/F");
  t1.Branch("pt_2b", &pt_2b, "pt_2b/F");
  t1.Branch("dM_A1_A2", &dM_A1_A2, "dM_A1_A2/F");
  t1.Branch("HT", &HT, "HT/F");
  t1.Branch("pt_Z", &pt_Z, "pt_Z/F");
  t1.Branch("dR_ll", &dR_ll, "dR_ll/F");
  t1.Branch("dPhi_ZH", &dPhi_ZH, "dPhi_ZH/F");
  t1.Branch("n_jets_after_cuts", &n_jets_after_cuts, "n_jets_after_cuts/I");
  t1.Branch("weight", &weight, "weight/F");
  t1.Branch("btag_1", &btag_1, "btag_1/F");
  t1.Branch("btag_2", &btag_2, "btag_2/F");
  t1.Branch("met_pt", &met_pt, "met_pt/F");
  

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "" << endl;
  cout << "Total number of events = " << nentries << endl;
  
  //N_expected

  float L_int     = 43.5E3;
  float N_expected(0.);
  

  //signal
  float sigma_signal = 0.8839*(3*0.0336)*0.75; //[pb]
  float N_expected_signal = sigma_signal*L_int;
  
  //DY + 4 Jets
  float sigma_DY_4  = 51.4; //[pb]
  float N_expected_DY_4  = sigma_DY_4*L_int;

  //DY + 3 Jets
  float sigma_DY_3 = 96.36; //[pb]
  float N_expected_DY_3 = sigma_DY_3*L_int;

  //DY + 2 Jets
  float sigma_DY_2 = 331.4; //[pb]
  float N_expected_DY_2 = sigma_DY_2*L_int;

  //DY + 1 Jet
  float sigma_DY_1 = 1016; //[pb]
  float N_expected_DY_1 = sigma_DY_1*L_int;

  //DY + 0 Jets
  float sigma_DY_0 = 4895; //[pb]
  float N_expected_DY_0 = sigma_DY_0*L_int;
  
  //TT Dileptonic
  float sigma_TT_Dileptonic = 88.29; //[pb]
  float N_expected_TT_Dileptonic = sigma_TT_Dileptonic*L_int;

  //TT Semileptonic
  float sigma_TT_Semileptonic = 365.34; //[pb]
  float N_expected_TT_Semileptonic = sigma_TT_Semileptonic*L_int;

  if(signal){

    cout << "" << endl;
    cout << "Number of expected events in signal:" << N_expected_signal << endl;
    weight = N_expected_signal/totalNumberofEvents;
    cout << "" << endl;
    cout << "Signal weight: " << weight << endl;
    N_expected = N_expected_signal;

  }
  
  if(fname == "histos_DY_4_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in DY + 4 jets background: " << N_expected_DY_4 << endl;
    weight = N_expected_DY_4/totalNumberofEvents;
    cout << "" << endl;
    cout << "Drell Yan + 4 jets weight: " << weight << endl;
    N_expected = N_expected_DY_4;

  }

  if(fname == "histos_DY_3_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in DY + 3 jets background: " << N_expected_DY_3 << endl;
    weight = N_expected_DY_3/totalNumberofEvents;
    cout << "" << endl;
    cout << "Drell Yan + 3 jets weight: " << weight << endl;
    N_expected = N_expected_DY_3;

  }

  if(fname == "histos_DY_2_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in DY + 2 jets background: " << N_expected_DY_2 << endl;
    weight = N_expected_DY_2/totalNumberofEvents;
    cout << "" << endl;
    cout << "Drell Yan + 2 jets weight: " << weight << endl;
    N_expected = N_expected_DY_2;

  }
  
  if(fname == "histos_DY_1_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in DY + 1 jet background: " << N_expected_DY_1 << endl;
    weight = N_expected_DY_1/totalNumberofEvents;
    cout << "" << endl;
    cout << "Drell Yan + 1 jet weight: " << weight << endl;
    N_expected = N_expected_DY_1;

  }
  

  if(fname == "histos_TT_Dileptonic_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in TT Dileptonic background: " << N_expected_TT_Dileptonic << endl;
    weight = N_expected_TT_Dileptonic/totalNumberofEvents;
    cout << "" << endl;
    cout << "TT Dileptonic weight: " << weight << endl;
    N_expected = N_expected_TT_Dileptonic;
    
  }

  if(fname == "histos_TT_Semileptonic_lowm.root"){

    cout << "" << endl;
    cout << "Number of expected events in TT Semileptonic background: " << N_expected_TT_Semileptonic << endl;
    weight = N_expected_TT_Semileptonic/totalNumberofEvents;
    cout << "" << endl;
    cout << "TT Semileptonic weight: " << weight << endl;
    N_expected = N_expected_TT_Semileptonic;
    
  }

  int count_step1(0);
  int count_step2(0);
  int count_step3(0);
  int count_step4(0);

  int matched_bb_pairs(0);
 
  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
       

      // if (jentry>140000) break;

    if(verbose && jentry<10) {
      cout << "Event # " << jentry << endl;
      cout << "" << endl;

    }
    //TLorentzVector pH;
    //TLorentzVector pZ;
    vector<TLorentzVector> vec_AA;
    vector<TLorentzVector> vec_bb1;                  //b quarks from A1 with mc_momidx=6
    vector<TLorentzVector> vec_bb2;                  //b quarks from A2 with mc_momidx=7
    vector<TLorentzVector> vec_ll;                   //antilepton (l1) and lepton (l2) from Z with mc_momidx=2
    //TLorentzVector pb1 , pb2, pb3, pb4;            //b quarks regardless of mc_momidx
    vector<TLorentzVector> pb;

    
    for (Int_t imc=0; imc<nmcparticles; imc++){

      if(verbose && jentry<10) {
	
	cout << "imcparticle " << imc << " : is a " << mc_id[imc] << ", and has a mother at: " << mc_momidx[imc] << "(" << mc_mom[imc] << " )" << "and has a 4-vector p = (" << mc_en[imc] << ", " << mc_px[imc] << ", " << mc_py[imc] << ", " << mc_pz[imc] << " ). Also, its status is:"<< mc_status[imc] << endl;
      }
      
      if(mc_id[imc]==25) { // found the Higgs boson
	
	TLorentzVector pH;
	pH.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);

	float ptH  = pH.Pt();
	float etaH = pH.Eta();
	float phiH = pH.Phi();
	float mH   = pH.M();
	
	h_ptH->Fill(ptH); h_etaH->Fill(etaH); h_phiH->Fill(phiH); h_mH->Fill(mH);
	 
	
      }

      if(mc_id[imc]==23) { //found the Z boson

	TLorentzVector pZ;
	pZ.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);

	
	float ptZ  = pZ.Pt();
	float etaZ = pZ.Eta();
	float phiZ = pZ.Phi();
	float mZ   = pZ.M();

	h_ptZ->Fill(ptZ); h_etaZ->Fill(etaZ); h_phiZ->Fill(phiZ); h_mZ->Fill(mZ);
       
      }
  
       
      if(mc_id[imc]==36 && imc==6){ //found the first A particle

	TLorentzVector pA1;
	pA1.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_AA.push_back(pA1);
	
      }

       if(mc_id[imc]==36 && imc==7){ //found the second A particle
	 
	TLorentzVector pA2;
	pA2.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_AA.push_back(pA2);

       }


       if(mc_id[imc]==5 && imc==8 && mc_mom[imc]==36 && mc_momidx[imc]==6){ //found the first b quark 

	TLorentzVector pb11;
	pb11.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_bb1.push_back(pb11);	
	
       }

       if(mc_id[imc]==-5 && imc==9 && mc_mom[imc]==36 &&  mc_momidx[imc]==6){ //found the second b quark
	 
	 TLorentzVector pb12;
	pb12.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_bb1.push_back(pb12);

       }

       if(mc_id[imc]==5 && imc==10 && mc_mom[imc]==36 && mc_momidx[imc]==7){ //found the third b quark
	 
	 TLorentzVector pb21;
	pb21.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_bb2.push_back(pb21);
	
       } 

       if(mc_id[imc]==-5 && imc==11 && mc_mom[imc]==36 && mc_momidx[imc]==7){ //found the fourth  b quark
	 
	 TLorentzVector pb22;
	pb22.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
	vec_bb2.push_back(pb22);
	
       }


       if(mc_id[imc]==-11 || mc_id[imc]==-13 || mc_id[imc]==-15) { // found the antilepton
	 
	 TLorentzVector pl1;
	pl1.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);
	
	if(pl1.Pt()>20 && fabs(pl1.Eta())<2.5){
	  vec_ll.push_back(pl1);
	}
	
       }

       if(mc_id[imc]==11 || mc_id[imc]==13 || mc_id[imc]==15) { // found the lepton

	 
	TLorentzVector pl2;
	pl2.SetPxPyPzE(mc_px[imc],mc_py[imc],mc_pz[imc],mc_en[imc]);

	if(pl2.Pt()>20 && fabs(pl2.Eta())<2.5){
	vec_ll.push_back(pl2);
	}

       }
 
       TLorentzVector itemp;
       itemp.SetPxPyPzE(mc_px[imc], mc_py[imc], mc_pz[imc], mc_en[imc]);
       /* 
       if(signal) {
	 if(abs(mc_id[imc])==5){
	   if (imc == 8  && itemp.Pt()>20 && fabs(itemp.Eta())<2.5)  pb.push_back(itemp);
	   if (imc == 9  && itemp.Pt()>20 && fabs(itemp.Eta())<2.5)  pb.push_back(itemp);
	   if (imc == 10 && itemp.Pt()>20 && fabs(itemp.Eta())<2.5)  pb.push_back(itemp);
	   if (imc == 11 && itemp.Pt()>20 && fabs(itemp.Eta())<2.5)  pb.push_back(itemp);
	 }
       }
       
       else { //DY background
       */	
       if( ( (abs(mc_id[imc])<=5 && abs(mc_id[imc])>=1) || mc_id[imc]==21) && (mc_status[imc]==1 || mc_status[imc]==23)) {
	   if(itemp.Pt()>20 && fabs(itemp.Eta())<2.5) pb.push_back(itemp);
	   //if(verbose) cout << "Drell-Yan: Found a q/g with id " << mc_id[imc] << " and Pt = " << mc_en[imc] << " and status " <<  mc_status[imc] << endl;
       }
       // }
	 
    } // end mc particle loop

    //q/g multiplicity
    h_nb->Fill(pb.size(), weight);

    //EVENT SELECTION CRITERIA
    
    //start the final state analysis
    if(vec_ll.size()>=2){
      if(pb.size()>=3){
	

	
	if(signal){

	  //kinematics
	  //A bosons:
	  h_ptA1->Fill(vec_AA[0].Pt(), weight); h_etaA1->Fill(vec_AA[0].Eta(), weight); h_phiA1->Fill(vec_AA[0].Phi(), weight); h_mA1->Fill(vec_AA[0].M(), weight);
	  h_ptA2->Fill(vec_AA[1].Pt(), weight); h_etaA2->Fill(vec_AA[1].Eta(), weight); h_phiA2->Fill(vec_AA[1].Phi(), weight); h_mA2->Fill(vec_AA[1].M(), weight); 
	  //leptons:
	  h_pt_l1->Fill(vec_ll[0].Pt(), weight);  h_eta_l1->Fill(vec_ll[0].Eta(), weight); h_phi_l1->Fill(vec_ll[0].Phi(), weight); 
	  h_pt_l2->Fill(vec_ll[1].Pt(), weight);  h_eta_l2->Fill(vec_ll[1].Eta(), weight); h_phi_l2->Fill(vec_ll[1].Phi(), weight);
	  //quarks:
	  h_ptb11->Fill(vec_bb1[0].Pt(), weight); h_etab11->Fill(vec_bb1[0].Eta(), weight); h_phib11->Fill(vec_bb1[0].Phi(), weight); h_mb11->Fill(vec_bb1[0].M(), weight);
	  h_ptb12->Fill(vec_bb1[1].Pt(), weight); h_etab12->Fill(vec_bb1[1].Eta(), weight); h_phib12->Fill(vec_bb1[1].Phi(), weight); h_mb12->Fill(vec_bb1[1].M(), weight);
	  h_ptb21->Fill(vec_bb2[0].Pt(), weight); h_etab21->Fill(vec_bb2[0].Eta(), weight); h_phib21->Fill(vec_bb2[0].Phi(), weight); h_mb21->Fill(vec_bb2[0].M(), weight);
	  h_ptb22->Fill(vec_bb2[1].Pt(), weight); h_etab22->Fill(vec_bb2[1].Eta(), weight); h_phib22->Fill(vec_bb2[1].Phi(), weight); h_mb22->Fill(vec_bb2[1].M(), weight);
	  
	}//end if signal
									
    
	// Define hadronic and leptonic system
	
	TLorentzVector pHad;
	
	for(std::vector<TLorentzVector>::size_type i = 0; i < pb.size(); i++){
	  if(i==4) break;
	  pHad += pb[i];
	}
	
    
	TLorentzVector pLep = vec_ll[0] + vec_ll[1];
    

	// SOME GENERATOR LEVEL KINEMATICS
	//--------------------------------
	
	//calculating deltaR for the pair of A particles
	
	if(signal){
	  float   detaA = fabs(vec_AA[0].Eta()   - vec_AA[1].Eta())  ;
	  float   dphiA = fabs(vec_AA[0].Phi()   - vec_AA[1].Phi())  ;
	  float   dRA   = sqrt(detaA*detaA + dphiA*dphiA);
	  
	  h_dRA->Fill(dRA, weight);
	  
	  
	  
	  //calculating deltaR for the first pair of b quarks: b11 and b12
	  
	  float   deta1 = fabs(vec_bb1[0].Eta()  - vec_bb1[1].Eta()) ;
	  float   dphi1 = fabs(vec_bb1[0].Phi()  - vec_bb1[1].Phi()) ;
	  float   dR1   = sqrt(deta1*deta1 + dphi1*dphi1);
	  
	  h_dR1->Fill(dR1, weight);
	  
	  //calculating deltaR for the second pair of b quarks: b21 and b22
	  
	  float   deta2 = fabs(vec_bb2[0].Eta()  - vec_bb2[1].Eta()) ;
	  float   dphi2 = fabs(vec_bb2[0].Phi()  - vec_bb2[1].Phi()) ;
	  float   dR2   = sqrt(deta2*deta2 + dphi2*dphi2);
	  
	  h_dR2->Fill(dR2, weight);
	  
	  
	  //calculating and storing transverse momenta for the 4 b quarks
	  
	  double ptB_truth[4]={vec_bb1[0].Pt(),vec_bb1[1].Pt(),vec_bb2[0].Pt(),vec_bb2[1].Pt()};
	  
	  //sorting the array in descending order
	  
	  std::sort(ptB_truth, ptB_truth + 4, std::greater<float>()); //now the maximum pT value is in pt[0] and the minimum pT in pt[3]
	  
	 	  
	  //fill histograms with the pt of each b quark in descending order
	  
	  h_pt_1_truth->Fill(ptB_truth[0], weight); 
	  h_pt_2_truth->Fill(ptB_truth[1], weight);
	  h_pt_3_truth->Fill(ptB_truth[2], weight); 
	  h_pt_4_truth->Fill(ptB_truth[3], weight);
	  
	}//end if signal
    

	//sorting the lepton vector based on pT in descending order
	
	std::sort(vec_ll.begin(), vec_ll.end(), [](const TLorentzVector& e, const TLorentzVector& f){
	  return e.Pt() > f.Pt();
	});

	//filling histogrmas with pT of the 2 most energetic leptons
	h_pt_leadL->Fill(vec_ll[0].Pt(), weight);
	h_pt_subleadL->Fill(vec_ll[1].Pt(), weight);
	
	//DR(ll)
	float   deta_ll_G = fabs(vec_ll[0].Eta()   - vec_ll[1].Eta())  ;
	float   dphi_ll_G = fabs(vec_ll[0].Phi()   - vec_ll[1].Phi())  ;
	float   dR_ll_G   = sqrt(deta_ll_G*deta_ll_G + dphi_ll_G*dphi_ll_G);
	h_dR_ll_G->Fill(dR_ll_G, weight);
	
	
	
	// CALCULATE GLOBAL EVENT VARIABLES
	//---------------------------------
	
	
	//calculating the scalar sum of of the pT of the q/g's (HT)
        
	float HT_G = 0;
	for (std::vector<TLorentzVector>::size_type i = 0; i < pb.size(); i++){
	  HT_G += pb[i].Pt();
	}
    
	h_HT_G->Fill(HT_G, weight);

	//calculating the vector sum of the pT of the q/g's
	
	float ptb_v_sum = pHad.Pt();
	h_ptb_v_sum->Fill(ptb_v_sum, weight);

    
	//sorting pb vector based on pT in descending order
	std::sort(pb.begin(), pb.end(), [](const TLorentzVector& g, const TLorentzVector& h) {
	  return g.Pt() > h.Pt();
	});

	//filling histograms with pT of the most energetic q/g's
	
	h_pt_1_global->Fill(pb[0].Pt(), weight);
	h_pt_2_global->Fill(pb[1].Pt(), weight);
	h_pt_3_global->Fill(pb[2].Pt(), weight);
	if(pb.size()>3)  h_pt_4_global->Fill(pb[3].Pt(), weight);
      
	//calculating DeltaPhi(hadronic,leptonic) and pHad.M(), pLep.M()

	float phi_H_G = pHad.Phi(); h_phi_H_G->Fill(phi_H_G);
	float phi_Z_G = pLep.Phi(); h_phi_Z_G->Fill(phi_Z_G);
	float dPhi_ZH_G = fabs(pHad.Phi()-pLep.Phi());
	if(dPhi_ZH_G >  M_PI) dPhi_ZH_G = 2*M_PI - dPhi_ZH_G;
	h_dPhi_ZH_G-> Fill(dPhi_ZH_G, weight);
	
	float m_H_G  = pHad.M();  h_m_H_G->Fill(m_H_G, weight);
	float m_Z_G  = pLep.M();  h_m_Z_G->Fill(m_Z_G, weight);
	
	float pt_Z_G = pLep.Pt(); h_pt_Z_G->Fill(pt_Z_G, weight);
	
	// if(pb.size()>=4){
	// //calculating the minimum and maximum difference between 2 b quark pair masses
	
	// float dM[3];
    
	// //6 combinations of b quarks: b1+b2, b1+b3, b1+b4, b2+b3, b2+b4, b3+b4 and 3 combinations of bb pairs
	 
	// dM[0] = fabs((pb[0]+pb[1]).M()-(pb[3]+pb[2]).M());
	// dM[1] = fabs((pb[0]+pb[2]).M()-(pb[1]+pb[3]).M());
	// dM[2] = fabs((pb[0]+pb[3]).M()-(pb[1]+pb[2]).M());  
 
	// //sorting the array in descending order

	// std::sort(dM, dM + 3, std::greater<float>()); //now the maximum dM value is in dM[0] and the minimum dM value is in dM[2]

	// h_dM_min->Fill(dM[2]);
	// h_dM_max->Fill(dM[0]);
	// }
    
	if(signal){
	  //pT_2b - m_2b plane truth
	  float m_2b_A1  = (vec_bb1[0]+vec_bb1[1]).M();
	  float m_2b_A2  = (vec_bb2[0]+vec_bb2[1]).M();
	  
	  float pT_2b_truth_vector = (vec_bb1[0] + vec_bb1[1]).Pt();
	  float pT_2b_truth_vector_other = (vec_bb2[0] + vec_bb2[1]).Pt();
	  
	  h_pt2b_m2b_truth_vector->Fill(m_2b_A1, pT_2b_truth_vector, weight);
	  
	  h_pt2b_m2b_truth_vector_other->Fill(m_2b_A2, pT_2b_truth_vector_other, weight);
	  
	  h_pt2b_truth_vector->Fill(pT_2b_truth_vector, weight);
	  
	  h_pt2b_truth_vector_other->Fill(pT_2b_truth_vector_other, weight);
	}
   
	//minimum dR analysis
    
	//6 possible pairs of b quarks: (1):b1b2, (2):b1b3, (3):b1b4, (4):b2b3, (5):b2b4, (6):b3b4
  
	float dR_min_1 = 999;
	TLorentzVector b1_min, b2_min;

	for (std::vector<TLorentzVector>::size_type i = 0; i < pb.size(); ++i) {
	  TLorentzVector b1 = pb[i];
	  
	  for (std::vector<TLorentzVector>::size_type j = i + 1; j < pb.size(); ++j) {
	    
	    TLorentzVector b2 = pb[j];
	    
	    float deta = fabs(b1.Eta() - b2.Eta());
	    float dphi = fabs(b1.Phi() - b2.Phi());
	    float dR = sqrt(deta * deta + dphi * dphi);
	    
	    
	    if (dR < dR_min_1) {
	      dR_min_1 = dR;
	      b1_min   = b1;
	      b2_min   = b2;
	      
	    }
	    
	  }
	}


	//MIN DR histograms
    
	float m_2b_min_dR = (b1_min + b2_min).M();
	float pT_2b_min_dR_vector = (b1_min + b2_min).Pt();

	float deta_min_G = fabs(b1_min.Eta() - b2_min.Eta());
	float dphi_min_G = fabs(b1_min.Phi() - b2_min.Phi());
	float dR_min_G   = sqrt(deta_min_G * deta_min_G + dphi_min_G * dphi_min_G);
	
	
        h_pt2b_m2b_mindR_vector->Fill(m_2b_min_dR, pT_2b_min_dR_vector, weight);
	
	h_dR_min_G->Fill(dR_min_G, weight);
	h_pT_min_vector->Fill(pT_2b_min_dR_vector, weight);
	h_m2b_mindR->Fill(m_2b_min_dR, weight);
	
	// Check if 2nd pair exists:
	vector<TLorentzVector> b_other;
     
	for (std::vector<TLorentzVector>::size_type i = 0; i<pb.size(); ++i){
	  
	  TLorentzVector b = pb[i];
	  
	  if(b!=b1_min && b!=b2_min){
	    b_other.push_back(pb[i]);
	  }
	}
	
	
	if(b_other.size()>=2) { // if there is a 2nd pair (ie total 4 quarks/gluons)
       
	  TLorentzVector  b1_other = b_other[0];
	  TLorentzVector  b2_other = b_other[1];
	  
	  //OTHER
	  float m_2b_other  = (b1_other + b2_other).M();
	  float pT_2b_other_vector = (b1_other + b2_other).Pt();
	  
	  h_pt2b_m2b_other_vector -> Fill(m_2b_other, pT_2b_other_vector, weight);
	  
	  //correlations
	  h_m1_m2->Fill(m_2b_min_dR, m_2b_other, weight);
	  h_pT1_pT2_vector->Fill(pT_2b_min_dR_vector, pT_2b_other_vector, weight);
	  
	  
	  float deta_other_G = fabs(b1_other.Eta() - b2_other.Eta());
	  float dphi_other_G = fabs(b1_other.Phi() - b2_other.Phi());
	  float dR_other_G   = sqrt(deta_other_G * deta_other_G + dphi_other_G * dphi_other_G);
	  
	  h_dRmin_dRother -> Fill(dR_min_G, dR_other_G, weight);
	  
	  h_dR_other_G->Fill(dR_other_G, weight);
	  
	  h_m2b_other->Fill(m_2b_other, weight);
	  
	  h_pT_other_vector->Fill(pT_2b_other_vector, weight);
	  
	  //new mass variables
	  
	  float dM_A1_A2_gen_level = m_2b_min_dR - m_2b_other;
	  h_dM_A1_A2_gen_level -> Fill(dM_A1_A2_gen_level, weight);
	  
 
	} // end if 2nd pair exists

      } //end q/g requirement
      
    } //end lepton requirement
     

    //START DETECTOR LEVEL ANALYSIS
    
     // Jets
    std::vector<TLorentzVector> vec_jet;
    std::vector<TLorentzVector> vec_bjet;
  
    // Leptons
    std::vector<TLorentzVector> vec_muons;
    std::vector<TLorentzVector> vec_ele;
    std::vector<TLorentzVector> vec_leptons;

    //muon configuration

    for (int i = 0; i < mn; i++)
    {
      TLorentzVector p_muon;
      p_muon.SetPxPyPzE(mn_px[i], mn_py[i], mn_pz[i], mn_en[i]);

      if (p_muon.Pt() < 20. || std::fabs(p_muon.Eta()) > 2.4)
        continue;

      // id + Isolation
      if (mn_passId[i] && mn_passIso[i])
      {
        vec_muons.push_back(p_muon);
	vec_leptons.push_back(p_muon);
      }
    }
    
    //electron configuration

    for (int i = 0; i < en; i++)
    {
      TLorentzVector p_electron;
      p_electron.SetPxPyPzE(en_px[i], en_py[i], en_pz[i], en_en[i]);

      if (p_electron.Pt() < 20. || std::fabs(p_electron.Eta()) > 2.4)
        continue;

      if (en_passId[i] && en_passIso[i])
      {
        vec_ele.push_back(p_electron);
        vec_leptons.push_back(p_electron);
      }
    }

    //vector that holds the TLorentzVector and the b tag disciminator
    std::vector<JetAndBtag> vec_jet_and_btag;

    //jets & cross cleaning
    
    for (int i = 0; i < jet; i++)
    {
      bool overlap = false;
      
      TLorentzVector p_jet;
      p_jet.SetPxPyPzE(jet_px[i], jet_py[i], jet_pz[i], jet_en[i]);

      if (p_jet.Pt() < 20. || std::fabs(p_jet.Eta()) > 2.5) continue;
	  
      for (std::vector<TLorentzVector>::size_type mn_count = 0; mn_count < vec_muons.size(); mn_count++)
      {
        float dR1_jet = getDeltaR(p_jet, vec_muons[mn_count]);
	//   h_dR_jet_mn_before->Fill(dR1_jet);
        if (dR1_jet < 0.4){
          overlap = true;
         }
      }

      for (std::vector<TLorentzVector>::size_type en_count = 0; en_count < vec_ele.size(); en_count++)
      {
        float dR1_jet = getDeltaR(p_jet, vec_ele[en_count]);
	//     h_dR_jet_en_before->Fill(dR1_jet);
        if (dR1_jet < 0.4){
         overlap = true;
        }
      }

     
      if (! overlap) {
      
      JetAndBtag jet_info; //a JetAndBtag structure
      jet_info.jet = p_jet; //assign the TLV p_jet to the jet member of the structure
      jet_info.btag = jet_btag1[i]; //assign the b tag discriminator jet_btag1 to the btag member of the structure
      vec_jet_and_btag.push_back(jet_info);
      
      h_btag_discriminator->Fill(jet_btag1[i], weight);
      
      if (jet_btag1[i] > 0.4941) {
	vec_bjet.push_back(p_jet); // add b-tagging	

      }

      }
        
    } //end jet loop


     
    //sorting b jet vector based on pT in descending order
    
    std::sort(vec_bjet.begin(), vec_bjet.end(), [](const TLorentzVector& a, const TLorentzVector& b){
      return a.Pt() > b.Pt();
    });

    //sorting the lepton vector based on pT in descending order

    std::sort(vec_leptons.begin(), vec_leptons.end(), [](const TLorentzVector& c, const TLorentzVector& d){
      return c.Pt() > d.Pt();
    });

    //sort the b-tagged jets by discriminator value   
    std::sort(vec_jet_and_btag.begin(), vec_jet_and_btag.end(), sortBtag);

    //fill the vec_jet
    for (std::vector<TLorentzVector>::size_type i=0; i<vec_jet_and_btag.size(); i++) {
      vec_jet.push_back(vec_jet_and_btag[i].jet);
    }

    
    //selection criteria
    if(vec_leptons.size()<2) continue; //at least 2 leptons
    count_step1++;
    
    float mll=(vec_leptons[0]+vec_leptons[1]).M();
    
    h_mll_before_cut->Fill(mll, weight); //dilepton mass before applying the Z mass window
    
    if (mll > 100. || mll < 80.) continue; //Z mass window
    count_step2++;
    
    h_n_jets->Fill(vec_jet.size(), weight);

    h_n_bjets->Fill(vec_bjet.size(), weight);
    
    if(vec_jet.size()<2) continue; //at least 2 jets
    count_step3++;
    
    if(vec_bjet.size()<2) continue; //at least 2 b jets
    count_step4++;
    
    n_jets_after_cuts = vec_jet.size(); h_n_jets_after_cuts->Fill(n_jets_after_cuts, weight);

    h_n_bjets_after_cuts->Fill(vec_bjet.size(), weight);

    pt_b1 = vec_bjet[0].Pt(); h_pt_b1->Fill(pt_b1, weight);
    
    //filling histograms with b tag discriminator
    btag_1 = vec_jet_and_btag[0].btag; h_btag1->Fill(btag_1, weight);
    btag_2 = vec_jet_and_btag[1].btag; h_btag2->Fill(btag_2, weight);
    
    float btag_3 = vec_jet_and_btag[2].btag; h_btag3->Fill(btag_3, weight);
    float btag_4 = vec_jet_and_btag[3].btag; h_btag4->Fill(btag_4, weight);

    TLorentzVector jet1 = vec_jet[0], jet2 = vec_jet[1];

    
    //filling histograms with pT of the highest b tagged jets
    h_pt_jet1->Fill(jet1.Pt(), weight); h_phi_jet1->Fill(jet1.Phi(), weight); h_eta_jet1->Fill(jet1.Eta(), weight); h_m_jet1->Fill(jet1.M(), weight);
    h_pt_jet2->Fill(jet2.Pt(), weight); h_phi_jet2->Fill(jet2.Phi(), weight); h_eta_jet2->Fill(jet2.Eta(), weight); h_m_jet2->Fill(jet2.M(), weight);
    if(vec_jet.size()>2) h_pt_jet3->Fill(vec_jet[2].Pt(), weight); h_phi_jet3->Fill(vec_jet[2].Phi(), weight); h_eta_jet3->Fill(vec_jet[2].Eta(), weight); h_m_jet3->Fill(vec_jet[2].M(), weight);
    if(vec_jet.size()>3) h_pt_jet4->Fill(vec_jet[3].Pt(), weight); h_phi_jet4->Fill(vec_jet[3].Phi(), weight); h_eta_jet4->Fill(vec_jet[3].Eta(), weight); h_m_jet4->Fill(vec_jet[3].M(), weight);
    
    //filling histograms with pT of the 2 most energetic leptons
    h_pt_lepton1->Fill(vec_leptons[0].Pt(), weight);
    h_pt_lepton2->Fill(vec_leptons[1].Pt(), weight);

    //deltaR between the leading and subleading lepton
    
    float   deta_ll = fabs(vec_leptons[0].Eta()   - vec_leptons[1].Eta())  ;
    float   dphi_ll = fabs(vec_leptons[0].Phi()   - vec_leptons[1].Phi())  ;

    dR_ll   = sqrt(deta_ll*deta_ll + dphi_ll*dphi_ll); h_dR_ll->Fill(dR_ll, weight);
    
    
    //Hadronic system
    
    TLorentzVector pHadronic;
    for(std::vector<TLorentzVector>::size_type i=0; i<vec_jet.size(); i++){
      if(i==4) break;
      pHadronic += vec_jet[i];
    }

    //Leptonic System
    
    TLorentzVector pLeptonic = vec_leptons[0] + vec_leptons[1];
   
    
    //Hadronic system pT (vector sum of the pT's of the jets) and mass
    
    pt_H = pHadronic.Pt();
    h_pt_H->Fill(pt_H, weight); 
    m_H = pHadronic.M();
    h_m_H->Fill(m_H, weight); 
    
    //Leptonic system pT and mass
    
    pt_Z = pLeptonic.Pt();
    h_pt_Z->Fill(pt_Z, weight);
    
    float m_Z = pLeptonic.M();
    h_m_Z->Fill(m_Z, weight); 

    float eta_H = pHadronic.Eta(); h_eta_H->Fill(eta_H);
    float eta_Z = pLeptonic.Eta(); h_eta_Z->Fill(eta_Z);

  
    float phi_H = pHadronic.Phi(); h_phi_H->Fill(phi_H); 
    float phi_Z = pLeptonic.Phi(); h_phi_Z->Fill(phi_Z); 

    dPhi_ZH = fabs(pHadronic.Phi()-pLeptonic.Phi());

    if(dPhi_ZH >  TMath::Pi()) dPhi_ZH = 2*TMath::Pi() - dPhi_ZH;

    h_dPhi_ZH-> Fill(dPhi_ZH, weight); 

  
    //scalar sum of the pT of the jets
    
    HT = 0;
    for(std::vector<TLorentzVector>::size_type i=0; i<vec_jet.size(); i++){
      HT += vec_jet[i].Pt();
    }
    h_HT->Fill(HT, weight);

    //pt(2b)
    pt_2b = jet1.Pt();
    h_pt_2b->Fill(pt_2b);
    // mass variable
    dM_A1_A2 = jet1.M() - jet2.M();
    h_dM_A1_A2->Fill(dM_A1_A2);
   
    // MET

    h_met_pt->Fill(met_pt,weight);
    met_pt_ = met_pt;

    // truth matched bb pairs


    if(signal){
      
      bool match_jet1_bb1 = jet_matched_lowm(jet1, vec_bb1);
      bool match_jet1_bb2 = jet_matched_lowm(jet1, vec_bb2);


      if(match_jet1_bb1 || match_jet1_bb2){
	matched_bb_pairs++;	
      }
      
      float jet1_pt = jet1.Pt();
      if(match_jet1_bb1){
        float pt_2b_truth = (vec_bb1[0] + vec_bb1[1]).Pt();
	h_ptjet_pt2b->Fill(jet1_pt, pt_2b_truth);	
      }
      else if(match_jet1_bb2){
	float pt_2b_truth = (vec_bb2[0] + vec_bb2[1]).Pt();
	h_ptjet_pt2b->Fill(jet1_pt, pt_2b_truth);
      }



    } // end if(signal)


    

        
    t1.Fill();
    
   }//end event loop

   cout << "TOTAL: " << totalNumberofEvents <<endl;
   cout << "" << endl;
   cout << "step1: " << count_step1 << endl;
   cout << "" << endl;
   cout << "step2: " << count_step2 << endl;
   cout << "" << endl;
   cout << "step3: " << count_step3 << endl;
   cout << "" << endl;
   cout << "step4: " << count_step4 << endl;
   cout << "" << endl;
   
   cout << "eff0 = " << float(totalNumberofEvents)/float(totalNumberofEvents) << endl;
   cout << "" << endl;
   cout << "eff1 = " << int(count_step1)/float(totalNumberofEvents) << endl; //step 1: at least 2 leptons
   cout << "" << endl;
   cout << "eff2 = " << int(count_step2)/float(totalNumberofEvents) << endl; //step 2: Z mass window
   cout << "" << endl;
   cout << "eff3 = " << int(count_step3)/float(totalNumberofEvents) << endl; //step 3: at least 3 jets
   cout << "" << endl;
   cout << "eff4 = " << int(count_step4)/float(totalNumberofEvents) << endl; //step 4: at least 3 b jets
   cout << "" << endl;
   cout << "N expected (final): " << float(N_expected)*(int(count_step4)/float(totalNumberofEvents)) << endl;
   cout << "" << endl;

   // fraction of truth-matched bb pairs

   double fraction_matched = static_cast<double>(matched_bb_pairs)/count_step4;
   cout << "fraction of correct truth - matched bb pairs:" << (fraction_matched)*100 << "%" << endl;

   t1.Write();
   
   // TFile f1("histos.root","RECREATE"); 
   h_ptH->Write();   h_etaH->Write();   h_phiH->Write();   h_mH->Write();
   h_ptZ->Write();   h_etaZ->Write();   h_phiZ->Write();   h_mZ->Write();
   h_ptA1->Write();  h_etaA1->Write();  h_phiA1->Write();  h_mA1->Write();
   h_ptA2->Write();  h_etaA2->Write();  h_phiA2->Write();  h_mA2->Write();
   h_ptb11->Write(); h_etab11->Write(); h_phib11->Write(); h_mb11->Write();
   h_ptb12->Write(); h_etab12->Write(); h_phib12->Write(); h_mb12->Write();
   h_ptb21->Write(); h_etab21->Write(); h_phib21->Write(); h_mb21->Write();
   h_ptb22->Write(); h_etab22->Write(); h_phib22->Write(); h_mb22->Write();
   h_pt_l1->Write(); h_eta_l1->Write(); h_phi_l1->Write(); //h_m_l1->Write();
   h_pt_l2->Write(); h_eta_l2->Write(); h_phi_l2->Write(); //h_m_l2->Write();
   h_dR1->Write();   h_dR2->Write();    h_dRA->Write();   
   
   
   h_HT_G->Write();h_ptb_v_sum->Write();
   h_dPhi_ZH_G->Write();
   h_dM_min->Write(); h_dM_max->Write(); h_pt_Z_G->Write();
   h_m_H_G->Write();  h_m_Z_G->Write();
   h_pt_1_truth->Write(); h_pt_2_truth->Write(); h_pt_3_truth->Write(); h_pt_4_truth->Write();
   h_pt_1_global->Write(); h_pt_2_global->Write(); h_pt_3_global->Write(); h_pt_4_global->Write();
   h_pt_leadL->Write(); h_pt_subleadL->Write(); h_dR_ll_G->Write();
    
   //pT(2b)-m(2b) plane
   h_pt2b_m2b_truth_vector->Write();
   h_pt2b_m2b_truth_vector_other->Write();
   h_pt2b_truth_vector->Write();
   h_pt2b_truth_vector_other->Write();
   
   //min dR analysis
   
   h_pt2b_m2b_mindR_vector->Write();
   h_pt2b_m2b_other_vector->Write();
   h_m1_m2->Write();
   h_pT1_pT2_vector->Write();
   h_dRmin_dRother->Write();
   h_dR_min_G->Write();
   h_dR_other_G->Write();
   h_m2b_mindR->Write(); h_m2b_other->Write();
   h_pT_min_vector->Write();
   h_pT_other_vector->Write();
   
   
   h_nb->Write();
   h_phi_H_G->Write();
   h_phi_Z_G->Write();
   
   h_dM_A1_A2_gen_level->Write();
   
   
   //detector level
   h_n_jets->Write(); h_n_bjets->Write(); h_n_jets_after_cuts->Write(); h_n_bjets_after_cuts->Write();
   h_btag_discriminator->Write();
   h_pt_H->Write(); h_m_H->Write(); h_phi_H->Write(); h_eta_H->Write();
   h_pt_Z->Write(); h_m_Z->Write(); h_phi_Z->Write(); h_eta_Z->Write();
   h_dPhi_ZH->Write();
   
   h_HT->Write();
   h_pt_jet1->Write(); h_eta_jet1->Write(); h_phi_jet1->Write(); h_m_jet1->Write();
   h_pt_jet2->Write(); h_eta_jet2->Write(); h_phi_jet2->Write(); h_m_jet2->Write();
   h_pt_jet3->Write(); h_eta_jet3->Write(); h_phi_jet3->Write(); h_m_jet3->Write();
   h_pt_jet4->Write(); h_eta_jet4->Write(); h_phi_jet4->Write(); h_m_jet4->Write();
   
   
   h_pt_lepton1->Write();
   h_pt_lepton2->Write();
   h_pt_2b->Write();

  

   h_dR_ll->Write();
   
   h_btag1->Write(); h_btag2->Write(); h_btag3->Write(); h_btag4->Write();

   h_dM_A1_A2->Write();
   

   h_pt_b1->Write();

   h_mll_before_cut->Write();

   h_met_pt->Write();

   h_best_matched_b_dR->Write();
   h_deltaR_correct_pair->Write();

   h_ptjet_pt2b->Write();

    fout.Close();

    
   
}//end loop method
