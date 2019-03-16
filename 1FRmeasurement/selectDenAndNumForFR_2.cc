#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include "TMath.h"
#include "TVector2.h"
#include "/u/user/joonblee/Headers/Object.h"
#include "/u/user/joonblee/Headers/NtupleHandle.h"

#include <iostream>

#include "./UsefulFunctions.C"

#include "/u/user/joonblee/Headers/RoCorr/RoccoR.cc"
#include "/u/user/joonblee/Headers/effSF_muon/DYAnalyzer.h"

#define trig "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
#define trig_ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"

using namespace std;

void dimuon_selection(vector<Muon> MuonCollection, vector<Muon> *add_dimuBox, bool *add_test_Zboson, int *add_a, int *add_b) {
  for(unsigned j=0; j<MuonCollection.size(); j++) {
    Muon mu = MuonCollection[j];
    if( mu.acceptance(20,2.4)
				&& mu.isGlobalMuon
        && mu.testMuonIDv2()
        && second_isolation( MuonCollection, mu ) < 0.1 ) {
      for(unsigned k=j+1; k<MuonCollection.size(); k++) {
        Muon mu_ = MuonCollection[k];
        if(mu.charge + mu_.charge != 0) continue;
        if( mu_.acceptance(10,2.4)
						&& mu_.isGlobalMuon
            && mu_.testMuonIDv2()
            && second_isolation( MuonCollection, mu_ ) < 0.1 ) {
          if( fabs((mu.Momentum + mu_.Momentum).M() - 91.2) < 7 ) {
            add_dimuBox->push_back(mu);
            add_dimuBox->push_back(mu_);
            *add_test_Zboson = true;
            *add_a = j; *add_b = k;
            break;
          }
        }
      }
      if(*add_test_Zboson) break;
    }
  }
}

void probe_muon_selection(vector<Muon> MuonCollection, vector<Muon> *add_probeBox, int a, int b ) {
  for(unsigned j=0; j!=MuonCollection.size(); j++) {
    if( j==a || j==b ) continue;
    Muon mu = MuonCollection[j];
    if( mu.acceptance(5,2.4)
        && mu.looseMuonID() ) {
      add_probeBox->push_back(mu);
      if(add_probeBox->size() >= 2) break;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////

TString dataset_name(TString dataSet, int input, int *add_index) {
  TString str_chain;
  TString str = to_string(input);

  if(dataSet == "dataB") {
    TString files;
    if(input == 0) files = "[01234]"; else if(input == 1) files = "*";
    else {
      files = "[56789]";
      str = "0";
    }
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunB/180227_010024/000"+str+"/ntuple_skim_*"+files+".root";
    *add_index += input;
  } //2
  else if(dataSet == "dataC") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunC/180227_010110/0000/ntuple_skim_*.root";
    *add_index += 10;
  }
  else if(dataSet == "dataD") {
    TString files;
    if(input == 0) files = "[01234]"; else files = "[56789]";
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunD/180227_010154/0000/ntuple_skim_*"+files+".root";
    *add_index += 20 + input;
  } //1
  else if(dataSet == "dataE") {
    TString files;
    if(input == 0) files = "[01234]"; else files = "[56789]";
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunE/180227_010605/0000/ntuple_skim_*"+files+".root";
    *add_index += 30 + input;
  } // 1
  else if(dataSet == "dataF") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunF/180227_010650/0000/ntuple_skim_*.root";
    *add_index += 40;
  }
  else if(dataSet == "dataG") {
    TString files;
    if(input == 0) files = "[01234]"; else if(input == 1) files = "*";
    else {
      files = "[56789]";
      str = "0";
    }
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunG/180227_010737/000"+str+"/ntuple_skim_*"+files+".root";
    *add_index += 50 + input;
  } //1
  else if(dataSet == "dataHv2") {
    TString files;
    if(input == 0) files = "[01234]"; else if(input == 1) files = "*";
    else {
      files = "[56789]";
      str = "0";
    }
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunHver2/180227_011449/000"+str+"/ntuple_skim_*"+files+".root";
    *add_index += 60 + input;
  }
  else if(dataSet == "dataHv3") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunHver3/180227_014825/0000/ntuple_skim_*.root";
    *add_index += 70;
  }

  /*
  // Signal //
  else if(dataSet=="H200A50") str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/HAA4Mu_Modified/HAA4Mu_Modified_H200A50_ntuple/170515_083815/0000/ntuple_skim_*.root";
  else if(dataSet=="H800A1") str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/HAA4Mu_Modified/HAA4Mu_Modified_H800A1_ntuple/170515_083652/0000/ntuple_skim_*.root";
  else if(dataSet=="H800A50") str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/HAA4Mu_Modified/HAA4Mu_Modified_H800A50_ntuple/170515_083722/0000/ntuple_skim_*.root";
  else if(dataSet=="H2000A1") str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/HAA4Mu_Modified/HAA4Mu_Modified_H2000A1_ntuple/170515_083559/0000/ntuple_skim_*.root";
  else if(dataSet=="H2000A50") str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/HAA4Mu_Modified/HAA4Mu_Modified_H2000A50_ntuple/170515_083625/0000/ntuple_skim_*.root";
  */

  // BKG MC //
  else if(dataSet == "qqToZZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ZZTo4L_13TeV_powheg_pythia8/crab_DYntuple_v20171022_EGMCorr_ZZTo4L/171230_140506/0000/ntuple_skim_*.root"; // qq->ZZ
    *add_index = 0;
  }
  else if(dataSet == "GGToZZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/crab_DYntuple_v20171022_EGMCorr_GGToZZTo4mu/171230_140643/0000/ntuple_skim_*.root"; // GG->ZZ
    *add_index = 1;
  }
  else if(dataSet == "DY") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v2p0_DYLL_M50toInf/171230_124754/0000/ntuple_skim_*"+str+".root"; // DY
    *add_index = 10+input;
  }


  else if(dataSet == "ZG") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v20171022_EGMCorr_ZGamma/180311_151437/0000/ntuple_skim_*.root"; // TT
    *add_index = 40+input;
  }


  else if(dataSet == "WZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/WZ_TuneCUETP8M1_13TeV-pythia8/crab_DYntuple_v20171022_EGMCorr_WZ/171230_140159/0000/ntuple_skim_*.root"; // WZ
    *add_index = 7;
  }
  else if(dataSet == "TT") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_DYntuple_v20171022_EGMCorr_TTbar/171230_135613/0000/ntuple_skim_*"+str+".root"; // TT
    *add_index = 20+input;
  }
  else if(dataSet == "ttZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_DYntuple_v20171022_EGMCorr_TTZToLLNuNu/171230_141418/0000/ntuple_skim_*.root"; // ttZ
    *add_index = 3;
  }
  else if(dataSet == "WWZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_DYntuple_v20171022_EGMCorr_WWZ/171230_141200/0000/ntuple_skim_*.root"; // WWZ
    *add_index = 4;
  }
  else if(dataSet == "WZZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_DYntuple_v20171022_EGMCorr_WZZ/171230_141244/0000/ntuple_skim_*.root"; // WZZ
    *add_index = 5;
  }
  else if(dataSet == "ZZZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_DYntuple_v20171022_EGMCorr_ZZZ/171230_141329/0000/ntuple_skim_*.root"; // ZZZ
    *add_index = 6;
  }
  else if(dataSet == "tW" && input == 0) { // tW
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_DYntuple_v20171022_EGMCorr_ST_tbarW/171230_140906/0000/ntuple_skim_*.root";
    *add_index = 8;
  }
  else if(dataSet == "tW" && input == 1) {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_DYntuple_v20171022_EGMCorr_ST_tW/171230_140812/0000/ntuple_skim_*.root";
    *add_index = 9;
  }
  else if(dataSet == "higgs") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/crab_DYntuple_v20171022_EGMCorr_GGHToZZTo4L/171230_140727/0000/ntuple*.root"; // higgs
    *add_index = 2;
  }
  else {
    cout<<"wrong sample"<<endl;
    str_chain = "ERROR";
  }

  return str_chain;
}



void selectDenAndNumForFR_2(TString dataSet, int input = 0) {
  clock_t begin_t, end_t;
  begin_t = clock();
  time();

  TString str = to_string(input);

  int index = 100;

  TChain *chain = new TChain("recoTree/DYTree");

  TString str_chain = dataset_name(dataSet, input, &index);
  chain->Add(str_chain);

  NtupleHandle *ntuple = new NtupleHandle( chain );
  ntuple->TurnOnBranches_Muon();
  ntuple->TurnOnBranches_MET();
  if(dataSet.Index("data")!=0)
    ntuple->TurnOnBranches_GenLepton();

  DYAnalyzer *analyzer = new DYAnalyzer();

  if( dataSet.Index("data")<0 ) {
    analyzer->SetupEfficiencyScaleFactor_BtoF();
		analyzer->SetupEfficiencyScaleFactor_GtoH();
    analyzer->SetupEfficiencyScaleFactorLoose_BtoF();
    analyzer->SetupEfficiencyScaleFactorLoose_GtoH();
	}

  TRandom3 *r1 = new TRandom3(0);
  RoccoR rc("/u/user/joonblee/Headers/RoCorr/rcdata.2016.v3");

  int entries = chain->GetEntries();
	//entries = 1000;

  // ---------- HISTOGRAM ---------- //
  TH1D* denominator_pt = new TH1D("denominator_pt","",500,0,500);
  TH1D* denominator_pt_barrel = new TH1D("denominator_pt_barrel","",500,0,500);
  TH1D* denominator_pt_endcap = new TH1D("denominator_pt_endcap","",500,0,500);

  TH1D* numerator_pt = new TH1D("numerator_pt","",500,0,500);
  TH1D* numerator_pt_barrel = new TH1D("numerator_pt_barrel","",500,0,500);
  TH1D* numerator_pt_endcap = new TH1D("numerator_pt_endcap","",500,0,500);

  denominator_pt->Sumw2();
  denominator_pt_barrel->Sumw2();
  denominator_pt_endcap->Sumw2();

  numerator_pt->Sumw2();
  numerator_pt_barrel->Sumw2();
  numerator_pt_endcap->Sumw2();

  TH1D* denominator_FourMu_pt = new TH1D("denominator_FourMu_pt","",500,0,500);
  TH1D* denominator_FourMu_pt_barrel = new TH1D("denominator_FourMu_pt_barrel","",500,0,500);
  TH1D* denominator_FourMu_pt_endcap = new TH1D("denominator_FourMu_pt_endcap","",500,0,500);

  TH1D* numerator_FourMu_pt = new TH1D("numerator_FourMu_pt","",500,0,500);
  TH1D* numerator_FourMu_pt_barrel = new TH1D("numerator_FourMu_pt_barrel","",500,0,500);
  TH1D* numerator_FourMu_pt_endcap = new TH1D("numerator_FourMu_pt_endcap","",500,0,500);

  denominator_FourMu_pt->Sumw2();
  denominator_FourMu_pt_barrel->Sumw2();
  denominator_FourMu_pt_endcap->Sumw2();

  numerator_FourMu_pt->Sumw2();
  numerator_FourMu_pt_barrel->Sumw2();
  numerator_FourMu_pt_endcap->Sumw2();

  TH1D* denominator_DiMu_pt = new TH1D("denominator_DiMu_pt","",500,0,500);
  TH1D* denominator_DiMu_pt_barrel = new TH1D("denominator_DiMu_pt_barrel","",500,0,500);
  TH1D* denominator_DiMu_pt_endcap = new TH1D("denominator_DiMu_pt_endcap","",500,0,500);

  TH1D* numerator_DiMu_pt = new TH1D("numerator_DiMu_pt","",500,0,500);
  TH1D* numerator_DiMu_pt_barrel = new TH1D("numerator_DiMu_pt_barrel","",500,0,500);
  TH1D* numerator_DiMu_pt_endcap = new TH1D("numerator_DiMu_pt_endcap","",500,0,500);

  denominator_DiMu_pt->Sumw2();
  denominator_DiMu_pt_barrel->Sumw2();
  denominator_DiMu_pt_endcap->Sumw2();

  numerator_DiMu_pt->Sumw2();
  numerator_DiMu_pt_barrel->Sumw2();
  numerator_DiMu_pt_endcap->Sumw2();

  TH1D* M_loose = new TH1D("M_loose","",1000,50,150); M_loose->Sumw2();

  TH1D* looseMu_pt               = new TH1D("looseMu_pt","",100,0,100);
  looseMu_pt->Sumw2();
  TH1D* looseMu_eta              = new TH1D("looseMu_eta","",480,-2.4,2.4);
  looseMu_eta->Sumw2();
  TH1D* looseMu_phi              = new TH1D("looseMu_phi","",640,-3.2,3.2);
  looseMu_phi->Sumw2();

  TH1D* M_tight = new TH1D("M_tight","",1000,50,150); M_tight->Sumw2();

  TH1D* tightMu_pt               = new TH1D("tightMu_pt","",100,0,100);
  tightMu_pt->Sumw2();
  TH1D* tightMu_eta              = new TH1D("tightMu_eta","",480,-2.4,2.4);
  tightMu_eta->Sumw2();
  TH1D* tightMu_phi              = new TH1D("tightMu_phi","",640,-3.2,3.2);

  TH1D* M_FourMu_loose = new TH1D("M_FourMu_loose","",1000,50,150); M_FourMu_loose->Sumw2();

  TH1D* looseMu_FourMu_pt               = new TH1D("looseMu_FourMu_pt","",100,0,100);
  looseMu_FourMu_pt->Sumw2();
  TH1D* looseMu_FourMu_eta              = new TH1D("looseMu_FourMu_eta","",480,-2.4,2.4);
  looseMu_FourMu_eta->Sumw2();
  TH1D* looseMu_FourMu_phi              = new TH1D("looseMu_FourMu_phi","",640,-3.2,3.2);
  looseMu_FourMu_phi->Sumw2();

  TH1D* M_FourMu_tight = new TH1D("M_FourMu_tight","",1000,50,150); M_FourMu_tight->Sumw2();

  TH1D* tightMu_FourMu_pt               = new TH1D("tightMu_FourMu_pt","",100,0,100);
  tightMu_FourMu_pt->Sumw2();
  TH1D* tightMu_FourMu_eta              = new TH1D("tightMu_FourMu_eta","",480,-2.4,2.4);
  tightMu_FourMu_eta->Sumw2();
  TH1D* tightMu_FourMu_phi              = new TH1D("tightMu_FourMu_phi","",640,-3.2,3.2);

  TH1D* M_DiMu_loose = new TH1D("M_DiMu_loose","",1000,50,150); M_DiMu_loose->Sumw2();

  TH1D* looseMu_DiMu_pt               = new TH1D("looseMu_DiMu_pt","",100,0,100);
  looseMu_DiMu_pt->Sumw2();
  TH1D* looseMu_DiMu_eta              = new TH1D("looseMu_DiMu_eta","",480,-2.4,2.4);
  looseMu_DiMu_eta->Sumw2();
  TH1D* looseMu_DiMu_phi              = new TH1D("looseMu_DiMu_phi","",640,-3.2,3.2);
  looseMu_DiMu_phi->Sumw2();

  TH1D* M_DiMu_tight = new TH1D("M_DiMu_tight","",1000,50,150); M_DiMu_tight->Sumw2();

  TH1D* tightMu_DiMu_pt               = new TH1D("tightMu_DiMu_pt","",100,0,100);
  tightMu_DiMu_pt->Sumw2();
  TH1D* tightMu_DiMu_eta              = new TH1D("tightMu_DiMu_eta","",480,-2.4,2.4);
  tightMu_DiMu_eta->Sumw2();
  TH1D* tightMu_DiMu_phi              = new TH1D("tightMu_DiMu_phi","",640,-3.2,3.2);

  cout<<"Start Iteration"<<endl;
  // ---------- Variables ---------- //
  double pt = 0; double eta = 0; double phi = 0; double iso = 0;
	double idSFofProbe = 1;
  double weight = 1.0;
  double weightedSum = 0;

  vector<Muon> probeBox, dimuBox, muons, MuonCollection, CorrMuonCollection;
	vector<GenLepton> GenLeptonCollection;
  // ---------- Start Iteration ---------- //
  for(int i=0; i<entries; i++) {
    ntuple->GetEvent(i);

    if(dataSet.Index("data") < 0) {
      if( ntuple->GENEvt_weight < 0 ) weight = -1;
			else weight = 1;
		}

    weightedSum += weight;

    if( ntuple->nMuon < 3 ) continue;

    MET pfMET; pfMET.FillFromNtuple(ntuple,0);
    if(pfMET.Pt > 25) continue;

    if( !(ntuple->isTriggered("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*")
         || ntuple->isTriggered("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")) )
      continue;


  MuonCollection.clear(); CorrMuonCollection.clear();

  for(int j=0; j<ntuple->nMuon; j++) {
    Muon mu; mu.FillFromNtuple( ntuple , j );

    if( mu.Pt < 4 || fabs(mu.eta) > 2.7 || !(mu.isGlobalMuon || mu.isTrackerMuon) ) continue;

    MuonCollection.push_back(mu);

    Double_t rndm[2]; r1->RndmArray(2, rndm);
    double SF = 0;
    Int_t s, m;

    if(dataSet.Index("data") >= 0)
      SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
    else
      SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);

    mu.Pt = SF*mu.Pt;

    TLorentzVector v;
    v.SetPtEtaPhiM( mu.Pt , mu.eta , mu.phi , M_Mu );

    mu.Momentum = v;

    CorrMuonCollection.push_back(mu);
  }



		if( MuonCollection.size() < 3 ) continue;

    // ---------- loose muon selection ---------- //
    dimuBox.clear();
    bool test_Zboson = false;
    int a, b;

    dimuon_selection(CorrMuonCollection, &dimuBox, &test_Zboson, &a, &b);

    if(!test_Zboson) continue;
    double mass = (dimuBox[0].Momentum + dimuBox[1].Momentum).M();

    probeBox.clear();
    probe_muon_selection(CorrMuonCollection, &probeBox, a, b);

    if(probeBox.size() != 1) continue;

    // QCD suppression
		muons.clear();
		muons.insert( muons.end(), dimuBox.begin(), dimuBox.end() );
		muons.push_back( probeBox[0] );
    if( !JPsi_suppression( muons ) ) continue;

    int nTrig = 0;
    if( probeBox[0].isTrigMatched( ntuple, trig ) || probeBox[0].isTrigMatched( ntuple, trig_ ) )
      nTrig++;
    if( dimuBox[0].isTrigMatched( ntuple, trig ) || dimuBox[0].isTrigMatched( ntuple, trig_ ) )
      nTrig++;
    if( dimuBox[1].isTrigMatched( ntuple, trig ) || dimuBox[1].isTrigMatched( ntuple, trig_ ) )
      nTrig++;
    if( nTrig < 2 ) continue;

    pt = probeBox[0].Pt; eta = probeBox[0].eta; phi = probeBox[0].phi;
    iso = second_isolation( CorrMuonCollection, probeBox[0] );

		bool Probe = false;
    if( probeBox[0].testMuonIDv2() && iso < 0.1 ) {
	    if( probeBox[0].isGlobalMuon )
				Probe = true;
			else if( dimuBox[0].Momentum.DeltaR( probeBox[0].Momentum ) < 0.3
	             || dimuBox[1].Momentum.DeltaR( probeBox[0].Momentum ) < 0.3 )
				Probe = true;
	  }

		if( dataSet.Index("data") < 0 ) {
			if( Probe ) 
				idSFofProbe = ( 19.72 * analyzer->MuonSF_BtoF_ID(probeBox[0]) + 16.146 * analyzer->MuonSF_GtoH_ID(probeBox[0]) ) / 35.866;
			else
				idSFofProbe = ( 19.72 * analyzer->LooseMuonSF_BtoF_ID(probeBox[0]) + 16.146 * analyzer->LooseMuonSF_GtoH_ID(probeBox[0]) ) / 35.866;
			weight = weight * idSFofProbe
							 * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_ID( dimuBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_ID( dimuBox ) ) / 35.866;
		}

    M_loose->Fill(mass,weight);
    looseMu_pt->Fill(probeBox[0].Pt,weight);
    looseMu_eta->Fill(probeBox[0].eta,weight);
    looseMu_phi->Fill(probeBox[0].phi,weight);

    denominator_pt->Fill(pt,weight);

    if( fabs(eta)<1.2 ) {
      denominator_pt_barrel->Fill(pt,weight);
    }
    else {
      denominator_pt_endcap->Fill(pt,weight);
    }

    if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 || dataSet.Index("TT") >= 0 ) {
      GenLeptonCollection.clear();
      GenLeptonCollection = DYspliter(ntuple);
      if( GenLeptonCollection.size() >= 4 ) {

        M_FourMu_loose->Fill( mass , weight );

					denominator_FourMu_pt->Fill( pt , weight );
					if( fabs(eta) < 1.2 ) 
						denominator_FourMu_pt_barrel->Fill( pt , weight );
					else
						denominator_FourMu_pt_endcap->Fill( pt , weight );

          looseMu_FourMu_pt->Fill( pt , weight );
          looseMu_FourMu_eta->Fill( eta , weight );
          looseMu_FourMu_phi->Fill( phi , weight );
      }
      else {
        M_DiMu_loose->Fill( mass , weight );

					denominator_DiMu_pt->Fill( pt , weight );
          if( fabs(eta) < 1.2 )
            denominator_DiMu_pt_barrel->Fill( pt , weight );
          else
            denominator_DiMu_pt_endcap->Fill( pt , weight );

          looseMu_DiMu_pt->Fill( pt , weight );
          looseMu_DiMu_eta->Fill( eta , weight );
          looseMu_DiMu_phi->Fill( phi , weight );
      }
    }

    if(!( Probe )) continue;

    M_tight->Fill(mass,weight);
    tightMu_pt->Fill(probeBox[0].Pt,weight);
    tightMu_eta->Fill(probeBox[0].eta,weight);
    tightMu_phi->Fill(probeBox[0].phi,weight);

    numerator_pt->Fill(pt,weight);

    if( fabs(eta)<1.2 ) {
      numerator_pt_barrel->Fill(pt,weight);
    }
    else {
      numerator_pt_endcap->Fill(pt,weight);
    }

    if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 || dataSet.Index("TT") >= 0 ) {
      if( GenLeptonCollection.size() >= 4 ) {

        M_FourMu_tight->Fill( mass , weight );

          numerator_FourMu_pt->Fill( pt , weight );
          if( fabs(eta) < 1.2 )
            numerator_FourMu_pt_barrel->Fill( pt , weight );
          else
            numerator_FourMu_pt_endcap->Fill( pt , weight );

          tightMu_FourMu_pt->Fill( pt , weight );
          tightMu_FourMu_eta->Fill( eta , weight );
          tightMu_FourMu_phi->Fill( phi , weight );
      }
      else {
        M_DiMu_tight->Fill( mass , weight );

          numerator_DiMu_pt->Fill( pt , weight );
          if( fabs(eta) < 1.2 )
            numerator_DiMu_pt_barrel->Fill( pt , weight );
          else
            numerator_DiMu_pt_endcap->Fill( pt , weight );

          tightMu_DiMu_pt->Fill( pt , weight );
          tightMu_DiMu_eta->Fill( eta , weight );
          tightMu_DiMu_phi->Fill( phi , weight );
      }
    }
  }

  TString str_index = to_string(index);

  TFile* f_ = new TFile("root_box/hist"+str_index+".root","RECREATE");

  denominator_pt->Write();
  denominator_pt_barrel->Write();
  denominator_pt_endcap->Write();

  numerator_pt->Write();
  numerator_pt_barrel->Write();
  numerator_pt_endcap->Write();

  denominator_FourMu_pt->Write();
  denominator_FourMu_pt_barrel->Write();
  denominator_FourMu_pt_endcap->Write();

  numerator_FourMu_pt->Write();
  numerator_FourMu_pt_barrel->Write();
  numerator_FourMu_pt_endcap->Write();

  denominator_DiMu_pt->Write();
  denominator_DiMu_pt_barrel->Write();
  denominator_DiMu_pt_endcap->Write();

  numerator_DiMu_pt->Write();
  numerator_DiMu_pt_barrel->Write();
  numerator_DiMu_pt_endcap->Write();

  f_->Close();


  TFile* f = new TFile("root_box/"+dataSet+"_"+str+".root","RECREATE");

  M_loose->Write();
  looseMu_pt->Write();
  looseMu_eta->Write();
  looseMu_phi->Write();
  M_tight->Write();
  tightMu_pt->Write();
  tightMu_eta->Write();
  tightMu_phi->Write();

  M_FourMu_loose->Write();
  looseMu_FourMu_pt->Write();
  looseMu_FourMu_eta->Write();
  looseMu_FourMu_phi->Write();
  M_FourMu_tight->Write();
  tightMu_FourMu_pt->Write();
  tightMu_FourMu_eta->Write();
  tightMu_FourMu_phi->Write();

  M_DiMu_loose->Write();
  looseMu_DiMu_pt->Write();
  looseMu_DiMu_eta->Write();
  looseMu_DiMu_phi->Write();
  M_DiMu_tight->Write();
  tightMu_DiMu_pt->Write();
  tightMu_DiMu_eta->Write();
  tightMu_DiMu_phi->Write();

  f->Close();

  cout<<"RESULT"<<endl;
  cout<<"weighted sum : "<<weightedSum<<endl;

  // ---------- time --------- //
  time();
  end_t = clock();
  cout<<"Running time : "<<((end_t-begin_t)/CLOCKS_PER_SEC)<<endl<<endl;

  cout<<"end"<<endl;
}









