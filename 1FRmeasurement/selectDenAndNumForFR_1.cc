#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TVector2.h>
#include "/u/user/joonblee/Headers/Object.h"
#include "/u/user/joonblee/Headers/NtupleHandle.h"
#include "./UsefulFunctions.C"

#include <iostream>

#define M_Z 91.2


#define trig "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
#define trig_ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"

using namespace std;

void FillHist( Muon mu , double wt , TH1D* hist , TH1D* hist_barrel , TH1D* hist_endcap );

void FakeDiMuon_selection(vector<Muon> pBox, vector<bool> bMuon, double *mass, bool *evt_selection, vector<Muon> *muons, vector<bool> *bPass) {
  muons->clear(); bPass->clear();
  double ZMass = 0;
  for(int j=0; j<pBox.size(); j++) {
    if( !bMuon[j] ) continue;
    if( !pBox[j].isGlobalMuon ) continue;
    for(int k=j+1; k<pBox.size(); k++) {
      if( !bMuon[k] ) continue;
      if( !pBox[k].isGlobalMuon ) continue;

      if( !(pBox[j].charge != pBox[k].charge) ) continue;
      ZMass = (pBox[j].Momentum + pBox[k].Momentum).M();
      if( !(fabs( ZMass - M_Z ) < 30) ) continue;

      for(int l=0; l<pBox.size(); l++) {
        if( l==j || l==k ) continue;
        for(int m=l+1; m<pBox.size(); m++) {
          if( m==j || m==k ) continue;
          if( !( pBox[l].Momentum.DeltaR( pBox[m].Momentum ) < 0.3 ) )
            continue;

          int Lead = 0; int Sub = 0;
          for(int ii=0; ii<4; ii++) {
            if( pBox[ii].Pt > 10 ) {
              Sub++;
              if( pBox[ii].Pt > 20 ) Lead++;
            }
          }
          if( Lead >= 1 && Sub >= 2 ) {
            *mass = (pBox[0].Momentum + pBox[1].Momentum + pBox[2].Momentum + pBox[3].Momentum).M();
            if( pBox[l].Pt > pBox[m].Pt ) {
              muons->push_back( pBox[l] ); muons->push_back( pBox[m] );
              bPass->push_back( bMuon[l] ); bPass->push_back( bMuon[m] );
            }
            else {
              muons->push_back( pBox[m] ); muons->push_back( pBox[l] );
              bPass->push_back( bMuon[m] ); bPass->push_back( bMuon[l] );
            }

            *evt_selection = true;
            break;
          }
        }
        if( *evt_selection ) break;
      }
      if( *evt_selection ) break;
    }
    if( *evt_selection ) break;
  }
}


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

  else if(dataSet == "datatest") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DoubleMuon/crab_DYntuple_v2p0_DoubleMuon_RunD/180227_010154/0000/ntuple_skim_*"+str+".root";
    *add_index += 20 + input;
  }

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
    //str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test02_GenOthers_DYLL_M50toInf/180320_052835/0000/ntuple_skim_*"+str+".root";
    *add_index = 10+input;
  }

  else if(dataSet == "DYlow") {
    TString files;
    if( input == 0 || input == 1 ) {
      if( input == 0 ) files = "[01234]"; else files = "[56789]";
      str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v2p0_DYLL_M10to50_v1/171230_124542/0000/ntuple_skim_*"+files+".root"; // DY
    }
    else if( 1 < input && input < 7 ) {
      if( input == 2 ) files = "[01]";
      else if( input == 3 ) files = "[23]";
      else if( input == 4 ) files = "[45]";
      else if( input == 5 ) files = "[67]";
      else if( input == 6 ) files = "[89]";
      str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v2p0_DYLL_M10to50_v2/171230_124625/0000/ntuple_skim_*"+files+".root";
    }
    else if( 6 < input ) {
      if( input == 7 ) files = "[123]";
      else if( input == 8 ) files = "[456]";
      else if( input == 9 ) files = "[7890]";
      str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v2p0_DYLL_M10to50_ext1v1/171230_124709/0000/ntuple_skim_*"+files+".root";
    }
    *add_index = 30+input;
  }

  else if(dataSet == "DY-M100to200") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test01_DYLL_M100to200/180309_185211/0000/ntuple_skim_*.root";
    *add_index = 1000+input;
  }
  else if(dataSet == "DY-M100to200_ext") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test01_DYLL_M100to200_ext/180309_185402/0000/ntuple_skim_*.root";
    *add_index = 1000+input;
  }
  else if(dataSet == "DY-M200to400") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test01_DYLL_M200to400/180309_185604/0000/ntuple_skim_*.root";
    *add_index = 1000+input;
  }
  else if(dataSet == "DY-M400to500") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test01_DYLL_M400to500/180309_185743/0000/ntuple_skim_*.root";
    *add_index = 1000+input;
  }

  else if(dataSet == "DY-M500to700") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_test01_DYLL_M500to700/180309_185922/0000/ntuple_skim_*.root";
    *add_index = 1000+input;
  }


  else if(dataSet == "WZ") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/WZ_TuneCUETP8M1_13TeV-pythia8/crab_DYntuple_v20171022_EGMCorr_WZ/171230_140159/0000/ntuple_skim_*.root"; // WZ
    *add_index = 7;
  }
  else if(dataSet == "TT") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_DYntuple_v20171022_EGMCorr_TTbar/171230_135613/0000/ntuple_skim_*"+str+".root"; // TT
    *add_index = 20+input;
  }
  else if(dataSet == "TTtest") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_DYntuple_v20171022_EGMCorr_TTbar/171230_135613/0000/ntuple_skim_*"+str+"0.root"; // TT
    *add_index = 20+input;
  }


  else if(dataSet == "ZG") {
    str_chain = "dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/joon/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYntuple_v20171022_EGMCorr_ZGamma/180311_151437/0000/ntuple_skim_*.root"; // TT
    *add_index = 40+input;
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

  else if(dataSet=="H200A50" || dataSet=="H300A50" || dataSet=="H400A50" || dataSet=="H500A50" || dataSet=="H750A50" || dataSet=="H1000A50" || dataSet=="H1250A50" || dataSet=="H1500A50" || dataSet=="H1750A50" || dataSet=="H2000A50") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/AutoProdForLimitCalc/ntuple_skim_"+dataSet+".root";
  }
  else if(dataSet=="H200A1" || dataSet=="H300A1" || dataSet=="H400A1" || dataSet=="H500A1" || dataSet=="H750A1" || dataSet=="H1000A1" || dataSet=="H1250A1" || dataSet=="H1500A1" || dataSet=="H1750A1" || dataSet=="H2000A1") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/A1/ntuple_skim_"+dataSet+".root";
  }
  else if(dataSet=="H200A10" || dataSet=="H300A10" || dataSet=="H400A10" || dataSet=="H500A10" || dataSet=="H750A10" || dataSet=="H1000A10" || dataSet=="H1250A10" || dataSet=="H1500A10" || dataSet=="H1750A10" || dataSet=="H2000A10") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/A10/ntuple_skim_"+dataSet+".root";
  }

  else {
    cout<<"wrong sample"<<endl;
    str_chain = "ERROR";
  }

  return str_chain;
}


void selectDenAndNumForFR_1(TString dataSet, int input = 0) {
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
  if(dataSet.First("data")!=0) {
    ntuple->TurnOnBranches_GenLepton();
    ntuple->TurnOnBranches_GenOthers();
  }

  int entries = chain->GetEntries();
  // ---------- HISTOGRAM ---------- //
	/*
  const int ptbinnum = 7;
  double ptbin[ptbinnum+1] = {5,7,10,15,22,40,100,500};
  const int ptbinnum_endcap = 7;
  double ptbin_endcap[ptbinnum_endcap+1] = {5,7,10,15,22,40,100,500};
	*/

  const int ptbinnum = 23;
  double ptbin[ptbinnum+1] = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,100,500};
  const int ptbinnum_endcap = 23;
  double ptbin_endcap[ptbinnum_endcap+1] = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,100,500};


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

  TH1D* ldenominator_pt = new TH1D("ldenominator_pt","",500,0,500);
  TH1D* ldenominator_pt_barrel = new TH1D("ldenominator_pt_barrel","",500,0,500);
  TH1D* ldenominator_pt_endcap = new TH1D("ldenominator_pt_endcap","",500,0,500);

  TH1D* lnumerator_pt = new TH1D("lnumerator_pt","",500,0,500);
  TH1D* lnumerator_pt_barrel = new TH1D("lnumerator_pt_barrel","",500,0,500);
  TH1D* lnumerator_pt_endcap = new TH1D("lnumerator_pt_endcap","",500,0,500);

  ldenominator_pt->Sumw2();
  ldenominator_pt_barrel->Sumw2();
  ldenominator_pt_endcap->Sumw2();

  lnumerator_pt->Sumw2();
  lnumerator_pt_barrel->Sumw2();
  lnumerator_pt_endcap->Sumw2();

  TH1D* sdenominator_pt = new TH1D("sdenominator_pt","",500,0,500);
  TH1D* sdenominator_pt_barrel = new TH1D("sdenominator_pt_barrel","",500,0,500);
  TH1D* sdenominator_pt_endcap = new TH1D("sdenominator_pt_endcap","",500,0,500);

  TH1D* snumerator_pt = new TH1D("snumerator_pt","",500,0,500);
  TH1D* snumerator_pt_barrel = new TH1D("snumerator_pt_barrel","",500,0,500);
  TH1D* snumerator_pt_endcap = new TH1D("snumerator_pt_endcap","",500,0,500);

  sdenominator_pt->Sumw2();
  sdenominator_pt_barrel->Sumw2();
  sdenominator_pt_endcap->Sumw2();

  snumerator_pt->Sumw2();
  snumerator_pt_barrel->Sumw2();
  snumerator_pt_endcap->Sumw2();

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

  TH1D* ldenominator_FourMu_pt = new TH1D("ldenominator_FourMu_pt","",500,0,500);
  TH1D* ldenominator_FourMu_pt_barrel = new TH1D("ldenominator_FourMu_pt_barrel","",500,0,500);
  TH1D* ldenominator_FourMu_pt_endcap = new TH1D("ldenominator_FourMu_pt_endcap","",500,0,500);

  TH1D* lnumerator_FourMu_pt = new TH1D("lnumerator_FourMu_pt","",500,0,500);
  TH1D* lnumerator_FourMu_pt_barrel = new TH1D("lnumerator_FourMu_pt_barrel","",500,0,500);
  TH1D* lnumerator_FourMu_pt_endcap = new TH1D("lnumerator_FourMu_pt_endcap","",500,0,500);

  ldenominator_FourMu_pt->Sumw2();
  ldenominator_FourMu_pt_barrel->Sumw2();
  ldenominator_FourMu_pt_endcap->Sumw2();

  lnumerator_FourMu_pt->Sumw2();
  lnumerator_FourMu_pt_barrel->Sumw2();
  lnumerator_FourMu_pt_endcap->Sumw2();

  TH1D* sdenominator_FourMu_pt = new TH1D("sdenominator_FourMu_pt","",500,0,500);
  TH1D* sdenominator_FourMu_pt_barrel = new TH1D("sdenominator_FourMu_pt_barrel","",500,0,500);
  TH1D* sdenominator_FourMu_pt_endcap = new TH1D("sdenominator_FourMu_pt_endcap","",500,0,500);

  TH1D* snumerator_FourMu_pt = new TH1D("snumerator_FourMu_pt","",500,0,500);
  TH1D* snumerator_FourMu_pt_barrel = new TH1D("snumerator_FourMu_pt_barrel","",500,0,500);
  TH1D* snumerator_FourMu_pt_endcap = new TH1D("snumerator_FourMu_pt_endcap","",500,0,500);

  sdenominator_FourMu_pt->Sumw2();
  sdenominator_FourMu_pt_barrel->Sumw2();
  sdenominator_FourMu_pt_endcap->Sumw2();

  snumerator_FourMu_pt->Sumw2();
  snumerator_FourMu_pt_barrel->Sumw2();
  snumerator_FourMu_pt_endcap->Sumw2();

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

  TH1D* ldenominator_DiMu_pt = new TH1D("ldenominator_DiMu_pt","",500,0,500);
  TH1D* ldenominator_DiMu_pt_barrel = new TH1D("ldenominator_DiMu_pt_barrel","",500,0,500);
  TH1D* ldenominator_DiMu_pt_endcap = new TH1D("ldenominator_DiMu_pt_endcap","",500,0,500);

  TH1D* lnumerator_DiMu_pt = new TH1D("lnumerator_DiMu_pt","",500,0,500);
  TH1D* lnumerator_DiMu_pt_barrel = new TH1D("lnumerator_DiMu_pt_barrel","",500,0,500);
  TH1D* lnumerator_DiMu_pt_endcap = new TH1D("lnumerator_DiMu_pt_endcap","",500,0,500);

  ldenominator_DiMu_pt->Sumw2();
  ldenominator_DiMu_pt_barrel->Sumw2();
  ldenominator_DiMu_pt_endcap->Sumw2();

  lnumerator_DiMu_pt->Sumw2();
  lnumerator_DiMu_pt_barrel->Sumw2();
  lnumerator_DiMu_pt_endcap->Sumw2();

  TH1D* sdenominator_DiMu_pt = new TH1D("sdenominator_DiMu_pt","",500,0,500);
  TH1D* sdenominator_DiMu_pt_barrel = new TH1D("sdenominator_DiMu_pt_barrel","",500,0,500);
  TH1D* sdenominator_DiMu_pt_endcap = new TH1D("sdenominator_DiMu_pt_endcap","",500,0,500);

  TH1D* snumerator_DiMu_pt = new TH1D("snumerator_DiMu_pt","",500,0,500);
  TH1D* snumerator_DiMu_pt_barrel = new TH1D("snumerator_DiMu_pt_barrel","",500,0,500);
  TH1D* snumerator_DiMu_pt_endcap = new TH1D("snumerator_DiMu_pt_endcap","",500,0,500);

  sdenominator_DiMu_pt->Sumw2();
  sdenominator_DiMu_pt_barrel->Sumw2();
  sdenominator_DiMu_pt_endcap->Sumw2();

  snumerator_DiMu_pt->Sumw2();
  snumerator_DiMu_pt_barrel->Sumw2();
  snumerator_DiMu_pt_endcap->Sumw2();


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


  // ---------- Setting ---------- //
  cout<<setprecision(5);
  cout<<right;
  //cout<<showpoint;

  // ---------- Variables ---------- //
  double pt = 0; double eta = 0; double iso = 0;
  double wt = 1.0;
  double wtSum = 0;
  double FR1_template, FR2_template;
  double FR;
  double mass;

  int passingEVT = 0;

  vector<Muon> glbBox, trkBox, fBox, pBox, lBox, muons, vMuon;
  vector<bool> bMuon, bPass;
  vector<GenLepton> GenLeptonCollection;
  // ---------- Start Iteration ---------- //
  for(int i=0; i<entries; i++) {
    ntuple->GetEvent(i);

    if(dataSet.First("data")!=0)
      ntuple->GENEvt_weight < 0 ? wt = -1 : wt = 1;

    wtSum += wt;

    if( ntuple->nMuon < 4 ) continue;

    if( !(ntuple->isTriggered( trig )
         || ntuple->isTriggered( trig_ )) )
      continue;

    pBox.clear(); fBox.clear(); glbBox.clear(); trkBox.clear(); lBox.clear();
    for(int j=0; j<ntuple->nMuon; j++) {
      Muon mu; mu.FillFromNtuple( ntuple , j );
      if( mu.acceptance( 5 , 2.4 )
          && mu.looseMuonID() ) {
        lBox.push_back( mu );
        if( mu.isGlobalMuon
            && mu.testMuonIDv2()
            && second_isolation( ntuple , mu ) < 0.1 ) 
          glbBox.push_back( mu );
        else if( !mu.isGlobalMuon
                 && mu.testMuonIDv2()
                 && second_isolation( ntuple , mu ) < 0.1 )
          trkBox.push_back( mu );
        else fBox.push_back( mu );
      }
    }

    if( lBox.size() != 4 ) continue;

    if( !JPsi_suppression( lBox ) ) continue;

    pBox.insert( pBox.end(), glbBox.begin(), glbBox.end() );

    for(int j=0; j<trkBox.size(); j++) {
      Muon mu = trkBox[j];

      bool TT = false;
      for(int k=0; k<glbBox.size(); k++) {
        Muon mu_ = glbBox[k];

        bool TBoost = false;
        for(int l=0; l<pBox.size(); l++) {
          Muon tmu = pBox[l];
          if( mu_.Pt == tmu.Pt ) continue;
          if( tmu.Momentum.DeltaR( mu_.Momentum ) < 0.3 ) {
            TBoost = true;
            break;
          }
        }
        if( TBoost ) continue;

        if( mu_.Momentum.DeltaR( mu.Momentum ) < 0.3 ) TT = true;
      }
      if( TT ) pBox.push_back( mu );
      else fBox.push_back( mu );    
    }

    bool evt_selection = false;
    mass = 0;
    if( pBox.size() >= 2 ) {
      vMuon.clear(); bMuon.clear();
      vMuon.insert( vMuon.end() , pBox.begin() , pBox.end() );
      vMuon.insert( vMuon.end() , fBox.begin() , fBox.end() );
      for(int j=0; j<pBox.size(); j++) {
        bMuon.push_back( true );
      }
      for(int j=0; j<fBox.size(); j++) {
        bMuon.push_back( false );
      }

      FakeDiMuon_selection(vMuon, bMuon, &mass, &evt_selection, &muons, &bPass);

      if( evt_selection ) {

        FillHist( muons[0] , wt , ldenominator_pt , ldenominator_pt_barrel , ldenominator_pt_endcap );
        FillHist( muons[1] , wt , sdenominator_pt , sdenominator_pt_barrel , sdenominator_pt_endcap );
        M_loose->Fill( mass , wt );
        if( bPass[0] ) 
          FillHist( muons[0] , wt , lnumerator_pt , lnumerator_pt_barrel , lnumerator_pt_endcap );
        if( bPass[1] )
          FillHist( muons[1] , wt , snumerator_pt , snumerator_pt_barrel , snumerator_pt_endcap );
        if( bPass[0] && bPass[1] )
          M_tight->Fill( mass ,wt );

        for(int j=0; j<2; j++) {
          FillHist( muons[j] , wt , denominator_pt , denominator_pt_barrel , denominator_pt_endcap );
          looseMu_pt->Fill( muons[j].Pt , wt );
          looseMu_eta->Fill( muons[j].eta , wt );
          looseMu_phi->Fill( muons[j].phi , wt );

          if( bPass[j] ) {
            FillHist( muons[j] , wt , numerator_pt , numerator_pt_barrel , numerator_pt_endcap );
            tightMu_pt->Fill( muons[j].Pt , wt );
            tightMu_eta->Fill( muons[j].eta , wt );
            tightMu_phi->Fill( muons[j].phi , wt );
          }
        }

        if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 ) {
          GenLeptonCollection.clear();
          GenLeptonCollection = DYspliter(ntuple);

          if( GenLeptonCollection.size() >= 4 ) {

	          FillHist( muons[0] , wt , ldenominator_FourMu_pt , ldenominator_FourMu_pt_barrel , ldenominator_FourMu_pt_endcap );
	          FillHist( muons[1] , wt , sdenominator_FourMu_pt , sdenominator_FourMu_pt_barrel , sdenominator_FourMu_pt_endcap );
	          M_FourMu_loose->Fill( mass , wt );

	          if( bPass[0] ) 
	            FillHist( muons[0] , wt , lnumerator_FourMu_pt , lnumerator_FourMu_pt_barrel , lnumerator_FourMu_pt_endcap );
	          if( bPass[1] ) 
	            FillHist( muons[1] , wt , snumerator_FourMu_pt , snumerator_FourMu_pt_barrel , snumerator_FourMu_pt_endcap );
	          if( bPass[0] && bPass[1] ) 
	            M_FourMu_tight->Fill( mass ,wt );

	          for(int j=0; j<2; j++) {
	            FillHist( muons[j] , wt , denominator_FourMu_pt , denominator_FourMu_pt_barrel , denominator_FourMu_pt_endcap );
	            looseMu_FourMu_pt->Fill( muons[j].Pt , wt );
	            looseMu_FourMu_eta->Fill( muons[j].eta , wt );
	            looseMu_FourMu_phi->Fill( muons[j].phi , wt );

	            if( bPass[j] ) {
	              FillHist( muons[j] , wt , numerator_FourMu_pt , numerator_FourMu_pt_barrel , numerator_FourMu_pt_endcap );
	              tightMu_FourMu_pt->Fill( muons[j].Pt , wt );
	              tightMu_FourMu_eta->Fill( muons[j].eta , wt );
	              tightMu_FourMu_phi->Fill( muons[j].phi , wt );
	            }
	          }
					}
					else {
            FillHist( muons[0] , wt , ldenominator_DiMu_pt , ldenominator_DiMu_pt_barrel , ldenominator_DiMu_pt_endcap );
            FillHist( muons[1] , wt , sdenominator_DiMu_pt , sdenominator_DiMu_pt_barrel , sdenominator_DiMu_pt_endcap );
            M_DiMu_loose->Fill( mass , wt );

            if( bPass[0] ) 
              FillHist( muons[0] , wt , lnumerator_DiMu_pt , lnumerator_DiMu_pt_barrel , lnumerator_DiMu_pt_endcap );
            if( bPass[1] ) 
              FillHist( muons[1] , wt , snumerator_DiMu_pt , snumerator_DiMu_pt_barrel , snumerator_DiMu_pt_endcap );
            if( bPass[0] && bPass[1] ) 
              M_DiMu_tight->Fill( mass ,wt );

            for(int j=0; j<2; j++) {
              FillHist( muons[j] , wt , denominator_DiMu_pt , denominator_DiMu_pt_barrel , denominator_DiMu_pt_endcap );
              looseMu_DiMu_pt->Fill( muons[j].Pt , wt );
              looseMu_DiMu_eta->Fill( muons[j].eta , wt );
              looseMu_DiMu_phi->Fill( muons[j].phi , wt );

              if( bPass[j] ) {
                FillHist( muons[j] , wt , numerator_DiMu_pt , numerator_DiMu_pt_barrel , numerator_DiMu_pt_endcap );
                tightMu_DiMu_pt->Fill( muons[j].Pt , wt );
                tightMu_DiMu_eta->Fill( muons[j].eta , wt );
                tightMu_DiMu_phi->Fill( muons[j].phi , wt );
              }
            }

					}
				}
      }
    }
  }

  TString str_index = to_string(index);

  TFile* f = new TFile("root_box/histB"+str_index+".root","RECREATE");

  denominator_pt->Write();
  denominator_pt_barrel->Write();
  denominator_pt_endcap->Write();

  numerator_pt->Write();
  numerator_pt_barrel->Write();
  numerator_pt_endcap->Write();

  ldenominator_pt->Write();
  ldenominator_pt_barrel->Write();
  ldenominator_pt_endcap->Write();

  lnumerator_pt->Write();
  lnumerator_pt_barrel->Write();
  lnumerator_pt_endcap->Write();

  sdenominator_pt->Write();
  sdenominator_pt_barrel->Write();
  sdenominator_pt_endcap->Write();

  snumerator_pt->Write();
  snumerator_pt_barrel->Write();
  snumerator_pt_endcap->Write();

  denominator_FourMu_pt->Write();
  denominator_FourMu_pt_barrel->Write();
  denominator_FourMu_pt_endcap->Write();

  numerator_FourMu_pt->Write();
  numerator_FourMu_pt_barrel->Write();
  numerator_FourMu_pt_endcap->Write();

  ldenominator_FourMu_pt->Write();
  ldenominator_FourMu_pt_barrel->Write();
  ldenominator_FourMu_pt_endcap->Write();

  lnumerator_FourMu_pt->Write();
  lnumerator_FourMu_pt_barrel->Write();
  lnumerator_FourMu_pt_endcap->Write();

  sdenominator_FourMu_pt->Write();
  sdenominator_FourMu_pt_barrel->Write();
  sdenominator_FourMu_pt_endcap->Write();

  snumerator_FourMu_pt->Write();
  snumerator_FourMu_pt_barrel->Write();
  snumerator_FourMu_pt_endcap->Write();

  denominator_DiMu_pt->Write();
  denominator_DiMu_pt_barrel->Write();
  denominator_DiMu_pt_endcap->Write();

  numerator_DiMu_pt->Write();
  numerator_DiMu_pt_barrel->Write();
  numerator_DiMu_pt_endcap->Write();

  ldenominator_DiMu_pt->Write();
  ldenominator_DiMu_pt_barrel->Write();
  ldenominator_DiMu_pt_endcap->Write();

  lnumerator_DiMu_pt->Write();
  lnumerator_DiMu_pt_barrel->Write();
  lnumerator_DiMu_pt_endcap->Write();

  sdenominator_DiMu_pt->Write();
  sdenominator_DiMu_pt_barrel->Write();
  sdenominator_DiMu_pt_endcap->Write();

  snumerator_DiMu_pt->Write();
  snumerator_DiMu_pt_barrel->Write();
  snumerator_DiMu_pt_endcap->Write();

  f->Close();

  TFile* f_ = new TFile("root_box/B_"+dataSet+"_"+str+".root","RECREATE");

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

  f_->Close();


  cout<<"- RESULT -"<<endl;
  cout<<"weighted sum : "<<wtSum<<endl;

  // ---------- time --------- //
  time();
  end_t = clock();
  cout<<"Running time : "<<((end_t-begin_t)/CLOCKS_PER_SEC)<<endl<<endl;

  cout<<"end"<<endl;
}


void FillHist( Muon mu , double wt , TH1D* hist , TH1D* hist_barrel , TH1D* hist_endcap ) {
  double pt = mu.Pt; double eta = fabs(mu.eta);
  hist->Fill( pt , wt );
  if( eta < 1.2 )
    hist_barrel->Fill( pt , wt );
  else
    hist_endcap->Fill( pt , wt );
}

