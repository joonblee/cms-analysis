#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TVector2.h>
#include "../Headers/Object.h"
#include "../Headers/NtupleHandle.h"

#include <iostream>

#include "./UsefulFunctions.C"

#include "../Headers/RoCorr/RoccoR.cc"
#include "../Headers/effSF_muon/DYAnalyzer.h"

#define M_Z 91.2

#define trig "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
#define trig_ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"

using namespace std;


void PPPP_selection(vector<Muon> pBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
      muons->clear();
      Muon mu[4];
      for(int j=0; j<pBox.size(); j++) {
        mu[0] = pBox[j];
        for(int k=j+1; k<pBox.size(); k++) {
          mu[1] = pBox[k];
          for(int l=k+1; l<pBox.size(); l++) {
            mu[2] = pBox[l];
            for(int m=l+1; m<pBox.size(); m++) {
              mu[3] = pBox[m];
              if( mu[0].charge + mu[1].charge + mu[2].charge + mu[3].charge == 0 ) {
                int Lead = 0; int Sub = 0;
                for(int ii=0; ii<4; ii++) {
                  if( mu[ii].Pt > 10 ) {
                    Sub++;
                    if( mu[ii].Pt > 20 ) Lead++;
                  }
                }
                if( Lead >= 1 && Sub >= 2 ) {
                  *mass = (mu[0].Momentum + mu[1].Momentum + mu[2].Momentum + mu[3].Momentum).M();
                  for(int iMu=0; iMu<4; iMu++) {
                    muons->push_back(mu[iMu]);
                  }
                  *evt_selection = true;
                  break;
                }
              }
            }
            if( *evt_selection ) break;
          }
          if( *evt_selection ) break;
        }
        if( *evt_selection ) break;
      }
}

void PPPF_selection(vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
      muons->clear();
      Muon mu[4];
      for(int j=0; j<pBox.size(); j++) {
        mu[0] = pBox[j];
        for(int k=j+1; k<pBox.size(); k++) {
          mu[1] = pBox[k];
          for(int l=k+1; l<pBox.size(); l++) {
            mu[2] = pBox[l];
            for(int m=0; m<fBox.size(); m++) {
              mu[3] = fBox[m];
              if( mu[0].charge + mu[1].charge + mu[2].charge + mu[3].charge == 0 ) {
                int Lead = 0; int Sub = 0;
                for(int ii=0; ii<4; ii++) {
                  if( mu[ii].Pt > 10 ) {
                    Sub++;
                    if( mu[ii].Pt > 20 ) Lead++;
                  }
                }
                if( Lead >= 1 && Sub >= 2 ) {

                  *mass = (mu[0].Momentum + mu[1].Momentum + mu[2].Momentum + mu[3].Momentum).M();

                  for(int iMu=0; iMu<4; iMu++) {
                    muons->push_back(mu[iMu]);
                  }

                  *evt_selection = true;
                  break;
                }
              }
            }
            if( *evt_selection ) break;
          }
          if( *evt_selection ) break;
        }
        if( *evt_selection ) break;
      }
}

void PPFF_selection(vector<Muon> pBox, vector<Muon> fBox,  double *mass, bool *evt_selection, vector<Muon> *muons) {
  muons->clear();
  Muon mu[4]; 
  for(int j=0; j<pBox.size(); j++) {
    mu[0] = pBox[j];
    for(int k=j+1; k<pBox.size(); k++) {
      mu[1] = pBox[k];
      for(int l=0; l<fBox.size(); l++) {
        mu[2] = fBox[l];
        for(int m=l+1; m<fBox.size(); m++) {
          mu[3] = fBox[m];
          if( mu[0].charge + mu[1].charge + mu[2].charge + mu[3].charge == 0 ) {
            int Lead = 0; int Sub = 0;
            for(int ii=0; ii<4; ii++) {
              if( mu[ii].Pt > 10 ) {
                Sub++;
                if( mu[ii].Pt > 20 ) Lead++;
              }
            }
            if( Lead >= 1 && Sub >= 2 ) {

              *mass = (mu[0].Momentum + mu[1].Momentum + mu[2].Momentum + mu[3].Momentum).M();

              for(int iMu=0; iMu<4; iMu++) {
                muons->push_back(mu[iMu]);
              }

              *evt_selection = true;
              break;
            }
          }
        }
        if( *evt_selection ) break;
      }
      if( *evt_selection ) break;
    }
    if( *evt_selection ) break;
  }
}

void PFFF_selection(vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
  muons->clear();
  Muon mu[4]; 
  for(int j=0; j<pBox.size(); j++) {
    mu[0] = pBox[j];
    for(int k=0; k<fBox.size(); k++) {
      mu[1] = fBox[k];
      for(int l=k+1; l<fBox.size(); l++) {
        mu[2] = fBox[l];
        for(int m=l+1; m<fBox.size(); m++) {
          mu[3] = fBox[m];
          if( mu[0].charge + mu[1].charge + mu[2].charge + mu[3].charge == 0 ) {
            int Lead = 0; int Sub = 0;
            for(int ii=0; ii<4; ii++) {
              if( mu[ii].Pt > 10 ) {
                Sub++;
                if( mu[ii].Pt > 20 ) Lead++;
              }
            }
            if( Lead >= 1 && Sub >= 2 ) {

              *mass = (mu[0].Momentum + mu[1].Momentum + mu[2].Momentum + mu[3].Momentum).M();

              for(int iMu=0; iMu<4; iMu++) {
                muons->push_back(mu[iMu]);
              }

              *evt_selection = true;
              break;
            }
          }
        }
        if( *evt_selection ) break;
      }
      if( *evt_selection ) break;
    }
    if( *evt_selection ) break;
  }
}


void FFFF_selection(vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
  muons->clear();
  Muon mu[4];
  for(int j=0; j<fBox.size(); j++) {
    mu[0] = fBox[j];
    for(int k=j+1; k<fBox.size(); k++) {
      mu[1] = fBox[k];
      for(int l=k+1; l<fBox.size(); l++) {
        mu[2] = fBox[l];
        for(int m=l+1; m<fBox.size(); m++) {
          mu[3] = fBox[m];
          if( mu[0].charge + mu[1].charge + mu[2].charge + mu[3].charge == 0 ) {
            int Lead = 0; int Sub = 0;
            for(int ii=0; ii<4; ii++) {
              if( mu[ii].Pt > 10 ) {
                Sub++;
                if( mu[ii].Pt > 20 ) Lead++;
              }
            }
            if( Lead >= 1 && Sub >= 2 ) {

              *mass = (mu[0].Momentum + mu[1].Momentum + mu[2].Momentum + mu[3].Momentum).M();

              for(int iMu=0; iMu<4; iMu++) {
                muons->push_back(mu[iMu]);
              }

              *evt_selection = true;
              break;
            }
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

  else if(dataSet=="H200A50" || dataSet=="H300A50" || dataSet=="H400A50" || dataSet=="H500A50" || dataSet=="H750A50" || dataSet=="H1000A50" || dataSet=="H1250A50" || dataSet=="H1500A50" || dataSet=="H1750A50" || dataSet=="H2000A50") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/AutoProdForLimitCalc/ntuple_skim_"+dataSet+".root";
  }
  else if(dataSet=="H200A1" || dataSet=="H300A1" || dataSet=="H400A1" || dataSet=="H500A1" || dataSet=="H750A1" || dataSet=="H1000A1" || dataSet=="H1250A1" || dataSet=="H1500A1" || dataSet=="H1750A1" || dataSet=="H2000A1") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/A1/ntuple_skim_"+dataSet+".root";
  }
  else if(dataSet=="zp200phi50" || dataSet=="zp300phi50" || dataSet=="zp400phi50" || dataSet=="zp500phi50" || dataSet=="zp750phi50" || dataSet=="zp1000phi50" || dataSet=="zp1250phi50" || dataSet=="zp1500phi50" || dataSet=="zp1750phi50" || dataSet=="zp2000phi50") {

    TString s = "";
    int i1, i2;
    i1 = dataSet.Index("zp") + 2;
    i2 = dataSet.Index("phi");

    for(int i=i1; i<i2; i++) {
      s += dataSet[i];
    }

    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/zp_phi50/ntuple_skim_"+s+"_50_v1.root";
    *add_index = 9999;
  }
  else if(dataSet=="zp200phi1" || dataSet=="zp300phi1" || dataSet=="zp400phi1" || dataSet=="zp500phi1" || dataSet=="zp750phi1" || dataSet=="zp1000phi1" || dataSet=="zp1250phi1" || dataSet=="zp1500phi1" || dataSet=="zp1750phi1" || dataSet=="zp2000phi1") {

    TString s = "";
    int i1, i2;
    i1 = dataSet.Index("zp") + 2;
    i2 = dataSet.Index("phi");

    for(int i=i1; i<i2; i++) {
      s += dataSet[i];
    }

    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/zp_phi1/ntuple_skim_"+s+"_1_v1.root";
    *add_index = 9999;
  }

  else {
    cout<<"wrong sample"<<endl;
    str_chain = "ERROR";
  }

  return str_chain;
}


void applyFR(TString dataSet, int input = 0) {
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
  if(dataSet.Index("data")!=0) {
    ntuple->TurnOnBranches_GenLepton();
		ntuple->TurnOnBranches_GenOthers();
	}

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
  // ---------- HISTOGRAM ---------- //
  double bins = 30000;
	if( index == 9999 ) bins = 12000;

  TH1D* PPPP_M = new TH1D("PPPP_M","",bins,0,3000); PPPP_M->Sumw2();
  TH1D* PPPP_FourMu = new TH1D("PPPP_FourMu","",bins,0,3000); PPPP_FourMu->Sumw2();
  TH1D* PPPP_DiMu = new TH1D("PPPP_DiMu","",bins,0,3000); PPPP_DiMu->Sumw2();

	///////////
  TH1D* PPPF_M = new TH1D("PPPF_M","",bins,0,3000); PPPF_M->Sumw2();
  TH1D* PPPF_M_ReWt = new TH1D("PPPF_M_ReWt","",bins,0,3000); PPPF_M_ReWt->Sumw2();
  TH1D* PPPF_FourMu = new TH1D("PPPF_FourMu","",bins,0,3000); PPPF_FourMu->Sumw2();
  TH1D* PPPF_DiMu = new TH1D("PPPF_DiMu","",bins,0,3000); PPPF_DiMu->Sumw2();
  TH1D* PPPF_FourMu_ReWt = new TH1D("PPPF_FourMu_ReWt","",bins,0,3000); PPPF_FourMu_ReWt->Sumw2();
  TH1D* PPPF_DiMu_ReWt = new TH1D("PPPF_DiMu_ReWt","",bins,0,3000); PPPF_DiMu_ReWt->Sumw2();

  TH1D* PPFF_M = new TH1D("PPFF_M","",bins,0,3000); PPFF_M->Sumw2();
  TH1D* PPFF_M_ReWt = new TH1D("PPFF_M_ReWt","",bins,0,3000); PPFF_M_ReWt->Sumw2();
  TH1D* PPFF_FourMu = new TH1D("PPFF_FourMu","",bins,0,3000); PPFF_FourMu->Sumw2();
  TH1D* PPFF_DiMu = new TH1D("PPFF_DiMu","",bins,0,3000); PPFF_DiMu->Sumw2();
  TH1D* PPFF_FourMu_ReWt = new TH1D("PPFF_FourMu_ReWt","",bins,0,3000); PPFF_FourMu_ReWt->Sumw2();
  TH1D* PPFF_DiMu_ReWt = new TH1D("PPFF_DiMu_ReWt","",bins,0,3000); PPFF_DiMu_ReWt->Sumw2();

  TH1D* PFFF_M = new TH1D("PFFF_M","",bins,0,3000); PFFF_M->Sumw2();
  TH1D* PFFF_M_ReWt = new TH1D("PFFF_M_ReWt","",bins,0,3000); PFFF_M_ReWt->Sumw2();

  TH1D* FFFF_M = new TH1D("FFFF_M","",bins,0,3000); FFFF_M->Sumw2();
  TH1D* FFFF_M_ReWt = new TH1D("FFFF_M_ReWt","",bins,0,3000); FFFF_M_ReWt->Sumw2();

	// R < 0.3 ////////////////////////////////////////////////////////////////////
  TH1D* PPPP_M_LR03 = new TH1D("PPPP_M_LR03","",bins,0,3000); PPPP_M_LR03->Sumw2();
  TH1D* PPPP_FourMu_LR03 = new TH1D("PPPP_FourMu_LR03","",bins,0,3000); PPPP_FourMu_LR03->Sumw2();
  TH1D* PPPP_DiMu_LR03 = new TH1D("PPPP_DiMu_LR03","",bins,0,3000); PPPP_DiMu_LR03->Sumw2();

  TH1D* PPPF_M_LR03 = new TH1D("PPPF_M_LR03","",bins,0,3000); PPPF_M_LR03->Sumw2();
  TH1D* PPPF_M_LR03_ReWt = new TH1D("PPPF_M_LR03_ReWt","",bins,0,3000); PPPF_M_LR03_ReWt->Sumw2();
  TH1D* PPPF_FourMu_LR03 = new TH1D("PPPF_FourMu_LR03","",bins,0,3000); PPPF_FourMu_LR03->Sumw2();
  TH1D* PPPF_DiMu_LR03 = new TH1D("PPPF_DiMu_LR03","",bins,0,3000); PPPF_DiMu_LR03->Sumw2();
  TH1D* PPPF_FourMu_LR03_ReWt = new TH1D("PPPF_FourMu_LR03_ReWt","",bins,0,3000); PPPF_FourMu_LR03_ReWt->Sumw2();
  TH1D* PPPF_DiMu_LR03_ReWt = new TH1D("PPPF_DiMu_LR03_ReWt","",bins,0,3000); PPPF_DiMu_LR03_ReWt->Sumw2();

  TH1D* PPFF_M_LR03 = new TH1D("PPFF_M_LR03","",bins,0,3000); PPFF_M_LR03->Sumw2();
  TH1D* PPFF_M_LR03_ReWt = new TH1D("PPFF_M_LR03_ReWt","",bins,0,3000); PPFF_M_LR03_ReWt->Sumw2();
  TH1D* PPFF_FourMu_LR03 = new TH1D("PPFF_FourMu_LR03","",bins,0,3000); PPFF_FourMu_LR03->Sumw2();
  TH1D* PPFF_DiMu_LR03 = new TH1D("PPFF_DiMu_LR03","",bins,0,3000); PPFF_DiMu_LR03->Sumw2();
  TH1D* PPFF_FourMu_LR03_ReWt = new TH1D("PPFF_FourMu_LR03_ReWt","",bins,0,3000); PPFF_FourMu_LR03_ReWt->Sumw2();
  TH1D* PPFF_DiMu_LR03_ReWt = new TH1D("PPFF_DiMu_LR03_ReWt","",bins,0,3000); PPFF_DiMu_LR03_ReWt->Sumw2();

  TH1D* PFFF_M_LR03 = new TH1D("PFFF_M_LR03","",bins,0,3000); PFFF_M_LR03->Sumw2();
  TH1D* PFFF_M_LR03_ReWt = new TH1D("PFFF_M_LR03_ReWt","",bins,0,3000); PFFF_M_LR03_ReWt->Sumw2();

  TH1D* FFFF_M_LR03 = new TH1D("FFFF_M_LR03","",bins,0,3000); FFFF_M_LR03->Sumw2();
  TH1D* FFFF_M_LR03_ReWt = new TH1D("FFFF_M_LR03_ReWt","",bins,0,3000); FFFF_M_LR03_ReWt->Sumw2();

	///////////////////////////////////////////////////////////////////////////////
  TH1D* PPPP_M_UR03 = new TH1D("PPPP_M_UR03","",bins,0,3000); PPPP_M_UR03->Sumw2();
  TH1D* PPPP_FourMu_UR03 = new TH1D("PPPP_FourMu_UR03","",bins,0,3000); PPPP_FourMu_UR03->Sumw2();
  TH1D* PPPP_DiMu_UR03 = new TH1D("PPPP_DiMu_UR03","",bins,0,3000); PPPP_DiMu_UR03->Sumw2();

  TH1D* PPPF_M_UR03 = new TH1D("PPPF_M_UR03","",bins,0,3000); PPPF_M_UR03->Sumw2();
  TH1D* PPPF_M_UR03_ReWt = new TH1D("PPPF_M_UR03_ReWt","",bins,0,3000); PPPF_M_UR03_ReWt->Sumw2();
  TH1D* PPPF_FourMu_UR03 = new TH1D("PPPF_FourMu_UR03","",bins,0,3000); PPPF_FourMu_UR03->Sumw2();
  TH1D* PPPF_DiMu_UR03 = new TH1D("PPPF_DiMu_UR03","",bins,0,3000); PPPF_DiMu_UR03->Sumw2();
  TH1D* PPPF_FourMu_UR03_ReWt = new TH1D("PPPF_FourMu_UR03_ReWt","",bins,0,3000); PPPF_FourMu_UR03_ReWt->Sumw2();
  TH1D* PPPF_DiMu_UR03_ReWt = new TH1D("PPPF_DiMu_UR03_ReWt","",bins,0,3000); PPPF_DiMu_UR03_ReWt->Sumw2();

  TH1D* PPFF_M_UR03 = new TH1D("PPFF_M_UR03","",bins,0,3000); PPFF_M_UR03->Sumw2();
  TH1D* PPFF_M_UR03_ReWt = new TH1D("PPFF_M_UR03_ReWt","",bins,0,3000); PPFF_M_UR03_ReWt->Sumw2();
  TH1D* PPFF_FourMu_UR03 = new TH1D("PPFF_FourMu_UR03","",bins,0,3000); PPFF_FourMu_UR03->Sumw2();
  TH1D* PPFF_DiMu_UR03 = new TH1D("PPFF_DiMu_UR03","",bins,0,3000); PPFF_DiMu_UR03->Sumw2();
  TH1D* PPFF_FourMu_UR03_ReWt = new TH1D("PPFF_FourMu_UR03_ReWt","",bins,0,3000); PPFF_FourMu_UR03_ReWt->Sumw2();
  TH1D* PPFF_DiMu_UR03_ReWt = new TH1D("PPFF_DiMu_UR03_ReWt","",bins,0,3000); PPFF_DiMu_UR03_ReWt->Sumw2();

  TH1D* PFFF_M_UR03 = new TH1D("PFFF_M_UR03","",bins,0,3000); PFFF_M_UR03->Sumw2();
  TH1D* PFFF_M_UR03_ReWt = new TH1D("PFFF_M_UR03_ReWt","",bins,0,3000); PFFF_M_UR03_ReWt->Sumw2();

  TH1D* FFFF_M_UR03 = new TH1D("FFFF_M_UR03","",bins,0,3000); FFFF_M_UR03->Sumw2();
  TH1D* FFFF_M_UR03_ReWt = new TH1D("FFFF_M_UR03_ReWt","",bins,0,3000); FFFF_M_UR03_ReWt->Sumw2();

	// ---------- Setting ---------- //
	cout<<setprecision(5);
	cout<<right;
	//cout<<showpoint;

  // ---------- Variables ---------- //
  double pt = 0; double eta = 0; double iso = 0;
  double wt = 1.0;
  double wtSum = 0;
  double FR, bFR;
  double mass, Z1, Z2, Z1_, Z2_;

  int passingEVT = 0;

  vector<Muon> glbBox, trkBox, fBox, pBox, lBox, muons, MuonCollection, CorrMuonCollection;
  vector<GenLepton> GenLeptonCollection;
	vector<Muon> Mu, Mu_;

  // ---------- Start Iteration ---------- //
  for(int i=0; i<entries; i++) {
    ntuple->GetEvent(i);

    if(dataSet.Index("data") < 0)
      ntuple->GENEvt_weight < 0 ? wt = -1 : wt = 1;

    wtSum += wt;

    if( ntuple->nMuon < 4 ) continue;

    if( !(ntuple->isTriggered( trig )
         || ntuple->isTriggered( trig_ )) )
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

	if( MuonCollection.size() < 4 ) continue;

    pBox.clear(); fBox.clear(); glbBox.clear(); trkBox.clear(); lBox.clear();
    for(int j=0; j<CorrMuonCollection.size(); j++) {
      Muon mu = CorrMuonCollection[j];
      if( mu.acceptance( 5 , 2.4 )
          && mu.looseMuonID() ) {
        lBox.push_back( mu );
        if( mu.isGlobalMuon
            && mu.testMuonIDv2()
            && second_isolation( CorrMuonCollection , mu ) < 0.1 ) 
          glbBox.push_back( mu );
        else if( !mu.isGlobalMuon
                 && mu.testMuonIDv2()
                 && second_isolation( CorrMuonCollection , mu ) < 0.1 )
          trkBox.push_back( mu );
        else fBox.push_back( mu );
      }
    }

    if( lBox.size() != 4 ) continue;

    if( !JPsi_suppression( lBox ) ) continue;

		Mu.clear(); Mu_.clear();
		for(int j=0; j<lBox.size(); j++) {
			if( lBox[j].charge < 0 ) Mu.push_back(lBox[j]);
			else Mu_.push_back(lBox[j]);
		}

    Z1 = (Mu[0].Momentum+Mu_[0].Momentum).M();
    Z2 = (Mu[1].Momentum+Mu_[1].Momentum).M();
    Z1_ = (Mu[0].Momentum+Mu_[1].Momentum).M();
    Z2_ = (Mu[1].Momentum+Mu_[0].Momentum).M();

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
    if( pBox.size() == 4 ) {
      PPPP_selection(pBox, &mass, &evt_selection, &muons);

      if( evt_selection ) {
        if( dataSet.Index("data") < 0 )
					wt = wt * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_ID( pBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_ID( pBox ) ) / 35.866;

        PPPP_M->Fill(mass,wt);

        bool TBoost = boost_selection(pBox, fBox);
        if( TBoost ) PPPP_M_LR03->Fill(mass,wt);
				else PPPP_M_UR03->Fill(mass,wt);

        if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 ) {
          GenLeptonCollection.clear();
          GenLeptonCollection = DYspliter(ntuple);

          if( GenLeptonCollection.size() >= 4 ) {
            PPPP_FourMu->Fill(mass,wt);
            if( TBoost ) PPPP_FourMu_LR03->Fill(mass,wt); // photon conversion
						else PPPP_FourMu_UR03->Fill(mass,wt);
          }
          else if( GenLeptonCollection.size() >= 2 ) {
            PPPP_DiMu->Fill(mass,wt);
            if( TBoost ) PPPP_DiMu_LR03->Fill(mass,wt);   // di-mu from a merged jet
						else PPPP_DiMu_UR03->Fill(mass,wt);
          }
        }
      }
    }
    else if( pBox.size() == 3 ) { // 3P1F //
      FR = 1; bFR = 1;
      PPPF_selection(pBox, fBox, &mass, &evt_selection, &muons);
      if( evt_selection ) {
		    if( dataSet.Index("data") < 0 ) 
	        wt = wt * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_ID( pBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_ID( pBox ) ) / 35.866 * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_LooseID( fBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_LooseID( fBox ) ) / 35.866;

				FR = FR_template( fBox );
        PPPF_M->Fill(mass,wt);
        PPPF_M_ReWt->Fill(mass,wt*FR);

        bool TBoost = boost_selection(pBox, fBox);

        if( TBoost ) {
					bFR = Boosted_FR_template( fBox );
          PPPF_M_LR03->Fill(mass,wt);
          PPPF_M_LR03_ReWt->Fill(mass,wt*bFR);
        }
				else {
          PPPF_M_UR03->Fill(mass,wt);
          PPPF_M_UR03_ReWt->Fill(mass,wt*FR);
        }

        if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 ) {
          GenLeptonCollection.clear();
          GenLeptonCollection = DYspliter(ntuple);

          if( GenLeptonCollection.size() >= 4 ) {
            PPPF_FourMu->Fill(mass,wt);
            PPPF_FourMu_ReWt->Fill(mass,wt*FR);
            if( TBoost ) {
              PPPF_FourMu_LR03->Fill(mass,wt);
              PPPF_FourMu_LR03_ReWt->Fill(mass,wt*bFR);
            }
            else {
              PPPF_FourMu_UR03->Fill(mass,wt);
              PPPF_FourMu_UR03_ReWt->Fill(mass,wt*FR);
            }
          }
          else if( GenLeptonCollection.size() >= 2 ) {
            PPPF_DiMu->Fill(mass,wt);
            PPPF_DiMu_ReWt->Fill(mass,wt*FR);
            if( TBoost ) {
              PPPF_DiMu_LR03->Fill(mass,wt);
              PPPF_DiMu_LR03_ReWt->Fill(mass,wt*bFR);
            }
						else {
              PPPF_DiMu_UR03->Fill(mass,wt);
              PPPF_DiMu_UR03_ReWt->Fill(mass,wt*FR);
            }
          }
        }

      }
    }
    else if( pBox.size() == 2 ) { // 2P2F //
      FR = 1; bFR = 1;
      PPFF_selection(pBox, fBox, &mass, &evt_selection, &muons);
      if( evt_selection ) {
        if( dataSet.Index("data") < 0 )
        	wt = wt * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_ID( pBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_ID( pBox ) ) / 35.866 * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_LooseID( fBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_LooseID( fBox ) ) / 35.866;

				FR = FR_template( fBox );
        PPFF_M->Fill(mass,wt);
        PPFF_M_ReWt->Fill(mass,wt*FR);

        bool TBoost = boost_selection(pBox, fBox);
        if( TBoost ) {
					bFR = Boosted_FR_template( fBox );
          PPFF_M_LR03->Fill(mass,wt);
          PPFF_M_LR03_ReWt->Fill(mass,wt*bFR);
        }
				else {
          PPFF_M_UR03->Fill(mass,wt);
          PPFF_M_UR03_ReWt->Fill(mass,wt*FR);
        }

        if( dataSet.Index("DY") >= 0 || dataSet.Index("ZG") >= 0 ) {
          GenLeptonCollection.clear();
          GenLeptonCollection = DYspliter(ntuple);

          if( GenLeptonCollection.size() >= 4 ) {
            PPFF_FourMu->Fill(mass,wt);
            PPFF_FourMu_ReWt->Fill(mass,wt*FR);
            if( TBoost ) {
              PPFF_FourMu_LR03->Fill(mass,wt);
              PPFF_FourMu_LR03_ReWt->Fill(mass,wt*bFR);
            }
						else {
              PPFF_FourMu_UR03->Fill(mass,wt);
              PPFF_FourMu_UR03_ReWt->Fill(mass,wt*FR);
            }
          }
          else if( GenLeptonCollection.size() >= 2 ) {
            PPFF_DiMu->Fill(mass,wt);
            PPFF_DiMu_ReWt->Fill(mass,wt*FR);
            if( TBoost ) {
              PPFF_DiMu_LR03->Fill(mass,wt);
              PPFF_DiMu_LR03_ReWt->Fill(mass,wt*bFR);
            }
						else {
              PPFF_DiMu_UR03->Fill(mass,wt);
              PPFF_DiMu_UR03_ReWt->Fill(mass,wt*FR);
            }
          }
        }
      }
    }
    else if( pBox.size() == 1 ) { // 1P3F //
      FR = 1; bFR = 1;
      PFFF_selection(pBox, fBox, &mass, &evt_selection, &muons);
      if( evt_selection ) {
        if( dataSet.Index("data") < 0 )
        	wt = wt * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_ID( pBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_ID( pBox ) ) / 35.866 * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_LooseID( fBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_LooseID( fBox ) ) / 35.866;

				FR = FR_template( fBox );
        PFFF_M->Fill(mass,wt);
        PFFF_M_ReWt->Fill(mass,wt*FR);

        bool TBoost = boost_selection(pBox, fBox);
        if( TBoost ) {
					bFR = Boosted_FR_template( fBox );
          PFFF_M_LR03->Fill(mass,wt);
          PFFF_M_LR03_ReWt->Fill(mass,wt*bFR);
        }
				else {
          PFFF_M_UR03->Fill(mass,wt);
          PFFF_M_UR03_ReWt->Fill(mass,wt*FR);
        }
      }
    }
    else if( pBox.size() == 0 ) {
      FR = 1; bFR = 1;
      FFFF_selection(pBox, fBox, &mass, &evt_selection, &muons);
      if( evt_selection ) {
        if( dataSet.Index("data") < 0 )
        	wt = wt * ( 19.72 * analyzer->EfficiencySF_EventWeight_HLT_BtoF_LooseID( fBox ) + 16.146 * analyzer->EfficiencySF_EventWeight_HLT_GtoH_LooseID( fBox ) ) / 35.866;

				FR = FR_template( fBox );
        FFFF_M->Fill(mass,wt);
        FFFF_M_ReWt->Fill(mass,wt*FR);

        bool TBoost = boost_selection(pBox, fBox);
        if( TBoost ) {
					bFR = Boosted_FR_template( fBox );
          FFFF_M_LR03->Fill(mass,wt);
          FFFF_M_LR03_ReWt->Fill(mass,wt*bFR);
        }
				else {
          FFFF_M_UR03->Fill(mass,wt);
          FFFF_M_UR03_ReWt->Fill(mass,wt*FR);
        }
      } 
    } 
    else cout<<"Fatal ERROR"<<endl;
  }

  TString str_index = to_string(index);

  TFile* f = new TFile("final_result/FR_"+dataSet+"_"+str+".root","RECREATE");

  PPPP_M->Write();
  PPPP_FourMu->Write();
  PPPP_DiMu->Write();

  PPPF_M->Write();
  PPPF_M_ReWt->Write();
  PPPF_FourMu->Write();
  PPPF_DiMu->Write();
  PPPF_FourMu_ReWt->Write();
  PPPF_DiMu_ReWt->Write();

  PPFF_M->Write();
  PPFF_M_ReWt->Write();
  PPFF_FourMu->Write();
  PPFF_DiMu->Write();
  PPFF_FourMu_ReWt->Write();
  PPFF_DiMu_ReWt->Write();

  PFFF_M->Write();
  PFFF_M_ReWt->Write();

  FFFF_M->Write();
  FFFF_M_ReWt->Write();

  PPPP_M_LR03->Write();
  PPPP_FourMu_LR03->Write();
  PPPP_DiMu_LR03->Write();

  PPPF_M_LR03->Write();
  PPPF_M_LR03_ReWt->Write();
  PPPF_FourMu_LR03->Write();
  PPPF_DiMu_LR03->Write();
  PPPF_FourMu_LR03_ReWt->Write();
  PPPF_DiMu_LR03_ReWt->Write();

  PPFF_M_LR03->Write();
  PPFF_M_LR03_ReWt->Write();
  PPFF_FourMu_LR03->Write();
  PPFF_DiMu_LR03->Write();
  PPFF_FourMu_LR03_ReWt->Write();
  PPFF_DiMu_LR03_ReWt->Write();

  PFFF_M_LR03->Write();
  PFFF_M_LR03_ReWt->Write();

  FFFF_M_LR03->Write();
  FFFF_M_LR03_ReWt->Write();

  PPPP_M_UR03->Write();
  PPPP_FourMu_UR03->Write();
  PPPP_DiMu_UR03->Write();

  PPPF_M_UR03->Write();
  PPPF_M_UR03_ReWt->Write();
  PPPF_FourMu_UR03->Write();
  PPPF_DiMu_UR03->Write();
  PPPF_FourMu_UR03_ReWt->Write();
  PPPF_DiMu_UR03_ReWt->Write();

  PPFF_M_UR03->Write();
  PPFF_M_UR03_ReWt->Write();
  PPFF_FourMu_UR03->Write();
  PPFF_DiMu_UR03->Write();
  PPFF_FourMu_UR03_ReWt->Write();
  PPFF_DiMu_UR03_ReWt->Write();

  PFFF_M_UR03->Write();
  PFFF_M_UR03_ReWt->Write();

  FFFF_M_UR03->Write();
  FFFF_M_UR03_ReWt->Write();

	f->Close();

  cout<<"- RESULT -"<<endl;
  cout<<"weighted sum : "<<wtSum<<endl;

  // ---------- time --------- //
  time();
  end_t = clock();
  cout<<"Running time : "<<((end_t-begin_t)/CLOCKS_PER_SEC)<<endl<<endl;

  cout<<"end"<<endl;
}

