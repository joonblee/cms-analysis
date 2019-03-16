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


void PPPP_selection(NtupleHandle *ntuple, vector<Muon> pBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
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

void PPPF_selection(NtupleHandle *ntuple, vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
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

void PPFF_selection(NtupleHandle *ntuple, vector<Muon> pBox, vector<Muon> fBox,  double *mass, bool *evt_selection, vector<Muon> *muons) {
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

void PFFF_selection(NtupleHandle *ntuple, vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
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


void FFFF_selection(NtupleHandle *ntuple, vector<Muon> pBox, vector<Muon> fBox, double *mass, bool *evt_selection, vector<Muon> *muons) {
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

  if(dataSet=="H200A50" || dataSet=="H300A50" || dataSet=="H400A50" || dataSet=="H500A50" || dataSet=="H750A50" || dataSet=="H1000A50" || dataSet=="H1250A50" || dataSet=="H1500A50" || dataSet=="H1750A50" || dataSet=="H2000A50") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/AutoProdForLimitCalc/ntuple_skim_"+dataSet+".root";
		*add_index = 9999;
  }
  else if(dataSet=="H200A1" || dataSet=="H300A1" || dataSet=="H400A1" || dataSet=="H500A1" || dataSet=="H750A1" || dataSet=="H1000A1" || dataSet=="H1250A1" || dataSet=="H1500A1" || dataSet=="H1750A1" || dataSet=="H2000A1") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/A1/ntuple_skim_"+dataSet+".root";
		*add_index = 9999;
  }
  else if(dataSet=="H200A10" || dataSet=="H300A10" || dataSet=="H400A10" || dataSet=="H500A10" || dataSet=="H750A10" || dataSet=="H1000A10" || dataSet=="H1250A10" || dataSet=="H1500A10" || dataSet=="H1750A10" || dataSet=="H2000A10") {
    str_chain = "/u/user/joonblee/fs/MET_ntuple/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/withEGMcorrection/A10/ntuple_skim_"+dataSet+".root";
		*add_index = 9999;
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


void signalRate(TString dataSet, int input = 0) {
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

	// ---------- Setting ---------- //
	cout<<setprecision(5);
	cout<<right;
	//cout<<showpoint;

  // ---------- Variables ---------- //
  double pt = 0; double eta = 0; double iso = 0;
  double wt = 1.0;
  double wtSum = 0;
  double FR, bFR;
  double mass;

  int passingEVT = 0;

  vector<Muon> glbBox, trkBox, fBox, pBox, lBox, muons;
  vector<GenLepton> GenLeptonCollection;

	double nFourMuEvt = 0;
	double nPassEvt = 0;
  // ---------- Start Iteration ---------- //
  for(int i=0; i<entries; i++) {
    ntuple->GetEvent(i);

    int nGenMu = 0;
    for(int j=0; j<ntuple->gnpair; j++) {
      GenLepton gn; gn.FillFromNtuple(ntuple,j);
      if( fabs(gn.ID) == 13 && gn.fromHardProcessFinalState ) nGenMu++;
    }
    if(nGenMu == 4) nFourMuEvt++;


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
          && mu.looseMuonIDv2() ) {
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
    if( pBox.size() == 4 ) {
      PPPP_selection(ntuple, pBox, &mass, &evt_selection, &muons);

      if( evt_selection ) {
        nPassEvt++;
      }
    }
  }

	TH1D* rate = new TH1D("rate","",3000,0,3000); rate->Sumw2();

	TString s = "";
	int i1, i2;
	if( dataSet.Index("H") >= 0 ) {
		i1 = dataSet.Index("H") + 1;
		i2 = dataSet.Index("A");
	}
	else {
		i1 = dataSet.Index("zp") + 2;
		i2 = dataSet.Index("phi");
	}

  for(int i=i1; i<i2; i++) {
    s += dataSet[i];
  }
  int M_H = s.Atoi();

	rate->Fill( M_H , nPassEvt/nFourMuEvt * 35.9 );
	rate->SetMarkerStyle(21);

  TFile* f = new TFile(dataSet+".root","RECREATE");
  rate->Write();
  f->Close();

  cout<<"- RESULT -"<<endl;
  cout<<"weighted sum : "<<wtSum<<endl;
	cout<<"# of 4 muons evt : "<<nFourMuEvt<<endl;
	cout<<"# of pass evt    : "<<nPassEvt<<endl;
	cout<<endl;

	cout<<"signal eff       : "<<nPassEvt/nFourMuEvt<<endl;
	cout<<"signal rate (L*e): "<<nPassEvt/nFourMuEvt*35.9<<endl;

  // ---------- time --------- //
  time();
  end_t = clock();
  cout<<"Running time : "<<((end_t-begin_t)/CLOCKS_PER_SEC)<<endl<<endl;

  cout<<"end"<<endl;
}

